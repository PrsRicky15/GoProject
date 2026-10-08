# Expression plotter and scientific calculator

[Documentation index](README.md)

## Responsibilities

The plotter sends formulas and domain settings to Go. Go parses once, evaluates
all requested points, and returns numeric arrays. Plotly renders those arrays in
the browser. The calculator evaluates a small expression in TypeScript without a
backend request. Both parsers implement a restricted arithmetic grammar; neither
runs user input with `eval` or `Function`.

## Modes

| UI / API mode | Meaning | Resolution per axis | Default |
| --- | --- | --- | --- |
| 1D curve / `line` | y=f(x) | 20–10,000 | 801, `exp(-x^2/2)`, x in [-5,5] |
| 2D field / `heatmap` | Color=f(x,y); heatmap or contours | 10–200 | 81, `sin(x)*cos(y)`, both axes [-5,5] |
| 3D surface / `surface` | Height z=f(x,y) | 10–150 | 81, `exp(-(x^2+y^2)/2)`, both axes [-4,4] |
| Isosurface / `iso` | Set f(x,y,z)=level | 10–60 | 31, `x^2+y^2+z^2`, axes [-2,2], level 1 |

For resolution n, sample count is n, n² or n³. At maximum isosurface resolution,
this is 216,000 values; raising each axis resolution has a cubic total cost.
The “3D surface” mode samples a function of two variables. The isosurface mode
samples a three-variable scalar field. Neither invokes a quantum solver.

## Arithmetic language

| Syntax | Meaning |
| --- | --- |
| `+ - * /` | Arithmetic; multiplication must be explicit |
| `^`, `**` | Power; right-associative |
| `%` | Remainder, not percentage |
| Parentheses | Group expressions |
| `pi`, `π`, `e` | Constants |
| `x`, `y`, `z` | Allowed only for active plot dimensions |
| `ans` | Calculator's previous answer only |
| Decimal/scientific literals | Examples: `.5`, `1.2`, `2e-3` |

Function names are case-sensitive. Supported one-argument functions are `sin`,
`cos`, `tan`, `asin`, `acos`, `atan`, `sinh`, `cosh`, `tanh`, `sqrt`, `abs`, `exp`,
`ln`, `log`, `log10`, `floor`, `ceil`, and `round`. Two-argument functions are
`pow(a,b)`, `min(a,b)`, and `max(a,b)`.

Angles are radians. `ln` is natural logarithm; `log` and `log10` are base 10.
Rounding ties go toward positive infinity, matching the calculator's convention:
`round(-1.5)=-1`. Power binds more strongly than unary sign: `-2^2=-4`, while
`(-2)^2=4`. `2^3^2=512`. No implicit multiplication, assignments, property access,
arrays, arbitrary function definitions, integration, symbolic algebra or units
conversion is supported.

Go caps the source at 1,024 bytes and 256 tokens. The TypeScript parser has its own
corresponding bounds and shares arithmetic fixtures with Go. Calculator input is
also constrained by its text field. Non-finite calculator answers are errors.

## Sampling and array orientation

Expression sampling includes both endpoints:

```text
coord[i] = min + i*(max-min)/(n-1), i=0,...,n-1
```

The last coordinate is explicitly set to max. Bounds must be finite/increasing,
and the step must be representable at the coordinate magnitude.
This is different from the endpoint-exclusive quantum grid.

The flat value array uses x as the fastest-changing index:

```text
1D: i = ix
2D: i = ix + n*iy
3D: i = ix + n*iy + n*n*iz
```

Heatmap/surface rows correspond to y; columns correspond to x. `figure.ts`
reshapes rows for Plotly. For isosurfaces it expands axis vectors into per-point
coordinates. `math.test.ts` checks that the renderer retains the Go orientation.

## Undefined values and limits of sampling

Go represents undefined/non-finite samples as JSON `null`, counts them in
`invalid`, and calculates min/max over finite values. The decoder converts null
to NaN in a typed array, never to zero. Curves/fields render gaps and do not connect
through missing values. If all samples are invalid, the request fails.

Isosurfaces require every sample to be finite and the level to be strictly between
sampled min and max. A constant field therefore cannot create an isosurface. A
valid level does not guarantee that a very small geometric feature is resolved.
Likewise, a pole between two curve samples can be missed. Sampling is not a
symbolic domain analysis.

## UI behavior

Use `/#tools`. Each mode keeps its draft settings while the workbench remains
mounted. Editing fields does not calculate until **Plot expression** is pressed.
Request cancellation and cleanup prevent a stale result from overwriting a later
request. The last successful result remains useful while controls are edited.

Plotly provides zoom, pan, reset, rotation for 3D, and PNG export. The Cartesian
bundle handles curves, heatmaps and contours; GL3D handles surfaces and isosurfaces.
WebGL availability is checked before spatial rendering.

The calculator supports inserting keys at the cursor, squaring or taking the
reciprocal of the entire expression, `ans`, and eight recent calculations. Clicking
a history entry restores its expression. Invalid input preserves the previous
answer. Clear clears expression/error, not the saved answer/history. State is
in-memory and is lost when the component is unmounted or the page is reloaded.

## Examples

| Goal | Expression / settings |
| --- | --- |
| Gaussian curve | `exp(-x^2/2)` |
| Harmonic potential surface | `0.5*(x^2+y^2)` |
| Oscillating field | `sin(x)*cos(y)` |
| Sphere | `x^2+y^2+z^2`, level 1, domain [-2,2]³ |
| Double-well curve | `0.5*x^4-x^2` |
| Calculator check | `sqrt(2)/2` |

The legacy `/calculator/` URL redirects to `/#tools` when served by Go. Other old
HTML/JavaScript demos are not the current TypeScript application.
