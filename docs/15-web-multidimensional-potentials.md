# Multidimensional potentials in the web calculator

[Documentation index](README.md) · [General nD library](14-multidimensional-quantum.md)

Added on 2026-10-10. Open **Calculations → N-dimensional states** at
`http://127.0.0.1:8080/#calculations` after rebuilding the frontend and restarting
the Go server. This mode calculates instantaneous energy eigenstates in 2, 3 or
4 Cartesian coordinates. The numerical work runs in Go; TypeScript controls
forms, validates responses and renders plots.

## Supported models

Define `q[d] = scale[d] * (x[d] - center[d])`, `r² = sum(q[d]²)`.
Strength/depth parameters in the web form are positive. Attractive potentials
explicitly include a negative sign.

| Choice | Definition | Interpretation |
| --- | --- | --- |
| Harmonic | `k/2 * sum(q[d]²)` | Independent oscillators, with per-axis mass/scale |
| Coupled harmonic | `k*(sum(q[d]²)/2 + c*sum(q[d]*q[d+1]))` | Open chain of coupled coordinates; `abs(c)<=0.45` keeps this form confining |
| Morse | `De*sum((1-exp(-alpha*q[d]))²)` | Additive extension of the original Morse code, not a radial bond potential |
| Quartic double wells | `k/4*sum((q[d]²-1)²)` | Two minima along every coordinate; the energy offset differs from the old 1D web double-well convention |
| Gaussian | `-D*exp(-r²/(2*sigma²))` | A localized attractive well; scaling makes it ellipsoidal in x |
| Super-Gaussian | `-D*exp(-(r/sigma)^p)` | A flatter-bottomed well; p is 2, 4, 6 or 8 |
| Soft-core | `-Z/sqrt(r²+a²)` | Smooth attractive Coulomb-like model; a is the softening length |
| Two Gaussians | `-D*sum_(±) exp(-abs(q ± s*e1/2)²/(2*sigma²))` | Two localized wells separated along q1; overlapping wells can be deeper than D |
| Sinusoidal lattice | `D*sum(sin(kappa*q[d])²)` | Periodic potential, subject to the selected box/boundary convention |
| Free particle | `0` | Finite-box or periodic free-particle spectrum |

The Gaussian, Morse, soft-core, super-Gaussian and polynomial implementations
reuse the existing potential types. Their `NDimEvaluateAt` methods are defined
in [nd_potentials.go](../gridData/nd_potentials.go). MultiGaussian also has a Go
library additive extension, but is not the same model as the web's two localized
radial wells. The old scalar methods and 1D interface are preserved.

For generic complex coordinates the radial extensions use analytic squares
`sum((z-center)²)`, not absolute squares. Super-Gaussian's radial extension agrees
with the 1D scalar formula for even order; the web exposes only even orders.
The web Hamiltonian uses real potentials and Cartesian kinetic derivatives.

## Controls and plots

### Enter your own potential

Select **Custom formula** in the potential menu. Type just the right-hand side
of `V(x1, …, xN)`, then run the calculation. For example, in 2D:

```text
0.5*(x1^2 + x2^2) + 0.2*x1*x2
```

This gives a coupled oscillator. With both masses equal to one, its continuum
ground-state energy is `(sqrt(1.2)+sqrt(0.8))/2 ≈ 0.994936153`; the grid result
approaches this as spatial resolution and domain converge.

Other examples:

```text
0.5*(x1^2 + x2^2 + x3^2) + 0.1*x1^4
-5*exp(-(x1^2 + x2^2)/(2*1.2^2))
0.5*((x1-1)^2 + 2*x2^2 + x3^2 + x4^2)
```

Choose three or four coordinates for formulas using `x3` or `x4`. Coordinate
names refer to the **actual grid coordinates**. Custom formulas do not apply
preset strength, center, or scale transformations. Write all coefficients,
offsets and widths explicitly. Per-axis masses still determine kinetic energy.
The **Hidden slice position** controls only the potential plot's fixed hidden
coordinates; it does not shift the Hamiltonian. The previous `center` API field
stores this position. `scale` remains in the axis schema but is unused by custom
formulas. Changing dimensions retains your formula; examples adapt to the current
dimension when selected. Switching presets also retains the custom formula.

Supported syntax:

- Variables `x1` through `xN`, constants `pi` (or `π`) and `e`, decimal and
  scientific-notation numbers.
- Explicit `+ - * / %`, powers `^` or `**`, parentheses and unary signs.
  Powers associate right to left; `-x1^2` means `-(x1^2)`.
- Unary functions `sin cos tan asin acos atan sinh cosh tanh sqrt abs exp ln
  log log10 floor ceil round`. Trigonometric arguments use radians. `ln` is
  natural logarithm; `log` and `log10` mean base ten.
- Two-argument functions `pow(a,b)`, `min(a,b)` and `max(a,b)`.

Go parses a bounded mathematical expression rather than executing source code.
Assignments, arbitrary function definitions, implicit multiplication and names
outside the selected dimensions are rejected. Limits are 1,024 UTF-8 bytes and
256 tokens. Every sampled value, including the displayed potential slice, must
be finite and satisfy `abs(V) <= 100000`. Errors identify invalid syntax or the
coordinates of a bad sample. Sampling cannot detect a singularity that lies
between all sampled points; check the formula and spatial convergence yourself.

The same full-grid Hamiltonian, eigensolver, marginal densities and isosurfaces
are used for custom and preset potentials. This feature currently applies to
**N-dimensional states (2D–4D)**, not the separate dynamics or 1D forms.

API: send `"potential":"custom"` and `nd.expression` containing the formula.
Preset-only coefficients can be omitted for a custom request. Required grid,
basis, mass, projection and state limits still apply. The response echoes the
formula and includes it in the calculation notes.

Implementation: [expression.go](../internal/plotmath/expression.go) shares the
existing plotter grammar through `CompileND`, with arbitrary positive
dimensionality and indexed coordinates; the web limits remain 2–4D.
Evaluation reuses the parsed tree without allocating at each sample. See
[parser tests](../internal/plotmath/nd_expression_test.go),
[solver tests](../internal/simulation/custom_potential_test.go) and
[API tests](../internal/webserver/simulations_test.go).

### Grid and result controls

Each axis has its own number of points, half-width, mass, center and coordinate
scale. Its grid spans `[-L,L)`. Changing the number of dimensions resets axis
settings to a grid that fits the calculation limit; changing a potential selects
its default strength and Gaussian width. Other edited parameters are retained.

Choose 1–4 states and which state supplies the 2D/3D density plots. State 0 is
the lowest state. Changing this selection or the projection axes requires another
run. Results keep the submitted parameters and show a previous-run badge when
the form differs from those parameters.

- The line plot shows each requested state's marginal density along the selected
  horizontal axis.
- The density heatmap integrates over all hidden coordinates. It can also be
  displayed as a 3D surface, whose vertical coordinate is probability density.
- The potential plot is a **slice** with hidden coordinates fixed at their
  centers. It is not a marginal or an integrated effective potential. Its
  colorbar/vertical axis is labeled V.
- For 3D and 4D calculations, an isosurface shows the selected state's density
  in x1,x2,x3. In 4D, x4 is integrated out. The slider chooses a percentage of
  peak density, not a percentage of enclosed probability. Nearly uniform states
  show an explanation instead of an empty isosurface.

Plot integration uses probability weights from normalized eigenvectors:
`rho_x[ix] += psi[i]²/dx`, `rho_xy[iy][ix] += psi[i]²/(dx*dy)`. The 3D volume
similarly divides by three displayed spacings and sums hidden indices. These
marginals integrate to one. Original last-coordinate-fastest grid ordering is
preserved; heatmaps use rows for y and columns for x.

WebGL is required for surface and isosurface views. Heatmaps and line plots stay
available when graphics acceleration is unavailable. Plotly's Cartesian and
3D bundles load on demand.

## Bounds and solver contract

| Input | Web limit |
| --- | --- |
| Coordinates | 2–4; the general Go library is not restricted to four |
| Points per axis | 4–64, with at most 8,192 total tensor points |
| Axis half-width | 2–12 |
| Axis mass | 0.1–10 |
| Axis coordinate scale | 0.25–2 |
| Axis center | Between `-L/2` and `L/2` |
| States | 1–4 |
| Strength/depth | 0.01–30 |
| Morse alpha | 0.05–1 |
| Gaussian width | 0.3–4 |
| Coupling | -0.45–0.45 |
| Softening | 0.2–3 |
| Separation | 0.2–6, in transformed q1 coordinates |
| Lattice wavenumber | 0.1–3 |

The solver uses `LowestStates`, a subspace of up to 96 vectors, at most 500
iterations, and residual tolerance `1e-7*max(1,abs(E))`. It checks cancellation
between operator applications. Nonconvergence returns an error; approximate
unconverged spectra are not labeled complete. Sampled potential magnitude must
not exceed 100,000. The shared HTTP limit remains one simulation at a time,
8,192-byte requests, and a ten-second calculation deadline.

Coarse grids, especially in 4D, are for model exploration. Numerical residuals
do not establish spatial convergence or prove that all lower levels were found.
Increase grid/domain and compare spectra. States above a dissociation threshold,
and states of free or lattice models, may describe the box rather than a bound
molecule. The mode does not add electronic structure, radial kinetic metrics,
chemical force fields or multidimensional time-evolution controls. General nD
time evolution remains available through the Go library.

## HTTP request

Use `POST /api/simulations` with `kind: "nd-spectrum"`. The `nd` object is mandatory
for this kind and rejected for other kinds. Required potential parameters must
be supplied even if the selected model does not use them; Go does not apply UI
defaults. Unknown fields at any nested level are rejected.

```json
{
  "kind": "nd-spectrum",
  "basis": "fourier",
  "potential": "coupled-harmonic",
  "strength": 1,
  "alpha": 0.25,
  "sigma": 1.2,
  "states": 4,
  "nd": {
    "axes": [
      { "points": 24, "halfWidth": 6, "mass": 1, "center": 0, "scale": 1 },
      { "points": 24, "halfWidth": 6, "mass": 1, "center": 0, "scale": 1 }
    ],
    "coupling": 0.2,
    "softening": 0.7,
    "separation": 3,
    "wavenumber": 1,
    "order": 4,
    "state": 0,
    "axisX": 0,
    "axisY": 1
  }
}
```

Axis and state indices are zero-based. The response extends the existing result
with optional `volume: {title,labels,x,y,z,values}`. These are paired flattened
coordinate/value arrays, at most 8,192 entries each. Charts can include
`valueLabel` for the colorbar/vertical axis; potential slices use `"V"`.

Implementation: [Go web adapter](../internal/simulation/nd_quantum.go),
[axis/model controls](../webpage/web-app/src/components/NDControls.tsx),
[3D density view](../webpage/web-app/src/components/DensityVolume.tsx).

## Verification

The preset regression runs all ten models in 2D, 3D and 4D, with both Fourier and
DVR: 60 combinations at the default grid resolutions. Additional tests check
analytic coupled-oscillator energy, marginal normalization with unequal axes,
hidden-coordinate integration, potential slices, malformed requests, complex
potential formulas and response validation. A 3D localized-well eigensystem is
also compared with a full dense reference on a smaller grid.

The full Go suite, static checks, targeted race checks, frontend tests, TypeScript,
lint and production build passed during development. Browser checks exercised a
2D coupled oscillator, a 3D two-Gaussian well, a 4D DVR model, projection/state
selection, heatmap/surface switching, isosurface threshold and oversized-grid
validation. These checks cover the defaults and listed behaviors, not arbitrary
physical parameter combinations.
