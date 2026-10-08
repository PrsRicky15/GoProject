# Grids, potentials and units

[Documentation index](README.md) · [Source: gridData](../gridData/grid.go)

## Spatial grid

`RadGrid` stores bounds, point count, spacing and conjugate-space parameters.
Despite the name, the active quantum page uses a Cartesian coordinate on a symmetric
interval, not a radial equation with angular-momentum terms.

For `NewRGrid(a, b, N)`:

```text
length = b - a
dx = (b - a) / N
x[j] = a + j*dx,   j = 0,...,N-1
dk = 2*pi / (b-a)
```

The interval is **[a,b)**. The final sample is `b-dx`, not `b`. This differs from
the expression plotter, which includes both domain endpoints. `NewFromLength(D,N)`
creates `[-D/2,D/2)`. The calculation page sets `a=-L`, `b=L`, so `dx=2L/N`.

| API | Contract |
| --- | --- |
| `NewRGrid(min,max,n)` | Checked grid construction; requires positive n and increasing bounds |
| `NewFromLength(length,n)` | Symmetric grid from total length |
| `RValues()` | Newly generated real-coordinate slice |
| `KValues()` | Newly generated FFT-ordered wavenumber slice |
| `RMin/RMax/NPoints/DeltaR/Length` | Real-grid metadata |
| `DeltaK/KMin/KMax/CutoffE` | Conjugate metadata; cutoff is a unit-mass estimate |
| `ReDefine`, `ReDefineMinMax`, `ReDefineLength` | Mutate the grid definition |
| `PotentialOnGrid`, `ForceOnGrid` | Evaluate a compatible potential on generated points |
| Display/Print methods | Console and text-file helpers; callers supply formatting |

Direct grid constructors do not comprehensively reject NaN/Infinity. The active
simulation validator and Hamiltonian constructor provide stronger checks. A
Hamiltonian snapshots its grid, so redefining the original grid does not resize
an existing Hamiltonian. Reconstruct the Hamiltonian after changing mass/grid.

## Fourier ordering

The sequence is zero, positive modes, then negative modes:

```text
N=5: [0, 1, 2, -2, -1] * dk
N=6: [0, 1, 2, -3, -2, -1] * dk
```

For even N, the Nyquist bin is negative. `KMin/KMax` describe nominal limits
`±pi/dx`; both are not necessarily present as samples. `CutoffE=KMax²/2` omits
mass and is therefore not a general mass-dependent kinetic-energy bound.

## Time grid

`NewTimeGrid(MacroDT, MacroSteps, MicroSteps)` defines:

```text
N = MacroSteps * MicroSteps
tMin = 0
tMax = MacroDT * MacroSteps
micro dt = MacroDT / MicroSteps
```

`TValues()` is endpoint-exclusive, and `WValues()` uses the same FFT ordering as
the spatial grid. The active web simulations do not use `TimeGrid`: they advance
explicitly from `i=0` to `steps` and include the final time `steps*dt` in output.

The legacy `TimeGrid.ReDefineLength` currently sets `tMax=length/2`, and general
redefinitions leave macro/micro metadata unchanged. Do not interpret this behavior
as a verified time-grid API; see the [limitations](11-status-and-limitations.md).

## File input/output

`NewRGridFromFile(directory)` reads `rgrid.inp`; `NewTGridFromFile(directory)` reads
`tgrid.inp`. The parser skips blank/comment lines and recognizes `key: value` pairs:

```text
# rgrid.inp
rMin: -8
rMax: 8
nPoints: 64
```

```text
# tgrid.inp
MacroDT: 0.1
MacroSteps: 100
MicroSteps: 10
```

These old readers ignore numeric parse errors and scanner errors, and do not close
their opened input files. Use validated constructors for new integration code.
Print helpers write coordinate/value columns; vector lengths must match the grid.
They are not a trajectory storage format or a browser import/export feature.

## Potential interface

`VarType` permits `float64` or `complex128`. The generic interface is:

```go
type PotentialOp[T VarType] interface {
    EvaluateAt(x T) T
    EvaluateOnGrid(x []T) []T
    ForceAt(x T) T
    ForceOnGrid(x []T) []T
}
```

The intended force convention is `F=-dV/dx`. Methods named `OnGrid` allocate a
result slice. Parameters are ordinary struct fields; not all potential types
validate parameter ranges internally. The web request layer validates its subset.

### Potentials selected by the current page

Let `q=x-center`:

| Potential | V(x) | F(x) |
| --- | --- | --- |
| Harmonic | `k*q²/2` | `-k*q` |
| Morse | `De*(1-exp(-alpha*q))²` | `-2*alpha*De*exp(-alpha*q)*(1-exp(-alpha*q))` |
| Symmetric double well | `s*(x⁴/2-x²)` | `2*s*x-2*s*x³` |
| Free | `0` | `0` |

The double well is a `Polynomial` with ascending-power coefficients
`[0,0,-s,0,s/2]`. It ignores `center`; its minima are at ±1 and the barrier relative
to a minimum is `s/2`. Morse has its minimum at zero energy and asymptote `De`.
Consequently, Morse energies at or above `De` are not bound levels on the full line,
even if a finite matrix returns discrete eigenvalues.

`Polynomial{Coeffs: []float64{c0,c1,...}}` means `sum(c[j]*x^j)`; its force uses
the negative derivative. Both real and complex evaluation paths exist.

### Other library potentials

| Type | Parameters / implemented energy | Status |
| --- | --- | --- |
| `SoftCore` | `Charge/sqrt((x-Centre)²+SoftParam²)` | Library only; charge sign controls attraction/repulsion |
| `Gaussian` | `Strength*exp(-(x-Cen)²/(2*Sigma²))` | Energy implementation exists; force sign defect |
| `MultiGaussian` | Sum of evenly spaced Gaussians, symmetric about zero | Library only; force behavior requires verification |
| `SuperGaussian` | Real energy `Strength*exp(-((x-Cen)/Sigma)^Order)` | Real/complex inconsistencies; not exposed |

Even `MultiGaussian.NumGauss` uses centers at `±Gap*(j+1/2)`; odd counts include a
center at zero plus pairs at `±Gap*(j+1)`. No page control currently selects these
additional potential types. `SuperGaussian` is usually used with even positive
order for symmetric localization, but its struct does not enforce that choice.

**Observed defects:** Gaussian's `ForceAt` returns `dV/dx` rather than `-dV/dx`.
SuperGaussian's real force has the analogous sign problem, and its complex energy
uses power 2 regardless of `Order`. Existing `TestGaussian_ForceAt` actually tests
a time grid, so it does not establish correctness of Gaussian forces.

## Time-dependent RHS interface

`TDPotentialOp` provides `EvaluateAt(x,t)`, `EvaluateOnRGrid(x,t)`, and
`EvaluateOnRGridInPlace(x,out,t)`. In ODE code, these represent the **right-hand
side of an ODE**, despite the name. A vector RHS may couple all components. It
must fill every output entry and must not retain the solver's temporary buffers.

For quantum integration, the vector is `[Re(psi), Im(psi)]`; for a general ODE,
its interpretation belongs to the caller. This interface is distinct from the
Hamiltonian's time-dependent potential callback described in [quantum](04-quantum.md).

## Units

The quantum and one-coordinate classical page uses atomic-unit conventions with
hbar=1. There is no automatic conversion from electron-volts, femtoseconds,
angstroms, or atomic mass units. In the Morse bond mode, mass is the reduced mass,
`mu=m1*m2/(m1+m2)`, in consistent units. The bond coordinate is a model coordinate;
the program does not prevent it from becoming negative.

The particle fluid instead uses reduced LJ units, explained in
[classical and MD](06-classical-and-md.md). Do not compare numbers across these
two unit systems without specifying a physical mapping.
