# Explicit and adaptive whole-vector propagation

[Index](README.md) · [Implicit vector solvers](16-vector-ode-and-hamiltonian-propagation.md)

These Go library APIs extend the existing solvers. Original scalar APIs, implicit
solvers, split operators and n-dimensional grids remain available. The browser
now exposes all 13 fixed RK methods, six adaptive pairs and four implicit methods
for quantum, classical and bond dynamics. Spatial wave-packet dynamics is still
1D in the browser; nD time propagation remains a Go library capability.

## Running from the webpage

Open `/#calculations`, choose **Quantum dynamics**, **Classical dynamics**, or
**Bond vibration**, then choose a propagator from its fixed, adaptive or implicit
group. Adaptive and implicit choices show absolute/relative tolerance controls.
For adaptive methods, `dt` is an output interval and maximum internal step, and
`steps` counts those intervals. Final time is always `dt*steps`. Results show RHS
evaluations and accepted/rejected solver steps, as well as physical diagnostics.
Classical and bond methods advance the complete `[position, velocity]` vector.
Quantum methods advance the full complex wavefunction directly. The web path
uses the new vector engines, never the deprecated scalar adaptive implementations.

The API accepts `absoluteTolerance` and `relativeTolerance` in `[1e-12,1e-3]`.
Omitted or zero values select `1e-9` and `1e-7`, respectively. Implicit quantum
methods use matrix-free GMRES for all four formulas; classical implicit methods
use coupled nonlinear Newton solves. The webserver retains its ten-second job
timeout and limits cumulative RHS evaluations across output intervals. Fixed
explicit quantum steps also have a spectral admission bound. Tight tolerances
or difficult models can exhaust the work limit; no partial job is returned.

The 2D particle fluid continues to use its separate velocity-Verlet integrator.
N-dimensional stationary-state calculations and original split/symplectic choices
remain available. Rebuild the frontend and restart the Go server after updating.

## Methods and compatibility

All methods below accept either a full `[]float64` or a full `[]complex128` state.
Each stage evaluates the entire coupled right-hand side, including nonlinear and
non-holomorphic complex equations. Dimensions are limited by workspace memory,
not by a fixed number of coordinates or spatial dimensions.

| `RKMethod` constant | Order | Corresponding original type |
| --- | --- | --- |
| `RKEuler` | 1 | `EulerExplicit` |
| `RKHeun2` | 2 | `HeunsExplicit` |
| `RKMidpoint2` | 2 | `MidPointExplicit` |
| `RKRalston2` | 2 | `Ralston2order` |
| `RKRalston3` | 3 | `Ralston3Order` |
| `RKHeun3` | 3 | `Huens3Explicit` |
| `RKWray3` | 3 | `VDHouwenExplicit` |
| `RKSSP3` | 3 | `SSPRungeKutta3` |
| `RKKutta3` | 3 | `RungeKutta3Explicit` |
| `RKClassic4` | 4 | `RungeKutta4Explicit` |
| `RKThreeEighths4` | 4 | `RungeKutta38` |
| `RKRalston4` | 4 | `Ralston4Order` |
| `RKNystrom5` | 5 | `Nystrom5Explicit` |
| `RKHeunEuler21` | 2, embedded 1 | `HuensEuler` |
| `RKFehlberg21` | 2, embedded 1 | `FehlbergRK12` |
| `RKBogackiShampine32` | 3, embedded 2 | `BogackiShampine` |
| `RKFehlberg54` | 5, embedded 4 | `RKFelberg` |
| `RKCashKarp54` | 5, embedded 4 | `CashKarp` |
| `RKDormandPrince54` | 5, embedded 4 | `DormandPrince` |

The six previously missing fixed-step `NextStepOnGrid` methods are also supplied
on their original types. The six original scalar adaptive types are **deprecated**
and retained unchanged apart from documentation; they still have the defects
listed in [the status guide](11-status-and-limitations.md). Use the new adaptive
API even for scalar calculations, representing a scalar as a length-one vector.
Do not assume a new implementation repaired calls to the old adaptive types.

## Fixed steps

Construct `NewVectorRK(rhs, method)` or `NewComplexRK(rhs, method)`, then call
`Step(y, t, dt)` or `StepContext(ctx, y, t, dt)`. All 19 selections are accepted;
an embedded pair uses its higher-order weights without adapting the fixed step.

The callback signatures are:

```go
type VectorRHS func(t float64, y, dy []float64) error
type ComplexRHS func(t float64, y, dy []complex128) error
```

Treat `y` as read-only, fill every component of `dy`, and retain neither slice.
Return an error for invalid model evaluations. The RHS must give the same result
for the same `(t,y)` during an integration: rejected stages and derivative reuse
are incompatible with callbacks that irreversibly advance external model state.

Stages use `t+c[j]*dt`. The caller advances its own time only after a successful
fixed step. A failure, non-finite stage or cancellation leaves its state unchanged.
Both positive and negative steps are supported. A zero step is a validated no-op.

## Adaptive intervals

Construct `NewVectorAdaptiveRK(rhs, method, options)` or
`NewComplexAdaptiveRK(rhs, method, options)`. Only the six embedded pairs are
accepted. Call `Integrate(y, tStart, tEnd)` or its `IntegrateContext` variant.

```go
import EquationSolver "GoProject/ODESolver"

rhs := func(t float64, psi, dpsi []complex128) error {
    // Two coupled quantum amplitudes, H=[[0,1],[1,0]], hbar=1.
    dpsi[0] = -1i * psi[1]
    dpsi[1] = -1i * psi[0]
    return nil
}
solver, err := EquationSolver.NewComplexAdaptiveRK(
    rhs, EquationSolver.RKDormandPrince54,
    EquationSolver.AdaptiveOptions{
        AbsoluteTolerance: 1e-12,
        RelativeTolerance: 1e-10,
    })
if err != nil { return err }
psi := []complex128{1, 0}
result, err := solver.Integrate(psi, 0, 1)
if err != nil { return err }
// result.Time==1; psi approximately [0.540302, -0.841471i].
// Inspect result.AcceptedSteps, RejectedSteps and RHSEvaluations.
```

`Integrate` traverses the entire requested interval, shortening its final step to
reach the endpoint. It handles rejection internally. On success `result.Time`
equals `tEnd`. On **any failure**, the caller's entire vector is unchanged and
`result.Time` remains `tStart`, even if some internal steps succeeded. Counts
describe attempted work; they do not indicate a partially committed state.
`NextStep` is a signed suggestion, not a hidden update to the caller's time.
Separate calls start afresh using the configured `InitialStep`.

| Option | Zero-value default | Meaning |
| --- | --- | --- |
| `AbsoluteTolerance` | `1e-9` | Per-component absolute local tolerance |
| `RelativeTolerance` | `1e-7` | Per-component relative local tolerance |
| `InitialStep` | One tenth of the interval | Initial positive step magnitude |
| `MinStep` | Any representable step | Smallest retry magnitude |
| `MaxStep` | Largest finite float | Largest internal step magnitude |
| `MaxAttempts` | `100000` | Total accepted plus rejected attempts per call |

All values must be finite and nonnegative; `MinStep <= MaxStep`. Initial steps
are clamped to the configured bounds. An endpoint step may be smaller than
`MinStep`; it must meet tolerance without subdividing further. Zero tolerances
select defaults rather than disabling a tolerance term.

For each component, the high/low difference is divided by
`atol + rtol*max(abs(yOld),abs(yHigh))`. Complex components use complex modulus.
The maximum component ratio must be at most one to accept a step. The controller
uses a safety factor of 0.9 and exponent `-1/(embeddedOrder+1)`, limiting each
factor to 0.2–5. Error estimates are local; they do not guarantee a specified global
error, conserved norm, or stable behavior on a stiff problem.

## Hamiltonians in one or n spatial dimensions

`Quantum.NewRKPropagator(h, method)` and `NewNDimRKPropagator(h, method)` provide
fixed steps. `NewAdaptiveRKPropagator(h, method, options)` and
`NewNDimAdaptiveRKPropagator(h, method, options)` integrate complete intervals.
They expose the same `Step`/`Integrate` and context variants described above.

The adapters apply `-i*H(t)*psi` to the whole complex wavefunction with `hbar=1`.
They use `ApplyHComplexAt`, so no dense Hamiltonian or real/imaginary packing is
needed. Static and time-dependent potentials work in DVR and Fourier bases;
time-dependent callbacks are evaluated at every stage's actual time. Potential
callbacks provide the **complete potential**, including any static part.
The nD state is the existing flattened tensor grid; all spatial coupling remains.

Explicit RK is not exactly unitary and no automatic normalization is performed.
Check phase, norm and step convergence for the chosen spectrum. Implicit midpoint
or the existing split methods remain available when their properties fit the
calculation. Abrupt changes to the RHS should be split into separate intervals.

## Cost, safety and validation

Workspace is proportional to `state length * stage count` and reused between
calls of the same size. A conservative 256 MiB workspace bound applies. No
matrices or Jacobians are allocated by these explicit engines. Fixed/adaptive
complex tests on 2048 components measure zero solver allocations after warm-up;
user callbacks can still allocate.

Rejected steps reuse their initial derivative. Bogacki–Shampine and
Dormand–Prince also reuse the final derivative as the next accepted step's first
derivative (FSAL). This cache is limited to one `Integrate` call, so changing the
state or model between calls does not reuse stale derivatives. Instances and
their Hamiltonians must not be copied or used concurrently.

Tests check all 19 higher-order formulas on coupled nonlinear real and complex
reference solutions, all six lower-order formulas, adaptive tolerance response,
forward/backward endpoints, rejection, rollback, cancellation, limits, allocation
and derivative counts, legacy grid compatibility, static DVR/Fourier phase
evolution, and all six adaptive pairs on a time-dependent 4D Fourier Hamiltonian.

Sources: [RK catalog and coefficients](../ODESolver/rk_tableau.go),
[engines](../ODESolver/vector_rk.go), [Hamiltonian adapters](../Quantum/rk_propagation.go).
Reference coefficient tables include the
[SUNDIALS ARKODE tables](https://sundials.readthedocs.io/en/v7.1.1/arkode/Butcher_link.html)
and [Fehlberg's NASA report](https://ntrs.nasa.gov/citations/19690021375).
