# Implementation status and known limitations

[Documentation index](README.md)

This is a source-based snapshot finalized on 2026-10-08. It distinguishes observed
defects from areas merely lacking validation. Documentation work has not silently
changed the numerical methods described here.

## Active, exercised paths

- Restricted Go plot expression parser and bounded 1D/2D/3D/isosurface sampling.
- TypeScript scientific calculator, typed request normalization and runtime decoders.
- DVR/Fourier Hamiltonian composition, real symmetric spectra and operator application.
- Midpoint split propagation and RK4/Ralston 3/Nyström 5 quantum grid integration.
- Classical velocity-Verlet/leapfrog/Yoshida trajectories and Morse bond model.
- 2D periodic shifted-force Lennard–Jones fluid and sampled playback.
- Same-origin HTTP service, strict top-level request fields, job bounds and cancellation.

These statements apply to the tested models and parameter ranges, not arbitrary
problems. See [testing](10-testing.md) for the evidence and coverage gaps.

## Completed repairs

| Area | Change now present |
| --- | --- |
| DVR exponentials | `ExpIdt`, `ExpIdtTo`, `ExpDtTo` perform their stated operations |
| DVR diagonalization | Uses a copy, preserving the cached kinetic matrix |
| Real matrix diagonalization | Replaced incompatible type assertions with symmetric eigensolver handling |
| Ralston 3 | Correct scalar third stage and actual third-order grid step |
| Ralston 2 / Nyström 5 | Grid steps implemented; Nyström scalar stage indices corrected |
| Yoshida | Middle drift factors corrected; proper name and coefficient initialization on redefinition |
| Störmer–Verlet | Previous-position initialization sign and repeat grid initialization corrected |
| Fourier operations | FFT ordering/normalization and kinetic application aligned; planner/destructor lock added |

## Concrete remaining problems

| Area / source | Observed issue | Consequence / next step |
| --- | --- | --- |
| [FiniteDifference.go](../OperatorAlgebra/FiniteDifference.go) | Constructor accepts 3/5/6, computation switches on 3/5/7/9; high-stencil coefficient signs conflict with the negative kinetic prefactor | Do not expose as a verified basis; reconcile stencil names, signs, boundary closure and public application API |
| [dvrBasis.go](../OperatorAlgebra/dvrBasis.go) | Half-line branch chosen by `abs(RMax+RMin)>=1` | Asymmetric domains can select an inappropriate formula; replace heuristic with explicit boundary/basis selection |
| [matrix_interface.go](../OperatorAlgebra/matrix_interface.go) | `zomplex64` backend is declared but never initialized | `ZGeev`/general complex diagonalization cannot be treated as working; configure a real implementation or replace the wrapper |
| [derivativeOp.go](../OperatorAlgebra/derivativeOp.go) | Exponential interfaces omit error returns present on concrete DVR in-place methods | Concrete DVR type does not satisfy these legacy interfaces as written |
| [exOperator.go](../ODESolver/exOperator.go) | Exponential interface also lacks concrete error-return signatures | Reconcile before generic use |
| [sympleticMethods.go](../ODESolver/sympleticMethods.go) | `MDodeSolver` signatures differ from implemented methods | Current adapters use concrete solvers; old interface is not a working abstraction |
| [adaptiveMethods.go](../ODESolver/adaptiveMethods.go) | Rejected steps can return unchanged state and nil error while changing hidden dt | Caller cannot track accepted time correctly; redesign acceptance/time contract |
| Same adaptive file | Bogacki–Shampine accesses index 3 of a three-element `b1Coefs` slice | Runtime panic in that path |
| Same adaptive file | `adaptiveRKBase` adds unscaled c fractions to time | Incorrect non-autonomous stage times; complete coefficient/time audit needed |
| [odeFixPoint.go](../ODESolver/odeFixPoint.go) | Heun predictor/iteration/output paths rebind local slice variables | Intended updates do not reliably reach caller/work buffers |
| Same fixed-point file | Midpoint copies initial state before allocating new buffers | First call after a size change loses the intended starting iterate |
| [odeNewtonRapson.go](../ODESolver/odeNewtonRapson.go) | Heun predictor rebinding; componentwise Newton derivatives | Initialization problem and no general coupled Jacobian solve |
| Same Newton file | Midpoint has `NextStepIm` but no complete standard scalar/grid interface | Not selectable as an `ODESolver` |
| [potential.go](../gridData/potential.go) | Gaussian force is the energy derivative rather than its negative | Incorrect force convention for direct classical use |
| Same potential file | SuperGaussian force sign; complex energy hardcodes order 2 | Inconsistent real/complex model; needs derivative and analytic-continuation tests |
| [grid.go](../gridData/grid.go), [timegrid.go](../gridData/timegrid.go) | Input readers ignore parse/scanner errors and leave files open | Malformed input can produce unintended values; repeated reads leak handles |
| [timegrid.go](../gridData/timegrid.go) | `ReDefineLength` halves length and redefinition leaves macro/micro metadata stale | Reported timing metadata can disagree with grid definition |
| [gridFunctionality.go](../gridData/gridFunctionality.go) | Base grid validation lacks comprehensive non-finite checks | Direct library callers need validation; web path validates separately |
| [numericalMethods.go](../gridData/NumericalRecipies/numericalMethods.go) | Chebyshev first-kind slice/matrix paths return input for order >1 | Only scalar real/complex recurrence is implemented beyond initial orders |

The listed issues do not all affect the current UI. In particular, the page uses
symmetric DVR grids, supported fixed-step methods and a limited potential set.
Isolation of those paths is not a repair of the unused library code.

## Partial declarations and unverified methods

- `PredictorCorrector`, `ChebyshevSecondKind`, `ChebyshevComplex`: declarations without algorithms.
- `NumericalRecipies.Integration`: interface only; trapezoidal/Simpson implementations absent.
- `NumericalRecipies/differentiation.go`, `Plots/basicPlot.go`, and `classical/moleulardyn.go`: package-only placeholders.
- `classical.ModelSystem` / `SystemInfo`: old data holders, with unexported fields and no complete evolution API.
- `InitialValueIntegrator`: interface with no concrete implementation in that old model path.
- Several named third/fourth-order explicit methods have scalar steps but no grid steps; see the complete [solver table](05-ode-solvers.md).
- MultiGaussian forces, most non-web potential variants, complex scaling, and legacy adaptive coefficient tables need dedicated validation.

The presence of complex-number support does not imply a working complex general
eigensolver. The existence of a solver name does not imply it implements the
project's complete solver interface.

## Current model limits

Quantum models are 1D, use real local potentials, finite domains, and have no
absorbing boundary. Fourier is periodic. Dense eigensolvers limit practical scale.
There is no automatic spatial refinement or time-step adaptation in active jobs.
Sampling a narrow packet on a coarse grid can be inaccurate even when validation passes.

Classical scalar-force methods do not supply arbitrary coupled particle forces.
The separate fluid is 2D, monatomic and reduced-unit, with no chemical topology,
thermostat, barostat, electrostatics, species, pressure or transport-property analysis.
The LJ fail-fast distance guard rejects bad steps rather than recovering them.

## Product and operational limits

No authentication, authorization, database, durable result storage, job queue,
cloud execution, GPU backend, team collaboration, molecular file import, general
data import/export or public-deployment configuration exists. The landing page's
GPU/cloud/data-management claims and some buttons are placeholders. “Sign In” is
not connected to a login flow.

The Go service binds to loopback by default. Cross-origin protection, limited JSON
and bounded jobs reduce specific risks but do not make the app a hosted multi-user
service. Dense factorization does not stop mid-operation when its context expires.
Library methods can panic on invalid dimensions. Public raw arrays can be misused.

## Suggested order for future work

1. Repair potential force definitions and add genuine force-gradient tests.
2. Make basis boundaries explicit; finish finite-difference stencils and convergence tests.
3. Reconcile ODE/operator interfaces and state ownership/error contracts.
4. Redesign adaptive accepted-step handling and implement coupled implicit solves.
5. Finish numerical-recipe placeholders only against concrete use cases and reference tests.
6. Optimize repeated RHS allocations, Fourier split propagation and MD neighbor search when scale requires it.
7. Add durable result export and automated browser checks before expanding the product surface.

This list is a proposed maintenance sequence, not a statement that those changes
have been implemented or scheduled.
