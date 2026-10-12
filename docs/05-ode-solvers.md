# ODE and symplectic solvers

[Documentation index](README.md) · [Source directory](../ODESolver/differentialEqnSolver.go)

The folder is `ODESolver`, while its Go package name is `EquationSolver`. Examples
may use `EquationSolver "GoProject/ODESolver"` to make that naming explicit.

## Contracts

The primary interface describes `dx/dt=f(x,t)`:

```go
type ODESolver interface {
    NextStep(x, time float64) (float64, error)
    NextStepOnGrid(x []float64, time float64) error
}
```

The scalar method returns a new state. The grid method modifies the caller's
vector. The supplied RHS is `gridData.TDPotentialOp`; its grid callback can couple
components. Time is passed in, not advanced for the caller. Fixed-step callers
must advance time by the configured dt after each successful step.

Constructor names vary (`NewDefine`, `NewDef`), as do redefinition names. Use the
concrete type's actual API rather than assuming a common factory. Many instances
reuse mutable buffers and should belong to one integration at a time.

## Explicit method inventory

“Present” means a method body exists, not that arbitrary equations have been validated.

| Concrete type | Intended order | Scalar | Grid | Current page |
| --- | --- | --- | --- | --- |
| `EulerExplicit` | 1 | Present | Present | No |
| `HeunsExplicit` | 2 | Present | Present | No |
| `MidPointExplicit` | 2 | Present | Present | No |
| `Ralston2order` | 2 | Present | Completed | No; convergence test exists |
| `Ralston3Order` | 3 | Corrected | Completed | Quantum |
| `Huens3Explicit` | 3 | Present | Completed | No |
| `VDHouwenExplicit` | 3 | Present | Completed | No |
| `SSPRungeKutta3` | 3 | Present | Completed | No |
| `RungeKutta3Explicit` | 3 | Present | Completed | No |
| `RungeKutta4Explicit` | 4 | Present | Present | Quantum |
| `RungeKutta38` | 4 | Present | Completed | No |
| `Ralston4Order` | 4 | Present | Completed | No |
| `Nystrom5Explicit` | 5 | Corrected | Completed | Quantum |

All listed fixed-step types now satisfy `ODESolver`. The order column reflects
the named algorithm, not a blanket certification of every legacy scalar table.
The new `NewVectorRK` and `NewComplexRK` constructors cover all 13 fixed formulas
and six embedded pairs with a shared, tested whole-vector engine and reusable
workspace. See [the explicit/adaptive guide](17-explicit-and-adaptive-vector-propagation.md).

## Explicit RK equations

For a tableau `(a,b,c)`:

```text
k[j] = f(x + dt*sum_{l<j} a[j,l]*k[l], t+c[j]*dt)
xNext = x + dt*sum_j b[j]*k[j]
```

`explicitGrid` calculates each stage from the complete vector. It rejects
non-finite stage/output values and copies the new state to the caller only after
success. It currently allocates stage arrays and a work vector per call.

The completed Ralston formulas are:

```text
Ralston 2:
  c = [0, 2/3]
  a21 = 2/3
  b = [1/4, 3/4]

Ralston 3:
  c = [0, 1/2, 3/4]
  a21 = 1/2, a32 = 3/4
  b = [2/9, 1/3, 4/9]
```

The earlier third-stage scalar Ralston error used the first stage instead of the
second. The old grid implementation was a second-order midpoint step; both were
corrected. Nyström 5 stage indices and coefficients were corrected, and its full
six-stage vector tableau is in `explicitGrid.go`. Tests establish scalar/grid
agreement and observed global order for a non-autonomous reference ODE.

RK4's familiar four slopes are evaluated at `t`, `t+dt/2`, `t+dt/2`, and `t+dt`.
It is useful for smooth nonstiff problems, but is not symplectic or unitary.
All explicit choices need a dt consistent with the equation's fastest timescale.

## Symplectic methods

These integrate `x'=v`, `v'=a(x)`. Their potential interface method is named
`ForceAt`, but the method must provide **acceleration**. The web adapter returns
`physicalForce/mass`, converts input momentum with `v=p/m`, and reports `p=m*v`.
Passing an unscaled physical force directly is only correct for unit mass.

| Type | Required state / behavior | Page use |
| --- | --- | --- |
| `VelocityVerlet` | Position and full-step velocity; second order | Classical and Morse |
| `LeapFrog` | Kick-drift-kick with internal half-step velocity; second order | Classical and Morse |
| `Yoshida` | Symmetric fourth-order composition with a negative intermediate stage | Classical and Morse |
| `StromerVerlet` | Position recurrence; initialize previous position from x and v | Library only |

Velocity-Verlet uses:

```text
xNext = x + dt*v + dt²*a(x)/2
vNext = v + dt*(a(x)+a(xNext))/2
```

The implemented leapfrog returns full-step velocities and is algebraically
equivalent to velocity-Verlet for these separable, position-dependent forces.
It is not an independent physical model.

Yoshida uses `w1=1/(2-cuberoot(2))`, `w0=-cuberoot(2)/(2-cuberoot(2))`.
The drift sequence is `[w1/2,(w0+w1)/2,(w0+w1)/2,w1/2]`; kicks are `[w1,w0,w1]`.
Its middle drift weights were previously halved twice. Tests now check fourth-order
convergence and scalar/grid agreement. Negative intermediate stages make it less
suitable for operations requiring strictly forward irreversible evolution.

Störmer–Verlet advances `xNext=2*x-xPrev+dt²*a(x)` with
`xPrev=x-dt*v+dt²*a(x)/2`. Initialization formerly reversed initial velocity.
`InitiateGrid` now reinitializes repeat runs, including changed vector length.

The older `MDodeSolver` interface has signatures that do not match these concrete
integrators. The web adapter selects concrete methods; it does not depend on that
interface. These scalar-force grid methods represent independent coordinates,
not arbitrary pairwise particle forces. The LJ model has its own coupled force loop.

## Legacy adaptive methods: deprecated

`adaptiveMethods.go` contains `HuensEuler`, `FehlbergRK12`, `BogackiShampine`,
`RKFelberg`, `CashKarp` and `DormandPrince`. They are not exposed in the current UI.
They lack a complete accepted-step/time contract: rejection can change internal
dt and return the unchanged state with nil error. A caller cannot safely assume
that a requested interval was traversed.

Concrete issues include an out-of-range `b1Coefs[3]` access in Bogacki–Shampine and
unscaled stage-time fractions in `adaptiveRKBase`. Embedded weights also need a
full audit. Do not treat the existence of these types as a usable adaptive solver.

Use `NewVectorAdaptiveRK` or `NewComplexAdaptiveRK` instead. These new APIs cover
all six pairs, integrate a complete interval, return committed time and step
counts, and support min/max step bounds, attempt limits and coupled-vector error
control. Failure leaves the caller's vector unchanged. The original scalar
implementations remain available with their known defects; they are not the
implementation used by the new API. See [usage and guarantees](17-explicit-and-adaptive-vector-propagation.md).

## Implicit methods (completed 2026-10-12)

For direct whole-vector callbacks, complex states, and matrix-free Hamiltonian
propagation, see [the vector ODE guide](16-vector-ode-and-hamiltonian-propagation.md).
These APIs extend all existing solvers without removing them.

All methods below implement both scalar `NextStep` and coupled-vector
`NextStepOnGrid`. The latter commits to the caller's buffer only after every
stage converges. Errors leave the input unchanged. `MidPointNewton.NextStepIm`
remains an alias for compatibility. Existing `TdFunc` and `DeltaT` fields remain;
all types also provide `NewDefine`, `ReDefine`, `Options`, `Stats` and `Name`.
`ReDefine` preserves options, including derivative callbacks; update those too
when replacing a model with different derivatives.

| Types | Order | Discrete equation / behavior |
| --- | --- | --- |
| `EulerImplicitFixPoint`, `EulerImplicitNewton` | 1 | `y = x + h*f(y,t+h)` |
| `HeunsImplicitFixPoint`, `HeunsImplicitNewton` | 2 | `y = x + h/2*(f(x,t)+f(y,t+h))`; implicit trapezoidal rule |
| `MidPointFixPoint`, `MidPointNewton` | 2 | `y = x + h*f((x+y)/2,t+h/2)` |
| `PredictorCorrector` | 2 | Euler prediction followed by trapezoidal fixed-point correction **to convergence**, not a single explicit correction |
| `SDIRK2` | 2 | Two-stage, stiffly accurate, L-stable method; damped modified Newton in each stage |

Fixed-point iteration still requires a contractive iteration map, which can fail
at a large time step on stiff equations even when the underlying implicit method
is stable. Use Newton for stiff systems. Backward Euler and SDIRK2 damp very fast
decaying modes; midpoint/trapezoidal are A-stable but do not strongly damp those
modes. Midpoint preserves quadratic invariants up to nonlinear-solve error;
SDIRK2 is intended for dissipative stiff equations, not unitary quantum dynamics.
These are Go library methods; the web method menus have not changed.

### Shared nonlinear solver

[implicit.go](../ODESolver/implicit.go) solves the stage residual

```text
R(y) = y - base - coefficient*f(z,time)
z    = (1-theta)*x + theta*y
M    = I - coefficient*theta*Jf(z,time)
M * correction = R(y)
yTrial = y - damping*correction
```

For a vector RHS the default finite-difference Jacobian perturbs each component
and evaluates the **whole vector**, so off-diagonal coupling is retained. Forward
differences use `sqrt(machine epsilon)*max(1,abs(z[j]))`, with the actual
representable increment used as denominator. Scalars also use scaled differences
unless a derivative is supplied. Finite differences assume a smooth RHS; for
badly scaled or restricted-domain models, supplying derivatives is preferable.

Newton solves use pivoted LU, never an explicit matrix inverse. The factors are
reused within a stage for up to four corrections by default, with earlier refresh
when residual reduction is poor. A backtracking search halves the correction
up to twelve times. If a stale Jacobian cannot produce an improving trial, the
solver retries with a fresh Jacobian. Singular/ill-conditioned solves, invalid
values, exhausted iterations and failed searches return errors.

Reuse of matrix data is a standard modified-Newton efficiency strategy; see
[SUNDIALS CVODE](https://computing.llnl.gov/projects/sundials/cvode) and its
[Jacobian reuse discussion](https://sundials.readthedocs.io/en/v6.6.1/cvode/Usage/index.html).
This implementation refreshes at each new step and SDIRK stage; it does **not**
assume the model or step size remains constant across calls.

### Configuration and derivative contracts

A complete, tested scalar example is in
[implicit_example_test.go](../ODESolver/implicit_example_test.go); run it with
`go test ./ODESolver -run ExampleEulerImplicitNewton -v`.

```go
solver := &EquationSolver.SDIRK2{
    TdFunc: rhs, // implements gridData.TDPotentialOp
    DeltaT: 0.01,
    Options: EquationSolver.ImplicitOptions{
        AbsoluteTolerance: 1e-10,
        RelativeTolerance: 1e-8,
        MaxIterations: 30,
        JacobianRefresh: 4,
    },
}
if err := solver.NextStepOnGrid(state, time); err != nil {
    // state and caller time are unchanged; retry with a smaller DeltaT or stop.
    return err
}
time += solver.DeltaT
stats := solver.Stats()
_ = stats // iterations, RHS calls, Jacobians, factorizations, backtracks, residual
```

Zero option values select the defaults shown. `JacobianRefresh: 1` recomputes the
Jacobian on every correction (full Newton). Convergence requires

```text
max_i |R_i| / (AbsoluteTolerance + RelativeTolerance*|xInitial_i|) <= 1.
```

The scale is fixed over the whole step, including both SDIRK stages. These are
nonlinear residual tolerances, **not** an estimate or guarantee of time-integration
error. Reduce `DeltaT` separately to check temporal convergence. There is no
automatic step-size adaptation. Zero dt is a validated no-op; negative finite dt
is allowed. Empty/non-finite states and non-advancing nonzero time steps fail.

Optional callbacks:

- `Options.Derivative func(x,t float64) float64`: scalar RHS derivative, used by
  `NextStep` only. Vector Jacobian callbacks are ignored by scalar steps.
- `Options.Jacobian func(x []float64,t float64,j []float64) error`: fill row-major
  `j[i*n+j] = df_i/dx_j` for a grid call. Supply the RHS Jacobian, not the Newton
  matrix. The destination is cleared first, so known zero entries can be omitted.
- `Options.DiagonalJacobian func(x []float64,t float64,d []float64) error`: fill
  only `df_i/dx_i` **when equations are actually independent**. This explicitly
  selects an O(n) diagonal solve and bypasses dense storage. It is incorrect for
  a coupled Hamiltonian or coupled particle equations. Supply either full or
  diagonal callbacks, not both.

Callbacks must treat inputs as read-only, write their outputs, and not retain
slices. Fixed-point methods do not use derivative callbacks. One solver instance
owns mutable reusable buffers: use it sequentially, do not share concurrently or
copy it after use. Buffer resizing is supported. Repeated fixed-size diagonal
steps do not allocate in the solver after initialization; RHS callbacks may still
allocate. Dense LU uses O(n²) memory and O(n³) factorization work; this path is
bounded to **1,024 components**. Large linear complex equations can now use
`ComplexLinearImplicit` and the Hamiltonian adapters described in the
[vector guide](16-vector-ode-and-hamiltonian-propagation.md). Large nonlinear
systems still need a sparse/banded or Newton-Krylov extension; preconditioning
is not yet implemented in the new linear GMRES path.

### SDIRK2 equations and stability

With `g = 1 - 1/sqrt(2)`:

```text
Y1 = x + h*g*f(Y1,t+g*h)
Y2 = x + h*(1-g)*f(Y1,t+g*h) + h*g*f(Y2,t+h)
xNext = Y2

 c | A
 g | g     0
 1 | 1-g   g
---+----------
 b | 1-g   g
```

This is the two-stage second-order L-stable SDIRK formula associated with
[Alexander (1977)](https://epubs.siam.org/doi/10.1137/0714068).
For `x'=lambda*x`, its stability function is
`R(z)=(1+(1-2*g)*z)/(1-g*z)^2`, `z=h*lambda`; it tends to zero as the magnitude
of a negative real z grows. Both stages solve n-dimensional systems, rather than
one coupled 2n-dimensional stage system. Failure in stage two still leaves the
original caller state untouched.

### Verification and performance

[implicit_test.go](../ODESolver/implicit_test.go) verifies all eight public types,
first/second-order convergence on a non-autonomous equation, scalar/grid agreement,
off-diagonal stiff coupling against analytical stability functions, supplied and
finite-difference Jacobians, damping, factor reuse, stiff decay, singular systems,
callback errors, atomic failure, zero solutions, backward steps, buffer resizing,
legacy aliases and the no-allocation diagonal path.

Run from the repository root in the configured Go/FFTW environment:

```bash
go test ./ODESolver
go test -race ./ODESolver
go test ./ODESolver -run '^$' -bench BenchmarkImplicitNewton -benchmem
```

On the current WSL machine, a 64-component **independent linear** model measured
about 56.7 µs/step with a numerical dense Jacobian, 41.2 µs with an analytic dense
Jacobian, and 1.67 µs with the analytic diagonal path. Dense paths used 72 bytes /
3 allocations per step; diagonal used zero. This compares configurations of the
new solver, not a general speedup over the old implementation. The diagonal result
cannot be extrapolated to coupled equations or an expensive RHS.
