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
| `Huens3Explicit` | 3 | Present | Missing | No |
| `VDHouwenExplicit` | 3 | Present | Missing | No |
| `SSPRungeKutta3` | 3 | Present | Missing | No |
| `RungeKutta3Explicit` | 3 | Present | Missing | No |
| `RungeKutta4Explicit` | 4 | Present | Present | Quantum |
| `RungeKutta38` | 4 | Present | Missing | No |
| `Ralston4Order` | 4 | Present | Missing | No |
| `Nystrom5Explicit` | 5 | Corrected | Completed | Quantum |

Scalar-only types do not satisfy `ODESolver`, because its grid method is required.
The order column reflects the named algorithm, not a blanket certification of
every legacy coefficient table.

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

## Adaptive methods: experimental

`adaptiveMethods.go` contains `HuensEuler`, `FehlbergRK12`, `BogackiShampine`,
`RKFelberg`, `CashKarp` and `DormandPrince`. They are not exposed in the current UI.
They lack a complete accepted-step/time contract: rejection can change internal
dt and return the unchanged state with nil error. A caller cannot safely assume
that a requested interval was traversed.

Concrete issues include an out-of-range `b1Coefs[3]` access in Bogacki–Shampine and
unscaled stage-time fractions in `adaptiveRKBase`. Embedded weights also need a
full audit. Do not treat the existence of these types as a usable adaptive solver.

A future API should explicitly return accepted state, accepted time, error estimate,
and suggested next step, or internally integrate a complete requested interval.
It also needs rejection limits, min/max step bounds and coupled-vector error norms.

## Implicit methods: library experiments

`odeFixPoint.go` includes Euler, Heun and midpoint fixed-point iterations.
`odeNewtonRapson.go` includes Euler/Heun Newton iterations and a scalar midpoint
method named `NextStepIm`. `PredictorCorrector` is only a struct declaration.

The shared constants are `maxIter=20`, `tolerance=1e-7`, finite-difference
`delta=1e-6`, and `adaptiveTolerance=1e-8`; no page exposes tolerance controls.
Several Heun grid helpers replace local slice variables rather than copying into
the caller's buffer. Newton grid paths use componentwise derivatives rather than
the full coupled Jacobian. They cannot be assumed to solve a coupled Schrödinger
system correctly. Midpoint interfaces are also inconsistent.

These files need contract and convergence work before promotion to active choices.
More details and repair priorities are in [status](11-status-and-limitations.md).
