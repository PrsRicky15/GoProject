# Whole-vector ODEs and implicit Hamiltonian propagation

[Documentation index](README.md) · [Implicit methods](05-ode-solvers.md)

For all 13 fixed explicit formulas and six adaptive pairs on real/complex vectors,
see [Explicit and adaptive vector propagation](17-explicit-and-adaptive-vector-propagation.md).

Added 2026-10-12. All existing scalar solvers, `NextStepOnGrid`, explicit methods
and split propagators are retained. These new Go library APIs accept entire real
or complex vectors without requiring a dummy scalar callback. Web menus are
unchanged.

## What “vector” means

For a coupled equation `dy/dt=F(t,y)`, every RHS call receives the **whole** vector
at the current stage/iterate. A component may depend on every other component.
It is not a loop that calls a scalar `f(y[i],t)` independently.

The existing `NextStepOnGrid` engine already handles this through
`EvaluateOnRGridInPlace`; its Newton path includes the full Jacobian. The new
interfaces remove the need to implement scalar methods or an allocating grid
method. They also add direct complex states and matrix-free linear solves.

| API | State | Equation | Solve strategy |
| --- | --- | --- | --- |
| `NewVectorImplicit` | `[]float64` | General coupled `F(t,y)` | Existing damped modified Newton, full real Jacobian |
| `NewComplexImplicit` | `[]complex128` | General coupled `F(t,psi)` | Real block representation, full real Jacobian |
| `NewComplexLinearImplicit` | `[]complex128` | **Linear** `A(t)psi` | Restarted GMRES with operator-vector products |
| `Quantum.NewImplicitPropagator` | `[]complex128` | 1D `-i H(t)psi` | Matrix-free implicit midpoint |
| `Quantum.NewNDimImplicitPropagator` | `[]complex128` | nD `-i H(t)psi` | Same matrix-free midpoint on the full tensor vector |

All `Step(state,t,dt)` calls commit to `state` only after a successful solve.
The caller advances time after success. Solvers own reusable buffers and must not
be copied after use or shared concurrently. RHS callbacks must fill all output
entries, treat input as read-only, and not retain slices. Callback errors are
wrapped and preserved for `errors.Is`. Callbacks should return the same result
for the same time and vector during a step.

## Real coupled equations

```go
import EquationSolver "GoProject/ODESolver"

rhs := func(t float64, y, dy []float64) error {
    // y=[position, velocity]; both components are read at the same iterate.
    dy[0] = y[1]
    dy[1] = -y[0]
    return nil
}
solver, err := EquationSolver.NewVectorImplicit(
    rhs, EquationSolver.ImplicitMidpointMethod,
    EquationSolver.ImplicitOptions{},
)
if err != nil { return err }
state := []float64{1, 0}
if err := solver.Step(state, 0, 0.01); err != nil { return err }
```

Available method constants are `ImplicitEulerMethod`,
`ImplicitTrapezoidalMethod`, `ImplicitMidpointMethod`, and `ImplicitSDIRK2Method`.
These new general vector APIs use Newton. The original fixed-point types remain
available through their unchanged grid APIs.

The full analytic Jacobian callback is optional. Without it, the solver perturbs
each component and reevaluates the complete vector. The original dense limit of
1,024 real components applies. The diagonal option is only valid for truly
independent equations and does not apply to the oscillator above.

## General complex, possibly nonlinear equations

```go
rhs := func(t float64, psi, dpsi []complex128) error {
    // Example: both derivatives depend on the entire two-entry state.
    dpsi[0] = -1i*psi[1]
    dpsi[1] = -1i*psi[0]
    return nil
}
solver, err := EquationSolver.NewComplexImplicit(
    rhs, EquationSolver.ImplicitMidpointMethod,
    EquationSolver.ImplicitOptions{},
)
if err != nil { return err }
psi := []complex128{1, 0}
if err := solver.Step(psi, 0, 0.01); err != nil { return err }
```

Internally, the vector is packed as `[Re(psi), Im(psi)]`. This supports nonlinear
coupling and non-holomorphic operations such as conjugation and `abs(psi)^2`.
A supplied `ImplicitOptions.Jacobian` is the **2n-by-2n real Jacobian** in that
ordering; a scalar derivative or diagonal callback is rejected. The dense limit
is 512 complex entries. These are nonlinear residual tolerances, not adaptive
time-step tolerances.

## Large linear vector equations

For `psi'=A(t)psi`, pass a callback that applies `A(t)` to the whole vector:

```go
solver, err := EquationSolver.NewComplexLinearImplicit(
    len(psi), rhs, EquationSolver.ImplicitMidpointMethod,
    EquationSolver.KrylovOptions{
        AbsoluteTolerance: 1e-12,
        RelativeTolerance: 1e-10,
        Restart: 30,
        MaxIterations: 300,
    },
)
if err != nil { return err }
if err := solver.Step(psi, time, dt); err != nil { return err }
time += dt
```

The RHS here **must be linear in the state**. Time dependence is allowed. A
state-dependent Hamiltonian or nonlinear Schrödinger equation must instead use
the nonlinear API (or a future Newton-Krylov solver); ordinary GMRES does not
solve nonlinear stage equations.

The same four method constants are supported. At each implicit stage, GMRES
solves `(I-c*A(tStage))*y=b` using only applications of `A`. Arnoldi uses modified
Gram-Schmidt with a second orthogonalization pass and complex Givens rotations.
Acceptance recomputes the true residual, including after restarts. The design
follows the [restarted GMRES algorithm in Netlib's Templates](https://www.netlib.org/linalg/old_html_templates/subsection2.6.3.4.html).

Defaults are shown above; restart is capped at 256 and workspace at 256 MiB.
The basis uses O(n*restart) memory, rather than a dense O(n²) Hamiltonian/Jacobian.
The vector dimension is fixed at construction. Buffers are reused: the test with
2,048 coupled entries performs zero solver allocations per warm step. Callback
and Hamiltonian allocations are outside that guarantee.

Convergence requires `norm(b-(I-c*A)y) <= atol + rtol*norm(b)` for each stage.
`Stats()` reports iterations, operator applications, restarts and the last
stage's true residual. `MaxIterations` applies separately to each stage.
`StepContext(ctx,...)` checks cancellation between operator applications and
before committing. A long individual callback cannot be interrupted internally.

There is currently no preconditioner, recycled Krylov basis, or automatic time-step
adaptation. Difficult spectra may need a smaller dt, a larger restart/iteration
limit, or a future problem-specific preconditioner. Failed solves return errors
without changing the original state. Do not infer temporal accuracy from a small
linear residual: also refine dt.

## Use your existing Hamiltonians directly

Both adapters use the Hamiltonian's existing `ApplyHComplexAt` and the complete
wavefunction vector, including n-dimensional coupling. Atomic units use hbar=1.

```go
// h is an existing *Quantum.NDimHamiltonian, DVR or Fourier.
propagator, err := Quantum.NewNDimImplicitPropagator(
    h, EquationSolver.KrylovOptions{},
)
if err != nil { return err }

for step := 0; step < steps; step++ {
    time := float64(step)*dt
    if err := propagator.StepContext(ctx, psi, time, dt); err != nil {
        return err
    }
}
```

For a 1D `*Quantum.HamiltonianOp`, use `Quantum.NewImplicitPropagator(h, options)`.
The Hamiltonian must stay open and must not be used concurrently with propagation.
For time dependence, use its existing `SetTimeDependentPotential` or nD
`SetTimeDependentFunction` setter; the callback retains its existing convention
of supplying the **complete** potential, not an additive correction.

The quantum adapter always uses implicit midpoint:

```text
Hmid = H(t + dt/2)
(I + i*dt*Hmid/2) psiNext = (I - i*dt*Hmid/2) psi
```

For static H this equals Crank-Nicolson. For time-varying H it differs from the
endpoint trapezoidal rule, which uses `H(t)` on the right and `H(t+dt)` on the left.
For Hermitian midpoint H, the Cayley transform preserves norm in exact arithmetic;
the implementation preserves it up to linear-solve and floating-point error.
It never rescales psi to hide drift. Midpoint is second order in time and has
phase error; it is not the exact matrix exponential. The dissipative Euler and
SDIRK methods are available in the general linear API but are not the quantum
adapter's default.

## Evidence and runnable examples

- [Vector examples](../ODESolver/vector_example_test.go): whole-vector oscillator
  and complex two-level system. Run
  `go test ./ODESolver -run 'ExampleNew(Vector|ComplexLinear)Implicit' -v`.
- [Vector tests](../ODESolver/vector_implicit_test.go): real coupled equations,
  nonlinear/non-holomorphic complex equations, all four linear methods, a
  noncommuting time-dependent two-level Hamiltonian with exact reference,
  GMRES restarts on 2,048 coupled entries, allocations and failure contracts.
- [Hamiltonian tests](../Quantum/implicit_propagation_test.go): DVR/Fourier
  superposition phase and norm, plus a 4D Fourier system of 2,048 entries with
  time-dependent potential, second-order convergence and cancellation.

```bash
go test ./ODESolver ./Quantum
go test -race ./ODESolver ./Quantum
```
