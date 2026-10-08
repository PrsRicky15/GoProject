# Testing, evidence and reproducibility

[Documentation index](README.md)

## Commands

Run from the repository root in the configured Go/FFTW environment:

```sh
go test ./...
go vet ./...
go test -race ./Quantum ./OperatorAlgebra ./internal/simulation ./internal/webserver
```

From `webpage/web-app`:

```sh
bun run typecheck
bun run lint
bun run test
bun run build
```

The latest implementation verification before these documentation-only edits
passed these checks on 2026-10-07 in WSL, using Go 1.26.4 and Bun 1.4.2. Bun ran
9 tests with 59 expectations. The production build passed with the expected large
Plotly chunk warning. This is recorded evidence, not a continuously updated CI badge.

The race detector covers Go-observable races in exercised paths; it is not a
proof of correctness of all C/FFTW internals or of safe shared-instance use.
No committed CI workflow, coverage target or browser end-to-end suite is present.

## Existing numerical tests

| Test source | What it checks |
| --- | --- |
| [grid_test.go](../gridData/grid_test.go) | Basic spatial/time constructors and file-input examples |
| [frequencies_test.go](../gridData/frequencies_test.go) | Odd/even FFT frequency ordering |
| [potential_test.go](../gridData/potential_test.go) | Despite its Gaussian test name, only a time-grid construction |
| [derivativeOp_test.go](../OperatorAlgebra/derivativeOp_test.go) | DVR matrix behavior covered by the existing test |
| [FourierBasis_test.go](../OperatorAlgebra/FourierBasis_test.go) | Momentum and kinetic operators on plane waves |
| [exponential_test.go](../OperatorAlgebra/exponential_test.go) | DVR diagonalization preserves K; real/complex eigenvector exponential, inverse and aliasing behavior |
| [hamiltonian_test.go](../Quantum/hamiltonian_test.go) | Matrix/application agreement, aliasing, harmonic spectrum, Fourier plane waves, time-dependent potential, snapshots and validation |
| [completed_methods_test.go](../ODESolver/completed_methods_test.go) | Ralston 2/3 and Nyström 5 scalar/grid agreement and global convergence order |
| [symplectic_test.go](../ODESolver/symplectic_test.go) | Yoshida order/grid agreement/redefinition; Störmer–Verlet initial velocity and repeat initialization |
| [lennard_jones_test.go](../classical/lennard_jones_test.go) | Seed repeatability, energy-gradient force agreement and dense-lattice rejection |
| [simulation_test.go](../internal/simulation/simulation_test.go) | End-to-end numerical adapters, quantum/classical references, MD diagnostics, limits and cancellation |

Tests are examples of verified properties. They do not cover all potentials,
asymmetric DVR domains, all RK coefficients, implicit/adaptive methods, every
parameter combination, long-time thermodynamic behavior or all failure paths.

## Analytic reference cases

### Harmonic spectrum

With m=k=1 and hbar=1, `E[n]=n+1/2`. The simulation test checks four low states in
both bases against 0.5, 1.5, 2.5 and 3.5 with 1e-6 tolerance, and residual below
1e-10. Agreement of low localized states does not validate box-sensitive excited states.

### Quantum motion

For a packet with x0=1 and p0=0 in that oscillator, the mean position at t=1 is
`cos(1)`. Tests run split, RK4, Ralston 3 and Nyström 5 with dt=0.002 and 500 steps;
mean-position tolerance is 3e-5 and norm drift must stay below 1e-6.

The driven test uses `Vdrive=A*x*sin(t)`, A=0.2, and zero initial mean position
and momentum. The mean obeys `x''+x=-A*sin(t)`, giving:

```text
<x>(t) = (A/2)*(t*cos(t)-sin(t))
```

At t=1 this is `0.1*(cos(1)-sin(1))`. The split test uses dt=0.005 and checks a
1e-5 position tolerance. These tests also exercise a genuinely time-dependent H.

### Scalar non-autonomous ODE

For `x'=x+t`, `x(0)=1`, the exact solution is `x(t)=2*exp(t)-t-1`.
The test compares steps of 0.2 and 0.1 at t=1. Observed global errors should decrease
by approximately `2^p` for order p in the asymptotic regime. It requires ratios of
at least 3, 6 and 24 for Ralston 2, Ralston 3 and Nyström 5 respectively, and
scalar/grid agreement within 1e-12.

### Classical trajectories and fluid

The classical test checks m=2 harmonic motion against `cos(t/sqrt(2))` with all
three exposed methods. Yoshida has a separate fourth-order convergence test.
The Morse test is a successful-run smoke test, not a general bond-spectrum validation.

The MD test uses 16 particles, box 6, initial temperature 0.2, seed 42, dt=0.001,
and 1,000 steps. It requires final momentum magnitude below 1e-10 and maximum
absolute energy drift below 1e-3. The force test differentiates the shifted energy
numerically and compares its negative gradient with the computed force.

## API and frontend tests

| File | Coverage |
| --- | --- |
| `internal/plotmath/plotmath_test.go` | Shared expression fixtures, rejected syntax, sampling, undefined values and bounds |
| `internal/webserver/server_test.go` | Plot success, null samples, no-store, malformed/oversized requests and origin protection |
| `internal/webserver/simulations_test.go` | Spectrum response and invalid content type, origin, fields, sizes and grid bounds |
| `src/lib/math.test.ts` | Calculator fixtures, invalid expressions, null-to-gap decoding, malformed result rejection, form normalization and plot orientation |
| `src/lib/simulation.test.ts` | Density dimensions, finite diagnostics and particle-frame dimensions |
| `testdata/math-expressions.json` | Shared arithmetic behavior between Go and TypeScript |

HTTP tests use an in-process test handler; they do not need a running server or
production build. Parser tests do not imply that arbitrary executable syntax is
supported—the rejection of such syntax is part of the expected behavior.

## Recorded browser checks

The built Go-served page was exercised for all five calculation kinds. Checks
observed harmonic levels and rendered densities, driven wave-packet results,
Yoshida classical results, Morse bond motion and MD play/pause advancing frame
time. Browser error/warning logs were empty for the inspected calculation runs.
These were manual automation checks, not committed regression scripts.

For future frontend changes, verify loading, successful rendering, parameter edits
marking old results, a representative error, cancel behavior, route cleanup, and
particle play/pause. Verify 1D/2D plots as well as WebGL-based plots when changing
shared `PlotCanvas`. Also inspect narrow layouts rather than relying on CSS alone.

## Reproducible scientific use

Record the complete request, code revision, runtime versions and relevant native
library environment. For quantum results record basis, L, N, potential and state
parameters, dt, steps and drive. For MD also record particle count, box, temperature
and seed. UI result parameters are available in a disclosure, but no durable run
export/database exists.

Perform three distinct checks:

1. **Algebraic:** eigenpair residual or force-gradient agreement.
2. **Numerical:** smaller dt, finer grid, larger domain and longer sampling.
3. **Physical:** whether the model, units and boundary conditions represent the problem.

Norm conservation alone does not establish correct phases. Energy conservation
alone does not establish an appropriate force field. Small eigenpair residuals
alone do not establish continuum convergence.
