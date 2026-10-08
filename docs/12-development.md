# Development and maintenance guide

[Documentation index](README.md)

## Preserve the current boundaries

Keep numerical models independent of HTTP and React. Go library packages define
operators/models; `internal/simulation` adapts bounded web jobs; `internal/webserver`
handles transport; TypeScript validates and presents results. Put only small
calculator arithmetic and visualization transformations in the frontend.

Use the existing implementation as a starting point, but read the status guide
before reusing legacy methods. Preserve user-owned changes when making edits.
Do not infer that a passing build certifies an unused numerical routine.

## Add a potential

1. Implement `gridData.PotentialOp` with a clear energy formula and `F=-dV/dx`.
2. Define parameter domains and behavior at singularities; do not rely on UI validation alone.
3. Test real energy/force against analytic values and a finite-difference gradient.
4. If supporting complex inputs, separately test analytic continuation and real-axis agreement.
5. Add its selection and construction in `internal/simulation`, with bounded parameter validation.
6. Add frontend selection, labels, formula and defaults; document units and meaning.
7. Test a representative trajectory/spectrum and update relevant docs.

No arbitrary Go source or browser JavaScript should be accepted as a potential.
For a future user-defined potential, design a bounded expression contract and
determine how its force will be derived and validated.

## Add a kinetic basis

`Quantum/kinetic.go` defines the internal `kineticEnergy` contract: real/complex
application, an independent dense matrix, and cleanup. A basis needs consistent
mass factors, signs, endpoint conventions and boundary conditions.

Add a checked Hamiltonian constructor. Test plane waves or suitable reference
states, matrix/application agreement, symmetry/Hermiticity, aliasing, convergence
and resource cleanup. Only then extend the request validator and basis selector.
Do not expose the current finite-difference file merely by adding a dropdown value.

Dense matrix construction must remain optional for operator-only use. If the new
basis needs nonlocal/time-dependent operators beyond T+local V, make that change
explicit rather than hiding it in a potential callback.

## Add a propagator

For a fixed-step vector method, specify whether success advances exactly dt,
whether failure leaves the input unchanged, and which buffers may alias. Calculate
coupled stages from complete intermediate states. Use a non-autonomous reference
problem so stage-time mistakes are visible, and measure the convergence ratio
under step halving. Test a coupled oscillator or Schrödinger state as well.

For an adaptive method, define accepted time and next-step advice before exposing
it. For an implicit method, supply a genuine coupled Jacobian/linear solve or
explicitly restrict the method to independent scalar equations. State any iteration
limits, tolerances and failure semantics.

Register a quantum method in `internal/simulation/quantum.go`, update validation
in `simulation.go`, expose a matching select option, choose a usable initial dt,
and document stability/accuracy limitations. Never silently renormalize a state
to make a conservation diagnostic look better.

## Add a calculation kind

The coordinated edit points are:

| Layer | Required changes |
| --- | --- |
| Go request | Fields, allowed values, numeric domains, job budget |
| Go dispatch | Route kind to a bounded context-aware numerical function |
| Result | Finite metrics and charts with explicitly defined units/axes |
| TypeScript | `CalculationKind`, request/result types, defaults and decoder bounds |
| UI | Kind button, relevant controls, results presentation and explanatory notes |
| Tests | Analytic/invariant reference, failure cases, API and decoder agreement |
| Docs | API, science guide, status and examples |

Avoid using a visible plot title as a stable machine identifier. If results need
export, comparison or long-term persistence, add explicit model/schema versioning.

## Extend plots and responses

Go plot modes and frontend `MODES` must agree on dimensions, axis ordering and
limits. Preserve x-fastest flattening unless the API contract deliberately changes.
New data fields need runtime validation, not just an interface declaration.
If extending isosurface support, keep Plotly-specific type workarounds localized.

Protect scientific meaning in the UI: distinguish density from amplitude, momentum
from velocity, initial temperature from controlled temperature, and sampled
playback time from physical time. Show units and submitted settings beside results.

## Performance work worth measuring

- Reuse stages/work buffers in generic explicit RK steps to reduce allocation pressure.
- Implement a Fourier-specialized split kinetic exponential for larger periodic grids.
- Consider partial/sparse eigensolvers when only a few states are needed at larger N.
- Add neighbor or cell lists before scaling MD far beyond its current particle cap.
- Reduce serialization and frontend copies if large plots dominate response time.
- Consider static result caching only with an explicit key, memory budget and invalidation policy.

Measure before and after using representative parameters and allocation profiles.
Preserve conservation/convergence checks; a faster implementation with altered
normalization, boundary behavior or unstable timesteps is not equivalent.

## Dependencies and generated files

Use Bun for frontend package changes and commit the resulting `bun.lock` with the
manifest. Use `bun install --frozen-lockfile` for reproducible installation. Do not
reintroduce a separate npm lockfile. `node_modules` and `dist` are generated assets,
not hand-edited source. TypeScript's `noEmit` check runs before Vite builds.

Go dependencies are recorded in `go.mod`/`go.sum`; FFTW is also a native runtime/build
dependency. Upgrading the wrapper or native library requires operator and lifecycle
tests, not only compilation. Avoid hardcoding workstation paths in application code.

## Validation and documentation workflow

Run checks relevant to the change, followed by required build checks. Numerical
changes need meaningful reference/convergence tests; wording/CSS changes generally
do not need implementation-mirroring unit tests. Inspect important UI behavior in
the running app and restart Go when backend code changes.

Update the documentation index if adding a guide. Keep API fields/defaults/limits
consistent with the validator and TypeScript modules. Update the status table
when a method becomes verified; retain the test that establishes the improvement.
The root README should remain a clear entry point rather than duplicating every
reference detail.

There is currently no committed CI pipeline or automated Markdown checker. For
documentation changes, verify relative links, JSON examples, code signatures and
whether examples use the correct working directory and runtime environment.
