# Repository and source-file map

[Documentation index](README.md)

Paths below are relative to the repository root. This inventory describes source
and configuration, not generated dependency trees. Some spellings are historical
(`NumericalRecipies`, `moleulardyn`, `sympleticMethods`, `Rapson`) and are preserved
here so the paths match the files.

## Top level

| Path | Responsibility |
| --- | --- |
| [README.md](../README.md) | Project entry point and documentation link |
| [go.mod](../go.mod), [go.sum](../go.sum) | Go module `GoProject` and dependency checksums |
| [main.go](../main.go) | Empty legacy program; not the web entry point |
| [cmd/web/main.go](../cmd/web/main.go) | Server flags, timeouts, listen address and graceful shutdown |
| [cmd/ndquantum/main.go](../cmd/ndquantum/main.go) | Runnable n-dimensional coupled oscillator, JSON energies and optional density CSV |
| [testdata/math-expressions.json](../testdata/math-expressions.json) | Shared arithmetic fixtures |
| `docs/` | These full-project notes |
| `.gitignore`, `.gitattributes` | Repository exclusion/text handling rules |
| `.idea/` | Existing IDE metadata; not application runtime logic |

## gridData package

| Source | Responsibility |
| --- | --- |
| [grid.go](../gridData/grid.go) | Spatial grid, metadata, input files, evaluation/output convenience methods |
| [gridHigherDim.go](../gridData/gridHigherDim.go) | Restored NDimGrid, original tensor layout, general potential interface and checked file/indexing helpers |
| [gridHigherDim_test.go](../gridData/gridHigherDim_test.go) | Multidimensional layout, validation, ownership and file regressions |
| [nd_potentials.go](../gridData/nd_potentials.go) | Explicit radial/additive nD extensions of existing potential types |
| [timegrid.go](../gridData/timegrid.go) | Time grid, macro/micro setup and input file reader |
| [gridFunctionality.go](../gridData/gridFunctionality.go) | Shared spacing, coordinate/frequency generation, display and file output |
| [potential.go](../gridData/potential.go) | Generic real/complex potentials and forces |
| [timedependentPoten.go](../gridData/timedependentPoten.go) | Scalar/vector time-dependent RHS interface |
| [grid_test.go](../gridData/grid_test.go) | Basic constructors/readers |
| [frequencies_test.go](../gridData/frequencies_test.go) | FFT ordering regression tests |
| [potential_test.go](../gridData/potential_test.go) | Misnamed Gaussian test; currently time-grid assertions |
| `testGrid/rgrid.inp`, `testGrid/tgrid.inp` | Example reader input fixtures |

### NumericalRecipies subpackage

| Source | Responsibility / status |
| --- | --- |
| [numericalMethods.go](../gridData/NumericalRecipies/numericalMethods.go) | Scalar Chebyshev first-kind recurrence and incomplete generic variants |
| [integration.go](../gridData/NumericalRecipies/integration.go) | Integration interface only |
| [differentiation.go](../gridData/NumericalRecipies/differentiation.go) | Package declaration only |

## OperatorAlgebra package

| Source | Responsibility |
| --- | --- |
| [dvrBasis.go](../OperatorAlgebra/dvrBasis.go) | DVR kinetic matrix/cache, complex scaling, real diagonalization and exponentials |
| [FourierBasis.go](../OperatorAlgebra/FourierBasis.go) | FFTW lifecycle and momentum/kinetic application |
| [fourier_multiplier.go](../OperatorAlgebra/fourier_multiplier.go) | Normalized complex spectral multiplier for kinetic propagation |
| [FiniteDifference.go](../OperatorAlgebra/FiniteDifference.go) | Incomplete/inconsistent finite-difference kinetic stencils |
| [derivativeOp.go](../OperatorAlgebra/derivativeOp.go) | Legacy operator interfaces; some signatures no longer match concrete types |
| [matrixoperation.go](../OperatorAlgebra/matrixoperation.go) | Real symmetric diagonalization and legacy complex wrapper |
| [matrix_interface.go](../OperatorAlgebra/matrix_interface.go) | Legacy LAPACK bridge; unconfigured complex backend |
| [derivativeOp_test.go](../OperatorAlgebra/derivativeOp_test.go) | DVR checks |
| [FourierBasis_test.go](../OperatorAlgebra/FourierBasis_test.go) | Fourier plane-wave checks |
| [exponential_test.go](../OperatorAlgebra/exponential_test.go) | Completed DVR exponential and cache-preservation checks |

## Quantum package

| Source | Responsibility |
| --- | --- |
| [hamiltonian.go](../Quantum/hamiltonian.go) | Checked constructors, T+V composition, potentials, matrix/eigen APIs and energy |
| [kinetic.go](../Quantum/kinetic.go) | Private DVR/Fourier adapters for Hamiltonian application |
| [propagation.go](../Quantum/propagation.go) | Cached dense kinetic eigensystem and midpoint Strang steps |
| [nd_hamiltonian.go](../Quantum/nd_hamiltonian.go) | Cartesian tensor DVR/Fourier Hamiltonian, coupled/time-dependent potentials and bounded dense operations |
| [nd_eigensolver.go](../Quantum/nd_eigensolver.go) | Restarted block Krylov/Rayleigh–Ritz low-energy solver, residual checks and cancellation |
| [nd_propagation.go](../Quantum/nd_propagation.go) | Axis-wise Strang propagation with Fourier phases or one-axis DVR eigensystems |
| [implicit_propagation.go](../Quantum/implicit_propagation.go) | Matrix-free complex midpoint adapters for 1D and nD Hamiltonians |
| [rk_propagation.go](../Quantum/rk_propagation.go) | Whole-complex-vector fixed/adaptive RK adapters for static/time-dependent 1D/nD Hamiltonians |
| [rk_tableau.go](../ODESolver/rk_tableau.go), [vector_rk.go](../ODESolver/vector_rk.go) | Catalog of 19 RK formulas, reusable real/complex engines and complete adaptive intervals |
| [vector_rk_test.go](../ODESolver/vector_rk_test.go), [rk_example_test.go](../ODESolver/rk_example_test.go), [rk_propagation_test.go](../Quantum/rk_propagation_test.go) | Convergence, adaptive control, rollback, reuse, examples and Hamiltonian references |
| [implicit_propagation_test.go](../Quantum/implicit_propagation_test.go) | Static superpositions, phase/norm checks and 4D time-dependent propagation |
| [nd_hamiltonian_test.go](../Quantum/nd_hamiltonian_test.go) | Coupled/anisotropic/degenerate spectra, plane waves, propagation convergence and compatibility |
| [nd_clustered_test.go](../Quantum/nd_clustered_test.go) | Localized 3D well spectrum compared with a full dense reference |
| [hamiltonian_test.go](../Quantum/hamiltonian_test.go) | Operators, spectra, time dependence, aliasing and validation |
| [README.md](../Quantum/README.md) | Focused package quick reference and executable example |

## ODESolver folder / EquationSolver package

| Source | Responsibility |
| --- | --- |
| [differentialEqnSolver.go](../ODESolver/differentialEqnSolver.go) | ODE interface, shared constants, explicit scalar/grid methods |
| [explicitGrid.go](../ODESolver/explicitGrid.go) | Generic full-vector RK stepping and completed methods |
| [adaptiveMethods.go](../ODESolver/adaptiveMethods.go) | Experimental adaptive embedded RK implementations |
| [odeFixPoint.go](../ODESolver/odeFixPoint.go) | Complete Euler/trapezoidal/midpoint fixed-point methods and converged predictor/corrector |
| [odeNewtonRapson.go](../ODESolver/odeNewtonRapson.go) | Euler/trapezoidal/midpoint Newton methods and two-stage SDIRK2 |
| [implicit.go](../ODESolver/implicit.go) | Shared residual solve, options/diagnostics, damped modified Newton, coupled Jacobians, LU/diagonal solve and reusable buffers |
| [implicit_test.go](../ODESolver/implicit_test.go) | Implicit convergence, coupled stiff systems, errors, allocations and benchmarks |
| [vector_implicit.go](../ODESolver/vector_implicit.go) | Direct real/complex whole-vector nonlinear APIs reusing the implicit engine |
| [complex_linear.go](../ODESolver/complex_linear.go) | Matrix-free linear complex implicit methods and restarted GMRES |
| [vector_implicit_test.go](../ODESolver/vector_implicit_test.go) | Coupling, complex/nonlinear equations, noncommuting Hamiltonians, restarts and failure contracts |
| [vector_example_test.go](../ODESolver/vector_example_test.go) | Executable vector-only usage examples |
| [sympleticMethods.go](../ODESolver/sympleticMethods.go) | Störmer/velocity Verlet, leapfrog and Yoshida |
| [exOperator.go](../ODESolver/exOperator.go) | Legacy exponential interface |
| [completed_methods_test.go](../ODESolver/completed_methods_test.go) | Ralston/Nyström consistency and convergence |
| [symplectic_test.go](../ODESolver/symplectic_test.go) | Yoshida and Störmer–Verlet regression checks |

## classical and Plots packages

| Source | Responsibility |
| --- | --- |
| [classical/lennard_jones.go](../classical/lennard_jones.go) | Working 2D periodic model fluid |
| [classical/lennard_jones_test.go](../classical/lennard_jones_test.go) | Pair forces, deterministic setup and invalid density checks |
| [classical/gridsystem.go](../classical/gridsystem.go) | Legacy SystemInfo/ModelSystem data holders |
| [classical/numericGridCalc.go](../classical/numericGridCalc.go) | Legacy InitialValueIntegrator interface |
| [classical/moleulardyn.go](../classical/moleulardyn.go) | Empty package file; new MD lives in lennard_jones.go |
| [Plots/basicPlot.go](../Plots/basicPlot.go) | Empty package file; active plotting is Go sampling plus Plotly |

## internal packages

Go's `internal` directory rules restrict these packages to this module tree.

| Source | Responsibility |
| --- | --- |
| [plotmath/expression.go](../internal/plotmath/expression.go) | Bounded arithmetic lexer/parser and reusable evaluator |
| [plotmath/sampling.go](../internal/plotmath/sampling.go) | Plot schemas, mode limits, coordinate generation and sampling |
| [plotmath/plotmath_test.go](../internal/plotmath/plotmath_test.go) | Parser/sampler tests |
| [simulation/simulation.go](../internal/simulation/simulation.go) | Common schemas, validation, potential construction, dispatch and finite-output checks |
| [simulation/quantum.go](../internal/simulation/quantum.go) | Spectrum and propagation jobs; direct whole-complex-vector evolution; retained legacy real/imaginary adapter |
| [simulation/propagators.go](../internal/simulation/propagators.go) | Web method catalog, adaptive/implicit dispatch, bounded RHS work and solver diagnostics |
| [simulation/propagators_test.go](../internal/simulation/propagators_test.go) | All 23 vector choices, coupled references, driven motion, rejection, limits and cancellation |
| [simulation/classical.go](../internal/simulation/classical.go) | Force/mass adapter and scalar/Morse trajectories |
| [simulation/md.go](../internal/simulation/md.go) | LJ job, sampling, metrics and result notes |
| [simulation/simulation_test.go](../internal/simulation/simulation_test.go) | Cross-package numerical reference tests |
| [webserver/server.go](../internal/webserver/server.go) | Plot endpoint, health, assets, headers and cross-origin protection |
| [webserver/simulations.go](../internal/webserver/simulations.go) | Simulation endpoint, capacity and request/deadline handling |
| [webserver/server_test.go](../internal/webserver/server_test.go) | Plot/transport tests |
| [webserver/simulations_test.go](../internal/webserver/simulations_test.go) | Simulation transport tests |

## Active frontend: webpage/web-app

| Path | Responsibility |
| --- | --- |
| `index.html` | Root element, metadata and TypeScript module entry |
| `package.json`, `bun.lock`, `bunfig.toml` | Scripts, frontend dependencies and Bun runtime configuration |
| `vite.config.ts` | React plugin and development/preview API proxies |
| `tsconfig.json` | Strict application type checks |
| `eslint.config.js` | TypeScript, React Hooks and refresh lint rules |
| `postcss.config.js`, `tailwind.config.js` | CSS tooling |
| `.gitignore` | Frontend generated-file exclusions |
| `README.md` | Frontend quick reference |
| `src/main.tsx`, `src/App.tsx` | React mount and hash routes |
| `src/index.css`, `src/App.css` | Global styles and retained starter stylesheet |
| `src/pages/HomePage.tsx` | Landing page |
| `src/pages/Workbench.tsx`, `src/pages/workbench.css` | Calculator/plotter workspace and shared lab appearance |
| `src/pages/Calculations.tsx`, `src/pages/calculations.css` | Numerical lab controls/results and responsive styling |
| `src/components/SolutionsSection.tsx` | Landing-page content |
| `src/components/Calculator.tsx` | Scientific calculator |
| `src/components/PlotCanvas.tsx` | Shared Plotly lifecycle |
| `src/components/SimulationChart.tsx` | Simulation chart presentation |
| `src/components/ParticlePlayback.tsx` | SVG particle playback |
| `src/lib/types.ts` | Plot data contracts |
| `src/lib/expression.ts` | Calculator parser |
| `src/lib/sampling.ts` | Plot defaults and request normalization |
| `src/lib/api.ts` | Plot HTTP client and decoder |
| `src/lib/figure.ts` | Plotly trace construction |
| `src/lib/simulation.ts` | Simulation defaults, contracts, decoder and client |
| `src/lib/math.test.ts`, `src/lib/simulation.test.ts` | Bun unit tests |
| `src/plotly.d.ts` | Partial-bundle module declarations |
| `node_modules/`, `dist/` | Generated dependency/build output; do not hand-edit |

`AppContent.jsx` and the former JSX entry/page files were replaced during the
TypeScript migration and are not the active source. The npm package lock and
standalone calculator JavaScript were removed in favor of Bun and the TS workspace.

## Legacy webpage files

`webpage/basic.html`, `basic.js`, `example2.html`, `example2.css`, `Harmonic.html`,
`style.css`, `light.css` and `dark.css` are older standalone demos/assets.
`webpage/calculator/index.html` and `style.css` retain the legacy calculator entry
and presentation. Go explicitly redirects `/calculator/` to the current tools page.
These siblings of `web-app` are not automatically copied into its `dist` build.
Their existence does not add endpoints or simulations to the active application.

For behavior rather than file locations, use the linked subject guides. For the
trustworthiness of a legacy implementation, consult the status guide before use.
