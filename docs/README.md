# Complete project notes

These notes describe **QuantumMLD**, whose Go module is named `GoProject`.
Prepared from the working source and finalized on **2026-10-08**, they cover the active app,
numerical libraries, tests, configuration, and older experimental files.

## Contents and reading order

| Document | What it defines |
| --- | --- |
| [01 — Architecture](01-architecture.md) | Purpose, component boundaries, data flow, ownership and performance |
| [02 — Setup and running](02-setup-and-running.md) | Dependencies, WSL paths, execution, build and troubleshooting |
| [03 — Grids and potentials](03-grids-and-potentials.md) | Grid conventions, FFT frequencies, potential formulas, units and input files |
| [04 — Quantum mechanics](04-quantum.md) | DVR/Fourier operators, Hamiltonian API, eigenstates and time propagation |
| [05 — ODE solvers](05-ode-solvers.md) | Method inventory, contracts, explicit RK and symplectic algorithms |
| [06 — Classical mechanics and MD](06-classical-and-md.md) | Trajectories, Morse bonds, Lennard–Jones forces, integration and diagnostics |
| [07 — Plotter and calculator](07-plotter-and-calculator.md) | Expression syntax, visualization modes, indexing and invalid samples |
| [08 — HTTP API](08-http-api.md) | Endpoints, full fields, bounds, examples, response formats and errors |
| [09 — Frontend](09-frontend.md) | Routes, components, state, request lifecycle, rendering and playback |
| [10 — Testing](10-testing.md) | Commands, reference solutions, coverage limits and reproducibility |
| [11 — Status and limitations](11-status-and-limitations.md) | Working features, partial methods, observed defects and missing features |
| [12 — Development](12-development.md) | Adding features, optimization opportunities and maintenance |
| [13 — Repository map](13-repository-map.md) | Source-file responsibilities, tests, configuration and legacy inventory |

For a first run, read 01 and 02, then the guide for your calculation. For numerical
development, read 03–06, 10 and 11. For web development, read 07–09 and 12.

## Vocabulary

| Term | Meaning here |
| --- | --- |
| Basis | Representation of the kinetic operator on a 1D spatial grid |
| DVR | Discrete variable representation; the page uses a symmetric-domain sinc matrix |
| Fourier | Periodic spectral kinetic operator using FFTW |
| Hamiltonian | `H = T + V`, applied to a real or complex grid vector |
| Propagator | Numerical method advancing a state by one time step |
| Grid step | A coupled-vector ODE update; not necessarily a spatial derivative |
| Eigenstate density | Squared wavefunction amplitude with unit grid integral |
| NVE | Evolution at fixed particle count and box without a thermostat |
| Reduced LJ units | Particle mass, sigma, epsilon and Boltzmann constant set to one |
| Atomic units | Consistent quantum/one-coordinate inputs with hbar set to one |

## Status conventions

- **Active:** reachable through the UI/API and exercised by the stated tests.
- **Library only:** callable from Go but not selectable in the UI; verification varies.
- **Partial/experimental:** missing implementations, contracts or numerical verification.
- **Known defect:** a concrete source-level issue described in the status guide.

A passing build proves compilation and only the tests actually present. It does
not prove convergence, stability or physical accuracy for arbitrary inputs.
Equations in these notes explain the code rather than a separate proposed design.
Field definitions come from the Go request structs and validators; frontend defaults
come from the TypeScript library files. Update notes alongside future changes.

[Project overview](../README.md) · [Hamiltonian quick reference](../Quantum/README.md)
· [Frontend quick reference](../webpage/web-app/README.md)
