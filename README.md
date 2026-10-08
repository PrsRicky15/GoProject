# QuantumMLD / GoProject

**Start with the [complete project documentation](docs/README.md).** The linked
guides describe the architecture, setup, mathematics, numerical packages, frontend,
HTTP API, testing, source files, and known incomplete methods as of 2026-10-08.

The current application is a local model-system workspace. It does not implement
electronic structure, a chemical force field, GPU computation, cloud jobs,
accounts, or shared projects. See the
[status guide](docs/11-status-and-limitations.md) for the actual supported scope.

## Quantum DVR Package

### Overview
A Go library for quantum mechanical calculations using Discrete Variable Representation (DVR) methods.
The numerical packages provide kinetic energy, momentum and Hamiltonian operators,
matrix diagonalization, time propagation, and classical model calculations.

### DVR Grid System:

- Operator Interfaces:- Unified interface for quantum operators
- Diagonalization :- Efficient eigenvalue/eigenvector computation
- Kinetic Energy  :- Multiple representations (DVR, canonical)
- Momentum        :- Momentum space calculations
- Hamiltonian     :- Build full quantum Hamiltonians

The [Hamiltonian guide](Quantum/README.md) describes the h2p-style operator
composition, DVR and Fourier constructors, and numerical verification.

### Web app

The Go server handles plot, quantum and classical calculations and serves the TypeScript/React app.
The frontend uses Bun 1.4.2. Build it, then start Go from the repository root:

```sh
cd webpage/web-app
bun install --frozen-lockfile
bun run build
cd ../..
go run ./cmd/web
```

Open `http://127.0.0.1:8080/#tools` for plotting, or
`http://127.0.0.1:8080/#calculations` for the quantum and classical lab.
The Go backend requires CGO and the FFTW3 development library for the Fourier basis.
If Bun and Go are installed in WSL, run these
commands in your WSL terminal. See the [web app README](webpage/web-app/README.md)
for development, API limits, and verification commands.

The root `main.go` is an empty experiment; use `./cmd/web` to start the application.
