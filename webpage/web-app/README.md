# Quantum chemistry web app

The [complete project notes](../../docs/README.md) cover the numerical backend,
[HTTP API](../../docs/08-http-api.md), [frontend architecture](../../docs/09-frontend.md),
and [unfinished methods](../../docs/11-status-and-limitations.md).

Go handles expression evaluation and grid sampling for every plot. The React
frontend is written in strict TypeScript and uses Vite and Bun for development,
production builds, tests, and linting. Simple calculator arithmetic stays in
TypeScript; Plotly renders the returned data in the browser.

## Requirements

Use Bun 1.4.2, the version recorded in `package.json`. If Bun is installed in WSL,
run the following commands from your WSL terminal.

`bunfig.toml` enables Bun for all package scripts, including tools with Node.js
shebangs. See [Bun's runtime configuration](https://bun.com/docs/runtime/bunfig).

## Development

Start the Go API from the repository root (Go 1.25.1 or newer, CGO enabled,
C compiler and FFTW3 development library available):

```sh
go run ./cmd/web
```

In a second terminal, from the repository root:

```sh
cd webpage/web-app
bun install --frozen-lockfile
bun run dev
```

Open the local URL printed by Vite, followed by `/#tools`. Vite proxies `/api`
to Go at `127.0.0.1:8080`, so the browser uses the same origin for both UI and API.
Changes to React components reload automatically. Restart Go after backend edits.

For a production build, run `bun run build` in this directory, then run
`go run ./cmd/web` from the repository root and open `http://127.0.0.1:8080/#tools`.
Go serves both the built app and its API; Bun is not a production server.
The command accepts `-addr` and `-dist` flags and binds to loopback by default.

## Build and check

Run these commands from `webpage/web-app`:

```sh
bun run lint
bun run test
bun run typecheck
bun run build
bun run preview
```

The build is written to `dist/`. The preview command serves that production build
locally. Use `bun run build` to run the Vite build script; `bun build` invokes Bun's
own bundler. `build` runs strict TypeScript checks before bundling. Keep Go running
when using Vite preview as well. From the repository root, run `go test ./...`
and `go vet ./...` for the backend. Go and TypeScript share expression fixtures in
`testdata/math-expressions.json` to check compatible arithmetic conventions.

## Dependencies

Commit `bun.lock` with dependency changes. Use `bun add <package>` or
`bun remove <package>` to change dependencies, and `bun install --frozen-lockfile`
for a reproducible installation from the committed lockfile.

## Plotter and calculator

Open `/#tools`, or select **Plotter & Calculator** / **Try Yourself** on the home
page. Choose a mode, edit the expression and domain, and press **Plot expression**.
Changes are explicit: editing a field does not resample until you press the button.
Each mode keeps its draft settings while the workspace remains open.

| Mode | Meaning | Example | Maximum samples per axis |
| --- | --- | --- | --- |
| 1D curve | y = f(x) | `exp(-x^2/2)` | 10,000 |
| 2D field | Heatmap or contour of f(x,y) | `sin(x)*cos(y)` | 200 |
| 3D surface | Height z = f(x,y) | `0.5*(x^2+y^2)` | 150 |
| Isosurface | Level set f(x,y,z) = c | `x^2+y^2+z^2`, level 1 | 60 |

Drag plots to pan or rotate. The plot toolbar provides zoom/reset and PNG export.
3D rendering requires WebGL. Isosurfaces require finite values throughout the
sampled domain and a level strictly inside the sampled range. Empty level sets
or features smaller than a grid cell may require changing the domain/resolution.
Undefined curve/field values appear as gaps; poles between sample points can
still be missed. All values are real numbers, and angles are radians.

Expressions support `+ - * / % ^` (`**` also works), parentheses, `pi`, `e`, and
the functions listed in the calculator's help. Use explicit multiplication.
`log` is base 10; `ln` is the natural logarithm. `%` is remainder, not percent.
The calculator supports `ans`, square/reciprocal of the complete expression, and
eight recent calculations. Invalid calculations keep the last valid answer.
Expression parsing uses a restricted arithmetic grammar, without `eval` or
`Function`. The old `/calculator/` URL redirects to the TypeScript workspace
when served by Go. The standalone JavaScript calculator has been removed.

Sampling runs in Go through `POST /api/plot`. Browser cancellation aborts the
request, and Go checks the request context while sampling. Go independently
enforces an 8 KiB request limit, a 1,024-byte/256-token expression limit, maximum
grid sizes, two concurrent calculations, and a five-second sampling deadline.
The handler rejects unknown JSON fields, malformed values and cross-origin
browser POSTs using Go's `http.CrossOriginProtection`. Undefined samples are
JSON nulls and become plot gaps. TypeScript checks response sizes, coordinate
arrays and numeric values at runtime; type annotations alone are not validation.
`GET /api/health` reports the Go engine status.

The main app and workspace load
separately, and the Plotly Cartesian and GL3D bundles load only when that family
of plots is used. Expect Vite's large-chunk warning for these vendor bundles
(roughly 508 KB and 571 KB gzipped); neither is part of the landing page bundle.
Plots update using `Plotly.react`; resize observers, requests, and graphics
resources are cleaned up when leaving the workspace. See Plotly's
[bundle documentation](https://github.com/plotly/plotly.js/blob/main/dist/README.md)
and [update API](https://plotly.com/javascript/plotlyjs-function-reference/).

Plot formulas are sent to the Go server; calculator expressions remain in the
browser. The default server is local. Hamiltonian and dynamics calculations use
the separate calculation API described below. This local application has
no user authentication; a public deployment needs its own access controls.

TypeScript compiles to JavaScript for browsers. Security here comes from the
restricted expression grammar, backend validation, bounded work, same-origin
requests and validated responses, rather than the source language alone.

## Quantum and classical lab

Open `/#calculations` or select **Calculations** in the navigation. Choose a model,
edit its parameters, and press **Run calculation**. Changing controls leaves the
previous result visible and marks it as a previous run. Results include parameters,
units, diagnostics, plots, and particle playback where applicable.

| Mode | Go implementation | Output |
| --- | --- | --- |
| Energy levels | Existing Hamiltonian composition in sinc DVR or periodic Fourier basis; symmetric eigensolver | Energies, eigenstate densities, eigenpair residuals |
| Quantum dynamics | Midpoint Strang splitting, RK4, completed Ralston 3 or Nyström 5 grid integration | Density versus position/time, norm, mean position, energy |
| Classical dynamics | Velocity-Verlet, leapfrog or corrected fourth-order Yoshida | Position, phase space, energy drift |
| Bond vibration | Morse potential and reduced mass, with the classical propagators | One-dimensional bond motion |
| Molecular dynamics | New 2D periodic Lennard–Jones model, velocity-Verlet | Particle playback, total energy and drift, temperature, momentum |

Quantum and one-particle models use atomic units with hbar=1. The quantum field
adds `A*x*sin(omega*t)` to the selected potential. Midpoint splitting caches a
dense kinetic eigensystem once per run, then applies its unitary exponential;
this small-grid implementation is not an FFT split-step method. Explicit RK
methods enforce a conservative step bound and report norm drift without hiding
it through renormalization. The grid is `[-L,L)` with `dx=2L/N`. Fourier boundaries
are periodic; the symmetric-domain sinc DVR matrix is truncated from the full-line
operator. There is no absorber. Check convergence with smaller dt and larger N/L.

The fluid uses reduced units `m=sigma=epsilon=kB=1`, a shifted-force cutoff at
`min(2.5, box/2)`, minimum-image forces and periodic wrapping. Seeded velocities
start with zero center-of-mass momentum and temperature defined with `2N-2`
degrees of freedom. Initial temperature is an initial condition, not a thermostat;
subsequent evolution is NVE. The model is monatomic and does not represent a
chemical force field. The Morse mode is a separate single bond-coordinate model.

`POST /api/simulations` accepts typed model parameters, with strict JSON validation,
an 8 KiB request limit, one active simulation per server handler, and a ten-second
deadline. Limits are 32–192 quantum grid points, up to 8 displayed eigenstates,
2,000 time steps, 40 million `N*N*steps` for quantum propagation, and 4–64 fluid
particles. Integration checks cancellation between steps; the bounded dense
eigensolver cannot be interrupted mid-factorization. The frontend validates
response dimensions and finite values before rendering.

Completed legacy methods include DVR `ExpIdt`, `ExpIdtTo`, `ExpDtTo`, non-mutating
real diagonalization, and Ralston/Nyström grid steps. Scalar Ralston/Nyström stage
errors, Yoshida drift coefficients and Störmer–Verlet initialization were corrected.
The DVR exponential API preserves its original `exp(+i*dt*K)` sign; the new
quantum propagator uses the Schrödinger `exp(-i*dt*H)` convention.
Legacy adaptive, implicit and finite-difference implementations are not exposed
on this page and still need separate completion and verification.

Tests compare harmonic spectra and motion with analytic values, driven quantum
motion with the resonant solution, RK/Yoshida convergence orders, scalar/grid
agreement, DVR eigenvector exponentials, Lennard–Jones force gradients, energy
and momentum conservation, and malformed API/frontend responses.
