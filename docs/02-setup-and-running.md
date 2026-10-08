# Setup, running and troubleshooting

[Documentation index](README.md)

## Requirements

| Dependency | Requirement / role |
| --- | --- |
| Go | `go.mod`: 1.25.1; recorded checks used 1.26.4 |
| Bun | `packageManager`: 1.4.2; installation, scripts, tests and builds |
| C compiler, CGO and FFTW3 development files | Native dependency of the Fourier package |
| Gonum | v0.16.0; linear algebra and BLAS helpers |
| go-fftw | Pinned pseudo-version in `go.mod`; FFTW bridge |
| Browser | Modern browser; WebGL for surfaces/isosurfaces |

The frontend uses React 19, Vite 7, TypeScript 6 and Plotly bundles 4.1.2. Exact
resolved versions are in `bun.lock`; package ranges are not an installed-version
inventory. `bunfig.toml` makes package scripts run using Bun.

On Debian/Ubuntu, FFTW development files are generally supplied by `libfftw3-dev`,
alongside compiler tooling. Install native dependencies in the same environment
where Go builds. The established workflow uses WSL; native Windows builds have not
been validated by the recorded checks. This guide does not install system packages.

## Workstation paths

These user-supplied paths are environment-specific, not application requirements:

```text
Windows repository: C:\Users\prash\PrashantProjects\GoWebserver\goProjects\GoProject
WSL repository: /mnt/c/Users/prash/PrashantProjects/GoWebserver/goProjects/GoProject
WSL Go: /mnt/c/Users/prash/Softwares/go/bin/go
WSL Bun: /home/prash/.bun/bin/bun
WSL distribution used: Ubuntu-26.04
```

In a WSL terminal:

```sh
export PATH=/home/prash/.bun/bin:/mnt/c/Users/prash/Softwares/go/bin:$PATH
cd /mnt/c/Users/prash/PrashantProjects/GoWebserver/goProjects/GoProject
go version
bun --version
go env CGO_ENABLED
```

## Build and serve

From the repository root:

```sh
cd webpage/web-app
bun install --frozen-lockfile
bun run build
cd ../..
go run ./cmd/web
```

| URL | Purpose |
| --- | --- |
| `http://127.0.0.1:8080/` | Landing page |
| `http://127.0.0.1:8080/#tools` | Calculator and plotter |
| `http://127.0.0.1:8080/#calculations` | Quantum/classical lab |
| `http://127.0.0.1:8080/api/health` | HTTP liveness |

Root `main.go` is empty. Use `./cmd/web`, not `go run .`. Stop the server with
Ctrl+C; it allows up to five seconds for graceful shutdown. Health confirms an
HTTP response, not successful execution of every solver.

An optional reusable binary:

```sh
go build -o /tmp/quantummld-web ./cmd/web
/tmp/quantummld-web -addr 127.0.0.1:8080 -dist webpage/web-app/dist
```

`-dist` defaults to a path relative to the process working directory. Run from the
root or pass an absolute asset path. A moved binary still needs FFTW runtime libraries.

## Vite development

Keep Go running in one terminal. In another:

```sh
cd webpage/web-app
bun run dev --host 127.0.0.1 --port 5177 --strictPort
```

Open `http://127.0.0.1:5177/#calculations`. Without a port flag, use Vite's printed
URL. The proxy targets 8080; change it deliberately if moving the backend port.
`bun run preview` serves the existing production build and also needs Go running.
Use `bun run build`, not `bun build`, to invoke this project's Vite build script.

## Verification commands

```sh
# Repository root
go test ./...
go vet ./...
go test -race ./Quantum ./OperatorAlgebra ./internal/simulation ./internal/webserver

# Frontend
cd webpage/web-app
bun run typecheck
bun run lint
bun run test
bun run build
```

## Troubleshooting

| Symptom | Action |
| --- | --- |
| Executable not found | Set PATH in the shell that starts the process, or use the supplied executable path |
| FFTW compile/link error | Check CGO, compiler, headers and matching library in the build environment |
| 503 requesting a frontend build | Build assets; verify working directory and `-dist` |
| Old controls appear | Rebuild production assets or restart Vite; reload and verify the port |
| Service unavailable | Check health, Go process, 8080 and the Vite proxy |
| Browser POST gets 403 | Check matching origin/proxy configuration rather than disabling protection |
| 429 busy | Finish/cancel the active job, then retry |
| Explicit dt rejected | Use the returned bound or choose split propagation |
| Morse potential too steep | Reduce domain, alpha or strength |
| Fluid lattice too dense | Increase box size or reduce particles; spacing must be at least 1.1 |
| Missing 3D plot | Check WebGL support; try 1D/2D to separate backend from graphics issues |
| Large-chunk warning | Expected for lazy Plotly bundles; not a build failure |

Serve only `dist` as static assets, not the repository root. The default server is
local and has no authentication. See [API protections](08-http-api.md) and
[validation evidence](10-testing.md).
