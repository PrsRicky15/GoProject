# Multidimensional grids and Hamiltonians

[Documentation index](README.md) · [Grid source](../gridData/gridHigherDim.go)
· [Hamiltonian source](../Quantum/nd_hamiltonian.go)

## Recovery and scope

The original `NDimGrid` was found in the Git history shared by
`origin/PrashantCommit` and `origin/saurabh_edits/webdev`; the version inspected
was `gridData/gridHigherDim.go` at `8cfc52f`. That file was absent from the current
`main` checkout. The earlier project notes described only that checkout and missed
this existing multidimensional design.

The public grid constructors, axis accessors, redefinition methods, indexing,
function/potential sampling, input reader and output writer are restored. The
original **last-coordinate-fastest** layout is preserved. Validation, ownership
and file-error handling are strengthened. The historical grid file did not by
itself supply an n-dimensional Hamiltonian: the Hamiltonian on that branch still
accepted `RadGrid`, and several `NDimEvaluateAt` methods were placeholders.

The new Go path connects the grid to Cartesian DVR and Fourier Hamiltonians,
coupled real potentials, low-energy eigenstates and time-dependent propagation.
Existing 1D constructors and web calculations continue to work.
**Calculations → N-dimensional states** and `/api/simulations` now expose 2D–4D
energy levels with model potentials and density projections. See the
[web guide](15-web-multidimensional-potentials.md). General nD propagation and
custom Go functions remain available through the Go API.
The expression plotter's 2D/3D/isosurface choices are separate visualization modes.

## Run a model

Run these from the project root in the WSL terminal where Go, a C compiler and
FFTW are installed. Bun is not needed for these Go calculations.

```sh
export PATH=/mnt/c/Users/prash/Softwares/go/bin:$PATH

# Three oscillators: E0=1.5; the next three levels are 2.5.
go run ./cmd/ndquantum -dims 3 -points 20 -states 4

# Coupled 2D potential V=(x²+y²)/2 + 0.3*x*y.
go run ./cmd/ndquantum -dims 2 -points 20 -coupling 0.3 -states 4

# Different grids, masses and frequencies per coordinate.
go run ./cmd/ndquantum -dims 2 -points 24,28 -lengths 12,14 \
  -masses 2,0.5 -frequencies 1,1.7 -basis dvr -states 3

# Optional coordinate/density export.
go run ./cmd/ndquantum -dims 2 -points 20 -states 4 -csv densities.csv
```

[Command source](../cmd/ndquantum/main.go). JSON output contains dimensions,
point counts, energies, residuals and iteration count. Failed convergence exits
with an error. `-tol`, `-iterations`, `-subspace` and `-timeout` control the solver;
defaults are `1e-9`, `500`, automatic, and `2m`. Ctrl+C cancels a run.

The example command supports 1–32 coordinates and this model:

```text
V(x) = sum_d m[d]*omega[d]²*x[d]²/2
       + coupling * sum_(d=1..D-1) x[d-1]*x[d]
```

The library has no fixed dimension limit and accepts general coupled functions.
Strong coupling can make the example's quadratic potential unconfined on infinite
space. A finite-box spectrum then does not describe bound levels of a confining
oscillator. Check the physical model and box convergence.

## Grid and original layout

`NewNDimGrid(mins,maxs,counts)` accepts one bound pair and count per axis.
`NewNDimGridUniform(d,min,max,n)` uses identical parameters on every axis.
`NewNDimGridFromLengths(lengths,counts)` centers each interval on zero.
Each axis samples `[min,max)` with `dx=(max-min)/count`. Only axis metadata is
stored; a table of every coordinate tuple is not allocated.

```text
N = product_d counts[d]
stride[d] = product_(j>d) counts[j]
flat = sum_d index[d]*stride[d]

shape [2,3,4] -> strides [12,4,1]
index [1,2,3] -> flat 23
```

| API | Contract |
| --- | --- |
| `NDim`, `TotalPoints`, `NPoints`, `Strides` | Shape, product and layout |
| `Min/Max/NumPoints/DeltaR/Length(d)` | Real-axis metadata |
| `DeltaK/KMin/KMax/CutoffE(d)` | Conjugate-axis metadata |
| `RValues(d)`, `KValues(d)` | Fresh slices; frequencies use the same FFT ordering as `RadGrid` |
| `GetPoint(flat)`, `GetFlatIndex(indices)` | Coordinate and index conversion |
| `PointTo(flat,out)` | Coordinate conversion without allocation |
| `VolumeElement()` | Product of all spacings |
| `Clone`, `Mins`, `Maxs`, `NPoints`, `Strides` | Independent copies |
| `ReDefine`, `ReDefineMinMax`, `ReDefineLengths` | Checked replacement; failure preserves the old grid |
| `FunctionOnGridInPlace(fn,out)` | Sample a function in flat-index order |
| `PotentialOnGrid`, `PotentialOnGridInPlace` | Sample an `NDimPotentialOp` |

Construction rejects inconsistent/empty inputs, invalid/nonfinite bounds, zero
counts, unrepresentable frequencies and total-point integer overflow. Invalid
point indices and wrong buffer lengths panic instead of silently wrapping or
leaving an output unfilled. Invalid `RValues` axes return errors; `KValues`
returns nil, preserving historical signatures. Scalar accessors require valid axes.

Callbacks must not retain or modify the reusable coordinate slice. Copy it when
coordinates need to be retained. Wavefunctions and potentials have length `N`,
not `sum(counts)`. Extremely scaled high-dimensional volume products may exceed
floating-point range; use consistent units and representable quadrature weights.

### File compatibility

`NewNDimGridFromFile(directory)` reads `ndgrid.inp`:

```text
nDim: 3
mins: -6,-6,-6
maxs: 6,6,6
nPoints: 20,20,20
```

Keys can occur in any order; blank lines and `#` comments are allowed. Duplicate
keys, malformed/missing numbers, unknown fields and mismatched list lengths return
errors. `PrintVectorToFile(vector,path,format)` requires exactly `N` values and
writes coordinates followed by the vector value, propagating I/O errors.

## General potentials

The nD interface requires only `NDimEvaluateAt([]float64) float64`.
`gridData.NDimPotentialFunc` adapts a function to that interface. Cross terms such
as `x[0]*x[1]` couple coordinates. Radial potentials can explicitly compute
`sqrt(sum(x[d]*x[d]))`. Separation into one-coordinate terms is not required.

`Harmonic[T].NDimEvaluateAt` sums the existing harmonic expression over coordinates.
`gridData.SeparablePotential{p0,p1,...}` constructs `sum_d p_d(x[d])` from existing
1D potential objects. Nil entries mean zero; a dimension mismatch panics.
Explicit additive Morse/polynomial/multi-Gaussian and radial Gaussian/soft-core/
super-Gaussian methods are now documented in the [potential guide](15-web-multidimensional-potentials.md).
The original
1D `PotentialOp` interface remains unchanged.

## Construct and solve in Go

```go
package main

import (
    "context"
    "fmt"
    "GoProject/Quantum"
    "GoProject/gridData"
)

func main() {
    grid, err := gridData.NewNDimGridFromLengths(
        []float64{12, 12}, []uint32{20, 20},
    )
    if err != nil { panic(err) }
    potential := gridData.NDimPotentialFunc(func(x []float64) float64 {
        return 0.5*(x[0]*x[0]+x[1]*x[1]) + 0.3*x[0]*x[1]
    })
    h, err := Quantum.NewNDimFourierHamiltonian(grid, []float64{1, 1}, potential)
    if err != nil { panic(err) }
    defer h.Close()
    result, err := h.LowestStates(context.Background(), 4, Quantum.EigenOptions{})
    if err != nil { panic(err) }
    fmt.Println(result.Values)
    fmt.Println(result.Residuals)
    // result.Vectors contains one eigenvector per column.
}
```

Use `NewNDimHamiltonian` for sinc DVR or `NewNDimFourierHamiltonian` for Fourier.
Supply one positive finite mass per axis; nil potential means a free particle.
Constructors snapshot the grid/potential; `Grid()` returns a clone. Reconstruct
the operator when changing grid, mass or basis.

```text
H = sum_d I ⊗ ... ⊗ T_d ⊗ ... ⊗ I + diag(V)
T_d = -1/(2*m[d]) * d²/dx[d]²            (hbar=1)
```

Large Kronecker matrices are not built during application. The code gathers a
coordinate line, applies the 1D kinetic operator and accumulates using the strides.
The potential has one value per full coordinate tuple.

For DVR, the kinetic-only grid is centered independently of physical coordinates
to avoid the legacy half-line heuristic. Moving a Cartesian box therefore does
not change the kinetic boundary convention. This is a truncated whole-line sinc
basis: wavefunctions should decay before the edges. Fourier axes are periodic.
Neither constructor supplies a radial/curvilinear metric, hard-wall sine basis,
mixed derivatives, vector potential, spin components or nonlocal potential.

## Eigenstates

`LowestStates(ctx,count,options)` uses restarted block Krylov/Rayleigh–Ritz
iteration. It starts with `count+2` deterministic random vectors (capped by `N`),
orthogonalizes twice, retains low Ritz vectors on restart and expands using
residuals. A block can represent degeneracies that a single starting Krylov
vector can miss. The potential is sampled once at `options.Time`: this solves
instantaneous `H(t)` and does not propagate a wavefunction.

```text
r[s] = ||H*v[s] - E[s]*v[s]||_2
success requires r[s] <= tolerance * max(1,abs(E[s])) for every requested state
```

Residuals are recomputed from the original operator. The result contains ascending
`Values`, column `Vectors`, `Residuals` and `Iterations`. Nonconvergence returns
the latest approximations AND an error. Callers must inspect the error. Small
residuals identify approximate eigenpairs; they do not prove every lower level
has been found. Increase the subspace/requested states, compare with small dense
solves, and check spatial convergence for difficult spectra.

`RealDiagonalize()` and `RealDiagonalizeAt(t)` return the full spectrum on a small
grid. `MatrixAt(t)` returns an independent dense matrix and an error. Unlike the
1D method, nD dense operations explicitly reject oversized allocations.

Eigenvectors have Euclidean norm one. For physical grid normalization use
`psiPhysical[i]=vector[i]/sqrt(VolumeElement())`. CSV densities apply this scaling.
`Energy` and `EnergyAt` divide by the state norm; uniform quadrature factors
cancel. A zero state has undefined (NaN) energy.

## Time dependence and propagation

Alongside the split propagator below, `NewNDimImplicitPropagator(h, options)` now
advances the whole complex wavefunction using implicit midpoint and matrix-free
GMRES. It handles both static and time-dependent Hamiltonians. See
[the vector propagation guide](16-vector-ode-and-hamiltonian-propagation.md)
for its API, accuracy, norm behavior and solver limits.

`ApplyHReal/ApplyHComplex` evaluate at zero; their `At` variants take time as the
last argument. They support overlapping input/output slices and can be used as
the RHS in coupled ODE solvers.

`SetStaticPotential(samples)` copies `N` finite samples and disables time dependence.
`SetTimeDependentPotential(func(t,out))` fills a full flattened potential;
`SetTimeDependentFunction(func(x,t) float64)` evaluates it by coordinates. Both
replace the complete static potential during application: include the static
term explicitly. Passing nil restores stored static samples.

```go
// Fragment: h is a two-coordinate Hamiltonian; math is imported.
h.SetTimeDependentFunction(func(x []float64, t float64) float64 {
    return 0.5*(x[0]*x[0]+x[1]*x[1]) + 0.3*x[0]*x[1] + 0.1*x[0]*math.Sin(t)
})
propagator, err := Quantum.NewNDimSplitPropagator(h)
if err != nil { panic(err) }
// psi is an initialized complex vector of length h.Dim().
if err := propagator.Step(psi, t, dt); err != nil { panic(err) }
```

The midpoint Strang step is

```text
psiNext = exp(-i*dt*Vmid/2) [product_d exp(-i*dt*T_d)] exp(-i*dt*Vmid/2) psi
Vmid = V(t+dt/2)
```

Cartesian kinetic terms commute across axes. Fourier uses normalized FFTs and
per-step phase arrays; DVR caches one-axis eigensystems. No full tensor kinetic
diagonalization is needed. Norm is not artificially rescaled. Invalid input or
nonfinite output returns an error without changing the input state. Choose `dt`
through convergence tests; unitarity alone does not prove accuracy.

An instance owns reusable work buffers and FFT plans. Do not use it concurrently
or retain callback buffers. The propagator cannot outlive its Hamiltonian.
`Close()` is idempotent and releases FFT plans.

## Scaling and limits

Let `N=product(n[d])` and `M` be the eigensolver subspace size.

| Operation | Work / storage |
| --- | --- |
| Grid metadata | O(D), no coordinate table |
| Potential/state buffers | O(N) |
| Fourier application/kinetic propagation | O(N sum log(n[d])) work; O(N+sum n[d]) storage |
| DVR application/kinetic propagation | O(N sum n[d]) work; O(N+sum n[d]²) storage |
| Iterative eigensolver | O(N M) basis/image storage plus projection/orthogonalization work |
| Full dense spectrum | O(N²) storage, O(N³) work |

Allocation guards:

- **1,048,576** total Hamiltonian points (`NDimMaxPoints`).
- **1,024** total points for full dense operations (`NDimDenseLimit`).
- **8,388,608** elements across DVR axis matrices.
- **1,024** points per axis for DVR propagation diagonalization.
- Iterative solver: at most **min(N,64)** states, **512** basis vectors, and
  `N*MaxSubspace <= 16,777,216`. Two stored bases alone can reach about 256 MiB;
  work and output arrays add memory.

These limits are not speed promises. Cancellation is checked between operator
applications/iterations, not inside one FFT, tensor application or factorization.
Dense solves do not take a context. This is a full tensor-grid method: state size
still grows exponentially with dimension. Sparse grids, tensor networks and
many-body symmetry reductions are not implemented.

## Verification

Tests cover indexing, unequal axis shapes, ownership, malformed files, overflow,
FFT ordering, dense/matrix-free agreement, aliasing and operator lifecycle.
Numerical references include plane waves in 1–4 dimensions; the 2D coupled
oscillator with frequencies `sqrt(1±0.3)`; 3D levels `1.5,2.5,2.5,2.5`; unequal
masses/frequencies; compatibility with the existing 1D path; translated boxes;
independent eigenpair residuals/orthogonality; cancellation and nonconvergence;
exact Fourier and time-dependent phase evolution; norm conservation, reversibility
and second-order propagation accuracy against dense diagonalization.

```sh
go test ./gridData ./Quantum ./OperatorAlgebra
go test ./...
go vet ./...
```

Sources: [grid tests](../gridData/gridHigherDim_test.go),
[Hamiltonian/solver/propagator tests](../Quantum/nd_hamiltonian_test.go).
