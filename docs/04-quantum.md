# Quantum operators, spectra and propagation

[Documentation index](README.md) · [Hamiltonian source](../Quantum/hamiltonian.go)

## Governing problem

The code represents the one-coordinate Hamiltonian, in hbar=1 units:

```text
H(t) = T + V(x,t)
T = -(1/(2m)) d²/dx²
i dpsi/dt = H(t) psi
```

The composition follows the operator strategy previously used in h2p: apply T,
then add pointwise multiplication by V. No runtime dependency on the h2p project
is needed. This is not a many-electron Hamiltonian or an orbital-basis chemistry solver.

## Basis implementations

### Sinc DVR

On the symmetric full-line grid, the truncated kinetic matrix is:

```text
T[i,i] = pi²/(6*m*dx²)
T[i,j] = (-1)^(i-j)/(m*dx²*(i-j)²),  i != j
```

`KeDvrBasis.GetMat()` lazily fills and caches this dense matrix, then returns the
cache itself. Treat the returned matrix as borrowed mutable data; changing it
changes subsequent low-level operator results. `RealDiagonalize()` copies it first.

The low-level code selects a different half-line formula when
`abs(RMax+RMin)>=1`. That heuristic does not reliably identify a physical half-line
problem. The web page uses `[-L,L)` specifically, selecting the full-line formula.
The asymmetric/half-line branch remains experimental.

`GetComplexMat(theta)` multiplies the real kinetic matrix by `exp(-2*i*theta)`.
This operation alone does not implement a complete complex-scaled resonance solver.
The legacy general complex eigensolver is not operational; see the status guide.

### Fourier

The Fourier basis uses periodic boundaries and FFT-ordered wavenumbers:

```text
T psi = inverseFFT( [k²/(2m)] * FFT(psi) )
p psi = inverseFFT( k * FFT(psi) )
```

FFTW transforms are unnormalized, so `operatorOp` includes the required `1/N`
factor. `LaplacianOp` is a legacy name: it applies **kinetic energy**, including
the mass and minus-one-half factors, not just `d²/dx²`.

`FFTInit` owns plans and a work array. `Clean()` releases plans; avoid reusing the
object afterward or calling raw cleanup repeatedly. Prefer Hamiltonian ownership
and its idempotent `Close()` for application code. Planning/destruction is locked
globally; mutable execution buffers remain specific to one caller.

## Hamiltonian API

| Method | Behavior |
| --- | --- |
| `NewHamiltonian(grid,mass,pot)` | Checked DVR constructor |
| `NewFourierHamiltonian(grid,mass,pot)` | Checked Fourier constructor |
| `NewHamil(...)` | Compatibility DVR wrapper; panics if construction fails |
| `Dim()` | State-vector length |
| `Grid()`, `VStatic()` | Return snapshots/copies rather than the internal mutable data |
| `SetStaticPotential(samples)` | Copies validated samples; disables time-dependent callback |
| `SetTimeDependentPotential(callback)` | Installs full V(t); nil restores stored static potential |
| `PotentialAt(t,out)` | Fills potential samples at t |
| `ApplyHReal`, `ApplyHComplex` | Apply at t=0 |
| `ApplyHRealAt`, `ApplyHComplexAt` | Apply at supplied t; support overlapping input/output |
| `MatrixAt(t)` | Independent dense H(t) |
| `EvaluateOp()` | Independent matrix at zero time |
| `Mat()` | Legacy method storing an internal dense matrix at zero time |
| `RealDiagonalize()` | Ascending eigenvalues and eigenvectors in columns, at t=0 |
| `Energy`, `EnergyAt` | Normalized real energy expectation; zero norm yields NaN |
| `Close()` | Releases resources; subsequent operator use is invalid |

Constructors snapshot grid and static potential, reject invalid mass/grid/potential
values, and permit nil potential for a free particle. Vector-size misuse of the
application methods is a programming error and can panic. Instances are not safe
for simultaneous calls because they reuse work buffers.

The callback must fill the **entire** potential, including its static term:

```go
base := h.VStatic()
x := h.Grid().RValues()
h.SetTimeDependentPotential(func(t float64, out []float64) {
    for i := range out {
        out[i] = base[i] + amplitude*x[i]*math.Sin(frequency*t)
    }
})
```

Do not retain `out`. A malformed custom callback can violate library assumptions;
the web page constructs its callback from bounded numeric parameters.

## Energy levels and normalization

The spectrum job forms H and uses a real symmetric eigensolver. Columns satisfy
`H*u[n]=E[n]*u[n]` and Euclidean normalization `sum(u[n]²)=1`. The plotted physical
density is `u[n]²/dx`, so `sum(density)*dx=1`.

The displayed residual is the largest `||H*u-E*u||₂` among the requested states.
A small residual verifies the matrix eigenproblem, not convergence to the
infinite-domain continuum problem. Increase N and L separately to check that.

For `m=k=1` harmonic motion, the reference spectrum is `E[n]=n+1/2`. In general,
`omega=sqrt(k/m)` and `E[n]=(n+1/2)*omega`. High finite-box Morse eigenstates may
represent discretized continuum states; compare energies with `De`.

## Initial wave packet and drive

The propagation job initializes and normalizes:

```text
psi(x,0) proportional to exp(-(x-x0)²/(2*sigma²)) * exp(i*p0*x)
sum |psi|² * dx = 1
V(x,t) = V0(x) + A*x*sin(omega*t)
```

`sigma` is the amplitude-envelope parameter. The untruncated density has standard
deviation `sigma/sqrt(2)`, not sigma. `p0` is momentum because hbar=1.
The validator requires `abs(x0)+3*sigma<L`, but this does not ensure that the packet
is adequately resolved by dx. Check convergence for narrow states and large p0.

### Midpoint split operator

`NewSplitPropagator(h)` diagonalizes the kinetic matrix once. A step applies:

```text
psi(t+dt) ≈ exp(-i*dt*V(t+dt/2)/2)
            exp(-i*dt*T)
            exp(-i*dt*V(t+dt/2)/2) psi(t)
```

Kinetic application uses cached real eigenvectors and phases `exp(-i*dt*Ekin)`.
The method is second order for smooth time-dependent potentials. It is unitary
up to numerical rounding for real potentials, but a small norm drift does not
prove small phase or trajectory error. It uses dense transforms for both bases.

### Explicit RK propagation

For `psi=a+i*b` and a real H, the adapter stores `[a,b]` and evaluates
`a'=H*b`, `b'=-H*a`. It delegates full coupled-vector steps to RK4, Ralston 3 or
Nyström 5. These methods are not exactly unitary.

The web layer estimates a spectral bound from the largest absolute row sum of
H(0), plus `abs(A)*L`. It requires `dt*bound<=2` for RK4 and `<=1` for the other
explicit choices. This is a conservative admission check, not an accuracy guarantee.
If norm leaves [0.5,1.5], the job fails and asks for a smaller dt. No renormalization
is applied to conceal integration error.

## Diagnostics and sampling

At each integration step Go checks context and norm. It records at stride
`ceil(steps/150)`, including initial and final states. At most 151 time samples
are generated by this scheme. Charts show density, normalized mean position,
norm, and energy expectation. Metrics include initial/final energy, maximum
absolute norm drift, and final mean position.

`EnergyAt` divides by the discrete norm; the common dx cancels in the quotient.
A driven Hamiltonian generally exchanges energy with the applied field. Energy
need not remain constant when A is nonzero. There is no absorber: Fourier states
can wrap periodically and finite DVR domains can contaminate boundary behavior.

## Low-level exponential sign convention

The completed DVR methods intentionally preserve their older API:

| Method | Mathematical operation |
| --- | --- |
| `ExpDt(dt,in)` / `ExpDtInPlace` / `ExpDtTo` | `exp(dt*T)` |
| `ExpIdt(dt,in)` / `ExpIdtInPlace` / `ExpIdtTo` | `exp(+i*dt*T)` |
| `SplitPropagator.Step(psi,t,dt)` | Schrödinger evolution using negative imaginary phases |

Therefore use a negative dt with the old `ExpIdt` if explicitly applying
`exp(-i*dt*T)`. Allocating and in-place DVR variants return errors; `...To`
wrappers retain their no-error signature and panic on invalid dimensions/errors.
Those convenience routines do not cache an eigensystem across successive calls;
use the split propagator for repeated web evolution.

## Minimal library example

```go
grid, err := gridData.NewFromLength(16, 64)
if err != nil { return err }
h, err := Quantum.NewHamiltonian(grid, 1,
    gridData.Harmonic[float64]{ForceConst: 1})
if err != nil { return err }
defer h.Close()
energies, vectors, err := h.RealDiagonalize()
if err != nil { return err }
// energies[0] is approximately 0.5; vectors stores eigenvectors in columns.
_ = energies
_ = vectors
```

This fragment belongs in a Go function returning error, with the two project
packages imported. A complete executable example is in the
[package README](../Quantum/README.md).
