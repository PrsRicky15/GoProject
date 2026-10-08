# One-dimensional Hamiltonian

For full-project context, see the [documentation index](../docs/README.md),
[detailed quantum notes](../docs/04-quantum.md), and
[known limitations](../docs/11-status-and-limitations.md).

This package follows the operator composition used in h2p's
`src/internal/physics/Hamiltonian.go`:

```text
H psi = T psi + V psi
T = -1/(2m) d²/dx²        (atomic units, hbar = 1)
```

The potential is stored as grid samples and multiplied point by point. Real and
complex application reuse work buffers and support in-place or overlapping input
and output slices. A single Hamiltonian instance is not safe for concurrent use.

## Choose a kinetic basis

| Constructor | Kinetic operator | Boundary convention |
| --- | --- | --- |
| `NewFourierHamiltonian(grid, mass, potential)` | FFT, as in h2p; O(N log N) application and O(N) storage | Periodic |
| `NewHamiltonian(grid, mass, potential)` | Existing DVR matrix; O(N²) application and storage | Existing `KeDvrBasis` domain convention |
| `NewHamil(grid, mass, potential)` | Compatibility wrapper for the DVR constructor | Same as `NewHamiltonian`; panics on invalid input |

Both checked constructors return `(*HamiltonianOp, error)`. A nil potential means
a free particle. They snapshot the grid and the sampled potential. Reconstruct
the operator when changing the grid or mass. Call `Close()` when finished,
particularly for Fourier operators, which own FFT plans.

For a localized state with the Fourier basis, choose a domain large enough that
the wavefunction is negligible at its periodic boundaries and check convergence
as the domain and grid resolution increase.

## Example

```go
package main

import (
    "fmt"
    "math"

    "GoProject/Quantum"
    "GoProject/gridData"
)

func main() {
    grid, err := gridData.NewFromLength(20, 128)
    if err != nil {
        panic(err)
    }
    h, err := Quantum.NewFourierHamiltonian(grid, 1,
        gridData.Harmonic[float64]{ForceConst: 1})
    if err != nil {
        panic(err)
    }
    defer h.Close()

    psi := make([]complex128, h.Dim())
    for i, x := range grid.RValues() {
        psi[i] = complex(math.Exp(-x*x/2), 0)
    }
    hpsi := make([]complex128, h.Dim())
    h.ApplyHComplex(psi, hpsi)
    fmt.Printf("Ground-state energy: %.6f\n", h.Energy(psi)) // 0.500000
}
```

`ApplyHReal` and `ApplyHComplex` evaluate at time zero. Their `At` variants accept
a time as the last argument. `EnergyAt(psi, t)` computes the normalized energy
expectation; a zero state has undefined (`NaN`) energy.

`SetStaticPotential(samples)` copies its input and disables any time-dependent
callback. `SetTimeDependentPotential(func(t float64, out []float64))` follows
h2p's convention: the callback fills the **complete** potential at time `t`,
including any static term. It must not retain `out`. Passing nil restores the
stored static potential.

`EvaluateOp()` builds an independent dense matrix at time zero; `MatrixAt(t)`
does so at a specified time. `RealDiagonalize()` returns ascending eigenvalues
and eigenvectors stored in columns at time zero. These dense operations are for
small grids; they are unnecessary when only applying H to a state.

## Verification

Use Go 1.25.1 or newer, with a C compiler and FFTW 3 development files available
for the existing FFTW dependency. From the repository root:

```sh
go test ./...
go vet ./...
```

Tests cover harmonic-oscillator energies and eigenvector residuals, free Fourier
modes with a constant potential, odd/even frequency ordering, matrix/application
agreement, overlapping slices, potential updates, and input validation.
