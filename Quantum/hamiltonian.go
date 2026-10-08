package Quantum

import (
	"GoProject/OperatorAlgebra"
	"GoProject/gridData"
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"
)

// HamiltonianOp applies H = T + V in atomic units (hbar = 1). Like h2p's
// Hamiltonian, it keeps the kinetic operator and diagonal potential separate.
// Work buffers are reused, so an instance must not be used concurrently.
type HamiltonianOp struct {
	grid    *gridData.RadGrid
	kinE    kineticEnergy
	vStatic []float64
	tdPot   func(t float64, out []float64)
	vbuf    []float64
	workR   []float64
	workC   []complex128
	hpsi    []complex128
	hmat    *mat.Dense
	closed  bool
}

// NewHamil preserves the original DVR constructor. Invalid inputs panic; use
// NewHamiltonian when the caller needs to handle a construction error.
func NewHamil(grid *gridData.RadGrid, mass float64, pot gridData.PotentialOp[float64]) *HamiltonianOp {
	op, err := NewHamiltonian(grid, mass, pot)
	if err != nil {
		panic(err)
	}
	return op
}

// NewHamiltonian constructs a DVR Hamiltonian. A nil potential means V = 0.
// The grid and sampled potential are snapshots; later caller mutations do not
// change the operator. Its kinetic boundary convention is that of KeDvrBasis.
func NewHamiltonian(grid *gridData.RadGrid, mass float64, pot gridData.PotentialOp[float64]) (*HamiltonianOp, error) {
	op, err := newHamiltonian(grid, mass, pot)
	if err != nil {
		return nil, err
	}
	op.kinE = &dvrKinetic{basis: OperatorAlgebra.NewKeDVR(op.grid, mass)}
	return op, nil
}

// NewFourierHamiltonian uses h2p's FFT kinetic-energy strategy with periodic
// boundaries. Applying H requires O(N log N) work and O(N) storage. Call Close
// when finished to release the FFT plans.
func NewFourierHamiltonian(grid *gridData.RadGrid, mass float64, pot gridData.PotentialOp[float64]) (*HamiltonianOp, error) {
	op, err := newHamiltonian(grid, mass, pot)
	if err != nil {
		return nil, err
	}
	op.kinE = &fourierKinetic{
		basis: OperatorAlgebra.FFTInit(op.grid, mass),
		work:  make([]complex128, op.Dim()),
	}
	return op, nil
}

func newHamiltonian(grid *gridData.RadGrid, mass float64, pot gridData.PotentialOp[float64]) (*HamiltonianOp, error) {
	if grid == nil || grid.NPoints() == 0 {
		return nil, fmt.Errorf("Hamiltonian: grid must contain at least one point")
	}
	if !finite(grid.RMin()) || !finite(grid.RMax()) || !finite(grid.Length()) ||
		!finite(grid.DeltaR()) || grid.DeltaR() <= 0 {
		return nil, fmt.Errorf("Hamiltonian: grid bounds and spacing must be finite, with positive spacing")
	}
	if !finite(mass) || mass <= 0 {
		return nil, fmt.Errorf("Hamiltonian: mass must be finite and positive")
	}
	gridCopy := *grid
	n := int(grid.NPoints())
	op := &HamiltonianOp{
		grid:    &gridCopy,
		vStatic: make([]float64, n),
		vbuf:    make([]float64, n),
		workR:   make([]float64, n),
		workC:   make([]complex128, n),
		hpsi:    make([]complex128, n),
	}
	if pot != nil {
		if err := op.SetStaticPotential(pot.EvaluateOnGrid(op.grid.RValues())); err != nil {
			return nil, err
		}
	}
	return op, nil
}

func finite(v float64) bool { return !math.IsNaN(v) && !math.IsInf(v, 0) }

func (op *HamiltonianOp) Dim() int { return len(op.vStatic) }

// Grid returns a copy of the grid used to construct this operator.
func (op *HamiltonianOp) Grid() *gridData.RadGrid {
	gridCopy := *op.grid
	return &gridCopy
}

func (op *HamiltonianOp) IsTimeDependent() bool { return op.tdPot != nil }

// VStatic returns a copy so callers cannot silently change the potential.
func (op *HamiltonianOp) VStatic() []float64 {
	return append([]float64(nil), op.vStatic...)
}

// SetStaticPotential copies the samples and disables a time-dependent potential.
func (op *HamiltonianOp) SetStaticPotential(v []float64) error {
	if op.closed {
		return fmt.Errorf("Hamiltonian: operator is closed")
	}
	if len(v) != op.Dim() {
		return fmt.Errorf("SetStaticPotential: len(v)=%d, want %d", len(v), op.Dim())
	}
	for i, vi := range v {
		if !finite(vi) {
			return fmt.Errorf("SetStaticPotential: non-finite value at index %d", i)
		}
	}
	copy(op.vStatic, v)
	op.tdPot = nil
	op.hmat = nil
	return nil
}

// SetTimeDependentPotential follows h2p's convention: eval supplies the complete
// potential, replacing the static samples during application. It must fill out
// without retaining that buffer. Passing nil restores the static potential.
func (op *HamiltonianOp) SetTimeDependentPotential(eval func(t float64, out []float64)) {
	op.checkReady()
	op.tdPot = eval
	op.hmat = nil
}

func (op *HamiltonianOp) currentV(t float64) []float64 {
	if op.tdPot == nil {
		return op.vStatic
	}
	clear(op.vbuf)
	op.tdPot(t, op.vbuf)
	return op.vbuf
}

func (op *HamiltonianOp) checkReady() {
	if op.closed {
		panic("Hamiltonian: operator is closed")
	}
}

func (op *HamiltonianOp) checkVectors(in, out int) {
	op.checkReady()
	if in != op.Dim() || out != op.Dim() {
		panic(fmt.Sprintf("Hamiltonian: vector lengths %d and %d, want %d", in, out, op.Dim()))
	}
}

func (op *HamiltonianOp) PotentialAt(t float64, out []float64) {
	op.checkVectors(len(out), len(out))
	copy(out, op.currentV(t))
}

func (op *HamiltonianOp) ApplyHReal(psi, out []float64) {
	op.ApplyHRealAt(psi, out, 0)
}

// ApplyHRealAt overwrites out with H(t) psi and supports overlapping slices.
func (op *HamiltonianOp) ApplyHRealAt(psi, out []float64, t float64) {
	op.checkVectors(len(psi), len(out))
	copy(op.workR, psi)
	op.kinE.applyReal(op.workR, out)
	for i, vi := range op.currentV(t) {
		out[i] += vi * op.workR[i]
	}
}

func (op *HamiltonianOp) ApplyHComplex(psi, out []complex128) {
	op.ApplyHComplexAt(psi, out, 0)
}

// ApplyHComplexAt overwrites out with H(t) psi and supports overlapping slices.
func (op *HamiltonianOp) ApplyHComplexAt(psi, out []complex128, t float64) {
	op.checkVectors(len(psi), len(out))
	copy(op.workC, psi)
	op.kinE.applyComplex(op.workC, out)
	for i, vi := range op.currentV(t) {
		out[i] += complex(vi, 0) * op.workC[i]
	}
}

// MatrixAt explicitly assembles H(t) for small-grid diagnostics or eigensolvers.
// It returns an independent matrix; changing it cannot change the operator.
func (op *HamiltonianOp) MatrixAt(t float64) *mat.Dense {
	op.checkReady()
	h := op.kinE.matrix()
	for i, vi := range op.currentV(t) {
		h.Set(i, i, h.At(i, i)+vi)
	}
	return h
}

// Mat rebuilds the dense matrix at t = 0, preserving the original API.
func (op *HamiltonianOp) Mat() { op.hmat = op.MatrixAt(0) }

func (op *HamiltonianOp) EvaluateOp() *mat.Dense {
	op.Mat()
	return op.hmat
}

// RealDiagonalize returns ascending eigenvalues and column eigenvectors at t=0.
// The stored kinetic operator and potential are not modified.
func (op *HamiltonianOp) RealDiagonalize() ([]float64, *mat.Dense, error) {
	h := op.MatrixAt(0)
	n := op.Dim()
	sym := mat.NewSymDense(n, nil)
	for i := 0; i < n; i++ {
		for j := i; j < n; j++ {
			sym.SetSym(i, j, h.At(i, j))
		}
	}
	var eig mat.EigenSym
	if !eig.Factorize(sym, true) {
		return nil, nil, fmt.Errorf("Hamiltonian: symmetric diagonalization failed")
	}
	var vectors mat.Dense
	eig.VectorsTo(&vectors)
	return eig.Values(nil), &vectors, nil
}

func (op *HamiltonianOp) Energy(psi []complex128) float64 { return op.EnergyAt(psi, 0) }

// EnergyAt returns Re(<psi|H(t)|psi>)/<psi|psi>; a zero state has NaN energy.
func (op *HamiltonianOp) EnergyAt(psi []complex128, t float64) float64 {
	op.ApplyHComplexAt(psi, op.hpsi, t)
	var numerator, norm float64
	for i, z := range psi {
		numerator += real(z)*real(op.hpsi[i]) + imag(z)*imag(op.hpsi[i])
		norm += real(z)*real(z) + imag(z)*imag(z)
	}
	return numerator / norm
}

// Close releases kinetic resources. Repeated calls are safe; application after
// closing is an error. DVR operators have no external resources to release.
func (op *HamiltonianOp) Close() {
	if !op.closed {
		op.kinE.close()
		op.closed = true
	}
}
