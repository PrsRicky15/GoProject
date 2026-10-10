package Quantum

import (
	"GoProject/OperatorAlgebra"
	"GoProject/gridData"
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"
)

// Limits bound allocations, not the number of dimensions. Dense work is optional.
const NDimMaxPoints = 1 << 20
const NDimDenseLimit = 1024
const ndMaxAxisMatrixElements = 8 << 20

// NDimHamiltonian applies H=sum_d[-1/(2*m_d) d²/dx_d²]+V(x,t), hbar=1,
// on the original last-coordinate-fastest NDimGrid. V may couple all coordinates.
// The axes must be Cartesian; this is not a curvilinear/radial Laplacian.
// An instance owns work buffers/FFT plans and must not be used concurrently.
type NDimHamiltonian struct {
	grid    *gridData.NDimGrid
	core    *HamiltonianOp
	kinetic *tensorKinetic
}

// NewNDimHamiltonian uses a Cartesian sinc DVR on each axis, independently of
// the location of the box. Wavefunctions should decay before the box boundaries.
func NewNDimHamiltonian(g *gridData.NDimGrid, masses []float64, pot gridData.NDimPotentialOp) (*NDimHamiltonian, error) {
	return newNDimHamiltonian(g, masses, pot, false)
}

// NewNDimFourierHamiltonian uses periodic Fourier kinetics on every axis.
// Call Close when finished, including after using an eigensolver/propagator.
func NewNDimFourierHamiltonian(g *gridData.NDimGrid, masses []float64, pot gridData.NDimPotentialOp) (*NDimHamiltonian, error) {
	return newNDimHamiltonian(g, masses, pot, true)
}

func newNDimHamiltonian(g *gridData.NDimGrid, masses []float64, pot gridData.NDimPotentialOp, fourier bool) (*NDimHamiltonian, error) {
	if g == nil || g.TotalPoints() == 0 || g.TotalPoints() > NDimMaxPoints {
		return nil, fmt.Errorf("NDimHamiltonian: grid must have 1..%d total points", NDimMaxPoints)
	}
	if len(masses) != int(g.NDim()) {
		return nil, fmt.Errorf("NDimHamiltonian: supply one mass per axis")
	}
	var axisElements uint64
	for d, m := range masses {
		if !finite(m) || m <= 0 {
			return nil, fmt.Errorf("NDimHamiltonian: invalid mass on axis %d", d)
		}
		dx := g.DeltaR(uint32(d))
		scale := 1 / (m * dx * dx)
		if !finite(scale) || scale <= 0 || !finite(math.Pi*math.Pi*scale) {
			return nil, fmt.Errorf("NDimHamiltonian: unrepresentable kinetic scale on axis %d", d)
		}
		n := uint64(g.NumPoints(uint32(d)))
		axisElements += n * n
	}
	if !fourier && axisElements > ndMaxAxisMatrixElements {
		return nil, fmt.Errorf("NDimHamiltonian: DVR axis matrices exceed storage limit; use Fourier or fewer points")
	}
	n := int(g.TotalPoints())
	core := &HamiltonianOp{vStatic: make([]float64, n), vbuf: make([]float64, n), workR: make([]float64, n), workC: make([]complex128, n), hpsi: make([]complex128, n)}
	g = g.Clone()
	if pot != nil {
		g.PotentialOnGridInPlace(pot, core.vStatic)
		for i, v := range core.vStatic {
			if !finite(v) {
				return nil, fmt.Errorf("NDimHamiltonian: non-finite potential at point %d", i)
			}
		}
	}
	k := &tensorKinetic{n: n}
	strides := g.Strides()
	for d, m := range masses {
		// Centering the kinetic-only grid avoids the legacy DVR radial-domain
		// heuristic. Cartesian derivatives are translation invariant; V uses g.
		a, err := gridData.NewFromLength(g.Length(uint32(d)), g.NumPoints(uint32(d)))
		if err != nil {
			k.close()
			return nil, err
		}
		na := int(a.NPoints())
		axis := tensorAxis{n: na, stride: int(strides[d]), realIn: make([]float64, na), realOut: make([]float64, na), complexIn: make([]complex128, na), complexOut: make([]complex128, na)}
		if fourier {
			axis.op = &fourierKinetic{basis: OperatorAlgebra.FFTInit(a, m), work: make([]complex128, na)}
			axis.energies = a.KValues()
			for j, waveNumber := range axis.energies {
				axis.energies[j] = waveNumber * waveNumber / (2 * m)
			}
		} else {
			axis.op = &dvrKinetic{basis: OperatorAlgebra.NewKeDVR(a, m)}
		}
		k.axes = append(k.axes, axis)
	}
	core.kinE = k
	return &NDimHamiltonian{grid: g, core: core, kinetic: k}, nil
}

func (h *NDimHamiltonian) Dim() int                             { return h.core.Dim() }
func (h *NDimHamiltonian) NDim() uint32                         { return h.grid.NDim() }
func (h *NDimHamiltonian) Grid() *gridData.NDimGrid             { return h.grid.Clone() }
func (h *NDimHamiltonian) Close()                               { h.core.Close() }
func (h *NDimHamiltonian) IsTimeDependent() bool                { return h.core.IsTimeDependent() }
func (h *NDimHamiltonian) VStatic() []float64                   { return h.core.VStatic() }
func (h *NDimHamiltonian) SetStaticPotential(v []float64) error { return h.core.SetStaticPotential(v) }
func (h *NDimHamiltonian) SetTimeDependentPotential(f func(t float64, out []float64)) {
	h.core.SetTimeDependentPotential(f)
}

// SetTimeDependentFunction evaluates the complete V(x,t), replacing static V.
// The callback must not retain or mutate x. Nil restores the static samples.
func (h *NDimHamiltonian) SetTimeDependentFunction(f func(x []float64, t float64) float64) {
	if f == nil {
		h.SetTimeDependentPotential(nil)
		return
	}
	x := make([]float64, h.NDim())
	h.SetTimeDependentPotential(func(t float64, out []float64) {
		for i := range out {
			h.grid.PointTo(uint64(i), x)
			out[i] = f(x, t)
		}
	})
}
func (h *NDimHamiltonian) PotentialAt(t float64, out []float64)      { h.core.PotentialAt(t, out) }
func (h *NDimHamiltonian) ApplyHReal(in, out []float64)              { h.core.ApplyHReal(in, out) }
func (h *NDimHamiltonian) ApplyHRealAt(in, out []float64, t float64) { h.core.ApplyHRealAt(in, out, t) }
func (h *NDimHamiltonian) ApplyHComplex(in, out []complex128)        { h.core.ApplyHComplex(in, out) }
func (h *NDimHamiltonian) ApplyHComplexAt(in, out []complex128, t float64) {
	h.core.ApplyHComplexAt(in, out, t)
}
func (h *NDimHamiltonian) Energy(psi []complex128) float64 { return h.core.Energy(psi) }
func (h *NDimHamiltonian) EnergyAt(psi []complex128, t float64) float64 {
	return h.core.EnergyAt(psi, t)
}

// MatrixAt is deliberately bounded: multidimensional application/LowestStates
// do not need this O(N²) allocation. Returned matrices are independent copies.
func (h *NDimHamiltonian) MatrixAt(t float64) (*mat.Dense, error) {
	if h.core.closed || !finite(t) {
		return nil, fmt.Errorf("NDimHamiltonian: closed operator or invalid time")
	}
	if h.Dim() > NDimDenseLimit {
		return nil, fmt.Errorf("dense Hamiltonian limited to %d points; use LowestStates or ApplyHReal/Complex", NDimDenseLimit)
	}
	v := h.core.currentV(t)
	for _, value := range v {
		if !finite(value) {
			return nil, fmt.Errorf("NDimHamiltonian: non-finite potential")
		}
	}
	m := h.kinetic.matrix()
	for i, value := range v {
		m.Set(i, i, m.At(i, i)+value)
	}
	return m, nil
}

func (h *NDimHamiltonian) RealDiagonalize() ([]float64, *mat.Dense, error) {
	return h.RealDiagonalizeAt(0)
}

// RealDiagonalizeAt returns the full instantaneous spectrum for a small grid.
func (h *NDimHamiltonian) RealDiagonalizeAt(t float64) ([]float64, *mat.Dense, error) {
	m, err := h.MatrixAt(t)
	if err != nil {
		return nil, nil, err
	}
	sym := mat.NewSymDense(h.Dim(), nil)
	for i := 0; i < h.Dim(); i++ {
		for j := i; j < h.Dim(); j++ {
			sym.SetSym(i, j, m.At(i, j))
		}
	}
	var eig mat.EigenSym
	if !eig.Factorize(sym, true) {
		return nil, nil, fmt.Errorf("NDimHamiltonian: diagonalization failed")
	}
	var vectors mat.Dense
	eig.VectorsTo(&vectors)
	return eig.Values(nil), &vectors, nil
}

type tensorAxis struct {
	n, stride             int
	op                    kineticEnergy
	energies              []float64 // Fourier kinetic eigenvalues, in FFT order
	realIn, realOut       []float64
	complexIn, complexOut []complex128
}
type tensorKinetic struct {
	n    int
	axes []tensorAxis
}

func (k *tensorKinetic) applyReal(in, out []float64) {
	clear(out)
	for d := range k.axes {
		a := &k.axes[d]
		for base := 0; base < k.n; base += a.n * a.stride {
			for offset := 0; offset < a.stride; offset++ {
				for j := range a.realIn {
					a.realIn[j] = in[base+offset+j*a.stride]
				}
				a.op.applyReal(a.realIn, a.realOut)
				for j, value := range a.realOut {
					out[base+offset+j*a.stride] += value
				}
			}
		}
	}
}
func (k *tensorKinetic) applyComplex(in, out []complex128) {
	clear(out)
	for d := range k.axes {
		a := &k.axes[d]
		for base := 0; base < k.n; base += a.n * a.stride {
			for offset := 0; offset < a.stride; offset++ {
				for j := range a.complexIn {
					a.complexIn[j] = in[base+offset+j*a.stride]
				}
				a.op.applyComplex(a.complexIn, a.complexOut)
				for j, value := range a.complexOut {
					out[base+offset+j*a.stride] += value
				}
			}
		}
	}
}
func (k *tensorKinetic) matrix() *mat.Dense {
	if k.n > NDimDenseLimit {
		panic("tensor kinetic: dense allocation exceeds limit")
	}
	m := mat.NewDense(k.n, k.n, nil)
	for _, a := range k.axes {
		axis := a.op.matrix()
		for base := 0; base < k.n; base += a.n * a.stride {
			for offset := 0; offset < a.stride; offset++ {
				for i := 0; i < a.n; i++ {
					for j := 0; j < a.n; j++ {
						row, col := base+offset+i*a.stride, base+offset+j*a.stride
						m.Set(row, col, m.At(row, col)+axis.At(i, j))
					}
				}
			}
		}
	}
	return m
}
func (k *tensorKinetic) close() {
	for _, a := range k.axes {
		a.op.close()
	}
}
