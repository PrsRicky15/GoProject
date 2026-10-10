package Quantum

import (
	"fmt"
	"math/cmplx"

	"gonum.org/v1/gonum/mat"
)

// NDimSplitPropagator uses midpoint Strang splitting. The commuting Cartesian
// kinetic terms are propagated axis by axis, with FFTs for Fourier or cached
// one-axis eigensystems for DVR. No full tensor matrix is formed.
type NDimSplitPropagator struct {
	h         *NDimHamiltonian
	axes      []ndAxisPropagator
	potential []float64
	work      []complex128
}
type ndAxisPropagator struct {
	values              []float64
	vectors             *mat.Dense // nil for a Fourier axis
	phase, coefficients []complex128
}

func NewNDimSplitPropagator(h *NDimHamiltonian) (*NDimSplitPropagator, error) {
	if h == nil || h.core.closed {
		return nil, fmt.Errorf("NDimSplitPropagator: nil or closed Hamiltonian")
	}
	p := &NDimSplitPropagator{h: h, potential: make([]float64, h.Dim()), work: make([]complex128, h.Dim())}
	for _, a := range h.kinetic.axes {
		axis := ndAxisPropagator{phase: make([]complex128, a.n), coefficients: make([]complex128, a.n)}
		if a.energies != nil {
			axis.values = append([]float64(nil), a.energies...)
		} else {
			if a.n > NDimDenseLimit {
				return nil, fmt.Errorf("DVR propagation limited to %d points per axis; use Fourier", NDimDenseLimit)
			}
			m := a.op.matrix()
			sym := mat.NewSymDense(a.n, nil)
			for i := 0; i < a.n; i++ {
				for j := i; j < a.n; j++ {
					sym.SetSym(i, j, m.At(i, j))
				}
			}
			var eig mat.EigenSym
			if !eig.Factorize(sym, true) {
				return nil, fmt.Errorf("NDimSplitPropagator: axis diagonalization failed")
			}
			axis.values, axis.vectors = eig.Values(nil), new(mat.Dense)
			eig.VectorsTo(axis.vectors)
		}
		p.axes = append(p.axes, axis)
	}
	return p, nil
}

// Step advances psi by exp(-i H dt). Errors leave psi unchanged. Norm is not
// rescaled. The propagator cannot outlive or run concurrently with its Hamiltonian.
func (p *NDimSplitPropagator) Step(psi []complex128, t, dt float64) error {
	if p.h.core.closed || len(psi) != p.h.Dim() || !finite(t) || !finite(dt) || !finite(t+dt) {
		return fmt.Errorf("NDimSplitPropagator: invalid state, time or closed Hamiltonian")
	}
	p.h.PotentialAt(t+dt/2, p.potential)
	for i, v := range p.potential {
		if !finite(v) || !finite(dt*v) || !finite(real(psi[i])) || !finite(imag(psi[i])) {
			return fmt.Errorf("NDimSplitPropagator: non-finite state or potential phase")
		}
		p.work[i] = psi[i] * cmplx.Exp(complex(0, -dt*v/2))
	}
	for d := range p.axes {
		a, prop := &p.h.kinetic.axes[d], &p.axes[d]
		for i, e := range prop.values {
			if !finite(dt * e) {
				return fmt.Errorf("NDimSplitPropagator: non-finite kinetic phase")
			}
			prop.phase[i] = cmplx.Exp(complex(0, -dt*e))
		}
		for base := 0; base < len(psi); base += a.n * a.stride {
			for offset := 0; offset < a.stride; offset++ {
				for j := range a.complexIn {
					a.complexIn[j] = p.work[base+offset+j*a.stride]
				}
				if f, ok := a.op.(*fourierKinetic); ok {
					f.basis.ApplyComplexSpectrum(a.complexIn, prop.phase)
				} else {
					for j := range prop.coefficients {
						var c complex128
						for i, v := range a.complexIn {
							c += complex(prop.vectors.At(i, j), 0) * v
						}
						prop.coefficients[j] = c * prop.phase[j]
					}
					for i := range a.complexIn {
						var v complex128
						for j, c := range prop.coefficients {
							v += complex(prop.vectors.At(i, j), 0) * c
						}
						a.complexIn[i] = v
					}
				}
				for j, v := range a.complexIn {
					p.work[base+offset+j*a.stride] = v
				}
			}
		}
	}
	for i, v := range p.potential {
		p.work[i] *= cmplx.Exp(complex(0, -dt*v/2))
		if !finite(real(p.work[i])) || !finite(imag(p.work[i])) {
			return fmt.Errorf("NDimSplitPropagator: non-finite propagated state")
		}
	}
	copy(psi, p.work)
	return nil
}
