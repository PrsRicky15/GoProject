package Quantum

import (
	"fmt"
	"gonum.org/v1/gonum/mat"
	"math/cmplx"
)

// SplitPropagator uses midpoint Strang splitting for H(t)=T+V(t).
// The small-grid kinetic eigensystem is cached for both DVR and Fourier bases.
// It shares its Hamiltonian's lifecycle and must not be used concurrently.
type SplitPropagator struct {
	h                  *HamiltonianOp
	values             []float64
	vectors            *mat.Dense
	potential          []float64
	work, coefficients []complex128
}

func NewSplitPropagator(h *HamiltonianOp) (*SplitPropagator, error) {
	if h == nil {
		return nil, fmt.Errorf("nil Hamiltonian")
	}
	h.checkReady()
	k := h.kinE.matrix()
	n := h.Dim()
	sym := mat.NewSymDense(n, nil)
	for i := 0; i < n; i++ {
		for j := i; j < n; j++ {
			sym.SetSym(i, j, k.At(i, j))
		}
	}
	var eig mat.EigenSym
	if !eig.Factorize(sym, true) {
		return nil, fmt.Errorf("kinetic diagonalization failed")
	}
	var vectors mat.Dense
	eig.VectorsTo(&vectors)
	return &SplitPropagator{h: h, values: eig.Values(nil), vectors: &vectors, potential: make([]float64, n), work: make([]complex128, n), coefficients: make([]complex128, n)}, nil
}

// Step advances exp(-i H dt), with hbar=1. Norm is not artificially rescaled.
func (p *SplitPropagator) Step(psi []complex128, t, dt float64) error {
	if len(psi) != p.h.Dim() || !finite(t) || !finite(dt) || !finite(t+dt) {
		return fmt.Errorf("invalid propagation input")
	}
	p.h.PotentialAt(t+dt/2, p.potential)
	for i, v := range p.potential {
		if !finite(v) || !finite(real(psi[i])) || !finite(imag(psi[i])) {
			return fmt.Errorf("non-finite state or potential")
		}
		p.work[i] = psi[i] * cmplx.Exp(complex(0, -dt*v/2))
	}
	for j, e := range p.values {
		var c complex128
		for i, w := range p.work {
			c += complex(p.vectors.At(i, j), 0) * w
		}
		p.coefficients[j] = c * cmplx.Exp(complex(0, -dt*e))
	}
	for i := range psi {
		var v complex128
		for j, c := range p.coefficients {
			v += complex(p.vectors.At(i, j), 0) * c
		}
		psi[i] = v * cmplx.Exp(complex(0, -dt*p.potential[i]/2))
	}
	return nil
}
