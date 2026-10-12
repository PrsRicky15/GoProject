package Quantum

import (
	EquationSolver "GoProject/ODESolver"
	"context"
	"fmt"
)

// RKPropagator evolves the whole wavefunction using dpsi/dt=-i H(t) psi
// (hbar=1). Explicit RK is not exactly unitary; no renormalization is applied.
// Hamiltonian and solver must be used sequentially, including shared adapters.
type RKPropagator struct{ solver *EquationSolver.ComplexRK }
type AdaptiveRKPropagator struct {
	solver *EquationSolver.ComplexAdaptiveRK
}

func hamiltonianRKRightHandSide(n int, apply func([]complex128, []complex128, float64), closed func() bool) EquationSolver.ComplexRHS {
	return func(t float64, psi, out []complex128) error {
		if closed() {
			return fmt.Errorf("Hamiltonian is closed")
		}
		if len(psi) != n || len(out) != n {
			return fmt.Errorf("wavefunction dimension must be %d", n)
		}
		apply(psi, out, t)
		for i, v := range out {
			out[i] = complex(0, -1) * v
		}
		return nil
	}
}
func NewRKPropagator(h *HamiltonianOp, method EquationSolver.RKMethod) (*RKPropagator, error) {
	if h == nil || h.closed {
		return nil, fmt.Errorf("nil or closed Hamiltonian")
	}
	s, err := EquationSolver.NewComplexRK(hamiltonianRKRightHandSide(h.Dim(), h.ApplyHComplexAt, func() bool { return h.closed }), method)
	if err != nil {
		return nil, err
	}
	return &RKPropagator{s}, nil
}
func NewNDimRKPropagator(h *NDimHamiltonian, method EquationSolver.RKMethod) (*RKPropagator, error) {
	if h == nil || h.core == nil || h.core.closed {
		return nil, fmt.Errorf("nil or closed n-dimensional Hamiltonian")
	}
	s, err := EquationSolver.NewComplexRK(hamiltonianRKRightHandSide(h.Dim(), h.ApplyHComplexAt, func() bool { return h.core.closed }), method)
	if err != nil {
		return nil, err
	}
	return &RKPropagator{s}, nil
}
func NewAdaptiveRKPropagator(h *HamiltonianOp, method EquationSolver.RKMethod, options EquationSolver.AdaptiveOptions) (*AdaptiveRKPropagator, error) {
	if h == nil || h.closed {
		return nil, fmt.Errorf("nil or closed Hamiltonian")
	}
	s, err := EquationSolver.NewComplexAdaptiveRK(hamiltonianRKRightHandSide(h.Dim(), h.ApplyHComplexAt, func() bool { return h.closed }), method, options)
	if err != nil {
		return nil, err
	}
	return &AdaptiveRKPropagator{s}, nil
}
func NewNDimAdaptiveRKPropagator(h *NDimHamiltonian, method EquationSolver.RKMethod, options EquationSolver.AdaptiveOptions) (*AdaptiveRKPropagator, error) {
	if h == nil || h.core == nil || h.core.closed {
		return nil, fmt.Errorf("nil or closed n-dimensional Hamiltonian")
	}
	s, err := EquationSolver.NewComplexAdaptiveRK(hamiltonianRKRightHandSide(h.Dim(), h.ApplyHComplexAt, func() bool { return h.core.closed }), method, options)
	if err != nil {
		return nil, err
	}
	return &AdaptiveRKPropagator{s}, nil
}
func (p *RKPropagator) Step(psi []complex128, t, dt float64) error { return p.solver.Step(psi, t, dt) }
func (p *RKPropagator) StepContext(ctx context.Context, psi []complex128, t, dt float64) error {
	return p.solver.StepContext(ctx, psi, t, dt)
}
func (p *AdaptiveRKPropagator) Integrate(psi []complex128, tStart, tEnd float64) (EquationSolver.AdaptiveResult, error) {
	return p.solver.Integrate(psi, tStart, tEnd)
}
func (p *AdaptiveRKPropagator) IntegrateContext(ctx context.Context, psi []complex128, tStart, tEnd float64) (EquationSolver.AdaptiveResult, error) {
	return p.solver.IntegrateContext(ctx, psi, tStart, tEnd)
}
