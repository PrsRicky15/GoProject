package Quantum

import (
	EquationSolver "GoProject/ODESolver"
	"context"
	"fmt"
)

// ImplicitPropagator advances the WHOLE complex wavefunction via implicit
// midpoint: (I+i*dt*H(t+dt/2)/2) psiNext = (I-i*dt*H(t+dt/2)/2) psi.
// This is second order and preserves norm for Hermitian H up to solve error.
// For static H it equals Crank-Nicolson; for varying H both sides use midpoint H.
// It applies H to vectors through GMRES; no dense Hamiltonian is constructed.
// The propagator and its Hamiltonian must be used sequentially and not copied.
type ImplicitPropagator struct {
	solver *EquationSolver.ComplexLinearImplicit
}

func newImplicitPropagator(n int, apply func([]complex128, []complex128, float64), closed func() bool, options EquationSolver.KrylovOptions) (*ImplicitPropagator, error) {
	rhs := func(t float64, psi, out []complex128) error {
		if closed() {
			return fmt.Errorf("Hamiltonian is closed")
		}
		apply(psi, out, t)
		for i, v := range out {
			out[i] = complex(0, -1) * v
		}
		return nil
	}
	solver, err := EquationSolver.NewComplexLinearImplicit(n, rhs, EquationSolver.ImplicitMidpointMethod, options)
	if err != nil {
		return nil, err
	}
	return &ImplicitPropagator{solver: solver}, nil
}

func NewImplicitPropagator(h *HamiltonianOp, options EquationSolver.KrylovOptions) (*ImplicitPropagator, error) {
	if h == nil || h.closed {
		return nil, fmt.Errorf("nil or closed Hamiltonian")
	}
	return newImplicitPropagator(h.Dim(), h.ApplyHComplexAt, func() bool { return h.closed }, options)
}
func NewNDimImplicitPropagator(h *NDimHamiltonian, options EquationSolver.KrylovOptions) (*ImplicitPropagator, error) {
	if h == nil || h.core == nil || h.core.closed {
		return nil, fmt.Errorf("nil or closed n-dimensional Hamiltonian")
	}
	return newImplicitPropagator(h.Dim(), h.ApplyHComplexAt, func() bool { return h.core.closed }, options)
}
func (p *ImplicitPropagator) Step(psi []complex128, t, dt float64) error {
	return p.solver.Step(psi, t, dt)
}
func (p *ImplicitPropagator) StepContext(ctx context.Context, psi []complex128, t, dt float64) error {
	return p.solver.StepContext(ctx, psi, t, dt)
}
func (p *ImplicitPropagator) Stats() EquationSolver.KrylovStats { return p.solver.Stats() }
