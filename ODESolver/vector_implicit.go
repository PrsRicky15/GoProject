package EquationSolver

import (
	"fmt"
	"math"
)

// VectorRHS evaluates dy/dt=F(t,y) using the ENTIRE vector at once.
// Treat y as read-only, fill all of dy, and do not retain either slice.
type VectorRHS func(t float64, y, dy []float64) error

// ComplexRHS is the corresponding whole-vector complex equation. It need not
// be complex analytic: conjugates and real/imaginary coupling are supported.
type ComplexRHS func(t float64, y, dy []complex128) error

type ImplicitMethod uint8

const (
	ImplicitEulerMethod ImplicitMethod = iota
	ImplicitTrapezoidalMethod
	ImplicitMidpointMethod
	ImplicitSDIRK2Method
)

func (m ImplicitMethod) valid() bool { return m <= ImplicitSDIRK2Method }

// VectorImplicit exposes the existing coupled nonlinear engine without requiring
// any scalar callbacks. Step takes dt explicitly and commits only on success.
// Instances own workspace and are not safe for concurrent use or copying.
type VectorImplicit struct {
	bridge  vectorRHSBridge
	method  ImplicitMethod
	options ImplicitOptions
	work    implicitWorkspace
}

func NewVectorImplicit(rhs VectorRHS, method ImplicitMethod, options ImplicitOptions) (*VectorImplicit, error) {
	if rhs == nil || !method.valid() {
		return nil, fmt.Errorf("invalid vector RHS or implicit method")
	}
	if _, err := options.normalized(); err != nil {
		return nil, err
	}
	return &VectorImplicit{bridge: vectorRHSBridge{rhs: rhs}, method: method, options: options}, nil
}

func (s *VectorImplicit) Step(y []float64, t, dt float64) error {
	s.bridge.err = nil
	err := s.work.step(&s.bridge, y, t, dt, s.options, implicitMethod(s.method), true, false)
	if s.bridge.err != nil {
		return fmt.Errorf("vector RHS: %w", s.bridge.err)
	}
	return err
}
func (s *VectorImplicit) Stats() ImplicitStats { return s.work.stats }

// This private adapter preserves the original TDPotentialOp-based API intact.
// Errors poison the residual so the legacy engine cannot commit a failed stage.
type vectorRHSBridge struct {
	rhs VectorRHS
	err error
}

func (b *vectorRHSBridge) EvaluateAt(float64, float64) float64 { return math.NaN() }
func (b *vectorRHSBridge) EvaluateOnRGrid(y []float64, t float64) []float64 {
	out := make([]float64, len(y))
	b.EvaluateOnRGridInPlace(y, out, t)
	return out
}
func (b *vectorRHSBridge) EvaluateOnRGridInPlace(y, out []float64, t float64) {
	if b.err == nil {
		b.err = b.rhs(t, y, out)
	}
	if b.err != nil {
		for i := range out {
			out[i] = math.NaN()
		}
	}
}

// ComplexImplicit uses the real block system [Re(y),Im(y)] internally. Thus the
// nonlinear Jacobian is genuinely coupled and supports non-holomorphic F.
// ImplicitOptions.Jacobian, if supplied, refers to this 2n-by-2n real system;
// DiagonalJacobian is deliberately rejected since it would drop complex coupling.
// Dense solves are limited to 512 complex entries. For large LINEAR equations,
// use ComplexLinearImplicit, which only needs operator-vector products.
type ComplexImplicit struct {
	rhs        ComplexRHS
	realSolver *VectorImplicit
	state      []float64
	y, dy      []complex128
}

func NewComplexImplicit(rhs ComplexRHS, method ImplicitMethod, options ImplicitOptions) (*ComplexImplicit, error) {
	if rhs == nil || options.DiagonalJacobian != nil || options.Derivative != nil {
		return nil, fmt.Errorf("complex equations require a whole-vector RHS and a full real-block Jacobian (or finite differences)")
	}
	s := &ComplexImplicit{rhs: rhs}
	v, err := NewVectorImplicit(s.evaluate, method, options)
	if err != nil {
		return nil, err
	}
	s.realSolver = v
	return s, nil
}

func (s *ComplexImplicit) evaluate(t float64, y, dy []float64) error {
	n := len(s.y)
	for i := 0; i < n; i++ {
		s.y[i] = complex(y[i], y[n+i])
	}
	if err := s.rhs(t, s.y, s.dy); err != nil {
		return err
	}
	for i, v := range s.dy {
		dy[i], dy[n+i] = real(v), imag(v)
	}
	return nil
}

func (s *ComplexImplicit) Step(y []complex128, t, dt float64) error {
	n := len(y)
	if n == 0 || n > 512 {
		return fmt.Errorf("nonlinear complex dense solve requires 1..512 entries; use ComplexLinearImplicit for large linear equations")
	}
	if len(s.y) != n {
		s.y = make([]complex128, n)
		s.dy = make([]complex128, n)
		s.state = make([]float64, 2*n)
	}
	for i, v := range y {
		s.state[i], s.state[n+i] = real(v), imag(v)
	}
	if err := s.realSolver.Step(s.state, t, dt); err != nil {
		return err
	}
	for i := range y {
		y[i] = complex(s.state[i], s.state[n+i])
	}
	return nil
}
func (s *ComplexImplicit) Stats() ImplicitStats { return s.realSolver.Stats() }
