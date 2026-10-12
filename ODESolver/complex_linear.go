package EquationSolver

import (
	"context"
	"fmt"
	"math"
	"math/cmplx"
)

// KrylovOptions controls restarted GMRES for linear implicit vector stages.
// Zero fields select defaults. Tolerances apply to the true Euclidean residual
// of each linear stage, not the time-integration error.
type KrylovOptions struct {
	AbsoluteTolerance float64 // default 1e-12
	RelativeTolerance float64 // default 1e-10
	Restart           int     // default 30, maximum 256
	MaxIterations     int     // default 300, per stage
}
type KrylovStats struct {
	Iterations, OperatorApplications, Restarts int
	Residual                                   float64 // true residual norm of last attempted stage
}

// ComplexLinearImplicit solves y'=A(t)y using only whole-vector products A(t)y.
// The supplied function MUST be linear in y; it may depend arbitrarily on time.
// No full matrix or finite-difference Jacobian is assembled. Storage is O(n*m)
// for restart length m. Instances reuse buffers and must not be shared concurrently.
type ComplexLinearImplicit struct {
	n, m                                   int
	apply                                  ComplexRHS
	method                                 ImplicitMethod
	options                                KrylovOptions
	initial, rhs, result, opWork, residual []complex128
	basis, h, g, coefficients, sines       []complex128
	cosines                                []float64
	stats                                  KrylovStats
}

func NewComplexLinearImplicit(n int, apply ComplexRHS, method ImplicitMethod, o KrylovOptions) (*ComplexLinearImplicit, error) {
	if n < 1 || apply == nil || !method.valid() {
		return nil, fmt.Errorf("invalid complex linear equation or method")
	}
	if o.AbsoluteTolerance == 0 {
		o.AbsoluteTolerance = 1e-12
	}
	if o.RelativeTolerance == 0 {
		o.RelativeTolerance = 1e-10
	}
	if o.Restart == 0 {
		o.Restart = 30
	}
	if o.MaxIterations == 0 {
		o.MaxIterations = 300
	}
	if !odeFinite(o.AbsoluteTolerance) || o.AbsoluteTolerance <= 0 || !odeFinite(o.RelativeTolerance) || o.RelativeTolerance < 0 || o.Restart < 1 || o.Restart > 256 || o.MaxIterations < 1 {
		return nil, fmt.Errorf("invalid Krylov tolerances or iteration limits")
	}
	m := min(n, o.Restart, o.MaxIterations)
	if 16*float64(n)*float64(m+6)+16*float64((m+1)*m+4*m+1) > 256*1024*1024 {
		return nil, fmt.Errorf("Krylov workspace exceeds 256 MiB; reduce restart length or vector size")
	}
	s := &ComplexLinearImplicit{n: n, m: m, apply: apply, method: method, options: o}
	data := make([]complex128, 5*n)
	s.initial = data[:n:n]
	s.rhs = data[n : 2*n : 2*n]
	s.result = data[2*n : 3*n : 3*n]
	s.opWork = data[3*n : 4*n : 4*n]
	s.residual = data[4*n : 5*n : 5*n]
	s.basis = make([]complex128, (m+1)*n)
	s.h = make([]complex128, (m+1)*m)
	s.g = make([]complex128, m+1)
	s.coefficients = make([]complex128, m)
	s.sines = make([]complex128, m)
	s.cosines = make([]float64, m)
	return s, nil
}

func complexFinite(z complex128) bool { return odeFinite(real(z)) && odeFinite(imag(z)) }
func complexVectorFinite(x []complex128) bool {
	for _, v := range x {
		if !complexFinite(v) {
			return false
		}
	}
	return true
}
func complexNorm(x []complex128) float64 {
	n := 0.
	for _, v := range x {
		n = math.Hypot(n, cmplx.Abs(v))
	}
	return n
}
func complexDot(x, y []complex128) complex128 {
	var dot complex128
	for i, v := range x {
		dot += cmplx.Conj(v) * y[i]
	}
	return dot
}

func (s *ComplexLinearImplicit) product(ctx context.Context, t float64, in, out []complex128) error {
	if err := ctx.Err(); err != nil {
		return err
	}
	if !complexVectorFinite(in) {
		return fmt.Errorf("non-finite operator input")
	}
	s.stats.OperatorApplications++
	if err := s.apply(t, in, out); err != nil {
		return fmt.Errorf("complex operator: %w", err)
	}
	if !complexVectorFinite(out) {
		return fmt.Errorf("non-finite operator output")
	}
	return nil
}
func (s *ComplexLinearImplicit) shifted(ctx context.Context, t, c float64, in, out []complex128) error {
	if err := s.product(ctx, t, in, out); err != nil {
		return err
	}
	for i := range in {
		out[i] = in[i] - complex(c, 0)*out[i]
	}
	if !complexVectorFinite(out) {
		return fmt.Errorf("non-finite shifted operator output")
	}
	return nil
}

// rotation generates [c s; -conj(s) c] with a real cosine, annihilating b.
func complexRotation(a, b complex128) (float64, complex128) {
	if b == 0 {
		return 1, 0
	}
	absA, absB := cmplx.Abs(a), cmplx.Abs(b)
	if absA == 0 {
		return 0, cmplx.Conj(b) / complex(absB, 0)
	}
	norm := math.Hypot(absA, absB)
	return absA / norm, (a / complex(absA, 0)) * (cmplx.Conj(b) / complex(norm, 0))
}

// solve computes (I-c*A(t))*result=rhs, starting from result. Modified
// Gram-Schmidt is repeated twice, followed by complex Givens QR rotations.
// Acceptance always checks the actual residual, including after a restart.
func (s *ComplexLinearImplicit) solve(ctx context.Context, t, c float64) error {
	s.stats.Residual = math.Inf(1)
	target := s.options.AbsoluteTolerance + s.options.RelativeTolerance*complexNorm(s.rhs)
	if !odeFinite(target) {
		return fmt.Errorf("linear tolerance scale overflow")
	}
	iterations := 0
	for {
		if err := s.shifted(ctx, t, c, s.result, s.opWork); err != nil {
			return err
		}
		for i := range s.rhs {
			s.residual[i] = s.rhs[i] - s.opWork[i]
		}
		beta := complexNorm(s.residual)
		s.stats.Residual = beta
		if !odeFinite(beta) {
			return fmt.Errorf("non-finite Krylov residual")
		}
		if beta <= target {
			return nil
		}
		if iterations >= s.options.MaxIterations {
			return fmt.Errorf("GMRES did not converge after %d iterations: residual %g, target %g; reduce dt or increase Krylov limits", iterations, beta, target)
		}
		if iterations > 0 {
			s.stats.Restarts++
		}
		clear(s.h)
		clear(s.g)
		s.g[0] = complex(beta, 0)
		for i := 0; i < s.n; i++ {
			s.basis[i] = s.residual[i] / complex(beta, 0)
		}
		used := 0
		for j := 0; j < s.m && iterations < s.options.MaxIterations; j++ {
			v := s.basis[j*s.n : (j+1)*s.n]
			if err := s.shifted(ctx, t, c, v, s.opWork); err != nil {
				return err
			}
			originalNorm := complexNorm(s.opWork)
			for pass := 0; pass < 2; pass++ {
				for i := 0; i <= j; i++ {
					vi := s.basis[i*s.n : (i+1)*s.n]
					d := complexDot(vi, s.opWork)
					s.h[i*s.m+j] += d
					for k := range vi {
						s.opWork[k] -= d * vi[k]
					}
				}
			}
			nextNorm := complexNorm(s.opWork)
			if !odeFinite(originalNorm) || !odeFinite(nextNorm) {
				return fmt.Errorf("non-finite Arnoldi basis")
			}
			breakdown := nextNorm <= 1e-14*originalNorm
			if !breakdown {
				s.h[(j+1)*s.m+j] = complex(nextNorm, 0)
				for k := 0; k < s.n; k++ {
					s.basis[(j+1)*s.n+k] = s.opWork[k] / complex(nextNorm, 0)
				}
			}
			for i := 0; i < j; i++ {
				a, b := s.h[i*s.m+j], s.h[(i+1)*s.m+j]
				cs, sn := complex(s.cosines[i], 0), s.sines[i]
				s.h[i*s.m+j] = cs*a + sn*b
				s.h[(i+1)*s.m+j] = -cmplx.Conj(sn)*a + cs*b
			}
			a, b := s.h[j*s.m+j], s.h[(j+1)*s.m+j]
			s.cosines[j], s.sines[j] = complexRotation(a, b)
			cs, sn := complex(s.cosines[j], 0), s.sines[j]
			s.h[j*s.m+j] = cs*a + sn*b
			s.h[(j+1)*s.m+j] = 0
			s.g[j+1] = -cmplx.Conj(sn) * s.g[j]
			s.g[j] *= cs
			iterations++
			s.stats.Iterations++
			used = j + 1
			if !complexFinite(s.g[j+1]) || !complexFinite(s.h[j*s.m+j]) {
				return fmt.Errorf("non-finite Krylov factorization")
			}
			if cmplx.Abs(s.g[j+1]) <= target || breakdown {
				break
			}
		}
		for i := used - 1; i >= 0; i-- {
			v := s.g[i]
			for j := i + 1; j < used; j++ {
				v -= s.h[i*s.m+j] * s.coefficients[j]
			}
			d := s.h[i*s.m+i]
			if d == 0 {
				return fmt.Errorf("singular Krylov system")
			}
			s.coefficients[i] = v / d
		}
		for j := 0; j < used; j++ {
			for i := 0; i < s.n; i++ {
				s.result[i] += s.basis[j*s.n+i] * s.coefficients[j]
			}
		}
	}
}

func (s *ComplexLinearImplicit) Stats() KrylovStats { return s.stats }
func (s *ComplexLinearImplicit) Step(y []complex128, t, dt float64) error {
	return s.StepContext(context.Background(), y, t, dt)
}

// StepContext supports cancellation between operator products. It never modifies
// y on failure, and never renormalizes the state to conceal propagation error.
func (s *ComplexLinearImplicit) StepContext(ctx context.Context, y []complex128, t, dt float64) error {
	s.stats = KrylovStats{Residual: math.Inf(1)}
	if ctx == nil {
		return fmt.Errorf("nil propagation context")
	}
	if err := ctx.Err(); err != nil {
		return err
	}
	if len(y) != s.n || !complexVectorFinite(y) || !odeFinite(t) || !odeFinite(dt) || !odeFinite(t+dt) {
		return fmt.Errorf("invalid complex state or time")
	}
	if dt == 0 {
		s.stats.Residual = 0
		return nil
	}
	if t+dt == t {
		return fmt.Errorf("time step is too small to advance time")
	}
	copy(s.initial, y)
	copy(s.rhs, y)
	copy(s.result, y)
	time, coefficient := t+dt, dt
	switch s.method {
	case ImplicitTrapezoidalMethod, ImplicitMidpointMethod:
		coefficient = dt / 2
		sourceTime := t
		if s.method == ImplicitMidpointMethod {
			time = t + dt/2
			sourceTime = time
		}
		if err := s.product(ctx, sourceTime, s.initial, s.opWork); err != nil {
			return err
		}
		for i := range y {
			s.rhs[i] += complex(dt/2, 0) * s.opWork[i]
		}
	case ImplicitSDIRK2Method:
		g := 1 - 1/math.Sqrt(2)
		time = t + g*dt
		coefficient = g * dt
	}
	if err := s.solve(ctx, time, coefficient); err != nil {
		return err
	}
	if s.method == ImplicitSDIRK2Method {
		s.stats.Residual = math.Inf(1)
		g := 1 - 1/math.Sqrt(2)
		if err := s.product(ctx, t+g*dt, s.result, s.opWork); err != nil {
			return err
		}
		for i := range y {
			s.rhs[i] = s.initial[i] + complex(dt*(1-g), 0)*s.opWork[i]
		}
		if err := s.solve(ctx, t+dt, g*dt); err != nil {
			return fmt.Errorf("SDIRK second stage: %w", err)
		}
	}
	if err := ctx.Err(); err != nil {
		return err
	}
	copy(y, s.result)
	return nil
}
