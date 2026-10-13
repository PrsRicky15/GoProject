package EquationSolver

import (
	"context"
	"fmt"
	"math"
	"math/cmplx"
)

type rkNumber interface{ float64 | complex128 }

// Each solver owns reusable workspace. Do not copy it or use it concurrently.
type VectorRK struct{ *rkEngine[float64] }
type ComplexRK struct{ *rkEngine[complex128] }

func NewVectorRK(rhs VectorRHS, method RKMethod) (*VectorRK, error) {
	e, err := newRK[float64](rhs, method, func(v, a float64) float64 { return v * a }, math.Abs)
	if err != nil {
		return nil, err
	}
	return &VectorRK{e}, nil
}
func NewComplexRK(rhs ComplexRHS, method RKMethod) (*ComplexRK, error) {
	e, err := newRK[complex128](rhs, method, func(v complex128, a float64) complex128 { return v * complex(a, 0) }, cmplx.Abs)
	if err != nil {
		return nil, err
	}
	return &ComplexRK{e}, nil
}

type rkEngine[T rkNumber] struct {
	rhs                func(float64, []T, []T) error
	scale              func(T, float64) T
	abs                func(T) float64
	q                  rkTableau
	k                  [][]T
	stage, next, state []T
}

func newRK[T rkNumber](rhs func(float64, []T, []T) error, method RKMethod, scale func(T, float64) T, abs func(T) float64) (*rkEngine[T], error) {
	q, ok := rkTable(method)
	if rhs == nil || !ok {
		return nil, fmt.Errorf("invalid RHS or Runge-Kutta method %q", method)
	}
	return &rkEngine[T]{rhs: rhs, scale: scale, abs: abs, q: q}, nil
}
func (e *rkEngine[T]) prepare(y []T, t, end float64) error {
	if len(y) == 0 || !odeFinite(t) || !odeFinite(end) || !odeFinite(end-t) {
		return fmt.Errorf("invalid state dimension or time interval")
	}
	for _, v := range y {
		if !odeFinite(e.abs(v)) {
			return fmt.Errorf("non-finite initial state")
		}
	}
	// Conservative bound using complex128 size for both scalar types.
	if len(y) > (256<<20)/(16*(len(e.q.b)+3)) {
		return fmt.Errorf("Runge-Kutta workspace exceeds 256 MiB")
	}
	if len(e.state) != len(y) {
		n := len(y)
		data := make([]T, n*(len(e.q.b)+3))
		e.k = make([][]T, len(e.q.b))
		for j := range e.k {
			e.k[j] = data[j*n : (j+1)*n]
		}
		o := len(e.k) * n
		e.stage = data[o : o+n]
		e.next = data[o+n : o+2*n]
		e.state = data[o+2*n:]
	}
	return nil
}

// Step advances the entire vector by dt, committing only on success.
func (e *rkEngine[T]) Step(y []T, t, dt float64) error {
	return e.StepContext(context.Background(), y, t, dt)
}
func (e *rkEngine[T]) StepContext(ctx context.Context, y []T, t, dt float64) error {
	if ctx == nil {
		return fmt.Errorf("nil context")
	}
	if err := ctx.Err(); err != nil {
		return err
	}
	if err := e.prepare(y, t, t+dt); err != nil {
		return err
	}
	if dt == 0 {
		return nil
	}
	if t+dt == t {
		return fmt.Errorf("step is too small to advance time")
	}
	_, err := e.trial(ctx, y, t, dt, false)
	if err != nil {
		return err
	}
	copy(y, e.next)
	return nil
}
func (e *rkEngine[T]) trial(ctx context.Context, y []T, t, h float64, reuseFirst bool) (calls int, err error) {
	for s := range e.k {
		if s == 0 && reuseFirst {
			continue
		}
		if err = ctx.Err(); err != nil {
			return
		}
		copy(e.stage, y)
		for j, a := range e.q.a[s] {
			if a != 0 {
				for i := range y {
					e.stage[i] += e.scale(e.k[j][i], h*a)
				}
			}
		}
		for _, v := range e.stage {
			if !odeFinite(e.abs(v)) {
				return calls, fmt.Errorf("non-finite RK stage %d state", s)
			}
		}
		st := t + e.q.c[s]*h
		if !odeFinite(st) {
			return calls, fmt.Errorf("non-finite RK stage time")
		}
		clear(e.k[s])
		calls++
		if err = e.rhs(st, e.stage, e.k[s]); err != nil {
			return calls, fmt.Errorf("RK stage %d RHS: %w", s, err)
		}
		for _, v := range e.k[s] {
			if !odeFinite(e.abs(v)) {
				return calls, fmt.Errorf("non-finite RK stage %d derivative", s)
			}
		}
	}
	copy(e.next, y)
	for j, b := range e.q.b {
		if b != 0 {
			for i := range y {
				e.next[i] += e.scale(e.k[j][i], h*b)
			}
		}
	}
	for _, v := range e.next {
		if !odeFinite(e.abs(v)) {
			return calls, fmt.Errorf("non-finite RK output")
		}
	}
	return calls, ctx.Err()
}

// AdaptiveOptions controls LOCAL componentwise error, not global error.
// Zero values select defaults. Step limits are positive magnitudes in either
// time direction; MinStep=0 permits any representable step.
type AdaptiveOptions struct {
	AbsoluteTolerance, RelativeTolerance float64
	InitialStep, MinStep, MaxStep        float64
	MaxAttempts                          int
}

func (o AdaptiveOptions) normalized() (AdaptiveOptions, error) {
	for _, v := range []float64{o.AbsoluteTolerance, o.RelativeTolerance, o.InitialStep, o.MinStep, o.MaxStep} {
		if !odeFinite(v) || v < 0 {
			return o, fmt.Errorf("invalid adaptive tolerance or step limit")
		}
	}
	if o.MaxAttempts < 0 {
		return o, fmt.Errorf("invalid adaptive attempt limit")
	}
	if o.AbsoluteTolerance == 0 {
		o.AbsoluteTolerance = 1e-9
	}
	if o.RelativeTolerance == 0 {
		o.RelativeTolerance = 1e-7
	}
	if o.MaxStep == 0 {
		o.MaxStep = math.MaxFloat64
	}
	if o.MaxAttempts == 0 {
		o.MaxAttempts = 100000
	}
	if o.MinStep > o.MaxStep {
		return o, fmt.Errorf("MinStep exceeds MaxStep")
	}
	return o, nil
}

// AdaptiveResult.Time is the time of the caller's committed vector: tEnd on
// success, tStart on failure. Counts include attempted work on failure.
// NextStep is a signed suggestion, not a hidden change to the caller's time.
type AdaptiveResult struct {
	Time, NextStep                               float64
	AcceptedSteps, RejectedSteps, RHSEvaluations int
}
type VectorAdaptiveRK struct{ *adaptiveRK[float64] }
type ComplexAdaptiveRK struct{ *adaptiveRK[complex128] }
type adaptiveRK[T rkNumber] struct {
	e       *rkEngine[T]
	options AdaptiveOptions
}

func NewVectorAdaptiveRK(rhs VectorRHS, method RKMethod, options AdaptiveOptions) (*VectorAdaptiveRK, error) {
	s, err := NewVectorRK(rhs, method)
	if err != nil {
		return nil, err
	}
	a, err := newAdaptive(s.rkEngine, options)
	if err != nil {
		return nil, err
	}
	return &VectorAdaptiveRK{a}, nil
}
func NewComplexAdaptiveRK(rhs ComplexRHS, method RKMethod, options AdaptiveOptions) (*ComplexAdaptiveRK, error) {
	s, err := NewComplexRK(rhs, method)
	if err != nil {
		return nil, err
	}
	a, err := newAdaptive(s.rkEngine, options)
	if err != nil {
		return nil, err
	}
	return &ComplexAdaptiveRK{a}, nil
}
func newAdaptive[T rkNumber](e *rkEngine[T], options AdaptiveOptions) (*adaptiveRK[T], error) {
	if len(e.q.low) == 0 {
		return nil, fmt.Errorf("method has no embedded error estimator")
	}
	o, err := options.normalized()
	if err != nil {
		return nil, err
	}
	return &adaptiveRK[T]{e: e, options: o}, nil
}

// Integrate advances a complete interval with internal accepted/rejected steps.
// On any failure y is unchanged, even if earlier internal steps succeeded.
func (a *adaptiveRK[T]) Integrate(y []T, tStart, tEnd float64) (AdaptiveResult, error) {
	return a.IntegrateContext(context.Background(), y, tStart, tEnd)
}
func (a *adaptiveRK[T]) IntegrateContext(ctx context.Context, y []T, tStart, tEnd float64) (r AdaptiveResult, err error) {
	r.Time = tStart
	if ctx == nil {
		return r, fmt.Errorf("nil context")
	}
	if err = ctx.Err(); err != nil {
		return
	}
	e, o := a.e, a.options
	if err = e.prepare(y, tStart, tEnd); err != nil {
		return
	}
	if tStart == tEnd {
		return
	}
	copy(e.state, y)
	direction := math.Copysign(1, tEnd-tStart)
	h := o.InitialStep
	if h == 0 {
		h = math.Abs(tEnd-tStart) / 10
		if h == 0 {
			h = math.Abs(tEnd - tStart)
		}
	}
	h = math.Min(o.MaxStep, math.Max(o.MinStep, h))
	t := tStart
	reuse := false
	for attempts := 0; attempts < o.MaxAttempts; attempts++ {
		if err = ctx.Err(); err != nil {
			return
		}
		remaining := math.Abs(tEnd - t)
		step := math.Min(h, remaining)
		signed := direction * step
		if step == 0 || t+signed == t {
			return r, fmt.Errorf("adaptive step cannot advance floating-point time")
		}
		calls, trialErr := e.trial(ctx, e.state, t, signed, reuse)
		r.RHSEvaluations += calls
		if trialErr != nil {
			return r, trialErr
		}
		reuse = true
		errorNorm := 0.0
		for i, v := range e.state {
			var delta T
			for j, b := range e.q.b {
				delta += e.scale(e.k[j][i], signed*(b-e.q.low[j]))
			}
			denom := o.AbsoluteTolerance + o.RelativeTolerance*math.Max(e.abs(v), e.abs(e.next[i]))
			if !odeFinite(denom) {
				return r, fmt.Errorf("adaptive error scale overflow")
			}
			errorNorm = math.Max(errorNorm, e.abs(delta)/denom)
		}
		if !odeFinite(errorNorm) {
			return r, fmt.Errorf("non-finite adaptive error estimate")
		}
		factor := 5.0
		if errorNorm > 0 {
			factor = math.Max(.2, math.Min(5, .9*math.Pow(errorNorm, -1/float64(e.q.lowOrder+1))))
		}
		if errorNorm <= 1 {
			r.AcceptedSteps++
			copy(e.state, e.next)
			h = math.Min(o.MaxStep, math.Max(o.MinStep, step*factor))
			r.NextStep = direction * h
			// A capped step can round directly to tEnd even when subtraction
			// made remaining a fraction of an ulp larger than step. Do not
			// attempt a zero-length extra step after reaching that endpoint.
			if step == remaining || t+signed == tEnd {
				if err = ctx.Err(); err != nil {
					return
				}
				copy(y, e.state)
				r.Time = tEnd
				return
			}
			t += signed
			reuse = e.q.fsal
			if reuse {
				copy(e.k[0], e.k[len(e.k)-1])
			}
		} else {
			r.RejectedSteps++
			if step <= o.MinStep {
				return r, fmt.Errorf("tolerance cannot be met at minimum step")
			}
			h = math.Max(o.MinStep, step*math.Min(1, factor))
		}
	}
	return r, fmt.Errorf("adaptive attempt limit exceeded")
}
