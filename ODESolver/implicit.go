package EquationSolver

import (
	"fmt"
	"math"

	"GoProject/gridData"
	"gonum.org/v1/gonum/mat"
)

// JacobianFunc fills the row-major n*n RHS Jacobian df_i/dx_j, not I-h*J.
// Input x is read-only. The destination is cleared before each call.
type JacobianFunc func(x []float64, t float64, jacobian []float64) error

// DiagonalJacobianFunc is an opt-in O(n) solve for independent equations ONLY.
// Coupled systems must use the full Jacobian (analytic or finite-difference).
type DiagonalJacobianFunc func(x []float64, t float64, diagonal []float64) error

// ImplicitOptions configures the nonlinear solve, not time-step error control.
// Zero values select defaults. Callbacks must not mutate or retain input slices.
type ImplicitOptions struct {
	AbsoluteTolerance float64 // default 1e-10
	RelativeTolerance float64 // default 1e-8
	MaxIterations     int     // default 30
	JacobianRefresh   int     // default 4; 1 gives full Newton
	Jacobian          JacobianFunc
	DiagonalJacobian  DiagonalJacobianFunc
	Derivative        func(x, t float64) float64 // scalar NextStep only
}

// ImplicitStats describes the last attempted step, including both SDIRK stages.
type ImplicitStats struct {
	Iterations, RHSEvaluations, JacobianEvaluations, Factorizations, Backtracks int
	Residual                                                                    float64 // maximum weighted residual of the last stage; <=1 succeeds
}

func (o ImplicitOptions) normalized() (ImplicitOptions, error) {
	if o.AbsoluteTolerance == 0 {
		o.AbsoluteTolerance = 1e-10
	}
	if o.RelativeTolerance == 0 {
		o.RelativeTolerance = 1e-8
	}
	if o.MaxIterations == 0 {
		o.MaxIterations = 30
	}
	if o.JacobianRefresh == 0 {
		o.JacobianRefresh = 4
	}
	if !odeFinite(o.AbsoluteTolerance) || o.AbsoluteTolerance <= 0 || !odeFinite(o.RelativeTolerance) || o.RelativeTolerance < 0 || o.MaxIterations < 1 || o.JacobianRefresh < 1 {
		return o, fmt.Errorf("invalid implicit tolerances or iteration limits")
	}
	if o.Jacobian != nil && o.DiagonalJacobian != nil {
		return o, fmt.Errorf("supply either a full or a diagonal Jacobian, not both")
	}
	return o, nil
}

type implicitMethod uint8

const (
	backwardEuler implicitMethod = iota
	trapezoidal
	implicitMidpoint
	sdirk2
)

// Solvers own their workspace: reuse sequentially, never concurrently or by
// copying a used solver. There is deliberately no factorization cache across
// steps: callers can change dt, the RHS or its parameters between calls.
type implicitWorkspace struct {
	initial, base, y, z, f, residual, trial, trialF, trialResidual, correction, perturbed, scale, diagonal []float64
	matrix                                                                                                 *mat.Dense
	lu                                                                                                     mat.LU
	rhsVec, correctionVec                                                                                  *mat.VecDense
	stats                                                                                                  ImplicitStats
	scalar                                                                                                 bool
	scalarState                                                                                            [1]float64
}

func odeFinite(x float64) bool { return !math.IsNaN(x) && !math.IsInf(x, 0) }
func odeFiniteSlice(x []float64) bool {
	for _, v := range x {
		if !odeFinite(v) {
			return false
		}
	}
	return true
}

func (w *implicitWorkspace) allocate(n int, dense bool) {
	if len(w.y) != n {
		data := make([]float64, 13*n)
		parts := []*[]float64{&w.initial, &w.base, &w.y, &w.z, &w.f, &w.residual, &w.trial, &w.trialF, &w.trialResidual, &w.correction, &w.perturbed, &w.scale, &w.diagonal}
		for i, p := range parts {
			*p = data[i*n : (i+1)*n]
		}
		w.matrix = nil
		w.lu = mat.LU{}
		w.rhsVec = mat.NewVecDense(n, w.residual)
		w.correctionVec = mat.NewVecDense(n, w.correction)
	}
	if dense && w.matrix == nil {
		w.matrix = mat.NewDense(n, n, nil)
	}
}

func (w *implicitWorkspace) evaluate(f gridData.TDPotentialOp, x, dst []float64, t float64) error {
	if !odeFiniteSlice(x) {
		return fmt.Errorf("non-finite implicit stage")
	}
	w.stats.RHSEvaluations++
	if w.scalar {
		dst[0] = f.EvaluateAt(x[0], t)
	} else {
		f.EvaluateOnRGridInPlace(x, dst, t)
	}
	if !odeFiniteSlice(dst) {
		return fmt.Errorf("non-finite ODE right-hand side at time %g", t)
	}
	return nil
}

// Each stage solves y-base-coefficient*f((1-theta)*initial+theta*y,time)=0.
func (w *implicitWorkspace) stageResidual(f gridData.TDPotentialOp, y, fy, residual []float64, time, coefficient, theta float64) (float64, error) {
	for i := range y {
		w.z[i] = (1-theta)*w.initial[i] + theta*y[i]
	}
	if err := w.evaluate(f, w.z, fy, time); err != nil {
		return 0, err
	}
	norm := 0.0
	for i := range y {
		residual[i] = y[i] - w.base[i] - coefficient*fy[i]
		value := math.Abs(residual[i]) / w.scale[i]
		if !odeFinite(value) {
			return 0, fmt.Errorf("non-finite implicit residual")
		}
		norm = math.Max(norm, value)
	}
	return norm, nil
}

func (w *implicitWorkspace) factor(f gridData.TDPotentialOp, o ImplicitOptions, time, coefficient, theta float64) error {
	n := len(w.y)
	for i := range w.y {
		w.z[i] = (1-theta)*w.initial[i] + theta*w.y[i]
	}
	w.stats.JacobianEvaluations++
	diagonal := w.scalar || o.DiagonalJacobian != nil
	if w.scalar {
		if o.Derivative != nil {
			w.diagonal[0] = o.Derivative(w.z[0], time)
		} else {
			x := w.z[0]
			step := math.Sqrt(math.Nextafter(1, 2)-1) * math.Max(1, math.Abs(x))
			w.z[0] = x + step
			step = w.z[0] - x
			err := w.evaluate(f, w.z, w.perturbed, time)
			w.z[0] = x
			if err != nil {
				return err
			}
			w.diagonal[0] = (w.perturbed[0] - w.f[0]) / step
		}
	} else if o.DiagonalJacobian != nil {
		clear(w.diagonal)
		if err := o.DiagonalJacobian(w.z, time, w.diagonal); err != nil {
			return fmt.Errorf("diagonal Jacobian: %w", err)
		}
	} else {
		data := w.matrix.RawMatrix().Data
		clear(data)
		if o.Jacobian != nil {
			if err := o.Jacobian(w.z, time, data); err != nil {
				return fmt.Errorf("Jacobian: %w", err)
			}
		} else {
			for j := 0; j < n; j++ {
				x := w.z[j]
				step := math.Sqrt(math.Nextafter(1, 2)-1) * math.Max(1, math.Abs(x))
				w.z[j] = x + step
				step = w.z[j] - x
				err := w.evaluate(f, w.z, w.perturbed, time)
				w.z[j] = x
				if err != nil {
					return err
				}
				for i := 0; i < n; i++ {
					data[i*n+j] = (w.perturbed[i] - w.f[i]) / step
				}
			}
		}
		for i := range data {
			data[i] *= -coefficient * theta
		}
		for i := 0; i < n; i++ {
			data[i*n+i] += 1
		}
		if !odeFiniteSlice(data) {
			return fmt.Errorf("non-finite Newton matrix")
		}
		w.lu.Factorize(w.matrix)
	}
	if diagonal {
		for i, derivative := range w.diagonal {
			v := 1 - coefficient*theta*derivative
			if !odeFinite(v) || math.Abs(v) <= 1e-14 {
				return fmt.Errorf("singular or non-finite Newton diagonal at component %d", i)
			}
			w.diagonal[i] = v
		}
	}
	w.stats.Factorizations++
	return nil
}

func (w *implicitWorkspace) solveStage(f gridData.TDPotentialOp, o ImplicitOptions, time, coefficient, theta float64, newton bool) error {
	w.stats.Residual = math.Inf(1)
	norm, err := w.stageResidual(f, w.y, w.f, w.residual, time, coefficient, theta)
	if err != nil {
		return err
	}
	w.stats.Residual = norm
	age, refresh := 0, true
	for iteration := 0; norm > 1 && iteration < o.MaxIterations; iteration++ {
		w.stats.Iterations++
		if !newton {
			for i := range w.y {
				w.trial[i] = w.y[i] - w.residual[i]
			}
			next, err := w.stageResidual(f, w.trial, w.trialF, w.trialResidual, time, coefficient, theta)
			if err != nil {
				return err
			}
			copy(w.y, w.trial)
			copy(w.f, w.trialF)
			copy(w.residual, w.trialResidual)
			norm = next
		} else {
			accepted := false
			// A failed stale-Jacobian search gets one retry with a fresh Jacobian.
			for attempt := 0; attempt < 2 && !accepted; attempt++ {
				fresh := refresh || age >= o.JacobianRefresh
				if fresh {
					if err := w.factor(f, o, time, coefficient, theta); err != nil {
						return err
					}
					age = 0
					refresh = false
				}
				if w.scalar || o.DiagonalJacobian != nil {
					for i := range w.y {
						w.correction[i] = w.residual[i] / w.diagonal[i]
					}
				} else if err := w.lu.SolveVecTo(w.correctionVec, false, w.rhsVec); err != nil {
					return fmt.Errorf("Newton linear solve: %w", err)
				}
				if !odeFiniteSlice(w.correction) {
					return fmt.Errorf("non-finite Newton correction")
				}
				alpha := 1.0
				for backtrack := 0; backtrack <= 12; backtrack++ {
					for i := range w.y {
						w.trial[i] = w.y[i] - alpha*w.correction[i]
					}
					next, trialErr := w.stageResidual(f, w.trial, w.trialF, w.trialResidual, time, coefficient, theta)
					if trialErr == nil && (next <= 1 || next <= (1-1e-4*alpha)*norm) {
						refresh = next > .7*norm
						copy(w.y, w.trial)
						copy(w.f, w.trialF)
						copy(w.residual, w.trialResidual)
						norm = next
						accepted = true
						age++
						break
					}
					if backtrack < 12 {
						w.stats.Backtracks++
						alpha *= .5
					}
				}
				if !accepted {
					if fresh {
						break
					}
					refresh = true
				}
			}
			if !accepted {
				return fmt.Errorf("Newton line search failed; reduce the time step (weighted residual %g)", norm)
			}
		}
		w.stats.Residual = norm
	}
	if norm > 1 {
		return fmt.Errorf("implicit step did not converge after %d iterations (weighted residual %g); reduce the time step", o.MaxIterations, norm)
	}
	return nil
}

func (w *implicitWorkspace) step(f gridData.TDPotentialOp, x []float64, t, dt float64, options ImplicitOptions, method implicitMethod, newton, scalar bool) error {
	w.stats = ImplicitStats{Residual: math.Inf(1)}
	o, err := options.normalized()
	if err != nil {
		return err
	}
	if f == nil || !odeFinite(t) || !odeFinite(dt) || !odeFinite(t+dt) || !odeFiniteSlice(x) || len(x) == 0 {
		return fmt.Errorf("invalid implicit initial state, RHS or time")
	}
	if dt == 0 {
		w.stats.Residual = 0
		return nil
	}
	if t+dt == t {
		return fmt.Errorf("time step is too small to advance time")
	}
	dense := newton && !scalar && o.DiagonalJacobian == nil
	if dense && len(x) > 1024 {
		return fmt.Errorf("dense implicit solves are limited to 1024 components; use a diagonal Jacobian only for independent equations")
	}
	w.scalar = scalar
	w.allocate(len(x), dense)
	copy(w.initial, x)
	copy(w.base, x)
	copy(w.y, x)
	for i, v := range x {
		w.scale[i] = o.AbsoluteTolerance + o.RelativeTolerance*math.Abs(v)
		if !odeFinite(w.scale[i]) {
			return fmt.Errorf("implicit error scale overflow")
		}
	}
	time, coefficient, theta := t+dt, dt, 1.0
	switch method {
	case trapezoidal:
		if err := w.evaluate(f, w.initial, w.f, t); err != nil {
			return err
		}
		for i := range x {
			w.base[i] += dt / 2 * w.f[i]
			w.y[i] += dt * w.f[i]
		}
		if !odeFiniteSlice(w.y) {
			copy(w.y, x)
		}
		coefficient = dt / 2
	case implicitMidpoint:
		time, theta = t+dt/2, .5
	case sdirk2:
		gamma := 1 - 1/math.Sqrt(2)
		time, coefficient = t+gamma*dt, gamma*dt
	}
	if err := w.solveStage(f, o, time, coefficient, theta, newton); err != nil {
		return err
	}
	if method == sdirk2 {
		gamma := 1 - 1/math.Sqrt(2)
		for i := range x {
			w.base[i] = w.initial[i] + dt*(1-gamma)*w.f[i]
		}
		if err := w.solveStage(f, o, t+dt, gamma*dt, 1, true); err != nil {
			return fmt.Errorf("SDIRK second stage: %w", err)
		}
	}
	copy(x, w.y) // Commit only after all stages have converged.
	return nil
}

func (w *implicitWorkspace) scalarStep(f gridData.TDPotentialOp, x, t, dt float64, options ImplicitOptions, method implicitMethod, newton bool) (float64, error) {
	w.scalarState[0] = x
	if err := w.step(f, w.scalarState[:], t, dt, options, method, newton, true); err != nil {
		return x, err
	}
	return w.scalarState[0], nil
}
