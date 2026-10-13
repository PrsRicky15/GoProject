package simulation

import (
	EquationSolver "GoProject/ODESolver"
	"context"
	"fmt"
)

var fixedRK = map[string]EquationSolver.RKMethod{
	"euler": EquationSolver.RKEuler, "heun2": EquationSolver.RKHeun2,
	"midpoint2": EquationSolver.RKMidpoint2, "ralston2": EquationSolver.RKRalston2,
	"ralston3": EquationSolver.RKRalston3, "heun3": EquationSolver.RKHeun3,
	"wray3": EquationSolver.RKWray3, "ssp3": EquationSolver.RKSSP3,
	"kutta3": EquationSolver.RKKutta3, "rk4": EquationSolver.RKClassic4,
	"rk38": EquationSolver.RKThreeEighths4, "ralston4": EquationSolver.RKRalston4,
	"nystrom5": EquationSolver.RKNystrom5,
}
var adaptiveRK = map[string]EquationSolver.RKMethod{
	"heun-euler21": EquationSolver.RKHeunEuler21, "fehlberg21": EquationSolver.RKFehlberg21,
	"bogacki-shampine32": EquationSolver.RKBogackiShampine32, "fehlberg54": EquationSolver.RKFehlberg54,
	"cash-karp54": EquationSolver.RKCashKarp54, "dormand-prince54": EquationSolver.RKDormandPrince54,
}
var implicitODE = map[string]EquationSolver.ImplicitMethod{
	"implicit-euler":    EquationSolver.ImplicitEulerMethod,
	"trapezoidal":       EquationSolver.ImplicitTrapezoidalMethod,
	"implicit-midpoint": EquationSolver.ImplicitMidpointMethod,
	"sdirk2":            EquationSolver.ImplicitSDIRK2Method,
}

func vectorMethod(method string) bool {
	_, f := fixedRK[method]
	_, a := adaptiveRK[method]
	_, i := implicitODE[method]
	return f || a || i
}

type solverWork struct{ evaluations, accepted, rejected, limit int }

func (w *solverWork) evaluation(ctx context.Context) error {
	if err := ctx.Err(); err != nil {
		return err
	}
	if w.evaluations >= w.limit {
		return fmt.Errorf("solver work limit reached; reduce the grid or duration, or relax tolerances")
	}
	w.evaluations++
	return nil
}
func (w *solverWork) adaptive(r EquationSolver.AdaptiveResult) {
	w.accepted += r.AcceptedSteps
	w.rejected += r.RejectedSteps
}
func (w *solverWork) metrics() []Metric {
	return []Metric{{"RHS evaluations", float64(w.evaluations)}, {"Accepted solver steps", float64(w.accepted)}, {"Rejected solver steps", float64(w.rejected)}}
}
func solverNotes(r Request) []string {
	if _, ok := adaptiveRK[r.Method]; ok {
		return []string{"Adaptive propagation: dt is the output interval and maximum internal step. Go chooses smaller internal steps from the absolute and relative tolerances. Counts include all output intervals; tolerance controls local error, not global error."}
	}
	if _, ok := implicitODE[r.Method]; ok {
		return []string{"Implicit propagation uses coupled solves. Tolerances control the equation solve, not time discretization error; reduce dt to check accuracy. Backward Euler and SDIRK2 can damp amplitudes and energy. Quantum midpoint and trapezoidal preserve norm for static Hermitian Hamiltonians up to solve error; driven trapezoidal need not preserve norm."}
	}
	return []string{"Explicit Runge–Kutta propagates the complete coupled state. It is not symplectic or exactly unitary; lower-order methods can accumulate substantial energy or norm drift. Compare smaller dt."}
}
func adaptiveOptions(r Request) EquationSolver.AdaptiveOptions {
	return EquationSolver.AdaptiveOptions{AbsoluteTolerance: r.AbsoluteTolerance, RelativeTolerance: r.RelativeTolerance, InitialStep: r.Dt, MaxStep: r.Dt, MaxAttempts: 10000}
}

func realVectorStepper(ctx context.Context, r Request, rhs EquationSolver.VectorRHS, state []float64, w *solverWork) (func(float64, float64) error, error) {
	boundedRHS := func(t float64, y, dy []float64) error {
		if err := w.evaluation(ctx); err != nil {
			return err
		}
		return rhs(t, y, dy)
	}
	if method, ok := adaptiveRK[r.Method]; ok {
		s, err := EquationSolver.NewVectorAdaptiveRK(boundedRHS, method, adaptiveOptions(r))
		if err != nil {
			return nil, err
		}
		return func(t, end float64) error {
			result, err := s.IntegrateContext(ctx, state, t, end)
			w.adaptive(result)
			return err
		}, nil
	}
	if method, ok := implicitODE[r.Method]; ok {
		s, err := EquationSolver.NewVectorImplicit(boundedRHS, method, EquationSolver.ImplicitOptions{AbsoluteTolerance: r.AbsoluteTolerance, RelativeTolerance: r.RelativeTolerance})
		if err != nil {
			return nil, err
		}
		return func(t, end float64) error {
			err := s.Step(state, t, end-t)
			if err == nil {
				w.accepted++
			}
			return err
		}, nil
	}
	s, err := EquationSolver.NewVectorRK(boundedRHS, fixedRK[r.Method])
	if err != nil {
		return nil, err
	}
	return func(t, end float64) error {
		err := s.StepContext(ctx, state, t, end-t)
		if err == nil {
			w.accepted++
		}
		return err
	}, nil
}

func complexVectorStepper(ctx context.Context, r Request, rhs EquationSolver.ComplexRHS, state []complex128, w *solverWork) (func(float64, float64) error, error) {
	boundedRHS := func(t float64, y, dy []complex128) error {
		if err := w.evaluation(ctx); err != nil {
			return err
		}
		return rhs(t, y, dy)
	}
	if method, ok := adaptiveRK[r.Method]; ok {
		s, err := EquationSolver.NewComplexAdaptiveRK(boundedRHS, method, adaptiveOptions(r))
		if err != nil {
			return nil, err
		}
		return func(t, end float64) error {
			result, err := s.IntegrateContext(ctx, state, t, end)
			w.adaptive(result)
			return err
		}, nil
	}
	if method, ok := implicitODE[r.Method]; ok {
		s, err := EquationSolver.NewComplexLinearImplicit(len(state), boundedRHS, method, EquationSolver.KrylovOptions{AbsoluteTolerance: r.AbsoluteTolerance, RelativeTolerance: r.RelativeTolerance})
		if err != nil {
			return nil, err
		}
		return func(t, end float64) error {
			err := s.StepContext(ctx, state, t, end-t)
			if err == nil {
				w.accepted++
			}
			return err
		}, nil
	}
	s, err := EquationSolver.NewComplexRK(boundedRHS, fixedRK[r.Method])
	if err != nil {
		return nil, err
	}
	return func(t, end float64) error {
		err := s.StepContext(ctx, state, t, end-t)
		if err == nil {
			w.accepted++
		}
		return err
	}, nil
}
