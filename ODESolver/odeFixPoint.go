package EquationSolver

import "GoProject/gridData"

// EulerImplicitFixPoint implements backward euler with fixed-point iteration.
// Reuse one instance per integration; it owns mutable workspace.
type EulerImplicitFixPoint struct {
	TdFunc  gridData.TDPotentialOp
	DeltaT  float64
	Options ImplicitOptions
	work    implicitWorkspace
}

func (s *EulerImplicitFixPoint) Name() string { return "Backward Euler with fixed-point iteration" }
func (s *EulerImplicitFixPoint) NewDefine(dt float64, f gridData.TDPotentialOp) *EulerImplicitFixPoint {
	return &EulerImplicitFixPoint{TdFunc: f, DeltaT: dt}
}

// ReDefine preserves configured options and reuses workspace. No matrix is cached across steps.
func (s *EulerImplicitFixPoint) ReDefine(dt float64, f gridData.TDPotentialOp) {
	s.TdFunc = f
	s.DeltaT = dt
}
func (s *EulerImplicitFixPoint) Stats() ImplicitStats { return s.work.stats }
func (s *EulerImplicitFixPoint) NextStep(x, t float64) (float64, error) {
	return s.work.scalarStep(s.TdFunc, x, t, s.DeltaT, s.Options, backwardEuler, false)
}
func (s *EulerImplicitFixPoint) NextStepOnGrid(x []float64, t float64) error {
	return s.work.step(s.TdFunc, x, t, s.DeltaT, s.Options, backwardEuler, false, false)
}

var _ ODESolver = (*EulerImplicitFixPoint)(nil)

// HeunsImplicitFixPoint implements implicit trapezoidal rule with fixed-point iteration.
// Reuse one instance per integration; it owns mutable workspace.
type HeunsImplicitFixPoint struct {
	TdFunc  gridData.TDPotentialOp
	DeltaT  float64
	Options ImplicitOptions
	work    implicitWorkspace
}

func (s *HeunsImplicitFixPoint) Name() string {
	return "Implicit trapezoidal rule with fixed-point iteration"
}
func (s *HeunsImplicitFixPoint) NewDefine(dt float64, f gridData.TDPotentialOp) *HeunsImplicitFixPoint {
	return &HeunsImplicitFixPoint{TdFunc: f, DeltaT: dt}
}

// ReDefine preserves configured options and reuses workspace. No matrix is cached across steps.
func (s *HeunsImplicitFixPoint) ReDefine(dt float64, f gridData.TDPotentialOp) {
	s.TdFunc = f
	s.DeltaT = dt
}
func (s *HeunsImplicitFixPoint) Stats() ImplicitStats { return s.work.stats }
func (s *HeunsImplicitFixPoint) NextStep(x, t float64) (float64, error) {
	return s.work.scalarStep(s.TdFunc, x, t, s.DeltaT, s.Options, trapezoidal, false)
}
func (s *HeunsImplicitFixPoint) NextStepOnGrid(x []float64, t float64) error {
	return s.work.step(s.TdFunc, x, t, s.DeltaT, s.Options, trapezoidal, false, false)
}

var _ ODESolver = (*HeunsImplicitFixPoint)(nil)

// PredictIni writes the explicit Euler predictor into the caller's buffer.
func (s *HeunsImplicitFixPoint) PredictIni(x, f, predictor []float64) {
	for i := range predictor {
		predictor[i] = x[i] + s.DeltaT*f[i]
	}
}

// MidPointFixPoint implements implicit midpoint with fixed-point iteration.
// Reuse one instance per integration; it owns mutable workspace.
type MidPointFixPoint struct {
	TdFunc  gridData.TDPotentialOp
	DeltaT  float64
	Options ImplicitOptions
	work    implicitWorkspace
}

func (s *MidPointFixPoint) Name() string { return "Implicit midpoint with fixed-point iteration" }
func (s *MidPointFixPoint) NewDefine(dt float64, f gridData.TDPotentialOp) *MidPointFixPoint {
	return &MidPointFixPoint{TdFunc: f, DeltaT: dt}
}

// ReDefine preserves configured options and reuses workspace. No matrix is cached across steps.
func (s *MidPointFixPoint) ReDefine(dt float64, f gridData.TDPotentialOp) {
	s.TdFunc = f
	s.DeltaT = dt
}
func (s *MidPointFixPoint) Stats() ImplicitStats { return s.work.stats }
func (s *MidPointFixPoint) NextStep(x, t float64) (float64, error) {
	return s.work.scalarStep(s.TdFunc, x, t, s.DeltaT, s.Options, implicitMidpoint, false)
}
func (s *MidPointFixPoint) NextStepOnGrid(x []float64, t float64) error {
	return s.work.step(s.TdFunc, x, t, s.DeltaT, s.Options, implicitMidpoint, false, false)
}

var _ ODESolver = (*MidPointFixPoint)(nil)

// PredictorCorrector implements euler predictor with converged trapezoidal corrector.
// Reuse one instance per integration; it owns mutable workspace.
type PredictorCorrector struct {
	TdFunc  gridData.TDPotentialOp
	DeltaT  float64
	Options ImplicitOptions
	work    implicitWorkspace
}

func (s *PredictorCorrector) Name() string {
	return "Euler predictor with converged trapezoidal corrector"
}
func (s *PredictorCorrector) NewDefine(dt float64, f gridData.TDPotentialOp) *PredictorCorrector {
	return &PredictorCorrector{TdFunc: f, DeltaT: dt}
}

// ReDefine preserves configured options and reuses workspace. No matrix is cached across steps.
func (s *PredictorCorrector) ReDefine(dt float64, f gridData.TDPotentialOp) {
	s.TdFunc = f
	s.DeltaT = dt
}
func (s *PredictorCorrector) Stats() ImplicitStats { return s.work.stats }
func (s *PredictorCorrector) NextStep(x, t float64) (float64, error) {
	return s.work.scalarStep(s.TdFunc, x, t, s.DeltaT, s.Options, trapezoidal, false)
}
func (s *PredictorCorrector) NextStepOnGrid(x []float64, t float64) error {
	return s.work.step(s.TdFunc, x, t, s.DeltaT, s.Options, trapezoidal, false, false)
}

var _ ODESolver = (*PredictorCorrector)(nil)
