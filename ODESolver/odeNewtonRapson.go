package EquationSolver

import "GoProject/gridData"

// EulerImplicitNewton implements backward euler with damped modified newton.
// Reuse one instance per integration; it owns mutable workspace.
type EulerImplicitNewton struct {
	TdFunc  gridData.TDPotentialOp
	DeltaT  float64
	Options ImplicitOptions
	work    implicitWorkspace
}

func (s *EulerImplicitNewton) Name() string { return "Backward Euler with damped modified Newton" }
func (s *EulerImplicitNewton) NewDefine(dt float64, f gridData.TDPotentialOp) *EulerImplicitNewton {
	return &EulerImplicitNewton{TdFunc: f, DeltaT: dt}
}

// ReDefine preserves configured options and reuses workspace. No matrix is cached across steps.
func (s *EulerImplicitNewton) ReDefine(dt float64, f gridData.TDPotentialOp) {
	s.TdFunc = f
	s.DeltaT = dt
}
func (s *EulerImplicitNewton) Stats() ImplicitStats { return s.work.stats }
func (s *EulerImplicitNewton) NextStep(x, t float64) (float64, error) {
	return s.work.scalarStep(s.TdFunc, x, t, s.DeltaT, s.Options, backwardEuler, true)
}
func (s *EulerImplicitNewton) NextStepOnGrid(x []float64, t float64) error {
	return s.work.step(s.TdFunc, x, t, s.DeltaT, s.Options, backwardEuler, true, false)
}

var _ ODESolver = (*EulerImplicitNewton)(nil)

// HeunsImplicitNewton implements implicit trapezoidal rule with damped modified newton.
// Reuse one instance per integration; it owns mutable workspace.
type HeunsImplicitNewton struct {
	TdFunc  gridData.TDPotentialOp
	DeltaT  float64
	Options ImplicitOptions
	work    implicitWorkspace
}

func (s *HeunsImplicitNewton) Name() string {
	return "Implicit trapezoidal rule with damped modified Newton"
}
func (s *HeunsImplicitNewton) NewDefine(dt float64, f gridData.TDPotentialOp) *HeunsImplicitNewton {
	return &HeunsImplicitNewton{TdFunc: f, DeltaT: dt}
}

// ReDefine preserves configured options and reuses workspace. No matrix is cached across steps.
func (s *HeunsImplicitNewton) ReDefine(dt float64, f gridData.TDPotentialOp) {
	s.TdFunc = f
	s.DeltaT = dt
}
func (s *HeunsImplicitNewton) Stats() ImplicitStats { return s.work.stats }
func (s *HeunsImplicitNewton) NextStep(x, t float64) (float64, error) {
	return s.work.scalarStep(s.TdFunc, x, t, s.DeltaT, s.Options, trapezoidal, true)
}
func (s *HeunsImplicitNewton) NextStepOnGrid(x []float64, t float64) error {
	return s.work.step(s.TdFunc, x, t, s.DeltaT, s.Options, trapezoidal, true, false)
}

var _ ODESolver = (*HeunsImplicitNewton)(nil)

// PredictIni writes the explicit Euler predictor into the caller's buffer.
func (s *HeunsImplicitNewton) PredictIni(x, f, predictor []float64) {
	for i := range predictor {
		predictor[i] = x[i] + s.DeltaT*f[i]
	}
}

// MidPointNewton implements implicit midpoint with damped modified newton.
// Reuse one instance per integration; it owns mutable workspace.
type MidPointNewton struct {
	TdFunc  gridData.TDPotentialOp
	DeltaT  float64
	Options ImplicitOptions
	work    implicitWorkspace
}

func (s *MidPointNewton) Name() string { return "Implicit midpoint with damped modified Newton" }
func (s *MidPointNewton) NewDefine(dt float64, f gridData.TDPotentialOp) *MidPointNewton {
	return &MidPointNewton{TdFunc: f, DeltaT: dt}
}

// ReDefine preserves configured options and reuses workspace. No matrix is cached across steps.
func (s *MidPointNewton) ReDefine(dt float64, f gridData.TDPotentialOp) { s.TdFunc = f; s.DeltaT = dt }
func (s *MidPointNewton) Stats() ImplicitStats                          { return s.work.stats }
func (s *MidPointNewton) NextStep(x, t float64) (float64, error) {
	return s.work.scalarStep(s.TdFunc, x, t, s.DeltaT, s.Options, implicitMidpoint, true)
}
func (s *MidPointNewton) NextStepOnGrid(x []float64, t float64) error {
	return s.work.step(s.TdFunc, x, t, s.DeltaT, s.Options, implicitMidpoint, true, false)
}

var _ ODESolver = (*MidPointNewton)(nil)

// NextStepIm preserves the legacy spelling; use NextStep in new code.
func (s *MidPointNewton) NextStepIm(x, t float64) (float64, error) { return s.NextStep(x, t) }

// SDIRK2 implements the two-stage, second-order L-stable SDIRK method
// with gamma=1-1/sqrt(2) and damped modified Newton stage solves.
// Reuse one instance per integration; it owns mutable workspace.
type SDIRK2 struct {
	TdFunc  gridData.TDPotentialOp
	DeltaT  float64
	Options ImplicitOptions
	work    implicitWorkspace
}

func (s *SDIRK2) Name() string { return "Second-order L-stable SDIRK with damped modified Newton" }
func (s *SDIRK2) NewDefine(dt float64, f gridData.TDPotentialOp) *SDIRK2 {
	return &SDIRK2{TdFunc: f, DeltaT: dt}
}

// ReDefine preserves configured options and reuses workspace. No matrix is cached across steps.
func (s *SDIRK2) ReDefine(dt float64, f gridData.TDPotentialOp) { s.TdFunc = f; s.DeltaT = dt }
func (s *SDIRK2) Stats() ImplicitStats                          { return s.work.stats }
func (s *SDIRK2) NextStep(x, t float64) (float64, error) {
	return s.work.scalarStep(s.TdFunc, x, t, s.DeltaT, s.Options, sdirk2, true)
}
func (s *SDIRK2) NextStepOnGrid(x []float64, t float64) error {
	return s.work.step(s.TdFunc, x, t, s.DeltaT, s.Options, sdirk2, true, false)
}

var _ ODESolver = (*SDIRK2)(nil)
