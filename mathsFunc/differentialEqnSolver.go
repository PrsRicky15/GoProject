package mathsFunc

import (
	"GoProject/gridData"
	"math"
	"slices"

	"gonum.org/v1/gonum/blas/blas64"
)

const (
	maxIter   = 20
	tolerance = 1e-7
)

// ODESolverExplicit interface for different solving methods
// dx(t)/dt = f(x,t)
type ODESolverExplicit interface {
	NextStepEx(funcAtT, x, time float64) float64
	NextStepExOnGrid(funcAtT, x []float64, time float64)
	Name() string
}

// ODESolverImplicit interface for different solving methods
// dx(t)/dt = f(x,t)
type ODESolverImplicit interface {
	NextStepIm(funcAtT, x, time float64) float64
	NextStepImOnGrid(funcAtT, x []float64, time float64)
	Name() string
}

// EulerSolver implements the basic Euler method
type EulerSolver struct {
	TdFunc gridData.TDPotentialOp
	DeltaT float64
	slope  []float64
}

func (e *EulerSolver) Name() string {
	return "Euler Method"
}

func (e *EulerSolver) NextStep(xt, t float64) float64 {
	slope := (e.TdFunc).EvaluateAt(xt, t)
	return xt + e.DeltaT*slope
}

func (e *EulerSolver) NextStepOnGrid(xt []float64, t float64) {
	nPoints := len(xt)
	if nPoints != len(e.slope) {
		e.slope = make([]float64, len(xt))
	}

	(e.TdFunc).EvaluateOnRGridInPlace(xt, e.slope, t)
	blas64.Axpy(e.DeltaT,
		blas64.Vector{N: nPoints, Data: e.slope, Inc: 1},
		blas64.Vector{N: nPoints, Data: xt, Inc: 1},
	)
}

// HeunsSolver implements the basic Euler method
type HeunsSolver struct {
	TdFunc    gridData.TDPotentialOp
	DeltaT    float64
	predictor []float64
	slope     []float64
}

func (h *HeunsSolver) Name() string {
	return "Heun's Method (Improved Euler method)"
}

func (h *HeunsSolver) NextStep(xt, t float64) float64 {
	fxt := (h.TdFunc).EvaluateAt(xt, t)
	predictor := xt + fxt*h.DeltaT
	slope2 := fxt + (h.TdFunc).EvaluateAt(predictor, t+h.DeltaT)
	return xt + 0.5*h.DeltaT*slope2
}

func (h *HeunsSolver) NextStepOnGrid(xt []float64, t float64) {

	nPoints := len(xt)
	if nPoints != len(h.slope) || nPoints != len(h.predictor) {
		h.slope = make([]float64, len(xt))
		h.predictor = make([]float64, len(xt))
	}

	(h.TdFunc).EvaluateOnRGridInPlace(xt, h.slope, t)
	h.predictor = slices.Clone(xt)

	blas64.Axpy(h.DeltaT,
		blas64.Vector{N: nPoints, Data: h.slope, Inc: 1},
		blas64.Vector{N: nPoints, Data: h.predictor, Inc: 1},
	)

	for i := range h.predictor {
		h.slope[i] += (h.TdFunc).EvaluateAt(h.predictor[i], t+h.DeltaT)
	}

	blas64.Axpy(0.5*h.DeltaT,
		blas64.Vector{N: nPoints, Data: h.slope, Inc: 1},
		blas64.Vector{N: nPoints, Data: xt, Inc: 1},
	)
}

type MidPointSolver struct {
	TdFunc gridData.TDPotentialOp
	DeltaT float64
	xtMid  []float64
	slope  []float64
}

func (m *MidPointSolver) Name() string {
	return "MidPoint Method (Implicit or Explicit)"
}
func (m *MidPointSolver) NextStep(xt, t float64) float64 {
	xtMid := xt + (m.TdFunc).EvaluateAt(xt, t)*m.DeltaT/2
	slope := (m.TdFunc).EvaluateAt(xtMid, t+m.DeltaT/2)
	return xt + m.DeltaT*slope
}

func (m *MidPointSolver) NextStepOnGrid(xt []float64, t float64) {
	nPoints := len(xt)

	(m.TdFunc).EvaluateOnRGridInPlace(xt, m.slope, t)

	m.xtMid = slices.Clone(xt)
	blas64.Axpy(m.DeltaT/2,
		blas64.Vector{N: nPoints, Data: m.slope, Inc: 1},
		blas64.Vector{N: nPoints, Data: m.xtMid, Inc: 1},
	)

	for i := range m.xtMid {
		m.slope[i] = (m.TdFunc).EvaluateOnRGrid(m.xtMid, t+m.DeltaT/2)[i]
	}

	blas64.Axpy(m.DeltaT,
		blas64.Vector{N: nPoints, Data: m.slope, Inc: 1},
		blas64.Vector{N: nPoints, Data: xt, Inc: 1},
	)
}

func (m *MidPointSolver) NextStepImplicit(xt, t float64) float64 {
	xtMid := xt + (m.TdFunc).EvaluateAt(xt, t)*m.DeltaT/2
	return xt + m.DeltaT*(m.TdFunc).EvaluateAt(xtMid, t+m.DeltaT/2)
}

func (m *MidPointSolver) NextStepOneGridImplicit(xt, t float64) float64 {
	return xt + t
}

// RungeKutta4 implements the basic Euler method
type RungeKutta4 struct {
	TdFunc gridData.TDPotentialOp
	DeltaT float64
}

func (rg *RungeKutta4) Name() string {
	return "Runge Kutta order 4"
}

func (rg *RungeKutta4) NextStep(xt, t float64) float64 {
	halfDt := 0.5 * rg.DeltaT
	k1 := (rg.TdFunc).EvaluateAt(xt, t)
	k2 := (rg.TdFunc).EvaluateAt(xt+k1*halfDt, t+halfDt)
	k3 := (rg.TdFunc).EvaluateAt(xt+k2*halfDt, t+halfDt)
	k4 := (rg.TdFunc).EvaluateAt(xt+k3*rg.DeltaT, t+rg.DeltaT)
	return xt + rg.DeltaT/6.0*(k1+2*k2+2*k3+k4)
}

func (rg *RungeKutta4) NextStepOnGrid(xt []float64, t float64) {
	halfDt := 0.5 * rg.DeltaT
	k1 := (rg.TdFunc).EvaluateOnRGrid(xt, t)

	xtPlusDt := slices.Clone(xt)
	blas64.Axpy(0.5*rg.DeltaT,
		blas64.Vector{N: len(k1), Data: k1, Inc: 1},
		blas64.Vector{N: len(xtPlusDt), Data: xtPlusDt, Inc: 1},
	)
	k2 := (rg.TdFunc).EvaluateOnRGrid(xtPlusDt, t+halfDt)
	xtPlusDt = slices.Clone(xt)
	blas64.Axpy(0.5*rg.DeltaT,
		blas64.Vector{N: len(k2), Data: k2, Inc: 1},
		blas64.Vector{N: len(xtPlusDt), Data: xtPlusDt, Inc: 1},
	)

	k3 := (rg.TdFunc).EvaluateOnRGrid(xtPlusDt, t+halfDt)

	xtPlusDt = slices.Clone(xt)
	blas64.Axpy(rg.DeltaT,
		blas64.Vector{N: len(k3), Data: k3, Inc: 1},
		blas64.Vector{N: len(xtPlusDt), Data: xtPlusDt, Inc: 1},
	)

	k4 := (rg.TdFunc).EvaluateOnRGrid(xtPlusDt, t+rg.DeltaT)

	Dtby6 := rg.DeltaT / 6
	for i := range xtPlusDt {
		xt[i] += Dtby6 * (k1[i] + 2*(k2[i]+k3[i]) + k4[i])
	}
}

// HeunsIterSolver implements the Heuns iterative method
type HeunsIterSolver struct {
	TdFunc    gridData.TDPotentialOp
	DeltaT    float64
	fatXt     []float64
	predictor []float64
	corrector []float64
	diff      []float64
}

func (hi *HeunsIterSolver) Name() string {
	return "Heun's Iterative Method"
}

func (hi *HeunsIterSolver) NextStep(xt, t float64) float64 {
	predictor := xt + hi.DeltaT*(hi.TdFunc).EvaluateAt(xt, t)
	corrector := predictor

	for i := 0; i < maxIter; i++ {
		correctorNew := xt + 0.5*hi.DeltaT*((hi.TdFunc).EvaluateAt(xt, t)+
			(hi.TdFunc).EvaluateAt(corrector, t+hi.DeltaT))

		var err float64
		if math.Abs(correctorNew) > 1e-10 {
			err = math.Abs(correctorNew-corrector) / math.Abs(correctorNew)
		} else {
			err = math.Abs(correctorNew - corrector)
		}

		if err < tolerance {
			return correctorNew
		}

		corrector = correctorNew
	}
	return corrector
}

func (hi *HeunsIterSolver) allocate(nPoints int) {

	if nPoints != len(hi.corrector) ||
		nPoints != len(hi.predictor) ||
		nPoints != len(hi.fatXt) ||
		nPoints != len(hi.diff) {

		hi.corrector = make([]float64, nPoints)
		hi.predictor = make([]float64, nPoints)
		hi.fatXt = make([]float64, nPoints)
		hi.diff = make([]float64, nPoints)
	}

}

func (hi *HeunsIterSolver) iterate(xt []float64, t float64) float64 {
	nPoints := len(xt)

	for i := range hi.predictor {
		hi.corrector[i] = hi.fatXt[i] + (hi.TdFunc).EvaluateAt(hi.predictor[i], t+hi.DeltaT)
	}

	hi.diff = slices.Clone(xt)
	blas64.Axpy(0.5*hi.DeltaT,
		blas64.Vector{N: nPoints, Data: hi.corrector, Inc: 1},
		blas64.Vector{N: nPoints, Data: hi.diff, Inc: 1},
	)

	hi.corrector = slices.Clone(hi.diff)

	for i := range hi.corrector {
		hi.diff[i] -= hi.predictor[i]
	}

	hi.predictor = slices.Clone(hi.corrector)

	return blas64.Nrm2(blas64.Vector{N: nPoints, Data: hi.diff, Inc: 1})
}

func (hi *HeunsIterSolver) NextStepOnGrid(xt []float64, t float64) {

	nPoints := len(xt)
	hi.allocate(nPoints)

	(hi.TdFunc).EvaluateOnRGridInPlace(xt, hi.fatXt, t)
	hi.predictor = slices.Clone(xt)

	blas64.Axpy(hi.DeltaT,
		blas64.Vector{N: nPoints, Data: hi.fatXt, Inc: 1},
		blas64.Vector{N: nPoints, Data: hi.predictor, Inc: 1},
	)

	var err float64
	for i := 0; i < maxIter; i++ {
		err = hi.iterate(xt, t)
		if err < tolerance {
			break
		}
	}
	xt = slices.Clone(hi.corrector)
}

// PredictorCorrector implements the basic Euler method
type PredictorCorrector struct {
	TdFunc gridData.TDPotentialOp
	DeltaT float64
}

type MDodeSolver interface {
	NextStep(forceAtT, x, v, time float64) float64
	NextStepOnGrid(forceAtT, x, v []float64, time float64)
	Name() string
}

type Verlet struct {
	ForceFunc gridData.TDPotentialOp
	Dt        float64
}

func (vi *Verlet) Name() string {
	return "Verlet Integrators"
}

func (vi *Verlet) NextStep(xt, t float64) float64 {
	return xt + t
}

type StormerVerlet struct {
	ForceFunc gridData.TDPotentialOp
	Dt        float64
}

type VelocityVerlet struct {
	ForceFunc gridData.TDPotentialOp
	Dt        float64
}

type LeapFrog struct {
	ForceFunc gridData.TDPotentialOp
	Dt        float64
}

type Beeman struct {
	ForceFunc gridData.TDPotentialOp
	Dt        float64
}
