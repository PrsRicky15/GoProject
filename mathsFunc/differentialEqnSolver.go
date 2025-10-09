package mathsFunc

import (
	"GoProject/gridData"
	"slices"

	"gonum.org/v1/gonum/blas/blas64"
)

// ODESolver interface for different solving methods
// dx(t)/dt = f(x,t)
type ODESolver interface {
	NextStep(funcAtT, x, time float64) float64
	NextStepOnGrid(funcAtT, x []float64, time float64)
	Name() string
}

// EulerSolver implements the basic Euler method
type EulerSolver struct {
	tdFunc *gridData.TDPotentialOp
	DeltaT float64
}

func (e *EulerSolver) Name() string {
	return "Euler Method"
}

func (e *EulerSolver) NextStep(xt, t float64) float64 {
	return xt + (*e.tdFunc).EvaluateAt(xt, t)*e.DeltaT
}

func (e *EulerSolver) NextStepOnGrid(xt []float64, t float64) {
	result := (*e.tdFunc).EvaluateOnRGrid(xt, t)
	blas64.Axpy(e.DeltaT,
		blas64.Vector{N: len(result), Data: result, Inc: 1},
		blas64.Vector{N: len(xt), Data: xt, Inc: 1},
	)
}

// HeunsSolver implements the basic Euler method
type HeunsSolver struct {
	tdFunc *gridData.TDPotentialOp
	DeltaT float64
}

func (h *HeunsSolver) Name() string {
	return "Heun's Method (Improved Euler method)"
}

func (h *HeunsSolver) NextStep(xt, t float64) float64 {
	k1 := (*h.tdFunc).EvaluateAt(xt, t)
	xtPlusDt := xt + k1*h.DeltaT
	return xt + 0.5*h.DeltaT*(k1+(*h.tdFunc).EvaluateAt(xtPlusDt, t+h.DeltaT))
}

func (h *HeunsSolver) NextStepOnGrid(xt []float64, t float64) {
	K1 := (*h.tdFunc).EvaluateOnRGrid(xt, t)
	xtPlusDt := slices.Clone(xt)

	blas64.Axpy(h.DeltaT,
		blas64.Vector{N: len(K1), Data: K1, Inc: 1},
		blas64.Vector{N: len(xtPlusDt), Data: xtPlusDt, Inc: 1},
	)

	for i := range xtPlusDt {
		K1[i] += (*h.tdFunc).EvaluateOnRGrid(xtPlusDt, t+h.DeltaT)[i]
	}

	blas64.Axpy(0.5*h.DeltaT,
		blas64.Vector{N: len(K1), Data: K1, Inc: 1},
		blas64.Vector{N: len(xt), Data: xt, Inc: 1},
	)
}

// RungeKutta4 implements the basic Euler method
type RungeKutta4 struct {
	tdFunc *gridData.TDPotentialOp
	DeltaT float64
}

func (rg *RungeKutta4) Name() string {
	return "Runge Kutta order 4"
}

func (rg *RungeKutta4) NextStep(xt, t float64) float64 {
	halfDt := 0.5 * rg.DeltaT
	k1 := (*rg.tdFunc).EvaluateAt(xt, t)
	k2 := (*rg.tdFunc).EvaluateAt(xt+k1*halfDt, t+halfDt)
	k3 := (*rg.tdFunc).EvaluateAt(xt+k2*halfDt, t+halfDt)
	k4 := (*rg.tdFunc).EvaluateAt(xt+k3*rg.DeltaT, t+rg.DeltaT)
	return xt + rg.DeltaT/6.0*(k1+2*k2+2*k3+k4)
}

func (rg *RungeKutta4) NextStepOnGrid(xt []float64, t float64) {
	halfDt := 0.5 * rg.DeltaT
	k1 := (*rg.tdFunc).EvaluateOnRGrid(xt, t)

	xtPlusDt := slices.Clone(xt)
	blas64.Axpy(0.5*rg.DeltaT,
		blas64.Vector{N: len(k1), Data: k1, Inc: 1},
		blas64.Vector{N: len(xtPlusDt), Data: xtPlusDt, Inc: 1},
	)
	k2 := (*rg.tdFunc).EvaluateOnRGrid(xtPlusDt, t+halfDt)
	xtPlusDt = slices.Clone(xt)
	blas64.Axpy(0.5*rg.DeltaT,
		blas64.Vector{N: len(k2), Data: k2, Inc: 1},
		blas64.Vector{N: len(xtPlusDt), Data: xtPlusDt, Inc: 1},
	)

	k3 := (*rg.tdFunc).EvaluateOnRGrid(xtPlusDt, t+halfDt)

	xtPlusDt = slices.Clone(xt)
	blas64.Axpy(rg.DeltaT,
		blas64.Vector{N: len(k3), Data: k3, Inc: 1},
		blas64.Vector{N: len(xtPlusDt), Data: xtPlusDt, Inc: 1},
	)

	k4 := (*rg.tdFunc).EvaluateOnRGrid(xtPlusDt, t+rg.DeltaT)

	Dtby6 := rg.DeltaT / 6
	for i := range xtPlusDt {
		xt[i] += Dtby6 * (k1[i] + 2*(k2[i]+k3[i]) + k4[i])
	}
}
