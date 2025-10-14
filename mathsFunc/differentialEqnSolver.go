package mathsFunc

import (
	"GoProject/gridData"
	"fmt"
	"math"

	"gonum.org/v1/gonum/blas/blas64"
)

const (
	maxIter   = 20
	tolerance = 1e-7
	delta     = 1e-6
)

// ODESolver interface for different solving methods
// dx(t)/dt = f(x,t)
type ODESolver interface {
	NextStep(x, time float64) (float64, error)
	NextStepOnGrid(x []float64, time float64) error
}

type EulerExplicit struct {
	TdFunc gridData.TDPotentialOp
	DeltaT float64
	fxt    blas64.Vector
}

func (eEx *EulerExplicit) Name() string {
	return "Euler Explicit Method"
}

func (eEx *EulerExplicit) NextStep(xt, t float64) (float64, error) {
	fxt := eEx.TdFunc.EvaluateAt(xt, t)
	if math.IsNaN(fxt) || math.IsInf(fxt, 0) {
		return xt, fmt.Errorf("the fxt is not valid")
	}
	return xt + eEx.DeltaT*fxt, nil
}

func (eEx *EulerExplicit) NextStepOnGrid(xt []float64, t float64) error {
	nPoints := len(xt)

	if eEx.fxt.N != nPoints {
		if cap(eEx.fxt.Data) >= nPoints {
			eEx.fxt.Data = eEx.fxt.Data[:nPoints]
		} else {
			eEx.fxt.Data = make([]float64, nPoints)
		}
		eEx.fxt.N = nPoints
		eEx.fxt.Inc = 1
	}

	eEx.TdFunc.EvaluateOnRGridInPlace(xt, eEx.fxt.Data, t)

	for _, v := range eEx.fxt.Data {
		if math.IsNaN(v) || math.IsInf(v, 0) {
			return fmt.Errorf("the fxt is not valid")
		}
	}

	blas64.Axpy(eEx.DeltaT, eEx.fxt, blas64.Vector{N: nPoints, Data: xt, Inc: 1})

	return nil
}

type HeunsExplicit struct {
	TdFunc    gridData.TDPotentialOp
	DeltaT    float64
	predictor blas64.Vector
	corrector blas64.Vector
	fxt       blas64.Vector
	halfDt    float64 // Pre-computed constant
}

func (hEx *HeunsExplicit) Name() string {
	return "Heun's Explicit Method (Improved Euler method)"
}

func (hEx *HeunsExplicit) NextStep(xt, t float64) (float64, error) {
	fxt := hEx.TdFunc.EvaluateAt(xt, t)
	predictor := xt + fxt*hEx.DeltaT
	slope2 := fxt + hEx.TdFunc.EvaluateAt(predictor, t+hEx.DeltaT)

	if math.IsNaN(slope2) || math.IsInf(slope2, 0) {
		return xt, fmt.Errorf("invalid slope value")
	}
	return xt + hEx.halfDt*slope2, nil
}

// PredictIni computes predictor: xP = xt + dt*fxt
func (hEx *HeunsExplicit) PredictIni(xt []float64, fxt, xP blas64.Vector) {
	copy(xP.Data, xt)
	blas64.Axpy(hEx.DeltaT, fxt, xP)
}

func (hEx *HeunsExplicit) allocate(nPoints int) {
	hEx.halfDt = hEx.DeltaT / 2
	if hEx.fxt.N != nPoints ||
		hEx.predictor.N != nPoints ||
		hEx.corrector.N != nPoints {
		hEx.fxt = blas64.Vector{N: nPoints, Inc: 1, Data: make([]float64, nPoints)}
		hEx.predictor = blas64.Vector{N: nPoints, Inc: 1, Data: make([]float64, nPoints)}
		hEx.corrector = blas64.Vector{N: nPoints, Inc: 1, Data: make([]float64, nPoints)}
	}
}

func (hEx *HeunsExplicit) NextStepOnGrid(xt []float64, t float64) error {
	nPoints := len(xt)
	hEx.allocate(nPoints)

	hEx.TdFunc.EvaluateOnRGridInPlace(xt, hEx.fxt.Data, t)
	hEx.PredictIni(xt, hEx.fxt, hEx.predictor)
	hEx.TdFunc.EvaluateOnRGridInPlace(hEx.predictor.Data, hEx.corrector.Data, t+hEx.DeltaT)

	for i := 0; i < nPoints; i++ {
		hEx.corrector.Data[i] += hEx.fxt.Data[i]
	}

	for _, v := range hEx.corrector.Data {
		if math.IsNaN(v) || math.IsInf(v, 0) {
			return fmt.Errorf("invalid corrector value")
		}
	}

	xNew := blas64.Vector{N: nPoints, Data: xt, Inc: 1}
	blas64.Axpy(hEx.halfDt, hEx.corrector, xNew)
	return nil
}

// SetDeltaT updates DeltaT and pre-computes halfDt
func (hEx *HeunsExplicit) SetDeltaT(dt float64) {
	hEx.DeltaT = dt
	hEx.halfDt = 0.5 * dt
}

type MidPointExplicit struct {
	TdFunc gridData.TDPotentialOp
	DeltaT float64
	xtMid  blas64.Vector
	fxt    blas64.Vector
}

func (mEx *MidPointExplicit) Name() string {
	return "MidPoint Explicit Method!"
}

func (mEx *MidPointExplicit) NextStep(xt, t float64) (float64, error) {
	halfDt := mEx.DeltaT * 0.5
	xtMid := xt + mEx.TdFunc.EvaluateAt(xt, t)*halfDt
	slope := mEx.TdFunc.EvaluateAt(xtMid, t+halfDt)

	if math.IsNaN(slope) || math.IsInf(slope, 0) {
		return xt, fmt.Errorf("the fxt is not valid")
	}
	return xt + mEx.DeltaT*slope, nil
}

func (mEx *MidPointExplicit) NextStepOnGrid(xt []float64, t float64) error {
	nPoints := len(xt)

	// Allocate buffers only if necessary
	if nPoints != mEx.fxt.N {
		mEx.fxt = blas64.Vector{N: nPoints, Data: make([]float64, nPoints), Inc: 1}
		mEx.xtMid = blas64.Vector{N: nPoints, Data: make([]float64, nPoints), Inc: 1}
	}

	halfDt := mEx.DeltaT * 0.5

	// Compute xtMid = xt + halfDt * f(xt, t)
	mEx.TdFunc.EvaluateOnRGridInPlace(xt, mEx.fxt.Data, t)
	copy(mEx.xtMid.Data, xt)
	blas64.Axpy(halfDt, mEx.fxt, mEx.xtMid)

	// Compute f(xtMid, t + halfDt)
	mEx.TdFunc.EvaluateOnRGridInPlace(mEx.xtMid.Data, mEx.fxt.Data, t+halfDt)

	// Validate result (single pass check)
	for i := 0; i < nPoints; i++ {
		if math.IsNaN(mEx.fxt.Data[i]) || math.IsInf(mEx.fxt.Data[i], 0) {
			return fmt.Errorf("the fxt is not valid")
		}
	}

	blas64.Axpy(mEx.DeltaT, mEx.fxt, blas64.Vector{N: nPoints, Data: xt, Inc: 1})

	return nil
}

// RungeKutta4Explicit implements the classic 4th-order Runge-Kutta method
type RungeKutta4Explicit struct {
	TdFunc gridData.TDPotentialOp
	DeltaT float64

	// Pre-allocated buffers for grid operations
	xPlusDt blas64.Vector
	k1      blas64.Vector
	k2      blas64.Vector
	k3      blas64.Vector
	k4      blas64.Vector
}

func (rgEx *RungeKutta4Explicit) Name() string {
	return "Runge-Kutta Order 4"
}

func (rgEx *RungeKutta4Explicit) NextStep(xt, t float64) (float64, error) {
	halfDt := 0.5 * rgEx.DeltaT

	k1 := rgEx.TdFunc.EvaluateAt(xt, t)
	k2 := rgEx.TdFunc.EvaluateAt(xt+k1*halfDt, t+halfDt)
	k3 := rgEx.TdFunc.EvaluateAt(xt+k2*halfDt, t+halfDt)
	k4 := rgEx.TdFunc.EvaluateAt(xt+k3*rgEx.DeltaT, t+rgEx.DeltaT)

	slope := (k1 + 2*(k2+k3) + k4) / 6.0

	if math.IsNaN(slope) || math.IsInf(slope, 0) {
		return xt, fmt.Errorf("invalid derivative: NaN or Inf encountered")
	}

	return xt + rgEx.DeltaT*slope, nil
}

func (rgEx *RungeKutta4Explicit) allocate(nPoints int) {
	// Only reallocate if size has changed
	if rgEx.k1.N == nPoints {
		return
	}

	rgEx.xPlusDt = blas64.Vector{N: nPoints, Inc: 1, Data: make([]float64, nPoints)}
	rgEx.k1 = blas64.Vector{N: nPoints, Inc: 1, Data: make([]float64, nPoints)}
	rgEx.k2 = blas64.Vector{N: nPoints, Inc: 1, Data: make([]float64, nPoints)}
	rgEx.k3 = blas64.Vector{N: nPoints, Inc: 1, Data: make([]float64, nPoints)}
	rgEx.k4 = blas64.Vector{N: nPoints, Inc: 1, Data: make([]float64, nPoints)}
}

func (rgEx *RungeKutta4Explicit) NextStepOnGrid(xt []float64, t float64) error {
	halfDt := 0.5 * rgEx.DeltaT
	nPoints := len(xt)

	rgEx.allocate(nPoints)

	// k1 = f(x, t)
	rgEx.TdFunc.EvaluateOnRGridInPlace(xt, rgEx.k1.Data, t)

	// k2 = f(x + dt/2 * k1, t + dt/2)
	copy(rgEx.xPlusDt.Data, xt)
	blas64.Axpy(halfDt, rgEx.k1, rgEx.xPlusDt)
	rgEx.TdFunc.EvaluateOnRGridInPlace(rgEx.xPlusDt.Data, rgEx.k2.Data, t+halfDt)

	// k3 = f(x + dt/2 * k2, t + dt/2)
	copy(rgEx.xPlusDt.Data, xt)
	blas64.Axpy(halfDt, rgEx.k2, rgEx.xPlusDt)
	rgEx.TdFunc.EvaluateOnRGridInPlace(rgEx.xPlusDt.Data, rgEx.k3.Data, t+halfDt)

	// k4 = f(x + dt * k3, t + dt)
	copy(rgEx.xPlusDt.Data, xt)
	blas64.Axpy(rgEx.DeltaT, rgEx.k3, rgEx.xPlusDt)
	rgEx.TdFunc.EvaluateOnRGridInPlace(rgEx.xPlusDt.Data, rgEx.k4.Data, t+rgEx.DeltaT)

	// x_new = x + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
	dtBy6 := rgEx.DeltaT / 6.0
	for i := 0; i < nPoints; i++ {
		xt[i] += dtBy6 * (rgEx.k1.Data[i] + 2*(rgEx.k2.Data[i]+rgEx.k3.Data[i]) + rgEx.k4.Data[i])
	}

	return nil
}
