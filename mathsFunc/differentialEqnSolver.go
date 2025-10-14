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
	timeFunc  gridData.TDPotentialOp
	deltaTime float64
	fxt       blas64.Vector
}

func (eEx *EulerExplicit) Name() string {
	return "Euler Explicit Method"
}

func (eEx *EulerExplicit) NewDefine(dt float64, tdFunc gridData.TDPotentialOp) *EulerExplicit {
	return &EulerExplicit{
		timeFunc:  tdFunc,
		deltaTime: dt,
	}
}

func (eEx *EulerExplicit) ReDefine(dt float64, tdFunc gridData.TDPotentialOp) {
	eEx.timeFunc = tdFunc
	eEx.deltaTime = dt
}

func (eEx *EulerExplicit) NextStep(xt, t float64) (float64, error) {
	fxt := eEx.timeFunc.EvaluateAt(xt, t)
	if math.IsNaN(fxt) || math.IsInf(fxt, 0) {
		return xt, fmt.Errorf("the fxt is not valid")
	}
	return xt + eEx.deltaTime*fxt, nil
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

	eEx.timeFunc.EvaluateOnRGridInPlace(xt, eEx.fxt.Data, t)

	for _, v := range eEx.fxt.Data {
		if math.IsNaN(v) || math.IsInf(v, 0) {
			return fmt.Errorf("the fxt is not valid")
		}
	}

	blas64.Axpy(eEx.deltaTime, eEx.fxt, blas64.Vector{N: nPoints, Data: xt, Inc: 1})

	return nil
}

type HeunsExplicit struct {
	timeFunc  gridData.TDPotentialOp
	deltaTime float64

	halfDt    float64
	predictor blas64.Vector
	corrector blas64.Vector
	fxt       blas64.Vector
}

func (hEx *HeunsExplicit) Name() string {
	return "Heun's Explicit Method (Improved Euler method)"
}

func (hEx *HeunsExplicit) NewDefine(dt float64, tdFunc gridData.TDPotentialOp) *HeunsExplicit {
	return &HeunsExplicit{
		deltaTime: dt,
		timeFunc:  tdFunc,
		halfDt:    dt / 2,
	}
}

func (hEx *HeunsExplicit) ReDefine(dt float64, tdFunc gridData.TDPotentialOp) {
	hEx.deltaTime = dt
	hEx.timeFunc = tdFunc
	hEx.halfDt = dt / 2
}

func (hEx *HeunsExplicit) NextStep(xt, t float64) (float64, error) {
	fxt := hEx.timeFunc.EvaluateAt(xt, t)
	predictor := xt + fxt*hEx.deltaTime
	slope2 := fxt + hEx.timeFunc.EvaluateAt(predictor, t+hEx.deltaTime)

	if math.IsNaN(slope2) || math.IsInf(slope2, 0) {
		return xt, fmt.Errorf("invalid slope value")
	}
	return xt + hEx.halfDt*slope2, nil
}

// PredictIni computes predictor: xP = xt + dt*fxt
func (hEx *HeunsExplicit) PredictIni(xt []float64, fxt, xP blas64.Vector) {
	copy(xP.Data, xt)
	blas64.Axpy(hEx.deltaTime, fxt, xP)
}

func (hEx *HeunsExplicit) allocate(nPoints int) {
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

	hEx.timeFunc.EvaluateOnRGridInPlace(xt, hEx.fxt.Data, t)
	hEx.PredictIni(xt, hEx.fxt, hEx.predictor)
	hEx.timeFunc.EvaluateOnRGridInPlace(hEx.predictor.Data, hEx.corrector.Data, t+hEx.deltaTime)

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

type MidPointExplicit struct {
	timeFunc gridData.TDPotentialOp
	delTime  float64

	halfDt float64
	xtMid  blas64.Vector
	fxt    blas64.Vector
}

func (mEx *MidPointExplicit) Name() string {
	return "MidPoint Explicit Method!"
}

func (mEx *MidPointExplicit) NewDefine(dt float64, tdFunc gridData.TDPotentialOp) *MidPointExplicit {
	return &MidPointExplicit{
		timeFunc: tdFunc,
		delTime:  dt,
		halfDt:   dt / 2,
	}
}

func (mEx *MidPointExplicit) ReDefine(dt float64, tdFunc gridData.TDPotentialOp) {
	mEx.delTime = dt
	mEx.timeFunc = tdFunc
	mEx.halfDt = dt / 2
}

func (mEx *MidPointExplicit) NextStep(xt, t float64) (float64, error) {
	xtMid := xt + mEx.timeFunc.EvaluateAt(xt, t)*mEx.halfDt
	slope := mEx.timeFunc.EvaluateAt(xtMid, t+mEx.halfDt)

	if math.IsNaN(slope) || math.IsInf(slope, 0) {
		return xt, fmt.Errorf("the fxt is not valid")
	}
	return xt + mEx.delTime*slope, nil
}

func (mEx *MidPointExplicit) NextStepOnGrid(xt []float64, t float64) error {
	nPoints := len(xt)

	// Allocate buffers only if necessary
	if nPoints != mEx.fxt.N {
		mEx.fxt = blas64.Vector{N: nPoints, Data: make([]float64, nPoints), Inc: 1}
		mEx.xtMid = blas64.Vector{N: nPoints, Data: make([]float64, nPoints), Inc: 1}
	}

	// Compute xtMid = xt + halfDt * f(xt, t)
	mEx.timeFunc.EvaluateOnRGridInPlace(xt, mEx.fxt.Data, t)
	copy(mEx.xtMid.Data, xt)
	blas64.Axpy(mEx.halfDt, mEx.fxt, mEx.xtMid)

	// Compute f(xtMid, t + halfDt)
	mEx.timeFunc.EvaluateOnRGridInPlace(mEx.xtMid.Data, mEx.fxt.Data, t+mEx.halfDt)

	for i := 0; i < nPoints; i++ {
		if math.IsNaN(mEx.fxt.Data[i]) || math.IsInf(mEx.fxt.Data[i], 0) {
			return fmt.Errorf("the fxt is not valid")
		}
	}

	blas64.Axpy(mEx.delTime, mEx.fxt, blas64.Vector{N: nPoints, Data: xt, Inc: 1})

	return nil
}

// RungeKutta4Explicit implements the classic 4th-order Runge-Kutta method
type RungeKutta4Explicit struct {
	timeFunc  gridData.TDPotentialOp
	deltaTime float64

	// Pre-allocated buffers for grid operations
	halfDt  float64
	dtBy6   float64
	xPlusDt blas64.Vector
	k1      blas64.Vector
	k2      blas64.Vector
	k3      blas64.Vector
	k4      blas64.Vector
}

func (rgEx *RungeKutta4Explicit) Name() string {
	return "Runge-Kutta Order 4"
}

func (rgEx *RungeKutta4Explicit) NewDef(dt float64, tdFunc gridData.TDPotentialOp) *RungeKutta4Explicit {
	return &RungeKutta4Explicit{
		timeFunc:  tdFunc,
		deltaTime: dt,
		halfDt:    dt / 2,
		dtBy6:     dt / 6,
	}
}

func (rgEx *RungeKutta4Explicit) ReDefine(dt float64, tdFunc gridData.TDPotentialOp) {
	rgEx.timeFunc = tdFunc
	rgEx.deltaTime = dt
	rgEx.dtBy6 = dt / 6
	rgEx.halfDt = dt / 2
}

func (rgEx *RungeKutta4Explicit) NextStep(xt, t float64) (float64, error) {
	k1 := rgEx.timeFunc.EvaluateAt(xt, t)
	k2 := rgEx.timeFunc.EvaluateAt(xt+k1*rgEx.halfDt, t+rgEx.halfDt)
	k3 := rgEx.timeFunc.EvaluateAt(xt+k2*rgEx.halfDt, t+rgEx.halfDt)
	k4 := rgEx.timeFunc.EvaluateAt(xt+k3*rgEx.deltaTime, t+rgEx.deltaTime)

	slope := (k1 + 2*(k2+k3) + k4) / 6.0

	if math.IsNaN(slope) || math.IsInf(slope, 0) {
		return xt, fmt.Errorf("invalid derivative: NaN or Inf encountered")
	}

	return xt + rgEx.deltaTime*slope, nil
}

func (rgEx *RungeKutta4Explicit) allocate(nPoints int) {
	// Only reallocate if size has changed or vectors are uninitialized
	if rgEx.k1.N != nPoints ||
		rgEx.k2.N != nPoints ||
		rgEx.k3.N != nPoints ||
		rgEx.k4.N != nPoints ||
		rgEx.xPlusDt.N != nPoints {

		rgEx.xPlusDt = blas64.Vector{N: nPoints, Inc: 1, Data: make([]float64, nPoints)}
		rgEx.k1 = blas64.Vector{N: nPoints, Inc: 1, Data: make([]float64, nPoints)}
		rgEx.k2 = blas64.Vector{N: nPoints, Inc: 1, Data: make([]float64, nPoints)}
		rgEx.k3 = blas64.Vector{N: nPoints, Inc: 1, Data: make([]float64, nPoints)}
		rgEx.k4 = blas64.Vector{N: nPoints, Inc: 1, Data: make([]float64, nPoints)}
	}
}

func (rgEx *RungeKutta4Explicit) NextStepOnGrid(xt []float64, t float64) error {
	nPoints := len(xt)

	rgEx.allocate(nPoints)

	// k1 = f(x, t)
	rgEx.timeFunc.EvaluateOnRGridInPlace(xt, rgEx.k1.Data, t)

	// k2 = f(x + dt/2 * k1, t + dt/2)
	copy(rgEx.xPlusDt.Data, xt)
	blas64.Axpy(rgEx.halfDt, rgEx.k1, rgEx.xPlusDt)
	rgEx.timeFunc.EvaluateOnRGridInPlace(rgEx.xPlusDt.Data, rgEx.k2.Data, t+rgEx.halfDt)

	// k3 = f(x + dt/2 * k2, t + dt/2)
	copy(rgEx.xPlusDt.Data, xt)
	blas64.Axpy(rgEx.halfDt, rgEx.k2, rgEx.xPlusDt)
	rgEx.timeFunc.EvaluateOnRGridInPlace(rgEx.xPlusDt.Data, rgEx.k3.Data, t+rgEx.halfDt)

	// k4 = f(x + dt * k3, t + dt)
	copy(rgEx.xPlusDt.Data, xt)
	blas64.Axpy(rgEx.deltaTime, rgEx.k3, rgEx.xPlusDt)
	rgEx.timeFunc.EvaluateOnRGridInPlace(rgEx.xPlusDt.Data, rgEx.k4.Data, t+rgEx.deltaTime)

	// x_new = x + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
	for i := 0; i < nPoints; i++ {
		xt[i] += rgEx.dtBy6 * (rgEx.k1.Data[i] + 2*(rgEx.k2.Data[i]+rgEx.k3.Data[i]) + rgEx.k4.Data[i])
	}

	return nil
}
