package EquationSolver

import (
	"GoProject/gridData"
	"fmt"
	"math"

	"gonum.org/v1/gonum/blas/blas64"
)

const (
	maxIter           = 20
	tolerance         = 1e-7
	delta             = 1e-6
	adaptiveTolerance = 1e-8
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
	fxt := eEx.timeFunc.EvaluateAtTime(xt, t)
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

	eEx.timeFunc.EvaluateOnRGridTimeInPlace(xt, eEx.fxt.Data, t)

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
	fxt := hEx.timeFunc.EvaluateAtTime(xt, t)
	predictor := xt + fxt*hEx.deltaTime
	slope2 := fxt + hEx.timeFunc.EvaluateAtTime(predictor, t+hEx.deltaTime)

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

	hEx.timeFunc.EvaluateOnRGridTimeInPlace(xt, hEx.fxt.Data, t)
	hEx.PredictIni(xt, hEx.fxt, hEx.predictor)
	hEx.timeFunc.EvaluateOnRGridTimeInPlace(hEx.predictor.Data, hEx.corrector.Data, t+hEx.deltaTime)

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
	xtMid := xt + mEx.timeFunc.EvaluateAtTime(xt, t)*mEx.halfDt
	slope := mEx.timeFunc.EvaluateAtTime(xtMid, t+mEx.halfDt)

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
	mEx.timeFunc.EvaluateOnRGridTimeInPlace(xt, mEx.fxt.Data, t)
	copy(mEx.xtMid.Data, xt)
	blas64.Axpy(mEx.halfDt, mEx.fxt, mEx.xtMid)

	// Compute f(xtMid, t + halfDt)
	mEx.timeFunc.EvaluateOnRGridTimeInPlace(mEx.xtMid.Data, mEx.fxt.Data, t+mEx.halfDt)

	for i := 0; i < nPoints; i++ {
		if math.IsNaN(mEx.fxt.Data[i]) || math.IsInf(mEx.fxt.Data[i], 0) {
			return fmt.Errorf("the fxt is not valid")
		}
	}

	blas64.Axpy(mEx.delTime, mEx.fxt, blas64.Vector{N: nPoints, Data: xt, Inc: 1})

	return nil
}

type Ralston2order struct {
	timeFunc gridData.TDPotentialOp
	delTime  float64

	dt2by3 float64
	dt1by4 float64
	dt3by4 float64

	xtMid blas64.Vector
	fxt   blas64.Vector
}

func (rEx2 *Ralston2order) Name() string {
	return "Second order Ralston Explicit Method!"
}

func (rEx2 *Ralston2order) NewDefine(dt float64, tdFunc gridData.TDPotentialOp) *Ralston2order {
	return &Ralston2order{
		timeFunc: tdFunc,
		delTime:  dt,
		dt2by3:   (2. * dt) / 3.,
		dt1by4:   dt / 4.,
		dt3by4:   (3. * dt) / 4.,
	}
}

func (rEx2 *Ralston2order) ReDefine(dt float64, tdFunc gridData.TDPotentialOp) {
	rEx2.delTime = dt
	rEx2.timeFunc = tdFunc
	rEx2.dt1by4 = dt / 4.
	rEx2.dt2by3 = (2. * dt) / 3.
	rEx2.dt3by4 = (3. * dt) / 4.
}

func (rEx2 *Ralston2order) NextStep(xt, t float64) (float64, error) {
	k1 := rEx2.timeFunc.EvaluateAtTime(xt, t)
	k2 := rEx2.timeFunc.EvaluateAtTime(xt+rEx2.dt2by3*k1, t+rEx2.dt2by3)
	return xt + k1*rEx2.dt1by4 + k2*rEx2.dt3by4, nil
}

// Ralston3Order need to be completed
type Ralston3Order struct {
	timeFunc gridData.TDPotentialOp
	delTime  float64

	halfDt     float64
	threeBy4Dt float64
	dtBy9      float64

	xtMid blas64.Vector
	fxt   blas64.Vector
}

func (rEx *Ralston3Order) Name() string {
	return "Third order Ralston Explicit Method!"
}

func (rEx *Ralston3Order) NewDefine(dt float64, tdFunc gridData.TDPotentialOp) *Ralston3Order {
	return &Ralston3Order{
		timeFunc:   tdFunc,
		delTime:    dt,
		halfDt:     dt / 2.,
		threeBy4Dt: 0.75 * dt,
		dtBy9:      dt / 9.,
	}
}

func (rEx *Ralston3Order) ReDefine(dt float64, tdFunc gridData.TDPotentialOp) {
	rEx.delTime = dt
	rEx.timeFunc = tdFunc
	rEx.halfDt = dt / 2
	rEx.threeBy4Dt = 0.75 * dt
	rEx.dtBy9 = dt / 9.
}

func (rEx *Ralston3Order) NextStep(xt, t float64) (float64, error) {
	k1 := rEx.timeFunc.EvaluateAtTime(xt, t)
	k2 := rEx.timeFunc.EvaluateAtTime(xt+rEx.halfDt*k1, t+rEx.halfDt)
	k3 := rEx.timeFunc.EvaluateAtTime(xt+rEx.threeBy4Dt*k1, t+rEx.threeBy4Dt)
	return xt + rEx.dtBy9*(2*k1+3*k2+4*k3), nil
}

func (rEx *Ralston3Order) NextStepOnGrid(xt []float64, t float64) error {
	nPoints := len(xt)

	// Allocate buffers only if necessary
	if nPoints != rEx.fxt.N {
		rEx.fxt = blas64.Vector{N: nPoints, Data: make([]float64, nPoints), Inc: 1}
		rEx.xtMid = blas64.Vector{N: nPoints, Data: make([]float64, nPoints), Inc: 1}
	}

	// Compute xtMid = xt + halfDt * f(xt, t)
	rEx.timeFunc.EvaluateOnRGridTimeInPlace(xt, rEx.fxt.Data, t)
	copy(rEx.xtMid.Data, xt)
	blas64.Axpy(rEx.halfDt, rEx.fxt, rEx.xtMid)

	// Compute f(xtMid, t + halfDt)
	rEx.timeFunc.EvaluateOnRGridTimeInPlace(rEx.xtMid.Data, rEx.fxt.Data, t+rEx.halfDt)

	for i := 0; i < nPoints; i++ {
		if math.IsNaN(rEx.fxt.Data[i]) || math.IsInf(rEx.fxt.Data[i], 0) {
			return fmt.Errorf("the fxt is not valid")
		}
	}

	blas64.Axpy(rEx.delTime, rEx.fxt, blas64.Vector{N: nPoints, Data: xt, Inc: 1})

	return nil
}

type Huens3Explicit struct {
	timeFunc gridData.TDPotentialOp
	delTime  float64

	dt2by3 float64
	dtBy3  float64
	dtBy4  float64
	dt3by4 float64
}

func (h3Ex *Huens3Explicit) Name() string {
	return "Huen's Order 3"
}

func (h3Ex *Huens3Explicit) NewDef(dt float64, tdFunc gridData.TDPotentialOp) *Huens3Explicit {
	return &Huens3Explicit{
		timeFunc: tdFunc,
		delTime:  dt,

		dtBy3:  dt / 3.,
		dt2by3: (2 * dt) / 3.,
		dtBy4:  dt / 4.,
		dt3by4: (3 * dt) / 4.,
	}
}

func (h3Ex *Huens3Explicit) ReDefine(dt float64, tdFunc gridData.TDPotentialOp) {
	h3Ex.timeFunc = tdFunc
	h3Ex.delTime = dt

	h3Ex.dtBy3 = dt / 3.
	h3Ex.dt2by3 = (2 * dt) / 3.
	h3Ex.dtBy4 = dt / 4.
	h3Ex.dt3by4 = (3 * dt) / 4.
}

func (h3Ex *Huens3Explicit) NextStep(xt, t float64) (float64, error) {
	k1 := h3Ex.timeFunc.EvaluateAtTime(xt, t)
	k2 := h3Ex.timeFunc.EvaluateAtTime(xt+h3Ex.dtBy3*k1, t+h3Ex.dtBy3)
	k3 := h3Ex.timeFunc.EvaluateAtTime(xt+h3Ex.dt2by3*k2, t+h3Ex.dt2by3)
	return xt + h3Ex.dtBy4*(k1+3.*k3), nil
}

type VDHouwenExplicit struct {
	timeFunc gridData.TDPotentialOp
	delTime  float64

	dt8by15 float64
	dt5by12 float64
	dt2by3  float64
	dtBy4   float64
}

func (vdh3Ex *VDHouwenExplicit) Name() string {
	return "Van der Houwen's/Wray's third-order method"
}

func (vdh3Ex *VDHouwenExplicit) NewDef(dt float64, tdFunc gridData.TDPotentialOp) *VDHouwenExplicit {
	return &VDHouwenExplicit{
		timeFunc: tdFunc,
		delTime:  dt,
		dtBy4:    dt / 4.,
		dt2by3:   (2 * dt) / 3.,
		dt5by12:  (5 * dt) / 12.,
		dt8by15:  (8 * dt) / 15.,
	}
}

func (vdh3Ex *VDHouwenExplicit) ReDefine(dt float64, tdFunc gridData.TDPotentialOp) {
	vdh3Ex.timeFunc = tdFunc
	vdh3Ex.delTime = dt
	vdh3Ex.dtBy4 = dt / 4.
	vdh3Ex.dt2by3 = (2 * dt) / 3.
	vdh3Ex.dt5by12 = (5 * dt) / 12.
	vdh3Ex.dt8by15 = (8 * dt) / 15.
}

func (vdh3Ex *VDHouwenExplicit) NextStep(xt, t float64) (float64, error) {
	k1 := vdh3Ex.timeFunc.EvaluateAtTime(xt, t)
	k2 := vdh3Ex.timeFunc.EvaluateAtTime(xt+vdh3Ex.dt8by15*k1, t+vdh3Ex.dt8by15)
	val := vdh3Ex.dtBy4*k1 + vdh3Ex.dt5by12*k2
	k3 := vdh3Ex.timeFunc.EvaluateAtTime(xt+val, t+vdh3Ex.dt2by3)
	return xt + vdh3Ex.dtBy4*(k1+3.*k3), nil
}

type SSPRungeKutta3 struct {
	timeFunc gridData.TDPotentialOp
	delTime  float64

	halfDt float64
	dtBy4  float64
	dtBy6  float64
	dt2by3 float64
}

func (ssprk3 *SSPRungeKutta3) Name() string {
	return "Third-order Strong Stability Preserving Runge-Kutta"
}

func (ssprk3 *SSPRungeKutta3) NewDef(dt float64, tdFunc gridData.TDPotentialOp) *SSPRungeKutta3 {
	return &SSPRungeKutta3{
		timeFunc: tdFunc,
		delTime:  dt,
		halfDt:   dt / 2,
		dtBy4:    dt / 4.,
		dtBy6:    dt / 6.,
		dt2by3:   (2 * dt) / 3.,
	}
}

func (ssprk3 *SSPRungeKutta3) ReDefine(dt float64, tdFunc gridData.TDPotentialOp) {
	ssprk3.timeFunc = tdFunc
	ssprk3.delTime = dt
	ssprk3.dtBy4 = dt / 4.
	ssprk3.dtBy6 = dt / 6.
	ssprk3.halfDt = dt / 2.
	ssprk3.dt2by3 = (2 * dt) / 3.
}

func (ssprk3 *SSPRungeKutta3) NextStep(xt, t float64) (float64, error) {
	k1 := ssprk3.timeFunc.EvaluateAtTime(xt, t)
	k2 := ssprk3.timeFunc.EvaluateAtTime(xt+ssprk3.delTime*k1, t+ssprk3.delTime)
	val := (k1 + k2) * ssprk3.dtBy4
	k3 := ssprk3.timeFunc.EvaluateAtTime(xt+val, t+ssprk3.halfDt)
	return xt + ssprk3.dtBy6*(k1+k2) + k3*ssprk3.dt2by3, nil
}

type RungeKutta3Explicit struct {
	timeFunc gridData.TDPotentialOp
	delTime  float64

	halfDt float64
	dt2by3 float64
	dtby6  float64
}

func (rg3Ex *RungeKutta3Explicit) Name() string {
	return "Runge-Kutta Order 3"
}

func (rg3Ex *RungeKutta3Explicit) NewDef(dt float64, tdFunc gridData.TDPotentialOp) *RungeKutta3Explicit {
	return &RungeKutta3Explicit{
		timeFunc: tdFunc,
		delTime:  dt,
		halfDt:   dt / 2,
		dt2by3:   (2 * dt) / 3.,
		dtby6:    dt / 6.,
	}
}

func (rg3Ex *RungeKutta3Explicit) ReDefine(dt float64, tdFunc gridData.TDPotentialOp) {
	rg3Ex.timeFunc = tdFunc
	rg3Ex.delTime = dt
	rg3Ex.dtby6 = dt / 6
	rg3Ex.halfDt = dt / 2
	rg3Ex.dt2by3 = (2 * dt) / 3.
}

func (rg3Ex *RungeKutta3Explicit) NextStep(xt, t float64) (float64, error) {
	k1 := rg3Ex.timeFunc.EvaluateAtTime(xt, t)
	k2 := rg3Ex.timeFunc.EvaluateAtTime(xt+rg3Ex.halfDt*k1, t+rg3Ex.halfDt)
	val := (-k1 + 2*k2) * rg3Ex.delTime
	k3 := rg3Ex.timeFunc.EvaluateAtTime(xt+val, t+rg3Ex.delTime)
	return xt + rg3Ex.dtby6*(k1+4.*k2+k3), nil
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
	k1 := rgEx.timeFunc.EvaluateAtTime(xt, t)
	k2 := rgEx.timeFunc.EvaluateAtTime(xt+k1*rgEx.halfDt, t+rgEx.halfDt)
	k3 := rgEx.timeFunc.EvaluateAtTime(xt+k2*rgEx.halfDt, t+rgEx.halfDt)
	k4 := rgEx.timeFunc.EvaluateAtTime(xt+k3*rgEx.deltaTime, t+rgEx.deltaTime)

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
	rgEx.timeFunc.EvaluateOnRGridTimeInPlace(xt, rgEx.k1.Data, t)

	// k2 = f(x + dt/2 * k1, t + dt/2)
	copy(rgEx.xPlusDt.Data, xt)
	blas64.Axpy(rgEx.halfDt, rgEx.k1, rgEx.xPlusDt)
	rgEx.timeFunc.EvaluateOnRGridTimeInPlace(rgEx.xPlusDt.Data, rgEx.k2.Data, t+rgEx.halfDt)

	// k3 = f(x + dt/2 * k2, t + dt/2)
	copy(rgEx.xPlusDt.Data, xt)
	blas64.Axpy(rgEx.halfDt, rgEx.k2, rgEx.xPlusDt)
	rgEx.timeFunc.EvaluateOnRGridTimeInPlace(rgEx.xPlusDt.Data, rgEx.k3.Data, t+rgEx.halfDt)

	// k4 = f(x + dt * k3, t + dt)
	copy(rgEx.xPlusDt.Data, xt)
	blas64.Axpy(rgEx.deltaTime, rgEx.k3, rgEx.xPlusDt)
	rgEx.timeFunc.EvaluateOnRGridTimeInPlace(rgEx.xPlusDt.Data, rgEx.k4.Data, t+rgEx.deltaTime)

	// x_new = x + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
	for i := 0; i < nPoints; i++ {
		xt[i] += rgEx.dtBy6 * (rgEx.k1.Data[i] + 2*(rgEx.k2.Data[i]+rgEx.k3.Data[i]) + rgEx.k4.Data[i])
	}

	return nil
}

// RungeKutta38 implements the classic 4th-order Runge-Kutta method
type RungeKutta38 struct {
	timeFunc  gridData.TDPotentialOp
	deltaTime float64

	// Pre-allocated buffers for grid operations
	dtBy3  float64
	dt2By3 float64
	dtBy8  float64

	xPlusDt blas64.Vector
	k1      blas64.Vector
	k2      blas64.Vector
	k3      blas64.Vector
	k4      blas64.Vector
}

func (rg38 *RungeKutta38) Name() string {
	return "Runge-Kutta-38 Order 4"
}

func (rg38 *RungeKutta38) NewDefine(dt float64, tdFunc gridData.TDPotentialOp) *RungeKutta38 {
	return &RungeKutta38{
		timeFunc:  tdFunc,
		deltaTime: dt,
		dtBy3:     dt / 3,
		dt2By3:    (2 * dt) / 3.,
		dtBy8:     dt / 8,
	}
}

func (rg38 *RungeKutta38) ReDefine(dt float64, tdFunc gridData.TDPotentialOp) {
	rg38.timeFunc = tdFunc
	rg38.deltaTime = dt
	rg38.dtBy3 = dt / 3
	rg38.dt2By3 = (2 * dt) / 3
	rg38.dtBy8 = dt / 8
}

func (rg38 *RungeKutta38) NextStep(xt, t float64) (float64, error) {
	k1 := rg38.timeFunc.EvaluateAtTime(xt, t)
	k2 := rg38.timeFunc.EvaluateAtTime(xt+k1*rg38.dtBy3, t+rg38.dtBy3)
	val := k2*rg38.deltaTime - k1*rg38.dtBy3
	k3 := rg38.timeFunc.EvaluateAtTime(xt+val, t+rg38.dt2By3)
	val = rg38.deltaTime * (k1 - k2 + k3)
	k4 := rg38.timeFunc.EvaluateAtTime(xt+val, t+rg38.deltaTime)
	slope := k1 + 3*(k2+k3) + k4

	if math.IsNaN(slope) || math.IsInf(slope, 0) {
		return xt, fmt.Errorf("invalid derivative: NaN or Inf encountered")
	}

	return xt + rg38.dtBy8*slope, nil
}

// Ralston4Order Ralston's fourth-order method
type Ralston4Order struct {
	timeFunc  gridData.TDPotentialOp
	deltaTime float64

	// Pre-allocated buffers for grid operations
	dtCoefs  []float64
	ksCoefs  []float64
	fnlCoefs []float64
}

func (ral4 *Ralston4Order) Name() string {
	return "Ralston's fourth-order method!"
}

func (ral4 *Ralston4Order) coefficients(dt float64, dts, ks, bs []float64) {
	dts[0] = (2 / 5.) * dt
	dts[1] = ((14. - 3*math.Sqrt(5.)) / 16.) * dt
	dts[2] = dt

	ks[0] = (2 / 5.) * dt
	ks[1] = ((-2889. + 1428*math.Sqrt(5.)) / 1024) * dt
	ks[2] = ((3785 - 1620*math.Sqrt(5.)) / 1024) * dt
	ks[3] = ((-3365 + 2094*math.Sqrt(5.)) / 6040) * dt
	ks[4] = ((-975 - 3046*math.Sqrt(5.)) / 2552) * dt
	ks[5] = ((467040 + 203968*math.Sqrt(5.)) / 240845) * dt

	bs[0] = ((263 + 24*math.Sqrt(5.)) / 1812) * dt
	bs[1] = ((125 - 1000*math.Sqrt(5.)) / 3828) * dt
	bs[2] = ((3426304 + 1661952*math.Sqrt(5.)) / 5924787) * dt
	bs[3] = ((30 - 4*math.Sqrt(5.)) / 123) * dt
}

func (ral4 *Ralston4Order) NewDefine(dt float64, tdFunc gridData.TDPotentialOp) *Ralston4Order {
	dts := make([]float64, 3)
	ks := make([]float64, 6)
	bs := make([]float64, 4)
	ral4.coefficients(dt, dts, ks, bs)

	return &Ralston4Order{
		timeFunc:  tdFunc,
		deltaTime: dt,
		dtCoefs:   dts,
		ksCoefs:   ks,
		fnlCoefs:  bs,
	}
}

func (ral4 *Ralston4Order) ReDefine(dt float64, tdFunc gridData.TDPotentialOp) {
	dts := make([]float64, 4)
	ks := make([]float64, 6)
	bs := make([]float64, 4)
	ral4.coefficients(dt, dts, ks, bs)

	ral4.timeFunc = tdFunc
	ral4.deltaTime = dt
	ral4.dtCoefs = dts
	ral4.ksCoefs = ks
	ral4.fnlCoefs = bs
}

func (ral4 *Ralston4Order) NextStep(xt, t float64) (float64, error) {
	k1 := ral4.timeFunc.EvaluateAtTime(xt, t)
	k2 := ral4.timeFunc.EvaluateAtTime(xt+k1*ral4.ksCoefs[0], t+ral4.dtCoefs[0])

	val := k1*ral4.ksCoefs[1] + k2*ral4.ksCoefs[2]
	k3 := ral4.timeFunc.EvaluateAtTime(xt+val, t+ral4.dtCoefs[1])

	val = k1*ral4.ksCoefs[3] + k2*ral4.ksCoefs[4] + k3*ral4.ksCoefs[5]
	k4 := ral4.timeFunc.EvaluateAtTime(xt+val, t+ral4.dtCoefs[2])

	integrant := k1*ral4.fnlCoefs[0] + k2*ral4.fnlCoefs[1] + k3*ral4.fnlCoefs[2] + k4*ral4.fnlCoefs[3]

	if math.IsNaN(integrant) || math.IsInf(integrant, 0) {
		return xt, fmt.Errorf("invalid derivative: NaN or Inf encountered")
	}

	return xt + integrant, nil
}

type Nystrom5Explicit struct {
	timeFunc  gridData.TDPotentialOp
	deltaTime float64

	// Pre-allocated buffers for grid operations
	dtCoefs  []float64
	ksCoefs  []float64
	fnlCoefs []float64
}

func (nRK5ex *Nystrom5Explicit) Name() string {
	return "Nystrom's fifth-order method"
}

func (nRK5ex *Nystrom5Explicit) coefficients(dt float64, dts, ks, bs []float64) {
	dts[0] = (1 / 3.) * dt
	dts[1] = (2. / 5.) * dt
	dts[2] = dt
	dts[3] = (2. / 3.) * dt
	dts[4] = (4. / 5.) * dt

	ks[0] = dts[0]
	ks[1] = (4. / 25.) * dt
	ks[2] = (6. / 25.) * dt
	ks[3] = (1 / 4.) * dt
	ks[4] = -3 * dt
	ks[5] = (15. / 4.) * dt
	ks[6] = (2. / 27.) * dt
	ks[7] = (10. / 9.) * dt
	ks[8] = (-50. / 81.) * dt
	ks[9] = (8. / 81.) * dt
	ks[10] = (2. / 25.) * dt
	ks[11] = (12. / 25.) * dt
	ks[12] = (2. / 15.) * dt
	ks[13] = (8. / 75.) * dt

	bs[0] = (23. / 192.) * dt
	bs[1] = (125. / 192.) * dt
	bs[2] = (-27. / 64) * dt
}

func (nRK5ex *Nystrom5Explicit) NewDefine(dt float64, tdFunc gridData.TDPotentialOp) *Nystrom5Explicit {
	dts := make([]float64, 5)
	ks := make([]float64, 14)
	bs := make([]float64, 3)
	nRK5ex.coefficients(dt, dts, ks, bs)

	return &Nystrom5Explicit{
		timeFunc:  tdFunc,
		deltaTime: dt,
		dtCoefs:   dts,
		ksCoefs:   ks,
		fnlCoefs:  bs,
	}
}

func (nRK5ex *Nystrom5Explicit) ReDefine(dt float64, tdFunc gridData.TDPotentialOp) {
	dts := make([]float64, 5)
	ks := make([]float64, 14)
	bs := make([]float64, 3)
	nRK5ex.coefficients(dt, dts, ks, bs)

	nRK5ex.timeFunc = tdFunc
	nRK5ex.deltaTime = dt
	nRK5ex.dtCoefs = dts
	nRK5ex.ksCoefs = ks
	nRK5ex.fnlCoefs = bs
}

func (nRK5ex *Nystrom5Explicit) NextStep(xt, t float64) (float64, error) {
	k1 := nRK5ex.timeFunc.EvaluateAtTime(xt, t)
	k2 := nRK5ex.timeFunc.EvaluateAtTime(xt+k1*nRK5ex.ksCoefs[0], t+nRK5ex.dtCoefs[0])

	val := k1*nRK5ex.ksCoefs[1] + k2*nRK5ex.ksCoefs[2]
	k3 := nRK5ex.timeFunc.EvaluateAtTime(xt+val, t+nRK5ex.dtCoefs[1])

	val = k1*nRK5ex.ksCoefs[3] + k2*nRK5ex.ksCoefs[4] + k3*nRK5ex.ksCoefs[5]
	k4 := nRK5ex.timeFunc.EvaluateAtTime(xt+val, t+nRK5ex.dtCoefs[2])

	val = k1*nRK5ex.ksCoefs[6] + k2*nRK5ex.ksCoefs[7] + k3*nRK5ex.ksCoefs[8] + k4*nRK5ex.ksCoefs[8]
	k5 := nRK5ex.timeFunc.EvaluateAtTime(xt+val, t+nRK5ex.dtCoefs[3])

	val = k1*nRK5ex.ksCoefs[9] + k2*nRK5ex.ksCoefs[10] + k3*nRK5ex.ksCoefs[11] + k4*nRK5ex.ksCoefs[12]
	k6 := nRK5ex.timeFunc.EvaluateAtTime(xt+val, t+nRK5ex.dtCoefs[5])

	integrant := k1*nRK5ex.fnlCoefs[0] + k3*nRK5ex.fnlCoefs[1] + k5*nRK5ex.fnlCoefs[2] + k6*nRK5ex.fnlCoefs[1]

	if math.IsNaN(integrant) || math.IsInf(integrant, 0) {
		return xt, fmt.Errorf("invalid derivative: NaN or Inf encountered")
	}

	return xt + integrant, nil
}
