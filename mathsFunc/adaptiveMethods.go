package mathsFunc

import (
	"GoProject/gridData"
	"math"
)

type HuensEuler struct {
	timeFunc  gridData.TDPotentialOp
	deltaTime float64
	halfDt    float64
}

func (HE *HuensEuler) Name() string {
	return "The simplest adaptive Runge–Kutta method !"
}

func (HE *HuensEuler) NewDefine(dt float64, tdFunc gridData.TDPotentialOp) *HuensEuler {
	return &HuensEuler{
		timeFunc:  tdFunc,
		deltaTime: dt,
		halfDt:    dt / 2.,
	}
}

func (HE *HuensEuler) ReDefine(dt float64, tdFunc gridData.TDPotentialOp) {
	HE.deltaTime = dt
	HE.timeFunc = tdFunc
	HE.halfDt = dt / 2
}

func (HE *HuensEuler) adaptiveDt(dt float64) {
	HE.deltaTime = dt
	HE.halfDt = dt / 2
}

func (HE *HuensEuler) NextStep(xt, t float64) (float64, error) {
	k1 := HE.timeFunc.EvaluateAt(xt, t)
	k2 := HE.timeFunc.EvaluateAt(xt+k1*HE.deltaTime, t+HE.deltaTime)

	xtPlusdt1 := xt + HE.halfDt*(k1+k2)
	xtPlusdt2 := xt + k1*HE.deltaTime

	Err := math.Abs(xtPlusdt1-xtPlusdt2) / HE.deltaTime

	if Err <= adaptiveTolerance {
		return xtPlusdt2, nil
	}

	HE.adaptiveDt(0.9 * HE.deltaTime * math.Sqrt(tolerance/Err))
	return xt, nil
}

type FehlbergRK12 struct {
	timeFunc  gridData.TDPotentialOp
	deltaTime float64
	halfDt    float64
	kCoefs    []float64
	b1Coefs   []float64
	b2Coefs   []float64
}

func (FRK12 *FehlbergRK12) Name() string {
	return "The simplest adaptive Runge–Kutta method !"
}

func (FRK12 *FehlbergRK12) NewDefine(dt float64, tdFunc gridData.TDPotentialOp) *FehlbergRK12 {

	return &FehlbergRK12{
		timeFunc:  tdFunc,
		deltaTime: dt,
		halfDt:    dt / 2.,
		kCoefs:    []float64{1. / 256., 255. / 256.},
		b1Coefs:   []float64{1. / 512., 255. / 256.},
		b2Coefs:   []float64{1. / 256., 255. / 256.},
	}
}

func (FRK12 *FehlbergRK12) ReDefine(dt float64, tdFunc gridData.TDPotentialOp) {
	FRK12.deltaTime = dt
	FRK12.timeFunc = tdFunc
	FRK12.halfDt = dt / 2
}

func (FRK12 *FehlbergRK12) adaptiveDt(dt float64) {
	FRK12.deltaTime = dt
	FRK12.halfDt = dt / 2
}

func (FRK12 *FehlbergRK12) NextStep(xt, t float64) (float64, error) {
	k1 := FRK12.timeFunc.EvaluateAt(xt, t)
	k2 := FRK12.timeFunc.EvaluateAt(xt+k1*FRK12.halfDt, t+FRK12.halfDt)

	val := k1*FRK12.kCoefs[0] + k2*FRK12.kCoefs[1]
	k3 := FRK12.timeFunc.EvaluateAt(xt+val, t+FRK12.deltaTime)

	xtPlusdt1 := xt + FRK12.deltaTime*(k1*FRK12.b1Coefs[0]+k2*FRK12.b1Coefs[1]+k3*FRK12.b1Coefs[0])
	xtPlusdt2 := xt + FRK12.deltaTime*(k1*FRK12.b2Coefs[0]+k2*FRK12.b1Coefs[1])

	Err := math.Abs(xtPlusdt1-xtPlusdt2) / FRK12.deltaTime

	if Err <= adaptiveTolerance {
		return xtPlusdt2, nil
	}

	FRK12.adaptiveDt(0.9 * FRK12.deltaTime * math.Sqrt(tolerance/Err))
	return xt, nil
}

type BogackiShampine struct {
	timeFunc  gridData.TDPotentialOp
	deltaTime float64
	halfDt    float64
	dt3By4    float64

	kCoefs  []float64
	b1Coefs []float64
	b2Coefs []float64
}

func (BS *BogackiShampine) Name() string {
	return "The simplest adaptive Runge–Kutta method !"
}

func (BS *BogackiShampine) NewDefine(dt float64, tdFunc gridData.TDPotentialOp) *BogackiShampine {

	return &BogackiShampine{
		timeFunc:  tdFunc,
		deltaTime: dt,
		halfDt:    dt / 2.,
		dt3By4:    (3 * dt) / 4.,
		kCoefs:    []float64{2. / 9., 1. / 3., 4. / 9.},
		b1Coefs:   []float64{2. / 9., 1. / 3., 4. / 9.},
		b2Coefs:   []float64{7. / 24., 0.25, 1. / 3., 1. / 8.},
	}
}

func (BS *BogackiShampine) ReDefine(dt float64, tdFunc gridData.TDPotentialOp) {
	BS.deltaTime = dt
	BS.timeFunc = tdFunc
	BS.halfDt = dt / 2
	BS.dt3By4 = (3 * dt) / 4.
}

func (BS *BogackiShampine) adaptiveDt(dt float64) {
	BS.deltaTime = dt
	BS.halfDt = dt / 2
	BS.dt3By4 = (3 * dt) / 4.
}

func (BS *BogackiShampine) NextStep(xt, t float64) (float64, error) {
	k1 := BS.timeFunc.EvaluateAt(xt, t)
	k2 := BS.timeFunc.EvaluateAt(xt+k1*BS.halfDt, t+BS.halfDt)
	k3 := BS.timeFunc.EvaluateAt(xt+k2*BS.dt3By4, t+BS.dt3By4)

	val := k1*BS.kCoefs[0] + k2*BS.kCoefs[1] + k3*BS.kCoefs[2]
	k4 := BS.timeFunc.EvaluateAt(xt+val, t+BS.deltaTime)

	xtPlusdt1 := xt + BS.deltaTime*(k1*BS.b1Coefs[0]+k2*BS.b1Coefs[1]+k3*BS.b1Coefs[2])
	xtPlusdt2 := xt + BS.deltaTime*(k1*BS.b1Coefs[0]+k2*BS.b1Coefs[1]+k3*BS.b1Coefs[2]+k4*BS.b1Coefs[3])

	Err := math.Abs(xtPlusdt1-xtPlusdt2) / BS.deltaTime

	if Err <= adaptiveTolerance {
		return xtPlusdt2, nil
	}

	BS.adaptiveDt(0.9 * BS.deltaTime * math.Sqrt(tolerance/Err))
	return xt, nil
}

type RKFelberg struct {
	timeFunc  gridData.TDPotentialOp
	deltaTime float64
}

type CashKarp struct {
	timeFunc  gridData.TDPotentialOp
	deltaTime float64
}

func (ckEx *CashKarp) NextStep(xt, t float64) (float64, error) {
	k1 := ckEx.timeFunc.EvaluateAt(xt, t)

	deltaTimeStep1 := ckEx.deltaTime / 5
	k2 := ckEx.timeFunc.EvaluateAt(xt+deltaTimeStep1*k1, t+deltaTimeStep1)

	deltaTimeStep2 := deltaTimeStep1 * 3
	k3 := ckEx.timeFunc.EvaluateAt(xt+deltaTimeStep2*(k1/8+3*k2/8), t+deltaTimeStep2/2)

	val := deltaTimeStep1 * (3*k1/2 - 9*k2/2 + 6*k3)
	k4 := ckEx.timeFunc.EvaluateAt(xt+val, t+3*deltaTimeStep1)

	val = 5*k2/2 + (-11*k1/2-35*(2*k3+k4))/27
	k5 := ckEx.timeFunc.EvaluateAt(xt+ckEx.deltaTime*val, t+ckEx.deltaTime)

	val = 1631*k1/55296 + 175*k2/512 + 575*k3/13824 + 44275*k4/110592 + 253*k5/4096
	k6 := ckEx.timeFunc.EvaluateAt(xt+ckEx.deltaTime*val, t+7*ckEx.deltaTime/8)
	return xt + k5*k6, nil
}

type DormandPrince struct {
	timeFunc  gridData.TDPotentialOp
	deltaTime float64
}
