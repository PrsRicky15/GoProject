package mathsFunc

import "GoProject/gridData"

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

type FehlbergRK12 struct {
	timeFunc  gridData.TDPotentialOp
	deltaTime float64
	halfDt    float64
	kCoefs    []float64
	b1Coefs   []float64
	b2Coefs   []float64
}

func (HE *FehlbergRK12) Name() string {
	return "The simplest adaptive Runge–Kutta method !"
}

func (HE *FehlbergRK12) NewDefine(dt float64, tdFunc gridData.TDPotentialOp) *FehlbergRK12 {

	return &FehlbergRK12{
		timeFunc:  tdFunc,
		deltaTime: dt,
		halfDt:    dt / 2.,
		kCoefs:    []float64{1. / 256., 255. / 256.},
		b1Coefs:   []float64{1. / 512., 255. / 256.},
		b2Coefs:   []float64{1. / 256., 255. / 256.},
	}
}

func (HE *FehlbergRK12) ReDefine(dt float64, tdFunc gridData.TDPotentialOp) {
	HE.deltaTime = dt
	HE.timeFunc = tdFunc
	HE.halfDt = dt / 2
}

type BogackiShampine struct {
	timeFunc  gridData.TDPotentialOp
	deltaTime float64
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
