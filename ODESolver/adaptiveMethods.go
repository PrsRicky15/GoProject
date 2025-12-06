package EquationSolver

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
	k1 := HE.timeFunc.EvaluateAtTime(xt, t)
	k2 := HE.timeFunc.EvaluateAtTime(xt+k1*HE.deltaTime, t+HE.deltaTime)

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
	k1 := FRK12.timeFunc.EvaluateAtTime(xt, t)
	k2 := FRK12.timeFunc.EvaluateAtTime(xt+k1*FRK12.halfDt, t+FRK12.halfDt)

	val := k1*FRK12.kCoefs[0] + k2*FRK12.kCoefs[1]
	k3 := FRK12.timeFunc.EvaluateAtTime(xt+val*FRK12.deltaTime, t+FRK12.deltaTime)

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
	k1 := BS.timeFunc.EvaluateAtTime(xt, t)
	k2 := BS.timeFunc.EvaluateAtTime(xt+k1*BS.halfDt, t+BS.halfDt)
	k3 := BS.timeFunc.EvaluateAtTime(xt+k2*BS.dt3By4, t+BS.dt3By4)

	val := k1*BS.kCoefs[0] + k2*BS.kCoefs[1] + k3*BS.kCoefs[2]
	k4 := BS.timeFunc.EvaluateAtTime(xt+val*BS.deltaTime, t+BS.deltaTime)

	xtPlusdt1 := xt + BS.deltaTime*(k1*BS.b1Coefs[0]+k2*BS.b1Coefs[1]+k3*BS.b1Coefs[2])
	xtPlusdt2 := xt + BS.deltaTime*(k1*BS.b1Coefs[0]+k2*BS.b1Coefs[1]+k3*BS.b1Coefs[2]+k4*BS.b1Coefs[3])

	Err := math.Abs(xtPlusdt1-xtPlusdt2) / BS.deltaTime

	if Err <= adaptiveTolerance {
		return xtPlusdt2, nil
	}

	BS.adaptiveDt(0.9 * BS.deltaTime * math.Sqrt(tolerance/Err))
	return xt, nil
}

// adaptiveRKBase contains the common implementation for adaptive RK methods
type adaptiveRKBase struct {
	timeFunc  gridData.TDPotentialOp
	deltaTime float64

	dtCoefs []float64
	kCoefs  []float64
	b1Coefs []float64
	b2Coefs []float64
}

func (ark *adaptiveRKBase) ReDefine(dt float64, tdFunc gridData.TDPotentialOp) {
	ark.deltaTime = dt
	ark.timeFunc = tdFunc
}

func (ark *adaptiveRKBase) adaptiveDt(dt float64) {
	ark.deltaTime = dt
}

func (ark *adaptiveRKBase) NextStep(xt, t float64) (float64, error) {
	k1 := ark.timeFunc.EvaluateAtTime(xt, t)
	k2 := ark.timeFunc.EvaluateAtTime(xt+k1*ark.kCoefs[0]*ark.deltaTime, t+ark.dtCoefs[0])

	val := k1*ark.kCoefs[1] + k2*ark.kCoefs[2]
	k3 := ark.timeFunc.EvaluateAtTime(xt+val*ark.deltaTime, t+ark.dtCoefs[1])

	val = k1*ark.kCoefs[3] + k2*ark.kCoefs[4] + k3*ark.kCoefs[5]
	k4 := ark.timeFunc.EvaluateAtTime(xt+val*ark.deltaTime, t+ark.dtCoefs[2])

	val = k1*ark.kCoefs[6] + k2*ark.kCoefs[7] + k3*ark.kCoefs[8] + k4*ark.kCoefs[9]
	k5 := ark.timeFunc.EvaluateAtTime(xt+val*ark.deltaTime, t+ark.dtCoefs[3])

	val = k1*ark.kCoefs[10] + k2*ark.kCoefs[11] + k3*ark.kCoefs[12] + k4*ark.kCoefs[13] + k5*ark.kCoefs[14]
	k6 := ark.timeFunc.EvaluateAtTime(xt+val*ark.deltaTime, t+ark.dtCoefs[4])

	integrant := k1*ark.b1Coefs[0] + k3*ark.b1Coefs[1] + k4*ark.b1Coefs[2] + k5*ark.b1Coefs[3] + k6*ark.b1Coefs[4]
	xtPlusdt1 := xt + ark.deltaTime*integrant

	integrant = k1*ark.b2Coefs[0] + k3*ark.b2Coefs[1] + k4*ark.b2Coefs[2] + k5*ark.b2Coefs[3]
	xtPlusdt2 := xt + ark.deltaTime*integrant

	err := math.Abs(xtPlusdt1-xtPlusdt2) / ark.deltaTime

	if err <= adaptiveTolerance {
		return xtPlusdt2, nil
	}

	ark.adaptiveDt(0.9 * ark.deltaTime * math.Sqrt(tolerance/err))
	return xt, nil
}

// RKFelberg - Runge-Kutta-Fehlberg method
type RKFelberg struct {
	adaptiveRKBase
}

func (RKF *RKFelberg) Name() string {
	return "The simplest adaptive Runge–Kutta method !"
}

func (RKF *RKFelberg) NewDefine(dt float64, tdFunc gridData.TDPotentialOp) *RKFelberg {
	return &RKFelberg{
		adaptiveRKBase: adaptiveRKBase{
			timeFunc:  tdFunc,
			deltaTime: dt,
			dtCoefs:   []float64{0.25, 3. / 8., 12. / 13., 1, 0.5},
			kCoefs: []float64{
				0.25,
				3. / 32., 9. / 32.,
				1932. / 2197., -7200. / 2197., 7296. / 2197.,
				439. / 216., -8, 3680. / 513., -845. / 4104.,
				-8. / 27., 2., -3544. / 2565., 1859. / 4104., -11. / 40.},
			b1Coefs: []float64{16. / 135., 6656. / 12825., 28561. / 56430., -9. / 50., 2. / 55.},
			b2Coefs: []float64{25. / 216., 1408. / 2565., 2197. / 4104., -1. / 5., 0},
		},
	}
}

// CashKarp - Cash-Karp method
type CashKarp struct {
	adaptiveRKBase
}

func (CK *CashKarp) Name() string {
	return "The Cash Karp method for adaptive ODE!"
}

func (CK *CashKarp) NewDefine(dt float64, tdFunc gridData.TDPotentialOp) *CashKarp {
	return &CashKarp{
		adaptiveRKBase: adaptiveRKBase{
			timeFunc:  tdFunc,
			deltaTime: dt,
			dtCoefs:   []float64{0.2, 3. / 10., 3. / 5., 1, 7. / 8.},
			kCoefs: []float64{
				0.2,
				3. / 40., 9. / 40.,
				3. / 10., -9. / 10., 6. / 5.,
				-11. / 54., 5. / 2., -70. / 27., 35. / 27.,
				1631. / 55296., 175. / 512., 575. / 13824., 44275. / 110592., 253. / 4096.},
			b1Coefs: []float64{37. / 378., 250. / 621., 125. / 594., 512. / 1771., 0},
			b2Coefs: []float64{2825. / 27648., 18575. / 48384., 13525. / 55296., 277. / 14336., 0.25},
		},
	}
}

type DormandPrince struct {
	timeFunc  gridData.TDPotentialOp
	deltaTime float64

	dtCoefs []float64
	kCoefs  []float64
	b1Coefs []float64
	b2Coefs []float64
}

func (DP *DormandPrince) Name() string {
	return "The simplest adaptive Runge–Kutta method !"
}

func (DP *DormandPrince) NewDefine(dt float64, tdFunc gridData.TDPotentialOp) *DormandPrince {

	return &DormandPrince{
		timeFunc:  tdFunc,
		deltaTime: dt,
		dtCoefs:   []float64{0.2, 3. / 10., 4. / 5., 8. / 9.},
		kCoefs: []float64{
			0.2,
			3. / 40., 9. / 40.,
			44. / 45., -56. / 15., 32. / 9.,
			19372. / 6561., 25360. / 2187., 64448. / 6561., -212. / 729.,
			9017. / 3168., 355. / 33., 46732. / 5247., 49. / 176., 5103. / 18656.,
		},
		b1Coefs: []float64{35. / 384., 500. / 1113., 125. / 192., -2187. / 6784., 11. / 84.},
		b2Coefs: []float64{5179. / 57600., 7571. / 16695., 393. / 640., -92097. / 339200.,
			187. / 2100., 1. / 40.},
	}
}

func (DP *DormandPrince) ReDefine(dt float64, tdFunc gridData.TDPotentialOp) {
	DP.deltaTime = dt
	DP.timeFunc = tdFunc
}

func (DP *DormandPrince) adaptiveDt(dt float64) {
	DP.deltaTime = dt
}

func (DP *DormandPrince) NextStep(xt, t float64) (float64, error) {
	k1 := DP.timeFunc.EvaluateAtTime(xt, t)
	k2 := DP.timeFunc.EvaluateAtTime(xt+k1*DP.kCoefs[0]*DP.deltaTime, t+DP.dtCoefs[0])

	val := k1*DP.kCoefs[1] + k2*DP.kCoefs[2]
	k3 := DP.timeFunc.EvaluateAtTime(xt+val*DP.deltaTime, t+DP.dtCoefs[1])

	val = k1*DP.kCoefs[3] + k2*DP.kCoefs[4] + k3*DP.kCoefs[5]
	k4 := DP.timeFunc.EvaluateAtTime(xt+val*DP.deltaTime, t+DP.dtCoefs[2])

	val = k1*DP.kCoefs[6] + k2*DP.kCoefs[7] + k3*DP.kCoefs[8] + k4*DP.kCoefs[9]
	k5 := DP.timeFunc.EvaluateAtTime(xt+val*DP.deltaTime, t+DP.dtCoefs[3])

	val = k1*DP.kCoefs[10] + k2*DP.kCoefs[11] + k3*DP.kCoefs[12] + k4*DP.kCoefs[13] + k5*DP.kCoefs[14]
	k6 := DP.timeFunc.EvaluateAtTime(xt+val*DP.deltaTime, t+DP.deltaTime)

	val = k1*DP.b1Coefs[0] + k3*DP.b1Coefs[1] + k4*DP.b1Coefs[2] + k5*DP.b1Coefs[3] + k6*DP.b1Coefs[4]
	k7 := DP.timeFunc.EvaluateAtTime(xt+val*DP.deltaTime, t+DP.deltaTime)

	integrant := k1*DP.b1Coefs[0] + k3*DP.b1Coefs[1] + k4*DP.b1Coefs[2] + k5*DP.b1Coefs[3] + k6*DP.b1Coefs[4]
	xtPlusdt1 := xt + DP.deltaTime*integrant

	integrant = k1*DP.b2Coefs[0] + k3*DP.b2Coefs[1] + k4*DP.b2Coefs[2] + k5*DP.b2Coefs[3] +
		k6*DP.b2Coefs[4] + k7*DP.b2Coefs[5]

	xtPlusdt2 := xt + DP.deltaTime*integrant

	err := math.Abs(xtPlusdt1-xtPlusdt2) / DP.deltaTime

	if err <= adaptiveTolerance {
		return xtPlusdt2, nil
	}

	DP.adaptiveDt(0.9 * DP.deltaTime * math.Sqrt(tolerance/err))
	return xt, nil
}
