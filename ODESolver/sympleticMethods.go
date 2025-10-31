package EquationSolver

import (
	"GoProject/gridData"
	"math"
)

type MDodeSolver interface {
	NextStep(forceAtT, x, v, time float64) float64
	NextStepOnGrid(forceAtT, x, v []float64, time float64)
	Name() string
}

// StromerVerlet method
type StromerVerlet struct {
	forceFunc gridData.PotentialOp[float64]
	dt        float64
	dt2       float64
	prevX     float64
	prevXGrid []float64
}

func (vi *StromerVerlet) Name() string {
	return "StromerVerlet Integrator"
}

func (vi *StromerVerlet) NewDef(poten gridData.PotentialOp[float64], dt float64) *StromerVerlet {
	return &StromerVerlet{
		forceFunc: poten,
		dt:        dt,
		dt2:       dt * dt,
	}
}

func (vi *StromerVerlet) Redefine(poten gridData.PotentialOp[float64], dt float64) {
	vi.forceFunc = poten
	vi.dt = dt
	vi.dt2 = dt * dt
}

func (vi *StromerVerlet) Initiate(x, v float64) {
	vi.prevX = x + v*vi.dt + 0.5*vi.dt2*vi.forceFunc.ForceAt(x)
}

func (vi *StromerVerlet) InitiateGrid(x, v []float64) {
	if vi.prevXGrid == nil {
		vi.prevXGrid = make([]float64, len(x))
		for i := range x {
			vi.prevXGrid[i] = x[i] + v[i]*vi.dt +
				0.5*vi.dt2*vi.forceFunc.ForceAt(x[i])
		}
	}
}

func (vi *StromerVerlet) NextStep(x float64) float64 {
	nextX := 2*x - vi.prevX + vi.forceFunc.ForceAt(x)*vi.dt2
	vi.prevX = x
	return nextX
}

func (vi *StromerVerlet) NextStepOnGrid(x []float64) {
	for i := range x {
		nextX := 2*x[i] - vi.prevXGrid[i] + vi.forceFunc.ForceAt(x[i])*vi.dt2
		vi.prevXGrid[i] = x[i]
		x[i] = nextX
	}
}

type VelocityVerlet struct {
	forceFunc gridData.PotentialOp[float64]
	dt        float64
	halfdt2   float64
	halfdt    float64

	acc []float64
}

func (vv *VelocityVerlet) Name() string {
	return "Velocity-Verlet Integrator"
}

func (vv *VelocityVerlet) NewDef(Poten gridData.PotentialOp[float64], dt float64) *VelocityVerlet {
	return &VelocityVerlet{
		forceFunc: Poten,
		dt:        dt,
		halfdt:    dt / 2.,
		halfdt2:   dt * dt / 2.,
	}
}

func (vv *VelocityVerlet) Redefine(Poten gridData.PotentialOp[float64], dt float64) {
	vv.forceFunc = Poten
	vv.dt = dt
	vv.halfdt = dt / 2.
	vv.halfdt2 = dt * dt / 2.
}

func (vv *VelocityVerlet) NextStep(x, v float64) (float64, float64) {
	accPre := vv.forceFunc.ForceAt(x)
	xNext := x + v*vv.dt + vv.halfdt2*accPre
	accNext := vv.forceFunc.ForceAt(xNext)
	vNext := v + vv.halfdt*(accPre+accNext)
	return xNext, vNext
}

func (vv *VelocityVerlet) NextStepOnGrid(x, v []float64) {

	if vv.acc == nil {
		vv.acc = make([]float64, len(v))
	}

	for i := range x {
		vv.acc[i] = vv.forceFunc.ForceAt(x[i])
		x[i] = x[i] + v[i]*vv.dt + vv.halfdt2*vv.acc[i]
		v[i] = v[i] + vv.halfdt*(vv.acc[i]+vv.forceFunc.ForceAt(x[i]))
	}
}

type LeapFrog struct {
	forceFunc gridData.PotentialOp[float64]
	dt        float64
	halfDt    float64

	vHalf []float64
}

func (lf *LeapFrog) Name() string {
	return "Leap-Frog Integrator"
}

func (lf *LeapFrog) NewDef(Poten gridData.PotentialOp[float64], dt float64) *LeapFrog {
	return &LeapFrog{
		forceFunc: Poten,
		dt:        dt,
		halfDt:    dt / 2.,
	}
}

func (lf *LeapFrog) Redefine(Poten gridData.PotentialOp[float64], dt float64) {
	lf.forceFunc = Poten
	lf.dt = dt
	lf.halfDt = dt / 2.
}

func (lf *LeapFrog) NextStep(x, v float64) (float64, float64) {
	vHalf := v + lf.forceFunc.ForceAt(x)*lf.halfDt
	xNext := x + vHalf*lf.dt
	force := lf.forceFunc.ForceAt(xNext)
	vNext := vHalf + lf.halfDt*force
	return xNext, vNext
}

func (lf *LeapFrog) NextStepOnGrid(x, v []float64, time float64) {
	if lf.vHalf == nil {
		lf.vHalf = make([]float64, len(v))
	}
	for i := range v {
		lf.vHalf[i] = v[i] + lf.forceFunc.ForceAt(x[i])*lf.halfDt
		x[i] = x[i] + lf.vHalf[i]*lf.dt
		v[i] = lf.vHalf[i] + lf.halfDt*lf.forceFunc.ForceAt(x[i])
	}
}

type Yoshida struct {
	fnc gridData.PotentialOp[float64]
	dt  float64

	w0       float64
	w1       float64
	wAvg     float64
	halfW1   float64
	halfWavg float64
}

func (yo *Yoshida) Name() string {
	return "Leap-Frog Integrator"
}

func (yo *Yoshida) NewDef(Poten gridData.PotentialOp[float64], dt float64) *Yoshida {
	w0 := -math.Cbrt(2.) / (2 - math.Cbrt(2.))
	w1 := 1. / (2 - math.Cbrt(2.))
	wAvg := (w0 + w1) / 2.

	return &Yoshida{
		fnc:      Poten,
		dt:       dt,
		w0:       w0,
		w1:       w1,
		wAvg:     wAvg,
		halfW1:   0.5 * w1,
		halfWavg: 0.5 * wAvg,
	}
}

func (yo *Yoshida) Redefine(Poten gridData.PotentialOp[float64], dt float64) {
	yo.fnc = Poten
	yo.dt = dt
}

func (yo *Yoshida) NextStep(x, v float64) (float64, float64) {
	xMid := x + yo.halfW1*v*yo.dt
	vMid := v + yo.w1*yo.fnc.ForceAt(xMid)*yo.dt

	xMid = xMid + yo.halfWavg*vMid*yo.dt
	vMid = vMid + yo.w0*yo.dt*yo.fnc.ForceAt(xMid)

	xMid = xMid + yo.halfWavg*vMid*yo.dt
	vMid = vMid + yo.w1*yo.dt*yo.fnc.ForceAt(xMid)

	return xMid + yo.halfW1*vMid*yo.dt, vMid
}

func (yo *Yoshida) NextStepOnGrid(x, v []float64) {
	for i := range x {
		x[i] = x[i] + yo.halfW1*yo.dt*v[i]
		v[i] = v[i] + yo.w1*yo.dt*yo.fnc.ForceAt(x[i])
	}

	for i := range x {
		x[i] = x[i] + yo.halfWavg*yo.dt*v[i]
		v[i] = v[i] + yo.w0*yo.dt*yo.fnc.ForceAt(x[i])
	}

	for i := range x {
		x[i] = x[i] + yo.halfWavg*yo.dt*v[i]
		v[i] = v[i] + yo.w1*yo.dt*yo.fnc.ForceAt(x[i])
		x[i] = x[i] + yo.halfW1*yo.dt*v[i]
	}
}
