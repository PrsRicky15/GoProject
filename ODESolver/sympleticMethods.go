package EquationSolver

import "GoProject/gridData"

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
