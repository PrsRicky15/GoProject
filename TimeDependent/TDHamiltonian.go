package TimeDependent

import (
	"GoProject/Quantum"
	"GoProject/gridData"
)

type GaugeType int

const (
	Dipole GaugeType = iota
	Velocity
	Acceleration
)

type CanTimeSolverOp interface {
	ExpDtTo(At float64, In []float64, Out []float64) error
	ExpDtInPlace(At float64, InOut []float64) error
}

type TimeSolverOp interface {
	ExpIdtTo(Dt float64, In []complex128, Out []complex128) error
	ExpIdtInPlace(Dt float64, InOut []complex128) error
}

type TDCalculation struct {
	Hamil    *Quantum.HamiltonianOp
	gauge    GaugeType
	extField gridData.TDPotentialOp
}

func NewTDCalculation(hamil *Quantum.HamiltonianOp, gauge GaugeType, extField gridData.TDPotentialOp) *TDCalculation {
	return &TDCalculation{
		Hamil:    hamil,
		gauge:    gauge,
		extField: extField,
	}
}

func (td *TDCalculation) SetGauge(gauge GaugeType) {
	td.gauge = gauge
}

func (td *TDCalculation) Calculate(psi []complex128, t float64) []complex128 {
	switch td.gauge {
	case Dipole:
		return td.calculateDipoleGauge(psi, t)
	case Velocity:
		return td.calculateVelocityGauge(psi, t)
	case Acceleration:
		return td.calculateAccelerationGauge(psi, t)
	default:
		panic("Unknown gauge type")
	}
}

func (td *TDCalculation) calculateDipoleGauge(psi []complex128, t float64) []complex128 {
	result := make([]complex128, len(psi))
	return result
}

func (td *TDCalculation) calculateVelocityGauge(psi []complex128, t float64) []complex128 {
	result := make([]complex128, len(psi))
	return result
}

func (td *TDCalculation) calculateAccelerationGauge(psi []complex128, t float64) []complex128 {
	result := make([]complex128, len(psi))
	return result
}
