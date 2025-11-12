package gridData

import (
	"fmt"
	"math"
)

// TDPotentialOp General interface for the evaluating the potential on a grid
// type TDPotentialOp interface {
// 	EvaluateAt(x float64, t float64) float64
// 	EvaluateOnRGrid(x []float64, t float64) []float64
// 	EvaluateOnRGridInPlace(x, res []float64, t float64)
// }

// TDPotentialOp General interface for evaluating the potential on a grid

type TDPotentialOp[T VarType] interface {
	PotentialOp[T]

	EvaluateAtT(x T, t float64) T
	EvaluateOnGridT(x []T, t float64) []T
	ForceAtT(x T, t float64) T
	ForceOnGridT(x []T, t float64) []T
}

// V(x,t) for array x, at time t

func onGridT[T VarType](f func(T, float64) T, x []T, t float64) []T {
	results := make([]T, len(x))
	for i, val := range x {
		results[i] = f(val, t)
	}
	return results
}

// Experimental: adapter so that any time independent argument

func switchPotential[T VarType](p any) TDPotentialOp[T] {
	switch potential := p.(type) {

	case TDPotentialOp[T]:
		return potential

	case PotentialOp[T]:
		return timeIndependentAdapter[T]{base: potential}

	default:
		panic("Make potential great again")
	}
}

// PotentialOp[T] will run on PotentialOpT[T] by ignoring the time argument

type timeIndependentAdapter[T VarType] struct {
	base PotentialOp[T]
}

func adaptTimeIndependent[T VarType](op PotentialOp[T]) TDPotentialOp[T] {
	return timeIndependentAdapter[T]{base: op}
}

func (a timeIndependentAdapter[T]) EvaluateAt(x T) T {
	return a.base.EvaluateAt(x)
}
func (a timeIndependentAdapter[T]) EvaluateOnGrid(x []T) []T {
	return a.base.EvaluateOnGrid(x)
}
func (a timeIndependentAdapter[T]) ForceAt(x T) T {
	return a.base.ForceAt(x)
}
func (a timeIndependentAdapter[T]) ForceOnGrid(x []T) []T {
	return a.base.ForceOnGrid(x)
}

func (a timeIndependentAdapter[T]) EvaluateAtT(x T, t float64) T { return a.base.EvaluateAt(x) }
func (a timeIndependentAdapter[T]) ForceAtT(x T, t float64) T    { return a.base.ForceAt(x) }
func (a timeIndependentAdapter[T]) EvaluateOnGridT(x []T, t float64) []T {
	return a.base.EvaluateOnGrid(x)
}
func (a timeIndependentAdapter[T]) ForceOnGridT(x []T, t float64) []T { return a.base.ForceOnGrid(x) }

// V(x,t) = Morse(x) + E_0*cos(omega t) * x
type DrivenMorse[T VarType] struct {
	Morse    Morse[T]
	FieldAmp float64
	Freq     float64
}

func (dm DrivenMorse[T]) String() string {
	return fmt.Sprintf(" Morse + E_0*cos(omega t)*x, where, FieldAmp: %g, Freq: %g", dm.FieldAmp, dm.Freq)
}

// Time-independent adaptations
func (dm DrivenMorse[T]) EvaluateAt(x T) T { return dm.Morse.EvaluateAt(x) }
func (dm DrivenMorse[T]) ForceAt(x T) T    { return dm.Morse.ForceAt(x) }
func (dm DrivenMorse[T]) EvaluateOnGrid(x []T) []T {
	return dm.Morse.EvaluateOnGrid(x)
}
func (dm DrivenMorse[T]) ForceOnGrid(x []T) []T {
	return dm.Morse.ForceOnGrid(x)
}

// Time-dependent implementations
func (dm DrivenMorse[T]) EvaluateAtT(x T, t float64) T {
	base := dm.Morse.EvaluateAt(x)
	oscill := dm.FieldAmp * math.Cos(dm.Freq*t)
	switch any(x).(type) {
	case float64:
		return any(any(base).(float64) + oscill*any(x).(float64)).(T)
	case complex128:
		xx := any(x).(complex128)
		add := complex(oscill*real(xx), oscill*imag(xx))
		return any(any(base).(complex128) + add).(T)
	default:
		panic("unsupported type")
	}
}

func (dm DrivenMorse[T]) ForceAtT(x T, t float64) T {
	base := dm.Morse.ForceAt(x)
	oscill := dm.FieldAmp * math.Cos(dm.Freq*t)
	switch any(x).(type) {
	case float64:
		return any(any(base).(float64) + oscill).(T)
	case complex128:
		return any(any(base).(complex128) + complex(oscill, 0)).(T)
	default:
		panic("unsupported type")
	}
}

func (dm DrivenMorse[T]) EvaluateOnGridT(x []T, t float64) []T { return onGridT(dm.EvaluateAtT, x, t) }
func (dm DrivenMorse[T]) ForceOnGridT(x []T, t float64) []T    { return onGridT(dm.ForceAtT, x, t) }
