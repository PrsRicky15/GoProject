package gridData

// TDPotentialOp General interface for the evaluating the potential on a grid
// type TDPotentialOp interface {
// 	EvaluateAt(x float64, t float64) float64
// 	EvaluateOnRGrid(x []float64, t float64) []float64
// 	EvaluateOnRGridInPlace(x, res []float64, t float64)
// }

import (
	"fmt"
	"math"
	"math/cmplx"
)

type VarType interface {
	float64 | complex128
}

// TDPotentialOp General interface for evaluating the potential on a grid

type TDPotentialOp interface {
	EvaluateAt(x T) T
	EvaluateOnGrid(x []T) []T
	ForceAt(x T) T
	ForceOnGrid(x []T) []T

	EvaluateAtT(x T, t float64) T
	EvaluateOnGrid(x []T, t float64) []T
	ForceAtT(x T, t float64) T
	ForceOnGrid(x []T, t float64) []T
}

// V(x,t) for array x, at time t

func onGridT[T VarType](f func(T, float64) T, x []T, t float64) []T {
	results :- make([]T, len(x))
	for i, val := range x {
		results[i] = f(val, t)
	}
	return results
}

// Experimental: adapter so that any time independent argument

func switchPotential[T VarType](v any) TDPotentialOp[T] {
	switch potential := v.(type) {

	case TDPotentialOp[T]:
		return v

	case PotentialOp[T]:
		return timeIndependentAdapter[T]{base: v}

	default:
		panic("Make potential great again")
	}
}

// PotentialOp[T] will run on PotentialOpT[T] by ignoring the time argument

type timeIndependentAdapter[T VarType] struct {
	base PontentialOp[T]
}

func adaptTimeIndependent[T VarType](op PotentialOp[T]) PotentialOpT[T] {
	return timeIndependentAdapter[T]{base: op}
}

func (a timeIndependentAdapter[T]) EvaluateAt(x T) T{
	return a.base.EvaluateAt(x)
}
func (a timeIndependentAdapter[T]) EvaluateOnGrid(x []T) []T{
	return a.base.EvaluateOnGrid(x)
}
func (a timeIndependentAdapter[T]) ForceAt(x T) T{
	return a.base.ForceAt(x)
}
func (a timeIndependentAdapter[T]) ForceOnGrid(x []T) []T {
	return a.base.ForceOnGrid(x)
}


func (a timeIndependentAdapter[T]) EvaluateAtT(x T, t float64) T { return a.base.EvaluateAt(x) }
func (a timeIndependentAdapter[T]) ForceAtT(x T, t float64) T { return a.base.ForceAt(x) }
func (a timeIndependentAdapter[T]) EvaluateOnGridT(x []T, t float64) []T { return a.base.EvaluateOnGrid(x) }
func (a timeIndependentAdapter[T]) ForceOnGridT(x []T, t float64) []T { return a.base.ForceOnGrid(x) }

