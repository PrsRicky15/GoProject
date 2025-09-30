package Quantum

import (
	"GoProject/gridData"

	"GoProject/matrix"

	"math"

	"golang.org/x/exp/constraints"

	"gonum.org/v1/gonum/mat"
)

type Number interface {
	constraints.Signed | constraints.Float
}

type MomDvrBasis struct {
	grid *gridData.RadGrid
	mass float64
}

type KeDvrBasis struct {
	grid *gridData.RadGrid
	mass float64
}

func (k *KeDvrBasis) CanonicalOpEvaluate(At float64) mat.Matrix {
	//TODO implement me
	panic("implement me")
}

func NewKeDVR(grid *gridData.RadGrid, mass float64) *KeDvrBasis {
	return &KeDvrBasis{
		grid: grid,
		mass: mass,
	}
}

func (k *KeDvrBasis) Evaluate() mat.Matrix {
	dim := int(k.grid.NPoints())
	dx2 := k.grid.DeltaR() * k.grid.DeltaR()
	massDx2 := k.mass * dx2
	diagTerm := math.Pi * math.Pi / (6. * massDx2)
	invMassDx2 := 1.0 / massDx2

	val := mat.NewDense(dim, dim, nil)

	for i := 0; i < dim; i++ {
		val.Set(i, i, diagTerm)
	}

	for i := 1; i < dim; i++ {
		for j := 0; j < i; j++ {
			diff := i - j
			sign := float64(1 - 2*(diff&1))
			diffSq := float64(diff * diff)
			kEval := sign * invMassDx2 / diffSq
			val.Set(i, j, kEval)
			val.Set(j, i, kEval)
		}
	}
	return val
}

func (k *KeDvrBasis) Diagonalize() ([]float64, *mat.Dense, error) {
	return matrix.RealDiagonalize(k, int(k.grid.NPoints()))
}
