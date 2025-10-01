package Quantum

import (
	"GoProject/gridData"
	"math"

	"golang.org/x/exp/constraints"

	"gonum.org/v1/gonum/mat"
)

type Number interface {
	constraints.Signed | constraints.Float
}

type MatrixOp interface {
	Mat()
	GetMat() mat.Matrix
	Diagonalize() ([]float64, *mat.Dense, error)
}

type CanTimeSolverOp interface {
	ExpDt(At float64, In []float64) []float64
	ExpDtTo(At float64, In []float64, Out []float64)
	ExpDtInPlace(At float64, InOut []float64)
}

type TimeSolverOp interface {
	ExpDt(In []float64) []float64
	ExpDtTo(In []float64, Out []float64)
	ExpDtInPlace(InOut []float64)
}

type MomentumOp interface {
	MatrixOp
	TimeSolverOp
}

type KineticOp interface {
	MatrixOp
	TimeSolverOp
}

type CanMomentumOp interface {
	MatrixOp
	CanTimeSolverOp
}

type CanKineticOp interface {
	MatrixOp
	CanTimeSolverOp
}

type CanonicalMomDvrBasis struct {
	grid *gridData.RadGrid
	mass float64
}

type MomDvrBasis struct {
	grid *gridData.RadGrid
	mass float64
}

type CanonicalKeDvrBasis struct {
	grid *gridData.RadGrid
	mass float64
	kMat *mat.Dense
}

type KeDvrBasis struct {
	grid *gridData.RadGrid
	mass float64
	kMat *mat.Dense
}

func NewKeDVR(grid *gridData.RadGrid, mass float64) *KeDvrBasis {
	return &KeDvrBasis{
		grid: grid,
		mass: mass,
	}
}

func (k *KeDvrBasis) Mat() {
	dim := int(k.grid.NPoints())
	dx2 := k.grid.DeltaR() * k.grid.DeltaR()
	massDx2 := k.mass * dx2
	diagTerm := math.Pi * math.Pi / (6. * massDx2)
	invMassDx2 := 1.0 / massDx2

	k.kMat = mat.NewDense(dim, dim, nil)

	for i := 0; i < dim; i++ {
		k.kMat.Set(i, i, diagTerm)
	}

	for i := 1; i < dim; i++ {
		for j := 0; j < i; j++ {
			diff := i - j
			sign := float64(1 - 2*(diff&1))
			diffSq := float64(diff * diff)
			kEval := sign * invMassDx2 / diffSq
			k.kMat.Set(i, j, kEval)
			k.kMat.Set(j, i, kEval)
		}
	}
}

func (k *KeDvrBasis) GetMat() mat.Matrix { return k.kMat }
func (k *KeDvrBasis) Diagonalize() ([]float64, *mat.Dense, error) {
	return RealDiagonalize(k, int(k.grid.NPoints()))
}

func (k *KeDvrBasis) ExpDt(In []float64) []float64 {
	//TODO implement me
	panic("implement me")
}

func (k *KeDvrBasis) ExpDtTo(In []float64, Out []float64) {
	//TODO implement me
	panic("implement me")
}

func (k *KeDvrBasis) ExpDtInPlace(InOut []float64) {
	//TODO implement me
	panic("implement me")
}
