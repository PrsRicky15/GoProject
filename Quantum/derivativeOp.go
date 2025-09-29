package Quantum

import (
	"GoProject/gridData"
	"fmt"

	"math"

	"golang.org/x/exp/constraints"

	"gonum.org/v1/gonum/mat"
)

type Number interface {
	constraints.Signed | constraints.Float
}

type EvaluateOp interface {
	Evaluate() mat.Matrix
	CanonicalEvaluate(At float64) mat.Matrix
	Diagonalize() ([]float64, *mat.Dense, error)
}

type MomentumOp interface {
	EvaluateOp
}

type KineticOp interface {
	EvaluateOp
}

func MatrixDiagonalize(Op EvaluateOp, n int) ([]float64, *mat.Dense, error) {
	matrix := Op.Evaluate()

	sym, ok := matrix.(*mat.SymDense)
	if !ok {
		sym = mat.NewSymDense(n, nil)
		for i := 0; i < n; i++ {
			for j := i; j < n; j++ {
				sym.SetSym(i, j, matrix.At(i, j))
			}
		}
	}

	var eig mat.EigenSym
	if ok := eig.Factorize(sym, true); !ok {
		return nil, nil, fmt.Errorf("failed to compute eigenvalues and eigenvectors")
	}

	eigenvalues := make([]float64, n)
	eig.Values(eigenvalues)

	eigenvectors := mat.NewDense(n, n, nil)
	eig.VectorsTo(eigenvectors)

	return eigenvalues, eigenvectors, nil
}

type KeDVR struct {
	grid *gridData.RadGrid
	mass float64
}

func NewKeDVR(grid *gridData.RadGrid, mass float64) *KeDVR {
	return &KeDVR{
		grid: grid,
		mass: mass,
	}
}

func (k *KeDVR) CanonicalEvaluate(float64) mat.Matrix {
	//TODO implement me
	panic("implement me")
}

func (k *KeDVR) Evaluate() mat.Matrix {
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

func (k *KeDVR) Diagonalize() ([]float64, *mat.Dense, error) {
	return MatrixDiagonalize(k, int(k.grid.NPoints()))
}
