package Quantum

import (
	"fmt"

	"math"

	"golang.org/x/exp/constraints"

	"gonum.org/v1/gonum/mat"

	"GoProject/gridData"
)

type Number interface {
	constraints.Signed | constraints.Float
}

func Abs[T Number](x T) T {
	if x < 0 {
		return -x
	}
	return x
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
	grid gridData.RadGrid
	mass float64
}

func (k *KeDVR) CanonicalEvaluate(At float64) mat.Matrix {
	//TODO implement me
	panic("implement me")
}

func keToeplitz(dim int, dx2 float64) []float64 {
	vals := make([]float64, dim)
	for i := 0; i < dim; i++ {
		vals[i] =
	}
	return vals
}

func (k *KeDVR) Evaluate() mat.Matrix {
	ngrid := int(k.grid.NPoints())
	dx2 := k.grid.DeltaR() * k.grid.DeltaR()
	diagTerm := math.Pow(math.Pi, 2) / (6. * dx2 * k.mass)

	keT := keToeplitz(ngrid, dx2)
	val := mat.NewDense(ngrid, ngrid, nil)
	for i := 0; i < int(k.grid.NPoints()); i++ {
		val.Set(i, i, diagTerm)
		for j := 0; j < i-1; j++ {
			val.Set(i, j, keT[Abs(i-j)])
			val.Set(j, i, keT[Abs(i-j)])
		}
	}

	return val
}

func (k *KeDVR) Diagonalize() ([]float64, *mat.Dense, error) {
	return MatrixDiagonalize(k, int(k.grid.NPoints()))
}