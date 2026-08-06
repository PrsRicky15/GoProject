package OperatorAlgebra

import "gonum.org/v1/gonum/mat"

type MatrixOp interface {
	ExpDtTo(Dt float64, In []float64, Out []float64) error
	ExpDtInPlace(Dt float64, InOut []float64) error
	RealDiagonalize() (eigenvalues []float64, eigenvectors *mat.Dense, err error)
}
