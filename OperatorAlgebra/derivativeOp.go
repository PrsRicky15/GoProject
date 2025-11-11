package OperatorAlgebra

import "gonum.org/v1/gonum/mat"

type MatrixOp interface {
	RealDiagonalize() (eigenvalues []float64, eigenvectors *mat.Dense, err error)
}

type CanTimeSolverOp interface {
	ExpDtTo(At float64, In []float64, Out []float64) error
	ExpDtInPlace(At float64, InOut []float64) error
}

type TimeSolverOp interface {
	ExpDtTo(Dt float64, In []float64, Out []float64) error
	ExpDtInPlace(Dt float64, InOut []float64) error

	ExpIdtTo(Dt float64, In []complex128, Out []complex128) error
	ExpIdtInPlace(Dt float64, InOut []complex128) error
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
