package OperatorAlgebra

import (
	"gonum.org/v1/gonum/mat"
)

type MatrixOp interface {
	Mat()
	GetMat() mat.Matrix
	GetCMat() mat.Matrix
	Diagonalize() ([]float64, mat.Matrix, error)
}

type CanTimeSolverOp interface {
	ExpDtTo(At float64, In []float64, Out []float64)
	ExpDtInPlace(At float64, InOut []float64)
}

type TimeSolverOp interface {
	ExpDtTo(Dt float64, In []float64, Out []float64)
	ExpDtInPlace(Dt float64, InOut []float64)

	ExpIdtTo(Dt float64, In []complex128, Out []complex128)
	ExpIdtInPlace(Dt float64, InOut []complex128)
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
