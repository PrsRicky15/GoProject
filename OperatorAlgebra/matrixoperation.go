package OperatorAlgebra

import (
	"gonum.org/v1/gonum/lapack"
	"gonum.org/v1/gonum/mat"
)

func ComplexGenDiagonalization(a mat.CMatrix, values []complex128, vectorsL mat.CMatrix, vectorsR mat.CMatrix) {
	ZGeev(lapack.LeftEVCompute, lapack.RightEVCompute, a, values, vectorsL.(*mat.CDense).RawCMatrix().Data,
		vectorsR.(*mat.CDense).RawCMatrix().Data)
}

func RealDiagonalizeLapack(eVecs mat.Matrix, evals []float64) error {
	err := dsyev(eVecs.(*mat.Dense).RawMatrix().Rows, eVecs.(*mat.SymDense).RawSymmetric(),
		evals, eVecs.(*mat.Dense).RawMatrix())
	if err != nil {
		return err
	}
	return nil
}
