package Quantum

import (
	"gonum.org/v1/gonum/lapack"
	"gonum.org/v1/gonum/mat"
)

func ComplexGenDiagonalization(a mat.CMatrix, values []complex128, vectorsL mat.CMatrix, vectorsR mat.CMatrix) {
	ZGeev(lapack.LeftEVCompute, lapack.RightEVCompute, a, values, vectorsL.(*mat.CDense).RawCMatrix().Data,
		vectorsR.(*mat.CDense).RawCMatrix().Data)
}

func RealDiagonalizeLapack(Op MatrixOp, evals []float64, eVecs mat.Matrix) error {
	matrix := Op.GetMat()
	err := dsyev(matrix.(*mat.Dense).RawMatrix().Rows, matrix.(*mat.SymDense).RawSymmetric(),
		evals, eVecs.(*mat.Dense).RawMatrix())
	if err != nil {
		return err
	}
	return nil
}
