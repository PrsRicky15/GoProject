package OperatorAlgebra

import (
	"fmt"
	"gonum.org/v1/gonum/lapack"
	"gonum.org/v1/gonum/mat"
)

func ComplexGenDiagonalization(a mat.CMatrix, values []complex128, vectorsL mat.CMatrix, vectorsR mat.CMatrix) {
	ZGeev(lapack.LeftEVCompute, lapack.RightEVCompute, a, values, vectorsL.(*mat.CDense).RawCMatrix().Data,
		vectorsR.(*mat.CDense).RawCMatrix().Data)
}

func RealDiagonalizeLapack(eVecs mat.Matrix, evals []float64) error {
	dense, ok := eVecs.(*mat.Dense)
	if !ok {
		return fmt.Errorf("eigenvector output must be a dense matrix")
	}
	n, m := dense.Dims()
	if n != m || len(evals) != n {
		return fmt.Errorf("eigensystem dimensions do not match")
	}
	if !mat.EqualApprox(dense, dense.T(), 1e-12) {
		return fmt.Errorf("matrix must be symmetric")
	}
	sym := mat.NewSymDense(n, nil)
	for i := 0; i < n; i++ {
		for j := i; j < n; j++ {
			sym.SetSym(i, j, dense.At(i, j))
		}
	}
	var eig mat.EigenSym
	if !eig.Factorize(sym, true) {
		return fmt.Errorf("symmetric eigensolver failed")
	}
	copy(evals, eig.Values(nil))
	eig.VectorsTo(dense)
	return nil
}
