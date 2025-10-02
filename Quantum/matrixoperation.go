package Quantum

import (
	"fmt"

	"gonum.org/v1/gonum/blas/blas64"
	"gonum.org/v1/gonum/lapack"
	"gonum.org/v1/gonum/lapack/lapack64"
	"gonum.org/v1/gonum/mat"
)

type zomplex interface {
	Zgeev(jobvl lapack.LeftEVJob, jobvr lapack.RightEVJob, n int, a []complex128, lda int, w []complex128,
		vl []complex128, ldvl int, vr []complex128, ldvr int, work []complex128, lwork int,
		rwork []float64) (first int)
}

var zomplex64 zomplex

func ZGeev(jobvl lapack.LeftEVJob, jobvr lapack.RightEVJob, a mat.CMatrix, w []complex128, vl,
	vr []complex128) (first int) {

	n := a.(*mat.CDense).RawCMatrix().Rows
	if a.(*mat.CDense).RawCMatrix().Cols != n {
		panic("lapack64: matrix not square")
	}

	lwork := 2 * n
	work := make([]complex128, lwork)
	rwork := make([]float64, lwork)

	return zomplex64.Zgeev(jobvl, jobvr, n, a.(*mat.CDense).RawCMatrix().Data,
		max(1, a.(*mat.CDense).RawCMatrix().Stride), w, vl, n, vr, n,
		work, lwork, rwork)
}

func ComplexGenDiagonalization(a mat.CMatrix, values []complex128, vectorsL mat.CMatrix, vectorsR mat.CMatrix) {
	ZGeev(lapack.LeftEVCompute, lapack.RightEVCompute, a, values, vectorsL.(*mat.CDense).RawCMatrix().Data,
		vectorsR.(*mat.CDense).RawCMatrix().Data)
}

func RealDiagonalizeLapack(Op MatrixOp, evals []float64, eVecs mat.Matrix) error {
	matrix := Op.GetMat()
	n := matrix.(*mat.Dense).RawMatrix().Rows

	w := make([]float64, n)
	work := make([]float64, 1)
	lwork := -1

	ok := lapack64.Syev(lapack.EVCompute, matrix.(*mat.SymDense).RawSymmetric(), w, work, lwork)
	if !ok {
		return fmt.Errorf("LAPACK Syev query failed")
	}

	optimalLwork := int(work[0])
	work = make([]float64, optimalLwork)
	lwork = optimalLwork

	ok = lapack64.Syev(lapack.EVCompute, matrix.(*mat.SymDense).RawSymmetric(), evals, work, lwork)
	if !ok {
		return fmt.Errorf("LAPACK Syev computation failed")
	}

	eVecs.(*mat.Dense).SetRawMatrix(blas64.General{
		Rows:   n,
		Cols:   n,
		Stride: n,
		Data:   matrix.(*mat.SymDense).RawSymmetric().Data,
	})

	return nil
}
