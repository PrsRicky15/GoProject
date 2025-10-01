package Quantum

import (
	"fmt"

	"gonum.org/v1/gonum/blas"
	"gonum.org/v1/gonum/blas/blas64"
	"gonum.org/v1/gonum/lapack"
	"gonum.org/v1/gonum/lapack/lapack64"
	"gonum.org/v1/gonum/mat"
)

func RealDiagonalize(Op MatrixOp, evals []float64, eVecs *mat.Dense) error {
	matrix := mat.Matrix(Op.GetMat())
	symMat, ok := matrix.(*mat.SymDense)
	if !ok {
		return fmt.Errorf("matrix must be *mat.SymDense")
	}

	var eig mat.EigenSym
	if ok := eig.Factorize(symMat, true); !ok {
		return fmt.Errorf("failed to compute eigenvalues and eigenvectors")
	}

	eig.Values(evals)
	eig.VectorsTo(eVecs)
	return nil
}

func RealDiagonalizeLapack(Op MatrixOp, evals []float64, eVecs *mat.Dense) error {
	matrix := Op.GetMat()
	_, n := matrix.Dims()

	data := matrix.RawMatrix().Data

	a := blas64.Symmetric{
		N:      n,
		Uplo:   blas.Upper,
		Data:   data,
		Stride: n,
	}

	w := make([]float64, n)
	work := make([]float64, 1)
	lwork := -1

	ok := lapack64.Syev(lapack.EVCompute, a, w, work, lwork)
	if !ok {
		return fmt.Errorf("LAPACK Syev query failed")
	}

	optimalLwork := int(work[0])
	work = make([]float64, optimalLwork)
	lwork = optimalLwork

	ok = lapack64.Syev(lapack.EVCompute, a, w, work, lwork)
	if !ok {
		return fmt.Errorf("LAPACK Syev computation failed")
	}

	// Copy results
	copy(evals, w)
	eVecs.SetRawMatrix(blas64.General{
		Rows:   n,
		Cols:   n,
		Stride: n,
		Data:   a.Data,
	})

	return nil
}
