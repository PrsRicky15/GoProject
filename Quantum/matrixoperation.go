package Quantum

import (
	"fmt"

	"gonum.org/v1/gonum/mat"
)

func RealDiagonalize(Op EvaluateOp, n int) ([]float64, *mat.Dense, error) {
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
