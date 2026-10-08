package Quantum

import (
	"GoProject/OperatorAlgebra"

	"gonum.org/v1/gonum/mat"
)

// Hamiltonian application does not require a time-propagation interface.
type kineticEnergy interface {
	applyReal(in, out []float64)
	applyComplex(in, out []complex128)
	matrix() *mat.Dense // returns an independent matrix
	close()
}

type dvrKinetic struct{ basis *OperatorAlgebra.KeDvrBasis }

func (k *dvrKinetic) applyReal(in, out []float64) {
	m := k.basis.GetMat().RawMatrix()
	for i := range out {
		row := m.Data[i*m.Stride : i*m.Stride+m.Cols]
		var sum float64
		for j, coefficient := range row {
			sum += coefficient * in[j]
		}
		out[i] = sum
	}
}

func (k *dvrKinetic) applyComplex(in, out []complex128) {
	m := k.basis.GetMat().RawMatrix()
	for i := range out {
		row := m.Data[i*m.Stride : i*m.Stride+m.Cols]
		var sum complex128
		for j, coefficient := range row {
			sum += complex(coefficient, 0) * in[j]
		}
		out[i] = sum
	}
}

func (k *dvrKinetic) matrix() *mat.Dense { return mat.DenseCopyOf(k.basis.GetMat()) }
func (k *dvrKinetic) close()             {}

type fourierKinetic struct {
	basis *OperatorAlgebra.FourierBasis
	work  []complex128
}

func (k *fourierKinetic) applyReal(in, out []float64) {
	for i, value := range in {
		k.work[i] = complex(value, 0)
	}
	k.basis.LaplacianOpInPlace(k.work)
	for i, value := range k.work {
		out[i] = real(value)
	}
}

func (k *fourierKinetic) applyComplex(in, out []complex128) {
	k.basis.LaplacianOp(in, out)
}

func (k *fourierKinetic) matrix() *mat.Dense {
	n := len(k.work)
	m := mat.NewDense(n, n, nil)
	unit := make([]float64, n)
	column := make([]float64, n)
	for j := range unit {
		unit[j] = 1
		k.applyReal(unit, column)
		m.SetCol(j, column)
		unit[j] = 0
	}
	return m
}

func (k *fourierKinetic) close() { k.basis.Clean() }
