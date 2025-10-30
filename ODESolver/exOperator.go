package EquationSolver

type ExponentialOperation interface {
	ExpDt(Dt float64, In []float64) []float64
	ExpIdt(Dt float64, In []complex128) []complex128
}
