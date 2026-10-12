package EquationSolver_test

import (
	EquationSolver "GoProject/ODESolver"
	"fmt"
)

type decay struct{}

func (decay) EvaluateAt(x, t float64) float64 { return -1000 * x }
func (f decay) EvaluateOnRGrid(x []float64, t float64) []float64 {
	out := make([]float64, len(x))
	f.EvaluateOnRGridInPlace(x, out, t)
	return out
}
func (decay) EvaluateOnRGridInPlace(x, out []float64, t float64) {
	for i, v := range x {
		out[i] = -1000 * v
	}
}

func ExampleEulerImplicitNewton() {
	solver := &EquationSolver.EulerImplicitNewton{
		TdFunc: decay{}, DeltaT: .1,
		Options: EquationSolver.ImplicitOptions{
			Derivative: func(x, t float64) float64 { return -1000 },
		},
	}
	x, err := solver.NextStep(1, 0)
	if err != nil {
		panic(err)
	}
	fmt.Printf("x = %.6f\n", x)
	// Output: x = 0.009901
}
