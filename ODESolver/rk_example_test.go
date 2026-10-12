package EquationSolver_test

import (
	EquationSolver "GoProject/ODESolver"
	"fmt"
)

func ExampleNewComplexAdaptiveRK() {
	// H = [[0,1],[1,0]], hbar=1: every component depends on another.
	rhs := func(t float64, psi, dpsi []complex128) error {
		dpsi[0] = -1i * psi[1]
		dpsi[1] = -1i * psi[0]
		return nil
	}
	solver, err := EquationSolver.NewComplexAdaptiveRK(rhs, EquationSolver.RKDormandPrince54,
		EquationSolver.AdaptiveOptions{AbsoluteTolerance: 1e-12, RelativeTolerance: 1e-10})
	if err != nil {
		panic(err)
	}
	psi := []complex128{1, 0}
	result, err := solver.Integrate(psi, 0, 1)
	if err != nil {
		panic(err)
	}
	fmt.Printf("time %.1f; psi = [%.6f, %.6fi]\n", result.Time, real(psi[0]), imag(psi[1]))
	// Output: time 1.0; psi = [0.540302, -0.841471i]
}
