package EquationSolver_test

import (
	EquationSolver "GoProject/ODESolver"
	"fmt"
)

func ExampleNewVectorImplicit() {
	// Every derivative evaluation receives both position and velocity.
	rhs := func(t float64, y, dy []float64) error {
		dy[0] = y[1]
		dy[1] = -y[0]
		return nil
	}
	solver, err := EquationSolver.NewVectorImplicit(rhs, EquationSolver.ImplicitMidpointMethod, EquationSolver.ImplicitOptions{})
	if err != nil {
		panic(err)
	}
	state := []float64{1, 0}
	if err = solver.Step(state, 0, .1); err != nil {
		panic(err)
	}
	fmt.Printf("%.6f %.6f\n", state[0], state[1])
	// Output: 0.995012 -0.099751
}

func ExampleNewComplexLinearImplicit() {
	// i*dpsi/dt = H*psi with H=[[0,1],[1,0]], in atomic units.
	apply := func(t float64, psi, dpsi []complex128) error {
		dpsi[0] = -1i * psi[1]
		dpsi[1] = -1i * psi[0]
		return nil
	}
	solver, err := EquationSolver.NewComplexLinearImplicit(2, apply, EquationSolver.ImplicitMidpointMethod, EquationSolver.KrylovOptions{})
	if err != nil {
		panic(err)
	}
	psi := []complex128{1, 0}
	if err = solver.Step(psi, 0, .1); err != nil {
		panic(err)
	}
	fmt.Printf("%.6f %.6f\n", real(psi[0]), imag(psi[1]))
	// Output: 0.995012 -0.099751
}
