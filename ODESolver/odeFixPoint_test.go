package EquationSolver

import (
	"math"
	"testing"
	"GoProject/gridData"
)

type ImplicitFixPointSolvers struct{}

func (s ImplicitFixPointSolvers) getSolvers() []SolverTestCase {
    return []SolverTestCase{
        {
            Name: "Euler Implicit",
            Factory: func(dt float64, op gridData.TDPotentialOp) ScalarODESolver {
                return (&EulerImplicitFixPoint{}).NewDefine(dt, op)
            },
            Order:        1,
            SupportsGrid: true,
        },
	}
}

func TestImplicitSolvers_SingleStepErrorBound(t *testing.T) {
    mockFunc := &MockPotential{}
    x0 := 1.0
    t0 := 0.0
    dt := 0.01
    exact := x0 * math.Exp(-dt)

    s := ImplicitFixPointSolvers{}
    for _, tc := range s.getSolvers() {
        t.Run(tc.Name+"_SingleStepError", func(t *testing.T) {
            solver := tc.Factory(dt, mockFunc)
            x1, err := solver.NextStep(x0, t0)
            if err != nil {
                t.Fatalf("NextStep failed: %v", err)
            }
            diff := math.Abs(x1 - exact)

            // loose bounds by nominal order to catch gross regressions
            bound := 5e-4 // order 1
            if tc.Order >= 2 {
                bound = 5e-5
            }
            if tc.Order >= 3 {
                bound = 5e-6
            }
            if tc.Order >= 4 {
                bound = 5e-7
            }
            if diff > bound {
                t.Errorf("single-step error too large: got=%e bound=%e", diff, bound)
            }
        })
    }
}