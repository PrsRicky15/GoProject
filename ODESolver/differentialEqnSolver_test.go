package EquationSolver

import (
	"fmt"
	"math"
	"testing"

	"GoProject/gridData"
)

// TestEulerExplicit_NextStep tests the scalar stepping function.
func TestEulerExplicit_NextStep(t *testing.T) {
	mockFunc := &MockPotential{}

	// Problem: dx/dt = -x, x(0) = 1.0
	// Euler explicit: x_new = x + dt * (-x) = x * (1 - dt)
	initialX := 1.0
	dt := 0.1
	expectedX := initialX * (1.0 - dt)

	solver := &EulerExplicit{}
	solver = solver.NewDefine(dt, mockFunc)

	// Step 1
	nextX, err := solver.NextStep(initialX, 0.0)
	if err != nil {
		t.Fatalf("NextStep failed: %v", err)
	}

	if math.Abs(nextX-expectedX) > 1e-9 {
		t.Errorf("NextStep wrong result. Got %f, expected %f", nextX, expectedX)
	}

	// Step 2
	expectedX2 := expectedX * (1.0 - dt)
	nextX2, err := solver.NextStep(nextX, dt)
	if err != nil {
		t.Fatalf("NextStep failed on second step: %v", err)
	}

	if math.Abs(nextX2-expectedX2) > 1e-9 {
		t.Errorf("NextStep wrong result on step 2. Got %f, expected %f", nextX2, expectedX2)
	}
}

// TestEulerExplicit_NextStepOnGrid tests the vector/grid stepping function.
func TestEulerExplicit_NextStepOnGrid(t *testing.T) {
	mockFunc := &MockPotential{}
	dt := 0.1
	solver := &EulerExplicit{}
	solver = solver.NewDefine(dt, mockFunc)

	// Grid inputs
	gridX := []float64{1.0, 2.0, 3.0}
	// Expected: x_new = x * (1 - dt)
	expectedGrid := []float64{
		1.0 * (1.0 - dt),
		2.0 * (1.0 - dt),
		3.0 * (1.0 - dt),
	}

	err := solver.NextStepOnGrid(gridX, 0.0)
	if err != nil {
		t.Fatalf("NextStepOnGrid failed: %v", err)
	}
	fmt.Println(gridX)

	for i, val := range gridX {
		if math.Abs(val-expectedGrid[i]) > 1e-9 {
			t.Errorf("Index %d: got %f, expected %f", i, val, expectedGrid[i])
		}
	}
}

// TestEulerExplicit_100StepsStable checks long-run behavior in a stable regime (dt < 2).
func TestEulerExplicit_100StepsStable(t *testing.T) {
	mockFunc := &MockPotential{}
	dt := 0.01
	steps := 100

	solver := (&EulerExplicit{}).NewDefine(dt, mockFunc)

	x := 1.0
	time := 0.0
	for i := 0; i < steps; i++ {
		var err error
		x, err = solver.NextStep(x, time)
		if err != nil {
			t.Fatalf("NextStep failed at step %d: %v", i, err)
		}
		time += dt
	}

	// exact Euler recurrence for this ODE: x_n = (1-dt)^n x0
	expectedEuler := math.Pow(1.0-dt, float64(steps))
	if math.Abs(x-expectedEuler) > 1e-12 {
		t.Errorf("100-step mismatch against Euler closed form. got=%g expected=%g", x, expectedEuler)
	}

	// optional sanity against analytical x(t)=exp(-t), allow looser tolerance (method error)
	expectedAnalytical := math.Exp(-float64(steps) * dt)
	if math.Abs(x-expectedAnalytical) > 2e-3 {
		t.Errorf("100-step mismatch against analytical solution. got=%g expected=%g", x, expectedAnalytical)
	}
}

// TestEulerExplicit_100StepsDivergence checks divergence in unstable regime (dt > 2).
func TestEulerExplicit_100StepsDivergence(t *testing.T) {
	mockFunc := &MockPotential{}
	dt := 2.1 // amplification factor = (1-dt) = -1.1, |factor| > 1 => divergence
	steps := 100

	solver := (&EulerExplicit{}).NewDefine(dt, mockFunc)

	x := 1.0
	time := 0.0
	for i := 0; i < steps; i++ {
		var err error
		x, err = solver.NextStep(x, time)
		if err != nil {
			t.Fatalf("NextStep failed at step %d: %v", i, err)
		}
		time += dt
	}

	expected := math.Pow(1.0-dt, float64(steps))
	if math.Abs((x-expected)/expected) > 1e-12 {
		t.Errorf("divergence recurrence mismatch. got=%g expected=%g", x, expected)
	}

	if math.Abs(x) <= 1e3 {
		t.Errorf("expected divergence after %d steps, got |x|=%g", steps, math.Abs(x))
	}
}

// TestEulerExplicit_ReDefine verifies that changing parameters works.
func TestEulerExplicit_ReDefine(t *testing.T) {
	mockFunc := &MockPotential{}
	dt1 := 0.1
	solver := &EulerExplicit{}
	solver = solver.NewDefine(dt1, mockFunc)

	// Check initial set up
	if solver.deltaTime != dt1 {
		t.Errorf("Initial dt mismatch. Got %f, expected %f", solver.deltaTime, dt1)
	}

	dt2 := 0.05
	solver.ReDefine(dt2, mockFunc)

	if solver.deltaTime != dt2 {
		t.Errorf("Redefined dt mismatch. Got %f, expected %f", solver.deltaTime, dt2)
	}
}

// TestEulerExplicit_Name verifies the Name method.
func TestEulerExplicit_Name(t *testing.T) {
	solver := &EulerExplicit{}
	expected := "Euler Explicit Method"
	if solver.Name() != expected {
		t.Errorf("Name mismatch. Got %s, expected %s", solver.Name(), expected)
	}
}

// A simple harmonic oscillator (approx) or just a linear growth for order testing
// Let's use dx/dt = t, x(0) = 0 => x(t) = t^2 / 2
// This helps check higher order exactness, but since our interface is float64 scalar for simple potential,
// we just stick to the decay one for stability and basic step correctness.

// The following interfaces are declared to understand the type assertion.
type ScalarODESolver interface {
    NextStep(x, time float64) (float64, error)
}

type GridODESolver interface {
    NextStepOnGrid(x []float64, time float64) error
}

type GetSolvers interface {
    getSolvers()
}

type ExplicitSolvers struct{}

type SolverTestCase struct {
    Name         string
    Factory      func(dt float64, op gridData.TDPotentialOp) ScalarODESolver
    Order        int
    SupportsGrid bool
}

func (s ExplicitSolvers) getSolvers() []SolverTestCase {
    return []SolverTestCase{
        {
            Name: "Euler Explicit",
            Factory: func(dt float64, op gridData.TDPotentialOp) ScalarODESolver {
                return (&EulerExplicit{}).NewDefine(dt, op)
            },
            Order:        1,
            SupportsGrid: true,
        },
        {
            Name: "Heun's Explicit",
            Factory: func(dt float64, op gridData.TDPotentialOp) ScalarODESolver {
                return (&HeunsExplicit{}).NewDefine(dt, op)
            },
            Order:        2,
            SupportsGrid: true,
        },
        {
            Name: "MidPoint Explicit",
            Factory: func(dt float64, op gridData.TDPotentialOp) ScalarODESolver {
                return (&MidPointExplicit{}).NewDefine(dt, op)
            },
            Order:        2,
            SupportsGrid: true,
        },
        {
            Name: "Ralston 2nd Order",
            Factory: func(dt float64, op gridData.TDPotentialOp) ScalarODESolver {
                return (&Ralston2order{}).NewDefine(dt, op)
            },
            Order:        2,
            SupportsGrid: false, // no NextStepOnGrid in implementation
        },
        {
            Name: "Ralston 3rd Order",
            Factory: func(dt float64, op gridData.TDPotentialOp) ScalarODESolver {
                return (&Ralston3Order{}).NewDefine(dt, op)
            },
            Order:        3,
            SupportsGrid: true,
        },
        {
            Name: "Heun's 3rd Order",
            Factory: func(dt float64, op gridData.TDPotentialOp) ScalarODESolver {
                return (&Huens3Explicit{}).NewDef(dt, op)
            },
            Order:        3,
            SupportsGrid: false,
        },
        {
            Name: "Van der Houwen's 3rd Order",
            Factory: func(dt float64, op gridData.TDPotentialOp) ScalarODESolver {
                return (&VDHouwenExplicit{}).NewDef(dt, op)
            },
            Order:        3,
            SupportsGrid: false,
        },
        {
            Name: "SSP Runge-Kutta 3",
            Factory: func(dt float64, op gridData.TDPotentialOp) ScalarODESolver {
                return (&SSPRungeKutta3{}).NewDef(dt, op)
            },
            Order:        3,
            SupportsGrid: false,
        },
        {
            Name: "Runge-Kutta 3 Explicit",
            Factory: func(dt float64, op gridData.TDPotentialOp) ScalarODESolver {
                return (&RungeKutta3Explicit{}).NewDef(dt, op)
            },
            Order:        3,
            SupportsGrid: false,
        },
        {
            Name: "Runge-Kutta 4 Explicit",
            Factory: func(dt float64, op gridData.TDPotentialOp) ScalarODESolver {
                return (&RungeKutta4Explicit{}).NewDef(dt, op)
            },
            Order:        4,
            SupportsGrid: true,
        },
        {
            Name: "Runge-Kutta 3/8 (Order 4)",
            Factory: func(dt float64, op gridData.TDPotentialOp) ScalarODESolver {
                return (&RungeKutta38{}).NewDefine(dt, op)
            },
            Order:        4,
            SupportsGrid: false,
        },
        {
            Name: "Ralston 4th Order",
            Factory: func(dt float64, op gridData.TDPotentialOp) ScalarODESolver {
                return (&Ralston4Order{}).NewDefine(dt, op)
            },
            Order:        4,
            SupportsGrid: false,
        },
        {
            Name: "Nystrom 5th Order",
            Factory: func(dt float64, op gridData.TDPotentialOp) ScalarODESolver {
                return (&Nystrom5Explicit{}).NewDefine(dt, op)
            },
            Order:        5,
            SupportsGrid: false,
        },
    }
}

func TestAllSolvers_GridDecay(t *testing.T) {
    mockFunc := &MockPotential{}
    s := ExplicitSolvers{}
    dt := 0.001
    steps := 10
    initialGrid := []float64{1.0, 2.0, -1.0, 0.5}

    for _, tc := range s.getSolvers() {
        t.Run(tc.Name+"_Grid", func(t *testing.T) {
            solver := tc.Factory(dt, mockFunc)
            gridSolver, ok := solver.(GridODESolver)
            if !ok {
                t.Skip("solver does not implement NextStepOnGrid")
            }

            currentGrid := make([]float64, len(initialGrid))
            copy(currentGrid, initialGrid)
            currentTime := 0.0

            for s := 0; s < steps; s++ {
                if err := gridSolver.NextStepOnGrid(currentGrid, currentTime); err != nil {
                    t.Fatalf("NextStepOnGrid failed at step %d: %v", s, err)
                }
                currentTime += dt
            }

            for i, val := range currentGrid {
                exact := initialGrid[i] * math.Exp(-currentTime)
                if math.Abs(val-exact) > 1e-4 {
                    t.Errorf("Grid point %d: got %f, expected %f", i, val, exact)
                }
            }
        })
    }
}

func TestAllSolvers_FactorySmoke(t *testing.T) {
    mockFunc := &MockPotential{}
    s := ExplicitSolvers{}
    for _, tc := range s.getSolvers() {
        t.Run(tc.Name+"_FactorySmoke", func(t *testing.T) {
            solver := tc.Factory(0.01, mockFunc)
            if solver == nil {
                t.Fatal("factory returned nil solver")
            }
            x, err := solver.NextStep(1.0, 0.0)
            if err != nil {
                t.Fatalf("NextStep failed: %v", err)
            }
            if math.IsNaN(x) || math.IsInf(x, 0) {
                t.Fatalf("invalid result from NextStep: %v", x)
            }
        })
    }
}

func TestAllSolvers_SingleStepErrorBound(t *testing.T) {
    mockFunc := &MockPotential{}
    x0 := 1.0
    t0 := 0.0
    dt := 0.01
    exact := x0 * math.Exp(-dt)

    s := ExplicitSolvers{}
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