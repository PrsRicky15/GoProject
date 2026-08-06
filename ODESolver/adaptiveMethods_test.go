package EquationSolver

import (
	"math"
	"testing"

	"GoProject/gridData"
)

type AdaptiveSolvers struct {
    Name         string
}

var adaptiveSolverNames = []AdaptiveSolvers{
    {Name: "HuensEuler"},
    {Name: "FehlbergRK12"},
    {Name: "BogackiShampine"},
    {Name: "RKFelberg"},
    {Name: "CashKarp"},
    {Name: "DormandPrince"},
}


type AdaptiveSolverTestCase struct{
    Solver ScalarODESolver
    GetDt func() float64
    Order        int
    SupportsGrid bool
}

func (s *AdaptiveSolvers) getAdaptiveSolvers(dt float64, op gridData.TDPotentialOp) AdaptiveSolverTestCase {
    var solver ScalarODESolver
    var getDt func() float64
    var order int
    var SupportsGrid bool

    switch s.Name {
        case "HuensEuler":
            solver = (&HuensEuler{}).NewDefine(dt, op)
            getDt = func() float64 { return solver.(*HuensEuler).deltaTime }
            order = 2
            SupportsGrid = false
        case "FehlbergRK12":
            solver = (&FehlbergRK12{}).NewDefine(dt, op)
            getDt = func() float64 { return solver.(*FehlbergRK12).deltaTime }
            order = 2
            SupportsGrid = false
        case "BogackiShampine":
            solver = (&BogackiShampine{}).NewDefine(dt, op)
            getDt = func() float64 { return solver.(*BogackiShampine).deltaTime }
            order = 3
            SupportsGrid = false
        case "RKFelberg":
            solver = (&RKFelberg{}).NewDefine(dt, op)
            getDt = func() float64 { return solver.(*RKFelberg).deltaTime }
            order = 4
            SupportsGrid = false
        case "CashKarp":
            solver = (&CashKarp{}).NewDefine(dt, op)
            getDt = func() float64 { return solver.(*CashKarp).deltaTime }
            order = 4
            SupportsGrid = false
        case "DormandPrince":
            solver = (&DormandPrince{}).NewDefine(dt, op)
            getDt = func() float64 { return solver.(*DormandPrince).deltaTime }
            order = 5
            SupportsGrid = false
    }

    return AdaptiveSolverTestCase{
        Solver: solver,
        GetDt:  getDt,
        Order:  order,
        SupportsGrid: SupportsGrid,
	}
}

func TestAllAdaptiveSolvers_SingleStepErrorBound(t *testing.T) {
    mockFunc := &MockGaussianPotential{} //&MockPotential{}
    x0 := 1.0
    t0 := 0.0
    dt := 0.001
    // exact := x0 * math.Exp(-dt)
    exact := x0 * math.Exp(-dt*dt) // for the Gaussian potential

    for _, tc := range adaptiveSolverNames {
        t.Run(tc.Name+"_SingleStepError", func(t *testing.T) {
            h := tc.getAdaptiveSolvers(dt, mockFunc)
            x1, err := h.Solver.NextStep(x0, t0)
            if err != nil {
                t.Fatalf("NextStep failed: %v", err)
            }

            // In adaptive methods, x1 == x0 means the error was > tolerance and step was rejected.
            // We cannot test accuracy on a rejected step.
            if x1 == x0 {
                t.Logf("Step was rejected (adaptive logic triggered). Skipping accuracy check for %s at dt=%v", tc.Name, dt)
                return
            }

            diff := math.Abs(x1 - exact) // accuracy of the single step

            // loose bounds by nominal order to catch gross regressions
            bound := 5e-4 // order 1
            if h.Order >= 2 {
                bound = 5e-5
            }
            if h.Order >= 3 {
                bound = 5e-6
            }
            if h.Order >= 4 {
                bound = 5e-7
            }
            if diff > bound {
                t.Errorf("single-step error too large: got=%e bound=%e", diff, bound)
            }
        })
    }
}

func TestAllAdaptiveSolvers_FullIntegration(t *testing.T) {
    mockFunc := &MockPotential{}
    x0 := 1.0
    t0 := 0.0
    dt := 0.5
    targetTime := 1.0
    // exact := x0 * math.Exp(-dt)

    for _, tc := range adaptiveSolverNames {
        t.Run(tc.Name+"_FullIntegration", func(t *testing.T) {
            h := tc.getAdaptiveSolvers(dt, mockFunc)
            dtCurrent := dt
            xCurrent := x0
            tCurrent := t0
            steps := 0
            maxSteps := 100000 // prevent infinite loops
            for tCurrent < targetTime{
                if steps >= maxSteps {
                    t.Fatalf("Exceeded maximum steps without reaching target time: Steps=%d, tCurrent=%.4f, dtCurrent=%.4f", steps, tCurrent, dtCurrent)
                }

                xNext, err := h.Solver.NextStep(xCurrent, tCurrent)
                if err != nil {
                    t.Fatalf("NextStep failed: %v", err)
                }
                if xNext == xCurrent {
                    // Step was rejected, reduce dt and try again
                    //dtOld := dtCurrent
                    dtCurrent = h.GetDt()
                    // t.Logf("[Step Rejected] Steps=%d, t=%.6f. dt reduced from %e to %e", steps, tCurrent, dtOld, dtCurrent)
                    steps++
                    continue
                }

                xCurrent = xNext
                tCurrent += dtCurrent
                dtCurrent = h.GetDt() // get the potentially updated dt from the solver
                steps++

                //t.Logf("Step %d: t=%.4f, x=%.6f, dt=%.6f\n", steps, tCurrent, xCurrent, dtCurrent)
            }

            
            exactAtFinal := math.Exp(-tCurrent) // to accomodate overshooting beyond target time
            errDiffFinal := math.Abs(xCurrent - exactAtFinal)

            bound := 5e-4 // order 1
            if h.Order >= 2 {
                bound = 5e-4
            }
            if h.Order >= 3 {
                bound = 5e-5
            }
            if h.Order >= 4 {
                bound = 5e-6
            }
            if errDiffFinal > bound {
                t.Errorf("Integration result inaccurate at t=%v. Got %v, expected %v (diff %e)", 
                    tCurrent, xCurrent, exactAtFinal, errDiffFinal)
            }
        })
    }
}