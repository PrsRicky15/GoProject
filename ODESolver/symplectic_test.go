package EquationSolver

import (
	"GoProject/gridData"
	"math"
	"testing"
)

func TestYoshidaConvergenceAndGrid(t *testing.T) {
	p := gridData.Harmonic[float64]{ForceConst: 1}
	errors := []float64{}
	for _, n := range []int{5, 10} {
		solver := (&Yoshida{}).NewDef(p, 1/float64(n))
		x, v := 1., 0.
		xs, vs := []float64{1}, []float64{0}
		for i := 0; i < n; i++ {
			x, v = solver.NextStep(x, v)
			solver.NextStepOnGrid(xs, vs)
		}
		if math.Abs(xs[0]-x) > 1e-14 || math.Abs(vs[0]-v) > 1e-14 {
			t.Fatal("scalar/grid disagreement")
		}
		errors = append(errors, math.Hypot(x-math.Cos(1), v+math.Sin(1)))
	}
	if errors[0]/errors[1] < 14 {
		t.Fatalf("not fourth order: %v", errors)
	}
	var redefined Yoshida
	redefined.Redefine(p, .1)
	x, _ := redefined.NextStep(1, 0)
	if math.Abs(x-math.Cos(.1)) > 1e-6 {
		t.Fatal("Redefine did not initialize coefficients")
	}
}

func TestStromerVerletInitialization(t *testing.T) {
	p := gridData.Polynomial[float64]{Coeffs: []float64{0}}
	solver := (&StromerVerlet{}).NewDef(p, .1)
	solver.Initiate(2, 3)
	if math.Abs(solver.NextStep(2)-2.3) > 1e-14 {
		t.Fatal("initial velocity reversed")
	}
	for _, start := range []float64{2, 5} {
		x := []float64{start}
		solver.InitiateGrid(x, []float64{3})
		solver.NextStepOnGrid(x)
		if math.Abs(x[0]-(start+.3)) > 1e-14 {
			t.Fatal("grid reinitialization failed")
		}
	}
}
