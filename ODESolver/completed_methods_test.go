package EquationSolver

import (
	"math"
	"testing"
)

type exponentialODE struct{}

func (exponentialODE) EvaluateAt(x, t float64) float64 { return x + t }
func (f exponentialODE) EvaluateOnRGrid(x []float64, t float64) []float64 {
	r := make([]float64, len(x))
	f.EvaluateOnRGridInPlace(x, r, t)
	return r
}
func (exponentialODE) EvaluateOnRGridInPlace(x, r []float64, t float64) {
	for i, v := range x {
		r[i] = v + t
	}
}

func TestCompletedExplicitMethods(t *testing.T) {
	for _, name := range []string{"ralston2", "ralston3", "nystrom5"} {
		errors := []float64{}
		for _, n := range []int{5, 10} {
			dt := 1 / float64(n)
			var solver ODESolver
			switch name {
			case "ralston2":
				solver = (&Ralston2order{}).NewDefine(dt, exponentialODE{})
			case "ralston3":
				solver = (&Ralston3Order{}).NewDefine(dt, exponentialODE{})
			case "nystrom5":
				solver = (&Nystrom5Explicit{}).NewDefine(dt, exponentialODE{})
			}
			x := 1.
			grid := []float64{1}
			for i := 0; i < n; i++ {
				var err error
				x, err = solver.NextStep(x, float64(i)*dt)
				if err != nil {
					t.Fatal(err)
				}
				if err = solver.NextStepOnGrid(grid, float64(i)*dt); err != nil {
					t.Fatal(err)
				}
			}
			if math.Abs(x-grid[0]) > 1e-12 {
				t.Fatalf("%s scalar/grid disagree", name)
			}
			errors = append(errors, math.Abs(x-(2*math.E-2)))
		}
		minRatio := 3.
		if name == "ralston3" {
			minRatio = 6
		}
		if name == "nystrom5" {
			minRatio = 24
		}
		if errors[0]/errors[1] < minRatio {
			t.Fatalf("%s convergence ratio %g", name, errors[0]/errors[1])
		}
	}
}
