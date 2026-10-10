package simulation

import (
	"context"
	"fmt"
	"math"
	"strings"
	"testing"
)

func TestCustomPotentialMatchesPreset(t *testing.T) {
	for _, dims := range []int{2, 3, 4} {
		for _, basis := range []string{"fourier", "dvr"} {
			t.Run(fmt.Sprintf("%dD/%s", dims, basis), func(t *testing.T) {
				r := testNDRequest(dims)
				r.Basis = basis
				r.States = 2
				r.Potential = "harmonic"
				terms := make([]string, dims)
				for d := range r.ND.Axes {
					r.ND.Axes[d].Points = 8
					terms[d] = fmt.Sprintf("x%d^2", d+1)
				}
				preset, err := Run(context.Background(), r)
				if err != nil {
					t.Fatal(err)
				}
				r.Potential = "custom"
				r.ND.Expression = "0.5*(" + strings.Join(terms, "+") + ")"
				// Custom formulas use raw coordinates, independent of preset transforms.
				for d := range r.ND.Axes {
					r.ND.Axes[d].Scale = 1.5
					r.ND.Axes[d].Center = .3
				}
				custom, err := Run(context.Background(), r)
				if err != nil {
					t.Fatal(err)
				}
				for s := 0; s < r.States; s++ {
					if math.Abs(custom.Metrics[s].Value-preset.Metrics[s].Value) > 1e-8 {
						t.Fatal("custom energies differ from preset")
					}
				}
				checkNDProjectionNorms(t, r, custom)
				if !strings.Contains(custom.Notes[0], r.ND.Expression) {
					t.Fatal("result must retain its formula")
				}
			})
		}
	}
}

func TestCustomCoupledOscillator(t *testing.T) {
	r := testNDRequest(2)
	r.Potential = "custom"
	r.States = 1
	r.ND.Expression = "0.5*(x1^2+x2^2)+0.2*x1*x2"
	result, err := Run(context.Background(), r)
	if err != nil {
		t.Fatal(err)
	}
	want := .5 * (math.Sqrt(1.2) + math.Sqrt(.8))
	if math.Abs(result.Metrics[0].Value-want) > 1e-5 {
		t.Fatalf("energy %g want %g", result.Metrics[0].Value, want)
	}
	cut := result.Charts[2]
	for iy, y := range cut.Y {
		for ix, x := range cut.X {
			if math.Abs(cut.Z[iy][ix]-(.5*(x*x+y*y)+.2*x*y)) > 1e-12 {
				t.Fatal("incorrect potential slice")
			}
		}
	}
}

func TestCustomPotentialErrors(t *testing.T) {
	for _, expression := range []string{"x3", "", "1/(x1-x1)", "sqrt(-1)", "exp(1000)", "100001"} {
		r := testNDRequest(2)
		r.Potential = "custom"
		r.ND.Expression = expression
		if _, err := Run(context.Background(), r); err == nil || !strings.Contains(err.Error(), "custom potential") {
			t.Fatalf("%q: %v", expression, err)
		}
	}
	// A hidden-coordinate slice may hit a singularity absent from the tensor grid.
	r := testNDRequest(3)
	r.States = 1
	r.Potential = "custom"
	r.ND.Expression = "1/(x3-0.3)^2"
	for d := range r.ND.Axes {
		r.ND.Axes[d].Points = 4
	}
	r.ND.Axes[2].Center = .3
	if _, err := Run(context.Background(), r); err == nil || !strings.Contains(err.Error(), "potential slice") {
		t.Fatalf("singular slice: %v", err)
	}
}
