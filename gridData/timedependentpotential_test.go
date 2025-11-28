package gridData

import (
	"math"
	"testing"
)

// A simple static potential for test: V(x) = x^2
type TestStatic struct{}

func (TestStatic) EvaluateAt(x float64) float64 {
	return x * x
}
func (TestStatic) EvaluateOnGrid(x []float64) []float64 {
	res := make([]float64, len(x))
	for i, xi := range x {
		res[i] = xi * xi
	}
	return res
}
func (TestStatic) ForceAt(x float64) float64 {
	return 2 * x
}
func (TestStatic) ForceOnGrid(x []float64) []float64 {
	res := make([]float64, len(x))
	for i, xi := range x {
		res[i] = 2 * xi
	}
	return res
}

func TestTimeDependentPotential_StaticOnly(t *testing.T) {
	static := TestStatic{}
	tp := NewTimeDependentPotential(static, nil)

	x := 2.0
	tval := 1.23

	got := tp.EvaluateAt(x, tval)
	want := x * x

	if got != want {
		t.Fatalf("StaticOnly: expected %.6f, got %.6f", want, got)
	}
}

func TestTimeDependentPotential_TDOnly(t *testing.T) {
	tdFunc := func(x, t float64) float64 {
		return math.Sin(t) * x
	}

	tp := NewTimeDependentPotential(nil, tdFunc)

	x := 3.0
	tval := 0.5

	got := tp.EvaluateAt(x, tval)
	want := math.Sin(tval) * x

	if math.Abs(got-want) > 1e-12 {
		t.Fatalf("TDOnly: expected %.6f, got %.6f", want, got)
	}
}

func TestTimeDependentPotential_StaticPlusTD(t *testing.T) {
	static := TestStatic{}              // V_static = x^2
	tdFunc := func(x, t float64) float64 { // V_td = x*t
		return x * t
	}

	tp := NewTimeDependentPotential(static, tdFunc)

	x := 4.0
	tval := 2.0

	got := tp.EvaluateAt(x, tval)
	want := x*x + x*tval // 16 + 8 = 24

	if got != want {
		t.Fatalf("Static+TD: expected %.6f, got %.6f", want, got)
	}
}

func TestTimeDependentPotential_GridEval(t *testing.T) {
	static := TestStatic{} // x^2
	tdFunc := func(x, t float64) float64 {
		return x + t
	}

	tp := NewTimeDependentPotential(static, tdFunc)

	xgrid := []float64{-1, 0, 1}
	tval := 1.0

	got := tp.EvaluateOnRGrid(xgrid, tval)
	want := []float64{
		(-1)*(-1) + (-1+tval), // 1 + 0 = 1
		0*0 + (0+tval),        // 0 + 1 = 1
		1*1 + (1+tval),        // 1 + 2 = 3
	}

	for i := range got {
		if math.Abs(got[i]-want[i]) > 1e-12 {
			t.Fatalf("GridEval: at index %d expected %.6f, got %.6f",
				i, want[i], got[i])
		}
	}
}

func TestTimeDependentPotential_GridEvalInPlace(t *testing.T) {
	static := TestStatic{} // x^2
	tdFunc := func(x, t float64) float64 {
		return x + t
	}

	tp := NewTimeDependentPotential(static, tdFunc)

	xgrid := []float64{-2, -1, 1, 2}
	res := make([]float64, len(xgrid))
	tval := 3.0

	tp.EvaluateOnRGridInPlace(xgrid, res, tval)

	want := make([]float64, len(xgrid))
	for i, x := range xgrid {
		want[i] = x*x + (x + tval)
	}

	for i := range res {
		if math.Abs(res[i]-want[i]) > 1e-12 {
			t.Fatalf("InPlace: at index %d expected %.6f, got %.6f",
				i, want[i], res[i])
		}
	}
}

func TestTimeDependentPotential_UsesRealStaticPotential(t *testing.T) {
	// Gaussian potential to ensure integration with potential.go works.
	g := Gaussian[float64]{Cen: 0, Sigma: 1, Strength: 2}

	tdFunc := func(x, t float64) float64 {
		return t // simple constant-in-x time shift
	}

	tp := NewTimeDependentPotential(g, tdFunc)

	x := 0.5
	tval := 0.3

	// expected = Gaussian(x) + t
	want := g.EvaluateAt(x) + tval
	got := tp.EvaluateAt(x, tval)

	if math.Abs(got-want) > 1e-12 {
		t.Fatalf("Gaussian+TD: expected %.6f, got %.6f", want, got)
	}
}

