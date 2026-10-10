package plotmath

import (
	"math"
	"strings"
	"testing"
)

func TestCompileND(t *testing.T) {
	for _, tc := range []struct {
		expression string
		want       float64
	}{
		{"0.5*(x1^2+x2^2+x3^2+x4^2)+0.2*x1*x2", 15.4},
		{"-x1^2 + 2^3^2", 511},
		{"exp(0)+ln(e)+log(100)+sqrt(x4)+min(x1,x2)+max(x3,x4)+pow(x2,3)", 19},
		{"sin(pi/2) + cos(0) + x4**2", 18},
	} {
		fn, err := CompileND(tc.expression, 4)
		if err != nil {
			t.Fatal(err)
		}
		if got := fn([]float64{1, 2, 3, 4}); math.Abs(got-tc.want) > 1e-12 {
			t.Errorf("%s = %g, want %g", tc.expression, got, tc.want)
		}
		if !math.IsNaN(fn([]float64{1})) {
			t.Fatal("incorrect coordinate count must return NaN")
		}
	}
	fn, err := CompileND("x12 + x1", 12)
	if err != nil {
		t.Fatal(err)
	}
	x := make([]float64, 12)
	x[0], x[11] = 2, 4
	if fn(x) != 6 {
		t.Fatal("indexed variables must support dimensions above nine")
	}
}

func TestCompileNDRejectsInvalidExpressions(t *testing.T) {
	for _, source := range []string{"", "x3", "x0", "x01", "x", "y", "V=x1", "2x1", "sin()", "min(x1)", "pow(x1,2,3)", "x1; exit()", "window.alert(1)", "1e999", strings.Repeat("1+", 128) + "1", strings.Repeat(" ", 1025)} {
		if _, err := CompileND(source, 2); err == nil {
			t.Errorf("accepted %q", source)
		}
	}
	if _, err := CompileND("1", 0); err == nil {
		t.Fatal("accepted zero dimensions")
	}
}

func TestExpressionSamplingDoesNotAllocate(t *testing.T) {
	legacy, err := Compile("sin(x)+y^2", 2)
	if err != nil {
		t.Fatal(err)
	}
	nd, err := CompileND("sin(x1)+x2^2", 2)
	if err != nil {
		t.Fatal(err)
	}
	x := []float64{1, 2}
	if n := testing.AllocsPerRun(100, func() { legacy(1, 2, 0); nd(x) }); n != 0 {
		t.Fatalf("sampling allocated %g times", n)
	}
}
