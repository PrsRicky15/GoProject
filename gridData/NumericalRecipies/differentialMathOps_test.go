package NumericalRecipies

import (
    "math"
    "testing"
)

func TestGradientScalar(t *testing.T) {
    // f(x) = x^2 + 2x + 1  -> f'(x) = 2x + 2
    f := func(x float64) float64 { return x*x + 2*x + 1 }
    x := 3.0
    g := (&GradientScalar{}).NewDefine(f, 1e-6)
	got := g.calculateGradientScalarAtX(x)
    want := 2*x + 2
    if math.Abs(got-want) > 1e-5 {
        t.Errorf("GradientScalar wrong: got %v, want %v", got, want)
    }
}

func TestGradientMultivariate(t *testing.T) {
    // F(x,y) = x^2 + 3y  -> ∇F = [2x, 3]
    F := func(v []float64) float64 {
        return v[0]*v[0] + 3*v[1]
    }
    pt := []float64{2.0, -1.5}
	g := (&Gradient{}).NewDefine(F, 1e-6)
    grad := g.calculateGradientAtX(pt)
    if len(grad) != 2 {
        t.Fatalf("expected gradient length 2, got %d", len(grad))
    }
    if math.Abs(grad[0]-(2*pt[0])) > 1e-5 || math.Abs(grad[1]-3) > 1e-5 {
        t.Errorf("Gradient multivariate incorrect: %v", grad)
    }
}

func TestJacobianSimple(t *testing.T) {
    // G(x,y) = [ x^2, x*y ]  -> J = [[2x,0],[y,x]]
    G := func(v []float64) []float64 {
        return []float64{v[0] * v[0], v[0] * v[1]}
    }
    pt := []float64{3.0, 4.0}
	g := (&Jacobian{}).NewDefine(G, 1e-6)
    J := g.calculateJacobianAtX(pt)
    if len(J) != 2 || len(J[0]) != 2 || len(J[1]) != 2 {
        t.Fatalf("Jacobian dimensions wrong: %v", J)
    }
    if math.Abs(J[0][0]-6.0) > 1e-4 || math.Abs(J[0][1]-0.0) > 1e-6 {
        t.Errorf("first row incorrect: %v", J[0])
    }
    if math.Abs(J[1][0]-pt[1]) > 1e-4 || math.Abs(J[1][1]-pt[0]) > 1e-4 {
        t.Errorf("second row incorrect: %v", J[1])
    }
}
