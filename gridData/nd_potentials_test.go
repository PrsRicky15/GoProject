package gridData

import (
	"math"
	"math/cmplx"
	"testing"
)

func TestNDimPotentialDefinitions(t *testing.T) {
	x := []float64{2, 3}
	for _, tc := range []struct {
		name string
		p    NDimPotentialOp
		want float64
	}{
		{"gaussian", Gaussian[float64]{Cen: 1, Sigma: 2, Strength: -3}, -3 * math.Exp(-5.0/8)},
		{"soft-core", SoftCore[float64]{Centre: 1, SoftParam: .5, Charge: -2}, -2 / math.Sqrt(5.25)},
		{"super-gaussian", SuperGaussian[float64]{Cen: 1, Sigma: 2, Strength: -3, Order: 4}, -3 * math.Exp(-25.0/16)},
		{"morse", Morse[float64]{Cen: 1, Alpha: .3, De: 2}, 2 * (math.Pow(1-math.Exp(-.3), 2) + math.Pow(1-math.Exp(-.6), 2))},
		{"polynomial", Polynomial[float64]{Coeffs: []float64{1, 2, 3}}, (1 + 4 + 12) + (1 + 6 + 27)},
	} {
		t.Run(tc.name, func(t *testing.T) {
			if got := tc.p.NDimEvaluateAt(x); math.Abs(got-tc.want) > 1e-13 {
				t.Fatalf("got %g want %g", got, tc.want)
			}
		})
	}
	g := Gaussian[float64]{Strength: 2, Sigma: 1}
	if math.Abs(g.NDimEvaluateAt([]float64{3, 4})-g.NDimEvaluateAt([]float64{0, 5})) > 1e-15 {
		t.Fatal("radial Gaussian is not rotation invariant")
	}
	m := MultiGaussian[float64]{Sigma: 1, Strength: -2, NumGauss: 3, Gap: 2}
	if math.Abs(m.NDimEvaluateAt(x)-m.EvaluateAt(2)-m.EvaluateAt(3)) > 1e-13 {
		t.Fatal("multi-Gaussian additive extension")
	}
	// The complex extension uses analytic squares, not |z|².
	z := []complex128{1 + 1i, 2 - .5i}
	sigma := 1.5
	gc := Gaussian[complex128]{Strength: -2, Sigma: sigma}
	want := -2 * cmplx.Exp(-(z[0]*z[0]+z[1]*z[1])/complex(2*sigma*sigma, 0))
	if cmplx.Abs(gc.NDimEvaluateAt(z)-want) > 1e-13 {
		t.Fatal("incorrect analytic complex extension")
	}
}
