package Quantum

import (
	"context"
	"math"
	"testing"

	"GoProject/gridData"
)

func TestNDLocalizedWellDegeneracyAgainstDense(t *testing.T) {
	g, _ := gridData.NewNDimGridUniform(3, -6, 6, 10)
	pot := gridData.NDimPotentialFunc(func(x []float64) float64 {
		return -4 * (math.Exp(-((x[0]-1.5)*(x[0]-1.5)+x[1]*x[1]+x[2]*x[2])/(2*1.2*1.2)) + math.Exp(-((x[0]+1.5)*(x[0]+1.5)+x[1]*x[1]+x[2]*x[2])/(2*1.2*1.2)))
	})
	h, err := NewNDimFourierHamiltonian(g, []float64{1, 1, 1}, pot)
	if err != nil {
		t.Fatal(err)
	}
	defer h.Close()
	want, _, err := h.RealDiagonalize()
	if err != nil {
		t.Fatal(err)
	}
	got, err := h.LowestStates(context.Background(), 4, EigenOptions{Tolerance: 1e-7, MaxSubspace: 96})
	if err != nil {
		t.Fatal(err)
	}
	for i := 0; i < 4; i++ {
		if math.Abs(got.Values[i]-want[i]) > 1e-6 {
			t.Fatalf("lost a low/degenerate level: got %v want %v", got.Values, want[:4])
		}
	}
}
