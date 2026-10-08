package gridData

import (
	"math"
	"slices"
	"testing"
)

func TestConjugatePointsFFTOrdering(t *testing.T) {
	for _, want := range [][]float64{{0}, {0, -1}, {0, 1, 2, -2, -1}, {0, 1, 2, -3, -2, -1}} {
		n := uint32(len(want))
		r, err := NewFromLength(2*math.Pi, n)
		if err != nil {
			t.Fatal(err)
		}
		times, err := NewTimeGrid(2*math.Pi, 1, n)
		if err != nil {
			t.Fatal(err)
		}
		if !slices.Equal(r.KValues(), want) || !slices.Equal(times.WValues(), want) {
			t.Fatalf("n=%d: k=%v, omega=%v, want %v", n, r.KValues(), times.WValues(), want)
		}
	}
}
