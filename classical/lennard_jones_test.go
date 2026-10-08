package classical

import (
	"math"
	"testing"
)

func TestLJForcesAndSeed(t *testing.T) {
	s, err := NewLJSystem(4, 5, .2, 42)
	if err != nil {
		t.Fatal(err)
	}
	other, _ := NewLJSystem(4, 5, .2, 42)
	for i := range s.X {
		if s.VX[i] != other.VX[i] {
			t.Fatal("seed is not deterministic")
		}
	}
	s.X[0] = 1.
	s.Y[0] = 1.
	s.X[1] = 2.2
	s.Y[1] = 1.
	if err = s.Forces(); err != nil {
		t.Fatal(err)
	}
	force := s.fx[0]
	epsilon := 1e-6
	s.X[0] += epsilon
	_ = s.Forces()
	plus := s.Potential
	s.X[0] -= 2 * epsilon
	_ = s.Forces()
	minus := s.Potential
	if math.Abs(force+(plus-minus)/(2*epsilon)) > 1e-6 {
		t.Fatal("LJ force is not minus energy gradient")
	}
	if _, err = NewLJSystem(64, 5, .2, 1); err == nil {
		t.Fatal("accepted overlapping initial lattice")
	}
}
