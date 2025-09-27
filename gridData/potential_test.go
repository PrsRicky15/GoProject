package gridData

import (
	"testing"
)

func TestGaussian_ForceAt(t *testing.T) {
	length := 10.0
	nPoints := uint32(100)

	grid, err := TimeGridFromLength(length, nPoints)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	if grid.tMin != 0.0 {
		t.Errorf("expected tMin = 0.0, got %v", grid.tMin)
	}
	if grid.tMax != length {
		t.Errorf("expected tMax = %v, got %v", length, grid.tMax)
	}
	if grid.nPoints != nPoints {
		t.Errorf("expected nPoints = %v, got %v", nPoints, grid.nPoints)
	}
}
