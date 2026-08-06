package gridData

import (
	"testing"
)

func TestGaussian_ForceAt(t *testing.T) {
	grid, err := NewTimeGrid(10, 20, 100)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	if grid.tMin != 0.0 {
		t.Errorf("expected tMin = 0.0, got %v", grid.tMin)
	}
	if grid.tMax != 200 {
		t.Errorf("expected tMax = %v, got %v", 200, grid.tMax)
	}
	if grid.nPoints != 2000 {
		t.Errorf("expected nPoints = %v, got %v", 2000, grid.nPoints)
	}
}
