package gridData

import (
	"testing"
)

func TestTimeGridFromLength(t *testing.T) {
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

func TestNewRGridFromFile(t *testing.T) {
	grid, err := NewRGridFromFile("GridInfo")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	if grid.rMin != -10 {
		t.Errorf("expected tMin = 0.0, got %v", grid.rMin)
	}
	if grid.rMax != 10. {
		t.Errorf("expected tMax = %v, got %v", 10, grid.rMax)
	}
	if grid.nPoints != 60 {
		t.Errorf("expected nPoints = %v, got %v", 60, grid.nPoints)
	}
}
