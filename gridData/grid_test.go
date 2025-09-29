package gridData

import (
	"testing"
)

func TestNewRGrid(t *testing.T) {
	grid, err := NewRGrid(0, 25, 100)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	if grid.rMin != 0. {
		t.Errorf("expected tMin = %v, got %v", 0., grid.rMin)
	}
	if grid.rMax != 25. {
		t.Errorf("expected tMax = %v, got %v", 25, grid.rMax)
	}
	if grid.nPoints != 100 {
		t.Errorf("expected nPoints = %v, got %v", 100, grid.nPoints)
	}
}

func TestNewFromLength(t *testing.T) {
	grid, err := NewFromLength(10, 100)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	if grid.rMin != -5. {
		t.Errorf("expected tMin = %v, got %v", -5, grid.rMin)
	}
	if grid.rMax != 5. {
		t.Errorf("expected tMax = %v, got %v", 5, grid.rMax)
	}
	if grid.nPoints != 100 {
		t.Errorf("expected nPoints = %v, got %v", 100, grid.nPoints)
	}
}

func TestTimeGridFromLength(t *testing.T) {
	grid, err := NewTimeGrid(10, 100, 100)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	if grid.tMin != 0.0 {
		t.Errorf("expected tMin = 0.0, got %v", grid.tMin)
	}
	if grid.tMax != 1000. {
		t.Errorf("expected tMax = %v, got %v", 1000, grid.tMax)
	}
	if grid.nPoints != 100*100 {
		t.Errorf("expected nPoints = %v, got %v", 100*100, grid.nPoints)
	}
}

func TestNewRGridFromFile(t *testing.T) {
	grid, err := NewRGridFromFile("testGrid")

	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if grid.rMin != -3. {
		t.Errorf("expected tMin = %v, got %v", -3, grid.rMin)
	}
	if grid.rMax != 3 {
		t.Errorf("expected tMax = %v, got %v", 3, grid.rMax)
	}
	if grid.nPoints != 18 {
		t.Errorf("expected nPoints = %v, got %v", 18, grid.nPoints)
	}

	grid.DisplayInfo()
}

func TestNewTGridFromFile(t *testing.T) {
	grid, err := NewTGridFromFile("testGrid")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	if grid.tMin != 0.0 {
		t.Errorf("expected tMin = 0.0, got %v", grid.tMin)
	}
	if grid.tMax != 100 {
		t.Errorf("expected tMax = %v, got %v", 100, grid.tMax)
	}
	if grid.nPoints != 1000 {
		t.Errorf("expected nPoints = %v, got %v", 1000, grid.nPoints)
	}

	grid.DisplayInfo()
}
