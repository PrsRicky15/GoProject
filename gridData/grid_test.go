package gridData

import (
	"fmt"
	"os"
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
	if len(os.Args) < 2 {
		fmt.Println("Usage: go run main.go <file-path>")
		os.Exit(1)
	}
	filePath := os.Args[1]
	grid, _ := NewRGridFromFile(filePath)
	grid.DisplayInfo()
	fmt.Println("kGrid points:")
	grid.DisplayRgrid()
}
