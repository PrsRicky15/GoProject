package gridData

import (
	"fmt"
	"os"
	"testing"
)

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
