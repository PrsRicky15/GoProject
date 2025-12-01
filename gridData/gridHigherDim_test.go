package gridData

import (
	"fmt"
	"testing"
)

func TestFromLength(t *testing.T) {
	grid, _ := NewNDimGridUniform(3, -10.0, 10.0, 128)
	fmt.Println(grid)

	grid2, _ := NewNDimGrid(
		[]float64{-5, -10},
		[]float64{5, 10},
		[]uint32{64, 128})

	fmt.Println(grid2)
}
