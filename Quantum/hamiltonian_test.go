package Quantum

import (
	"GoProject/gridData"
	"testing"
)

func TestHamiltonianOp(t *testing.T) {
	grid, _ := gridData.NewFromLength(16., 48)
	PotE := gridData.Harmonic[float64]{ForceConst: 1.}
	Harmonic := NewHamil(grid, 1., PotE)
	Harmonic.Mat()
}
