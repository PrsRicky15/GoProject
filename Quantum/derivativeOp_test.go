package Quantum

import (
	"GoProject/gridData"
	"testing"
)

func TestKeDVR_Evaluate(t *testing.T) {
	rgrid, err := gridData.NewFromLength(10., 30)
	if err != nil {
		panic(err)
	}
	kinE := NewKeDVR(rgrid, 1.)
	kinE.Evaluate()
}

func TestHamiltonianOp(t *testing.T) {
	grid, _ := gridData.NewFromLength(16., 48)
	PotE := gridData.Harmonic{ForceConst: 1.}
	Harmonic := NewHamil(grid, 1., PotE)
	Harmonic.Mat()
}
