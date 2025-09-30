package Quantum

import (
	"GoProject/gridData"
	"GoProject/matrix"
)

type HamiltonianOp struct {
	grid gridData.RadGrid
	kinE matrix.KineticOp
	potE gridData.PotentialOp
}

func NewHamil(grid gridData.RadGrid, mass float64, Pot gridData.PotentialOp) *HamiltonianOp {
	kinE := NewKeDVR(&grid, mass)
	return &HamiltonianOp{
		grid: grid,
		kinE: kinE,
		potE: Pot,
	}
}
