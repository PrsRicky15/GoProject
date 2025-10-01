package Quantum

import (
	"GoProject/gridData"

	"gonum.org/v1/gonum/mat"
)

type HamiltonianOp struct {
	grid *gridData.RadGrid
	kinE KineticOp
	potE *gridData.PotentialOp
	hmat *mat.Dense
}

func NewHamil(grid *gridData.RadGrid, mass float64, Pot gridData.PotentialOp) *HamiltonianOp {
	kinE := NewKeDVR(grid, mass)
	return &HamiltonianOp{
		grid: grid,
		kinE: kinE,
		potE: &Pot,
	}
}

func (op *HamiltonianOp) Mat() {
	vPot := op.grid.PotentialOnGrid(*op.potE)
	err := op.grid.PrintVectorToFile(vPot, "potent.dat", "%21.14e")
	if err != nil {
		return
	}
	op.hmat = op.kinE.(*KeDvrBasis).kMat

	for i := 0; i < int(op.grid.NPoints()); i++ {
		op.hmat.Set(i, i, op.hmat.At(i, i)+vPot[i])
	}
}

func (op *HamiltonianOp) EvaluateOp() *mat.Dense {
	return nil
}
