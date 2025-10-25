package Quantum

import (
	"GoProject/OperatorAlgebra"
	"GoProject/gridData"

	"gonum.org/v1/gonum/mat"
)

type HamiltonianOp struct {
	grid *gridData.RadGrid
	kinE OperatorAlgebra.KineticOp
	potE gridData.PotentialOp[float64]
	hmat mat.Matrix
}

func NewHamil(grid *gridData.RadGrid, mass float64, Pot gridData.PotentialOp[float64]) *HamiltonianOp {
	kinE := OperatorAlgebra.NewKeDVR(grid, mass)
	Hmat := mat.Matrix(mat.NewDense(int(grid.NPoints()), int(grid.NPoints()), nil))
	return &HamiltonianOp{
		grid: grid,
		kinE: kinE,
		potE: Pot,
		hmat: Hmat,
	}
}

func (op *HamiltonianOp) Mat() {
	vPot := op.grid.PotentialOnGrid(op.potE)
	err := op.grid.PrintVectorToFileRe(vPot, "potent.dat", "%21.14e")
	if err != nil {
		return
	}
	// need a change
	op.hmat = op.kinE.KMat

	for i := 0; i < int(op.grid.NPoints()); i++ {
		op.hmat.(*mat.Dense).Set(i, i, op.hmat.At(i, i)+vPot[i])
	}
}

func (op *HamiltonianOp) EvaluateOp() *mat.Dense {
	return nil
}
