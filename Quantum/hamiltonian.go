package Quantum

import (
	"GoProject/OperatorAlgebra"
	"GoProject/gridData"

	"gonum.org/v1/gonum/mat"
)

type HamiltonianOp[T gridData.VarType] struct {
	grid *gridData.RadGrid
	kinE OperatorAlgebra.KineticOp
	potE *gridData.PotentialOp[T]
	hmat mat.Matrix
}

func NewHamil(grid *gridData.RadGrid, mass float64, Pot gridData.PotentialOp[float64]) *HamiltonianOp[float64] {
	kinE := OperatorAlgebra.NewKeDVR(grid, mass)
	Hmat := mat.Matrix(mat.NewDense(int(grid.NPoints()), int(grid.NPoints()), nil))
	return &HamiltonianOp[float64]{
		grid: grid,
		kinE: kinE,
		potE: &Pot,
		hmat: Hmat,
	}
}

func (op *HamiltonianOp[float64]) Mat() {
	vPot := op.grid.PotentialOnGrid(*op.potE)
	err := op.grid.PrintVectorToFile(vPot, "potent.dat", "%21.14e")
	if err != nil {
		return
	}
	op.hmat = op.kinE.(*OperatorAlgebra.KeDvrBasis).KMat

	for i := 0; i < int(op.grid.NPoints()); i++ {
		op.hmat.(*mat.Dense).Set(i, i, op.hmat.At(i, i)+vPot[i])
	}
}

func (op *HamiltonianOp[float64]) EvaluateOp() *mat.Dense {
	return nil
}
