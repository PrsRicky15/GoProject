package BasicOneD

import (
	"GoProject/gridData"
	"errors"
	"fmt"
	"math"
)

type PoissonSolver1D struct {
	grid *gridData.RadGrid
	fx   gridData.Rfunc

	maxIter int
	tol     float64

	nMin1  int
	dx2    float64
	rho    []float64
	phiPre []float64
	phiNew []float64
}

func NewOneDimPoissonSolver(grid *gridData.RadGrid, fx gridData.Rfunc) (*PoissonSolver1D, error) {

	if grid == nil || fx == nil {
		return nil, errors.New("grid or function is not defined")
	}

	return &PoissonSolver1D{
		grid: grid,
		fx:   fx,
	}, nil
}

func (PS1D *PoissonSolver1D) Redefine(grid *gridData.RadGrid, fx gridData.Rfunc) error {
	if grid == nil {
		return errors.New("grid cannot be nil")
	}
	if fx == nil {
		return errors.New("function cannot be nil")
	}

	PS1D.grid = grid
	PS1D.fx = fx
	return nil
}

func (PS1D *PoissonSolver1D) RedefineFunc(fx gridData.Rfunc) error {
	if fx == nil {
		return errors.New("function cannot be nil")
	}

	PS1D.fx = fx
	return nil
}

func (PS1D *PoissonSolver1D) Print() {
	fmt.Println(PS1D.fx)
	fmt.Println(PS1D.grid)
}

func (PS1D *PoissonSolver1D) Initialize(maxIter int, tolerance, phiA, phiB float64) error {
	if maxIter < 1 {
		return errors.New("maxIter must be greater than zero")
	}
	if tolerance <= 0 || tolerance > 1 {
		return errors.New("tolerance must be between 0 and 1")
	}
	if math.IsNaN(phiA) || math.IsInf(phiA, 0) {
		return errors.New("phiA must be a finite number")
	}
	if math.IsNaN(phiB) || math.IsInf(phiB, 0) {
		return errors.New("phiB must be a finite number")
	}

	PS1D.maxIter = maxIter
	PS1D.tol = tolerance
	PS1D.dx2 = PS1D.grid.DeltaR() * PS1D.grid.DeltaR()
	PS1D.nMin1 = int(PS1D.grid.NPoints() - 1)

	// Boundary condition
	PS1D.phiPre = make([]float64, PS1D.grid.NPoints())
	PS1D.phiNew = make([]float64, PS1D.grid.NPoints())
	PS1D.phiPre[0] = phiA
	PS1D.phiPre[len(PS1D.phiPre)-1] = phiB

	copy(PS1D.phiNew, PS1D.phiPre)

	PS1D.rho = make([]float64, PS1D.grid.NPoints())
	PS1D.grid.FunctionOnGridInPlace(PS1D.fx, PS1D.rho)

	return nil
}

func (PS1D *PoissonSolver1D) IteratingStepMethod() ([]float64, int, error) {

	for iter := 0; iter < PS1D.maxIter; iter++ {

		for i := 1; i < PS1D.nMin1; i++ {
			PS1D.phiNew[i] = 0.5 * (PS1D.phiPre[i+1] + PS1D.phiPre[i-1] - PS1D.rho[i]*PS1D.dx2)
		}

		errL2Norm := 0.0
		for i := 1; i < PS1D.nMin1; i++ {
			errL2Norm += (PS1D.phiNew[i] - PS1D.phiPre[i]) * (PS1D.phiNew[i] - PS1D.phiPre[i])
		}

		copy(PS1D.phiPre, PS1D.phiNew)

		if math.Sqrt(errL2Norm) < PS1D.tol {
			return PS1D.phiNew, iter + 1, nil
		}
	}

	return nil, PS1D.maxIter, fmt.Errorf("did not converge in %d iterations", PS1D.maxIter)
}

func (PS1D *PoissonSolver1D) IteratingStepWithOmega(omega float64) ([]float64, int, error) {
	for iter := 0; iter < PS1D.maxIter; iter++ {
		for i := 1; i < PS1D.nMin1; i++ {
			PS1D.phiNew[i] = 1 / (2 + omega) * (PS1D.phiPre[i+1] - omega*PS1D.phiPre[i] +
				PS1D.phiPre[i-1] - PS1D.rho[i]/PS1D.dx2)
		}

		errL2Norm := 0.0
		for i := 1; i < PS1D.nMin1; i++ {
			errL2Norm += (PS1D.phiNew[i] - PS1D.phiPre[i]) * (PS1D.phiNew[i] - PS1D.phiPre[i])
		}

		copy(PS1D.phiPre, PS1D.phiNew)

		if math.Sqrt(errL2Norm) < PS1D.tol {
			return PS1D.phiNew, iter + 1, nil
		}
	}

	return nil, PS1D.maxIter, fmt.Errorf("did not converge in %d iterations", PS1D.maxIter)
}

func (PS1D *PoissonSolver1D) CheckError() error {

	err := PS1D.Initialize(1000, 5.0e-3, 0., 0.)
	if err != nil {
		return err
	}
	phi, iter, err := PS1D.IteratingStepMethod()

	err = PS1D.Initialize(1000, 5.0e-3, 0., 0.)
	if err != nil {
		return err
	}

	phiOmg, iterOmg, errOmg := PS1D.IteratingStepWithOmega(1)

	fmt.Printf("Error Without Omega:\t %v %v\n", iter, err)
	fmt.Printf("Error With Omega:\t %v %v\n", iterOmg, errOmg)

	err = PS1D.grid.PrintVectorToFileRe(phi, "phiWO_omega.dat", "%21.14e")
	if err != nil {
		return err
	}

	err = PS1D.grid.PrintVectorToFileRe(phiOmg, "phi_with_omega.dat", "%21.14e")
	if err != nil {
		return err
	}

	return nil
}
