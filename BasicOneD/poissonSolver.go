package BasicOneD

import (
	"GoProject/gridData"
	"errors"
	"fmt"
)

type PoissonSolver1D struct {
	grid *gridData.RadGrid
	fx   gridData.Rfunc
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
