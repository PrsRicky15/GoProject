package BasicOneD

import (
	"GoProject/gridData"
	"fmt"
)

func OneDimPoissonSolver() error {
	grid, err := gridData.NewRGrid(-10., 10, 200)
	if err != nil {
		return err
	}

	gaus := gridData.Gaussian[float64]{Sigma: 4., Strength: 1}
	err = grid.PrintPotentToFileRe(gaus, "GaussFunc.dat", "%21.14e")
	if err != nil {
		return err
	}

	Solver1D, err := NewOneDimPoissonSolver(grid, gaus)
	if err != nil {
		return err
	}
	Solver1D.Print()
	return nil
}

func BarrierPotential(v0x, aPoint []float64) error {
	grid, err := gridData.NewRGrid(aPoint[0]-5, aPoint[len(aPoint)-1], 100)
	fmt.Println(grid, err)
	if err != nil {
		return err
	}

	myfunc := gridData.NewRectBarrier(v0x, aPoint)

	fgrid := make([]float64, grid.NPoints())
	grid.FunctionOnGridInPlace(myfunc, fgrid)

	err = grid.PrintVectorToFileRe(fgrid, "barrier.dat", "%21.14e")
	if err != nil {
		return err
	}

	NewFiniteBarrier(grid, myfunc, 1.)

	return nil
}

func GaussianBarrier() error {
	grid, err := gridData.NewRGrid(-20, 20, 200)
	if err != nil {
		return err
	}
	gauss := gridData.Gaussian[float64]{Strength: 1., Sigma: 2}
	err = grid.PrintPotentToFileRe(gauss, "GaussBarrier.dat", "%21.14e")
	if err != nil {
		return err
	}
	fmt.Println("gauss:", gauss)

	gaussBar := NewFiniteBarrier(grid, gauss, 1.)
	gaussBar.ConvertFuncToBarrier(20)
	gaussBar.DisplayMinMaxBarrier()
	err = gaussBar.PrintFuncToDeltaToFile("GaussBarrier.dat", "%21.14e")
	if err != nil {
		return err
	}
	return nil
}

func SuperGaussianBarrier() {
	grid, err := gridData.NewRGrid(-20, 20, 200)
	if err != nil {
		panic(err)
	}
	gauss := gridData.SuperGaussian[float64]{Strength: 1., Sigma: 2, Order: 6}
	err = grid.PrintPotentToFileRe(gauss, "SuperGaussBarrier.dat", "%21.14e")
	if err != nil {
		panic(err)
	}
	fmt.Println("gauss:", gauss)
}

func CompositeFunction() {
	grid, err := gridData.NewRGrid(-10, 10, 200)
	if err != nil {
		panic(err)
	}

	funcSin := gridData.NewSin(1., 3, 0.)
	funcGauss := gridData.Gaussian[float64]{Sigma: 2, Strength: 1.}
	funcMult := gridData.NewProductFunc(funcSin, funcGauss)

	fgrid := make([]float64, grid.NPoints())
	grid.FunctionOnGridInPlace(funcMult, fgrid)

	err = grid.PrintVectorToFileRe(fgrid, "barrier.dat", "%21.14e")
	if err != nil {
		return
	}
}
