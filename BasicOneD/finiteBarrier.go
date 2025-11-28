package BasicOneD

import (
	"GoProject/gridData"
	"fmt"
)

func Barrier() {
	grid, err := gridData.NewRGrid(-3, 8, 100)
	if err != nil {
		panic(err)
	}

	v0x := make([]float64, 3)
	v0x[0] = 0
	v0x[1] = 1
	v0x[2] = -0.5

	aPoint := make([]float64, 3)
	aPoint[0] = 0
	aPoint[1] = 5
	aPoint[2] = 8

	myfunc := gridData.NewRectBarrier(v0x, aPoint)

	fgrid := make([]float64, grid.NPoints())
	grid.FunctionOnGridInPlace(myfunc, fgrid)

	err = grid.PrintVectorToFileRe(fgrid, "barrier.dat", "%21.14e")
	if err != nil {
		return
	}
}

func GaussianBarrier() {
	grid, err := gridData.NewRGrid(-20, 20, 200)
	if err != nil {
		panic(err)
	}
	gauss := gridData.Gaussian[float64]{Strength: 1., Sigma: 2}
	err = grid.PrintPotentToFileRe(gauss, "GaussBarrier.dat", "%21.14e")
	if err != nil {
		panic(err)
	}
	fmt.Println("gauss:", gauss)
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
