package model

import "GoProject/gridData"

// Rectangle Numerical integration
func Rectangle(f func(x float64) float64, a, b float64, nPoints int) float64 {
	dx := (b - a) / float64(nPoints)
	sum := 0.0
	for i := 0; i < nPoints; i++ {
		x := a + float64(i)*dx
		sum += f(x)
	}
	return sum * dx
}

// Trapezoidal Numerical integration
func Trapezoidal(f func(x float64) float64, a, b float64, nPoints int) float64 {
	dx := (b - a) / float64(nPoints)
	sum := 0.0
	for i := 1; i < nPoints-1; i++ {
		x := a + float64(i)*dx
		sum += f(x)
	}
	return sum*dx + 0.5*(f(a)+f(b))
}

// Simpson Numerical integration
func Simpson(f func(x float64) float64, a, b float64, nPoints int) float64 {
	dx := (b - a) / float64(nPoints)
	sum := 0.0
	for i := 1; i < nPoints-1; i++ {
		x := a + float64(i)*dx
		sum += f(x)
	}
	return sum*dx + 0.5*(f(a)+f(b))
}

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
	wrapped := gridData.FuncWrapper{Fn: myfunc.EvaluateAt}
	grid.FunctionOnGridInPlace(wrapped, fgrid)

	err = grid.PrintVectorToFileRe(fgrid, "barrier.dat", "%21.14e")
	if err != nil {
		return
	}
}

func CompositeFunction() {
	grid, err := gridData.NewRGrid(-10, 10, 200)
	if err != nil {
		panic(err)
	}

	funcSin := gridData.NewSin(1., 3, 0.)
	funcGauss := gridData.Gaussianf64{Sigma: 2, Strength: 1.}
	funcMult := gridData.NewProductFunc(funcSin, funcGauss)

	fgrid := make([]float64, grid.NPoints())
	wrapped := gridData.FuncWrapper{Fn: funcMult.EvaluateAt}
	grid.FunctionOnGridInPlace(wrapped, fgrid)

	err = grid.PrintVectorToFileRe(fgrid, "barrier.dat", "%21.14e")
	if err != nil {
		return
	}
}
