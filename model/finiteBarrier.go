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

func Barrier(grid gridData.RadGrid, f func(x float64) float64) {
	fgrid := make([]float64, grid.NPoints())
	wrapped := gridData.FuncWrapper{Fn: f}
	grid.FunctionOnGridInPlace(wrapped, fgrid)

	err := grid.PrintVectorToFileRe(fgrid, "barrier.dat", "%21.14e")
	if err != nil {
		return
	}
}
