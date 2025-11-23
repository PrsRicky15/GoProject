package model

import (
	"GoProject/gridData"
	"fmt"
	"math"
	"strings"
)

// Rectangle Numerical integration
func Rectangle(f func(x float64) float64, a, b float64, nPoints int) float64 {
	dx := (b - a) / float64(nPoints)
	sum := 0.0
	x := a
	for i := 0; i < nPoints; i++ {
		sum += f(x)
		x += dx
	}
	return sum * dx
}

// Trapezoidal Numerical integration
func Trapezoidal(f func(x float64) float64, a, b float64, nPoints int) float64 {
	dx := (b - a) / float64(nPoints)
	sum := 0.0
	for i := 1; i < nPoints; i++ {
		x := a + float64(i)*dx // Calculate x based on index
		sum += f(x)
	}
	return (sum + 0.5*(f(a)+f(b))) * dx
}

// Simpson Numerical integration
func Simpson(f func(x float64) float64, a, b float64, nPoints int) float64 {
	dx := (b - a) / float64(nPoints)
	sumOdd := 0.0
	sumEven := 0.0

	for i := 1; i < nPoints; i += 2 {
		x := a + float64(i)*dx
		sumOdd += f(x)
	}

	for i := 2; i < nPoints; i += 2 {
		x := a + float64(i)*dx
		sumEven += f(x)
	}

	return dx / 3.0 * (f(a) + 4*sumOdd + 2*sumEven + f(b))
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
	grid.FunctionOnGridInPlace(myfunc, fgrid)

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
	grid.FunctionOnGridInPlace(funcMult, fgrid)

	err = grid.PrintVectorToFileRe(fgrid, "barrier.dat", "%21.14e")
	if err != nil {
		return
	}
}

// Test functions with known analytical solutions
func testFunc1(x float64) float64 {
	return x * x // integral from 0 to 1 = 1/3
}

func AnalyzeConvergence() {
	f := testFunc1
	a, b := 0.0, 1.0
	expected := 1.0 / 3.0

	gridPoints := []int{10, 20, 50, 100, 200, 500, 1000, 2000, 5000}

	fmt.Println("=== Error Analysis: Integrating x^2 from 0 to 1 (exact = 0.333333) ===")
	fmt.Printf("%-10s %-20s %-20s %-20s\n", "N Points", "Rectangle Error", "Trapezoidal Error", "Simpson Error")
	fmt.Println(strings.Repeat("-", 75))

	for _, n := range gridPoints {
		rectResult := Rectangle(f, a, b, n)
		trapResult := Trapezoidal(f, a, b, n)
		simpResult := Simpson(f, a, b, n)

		rectError := math.Abs(rectResult - expected)
		trapError := math.Abs(trapResult - expected)
		simpError := math.Abs(simpResult - expected)

		fmt.Printf("%-10d %-20.10e %-20.10e %-20.10e\n", n, rectError, trapError, simpError)
		fmt.Printf("%-10d %-20.10e %-20.10e %-20.10e\n", n, rectResult, trapResult, simpResult)
	}

	// Convergence rate analysis
	fmt.Println("=== Convergence Rate Analysis ===")
	fmt.Println("Rectangle Method: Expected O(h) = O(1/N)")
	fmt.Println("Trapezoidal Method: Expected O(h^2) = O(1/N^2)")
	fmt.Println("Simpson Method: Expected O(h^4) = O(1/N^4)")

	// Calculate observed convergence rates
	n1, n2 := 100, 1000
	rect1 := math.Abs(Rectangle(f, a, b, n1) - expected)
	rect2 := math.Abs(Rectangle(f, a, b, n2) - expected)
	trap1 := math.Abs(Trapezoidal(f, a, b, n1) - expected)
	trap2 := math.Abs(Trapezoidal(f, a, b, n2) - expected)
	simp1 := math.Abs(Simpson(f, a, b, n1) - expected)
	simp2 := math.Abs(Simpson(f, a, b, n2) - expected)

	rectRate := math.Log(rect1/rect2) / math.Log(float64(n2)/float64(n1))
	trapRate := math.Log(trap1/trap2) / math.Log(float64(n2)/float64(n1))
	simpRate := math.Log(simp1/simp2) / math.Log(float64(n2)/float64(n1))

	fmt.Printf("\nObserved convergence rates (from N=%d to N=%d):\n", n1, n2)
	fmt.Printf("Rectangle: %.2f (expected ~1.0)\n", rectRate)
	fmt.Printf("Trapezoidal: %.2f (expected ~2.0)\n", trapRate)
	fmt.Printf("Simpson: %.2f (expected ~4.0)\n", simpRate)
}
