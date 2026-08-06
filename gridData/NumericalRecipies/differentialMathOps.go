package NumericalRecipies

import (
    "fmt"
	"GoProject/gridData"
	"reflect"
)

const DefaultDelta = 1e-6

type Differentiator1D interface {
	diffrentiateAtX(f gridData.RadGrid, x float64, delta float64, order int) float64
	Name() string
}

type Differentiator2D interface {
	differentiateAtX(f gridData.RadGrid, x float64, delta float64, order int) float64
	Name() string
}

type FiniteDifference struct{
	pointFormula        int 
}

type CentralDifferenceStencil struct {
    firstDerivativeCoeffs map[int][]float64
    secondDervativeCoeffs map[int][]float64
}

func (cds *CentralDifferenceStencil) newCentranlDiffefferenceStencil() *CentralDifferenceStencil {
    return &CentralDifferenceStencil{
        firstDerivativeCoeffs : map[int][]float64{
	        3: {1.0 / 2.0},                                                         // 3-point (only positive side)
	        5: {2.0 / 3.0, -1.0 / 12.0},                                            // 5-point
	        7: {3.0 / 4.0, -3.0 / 20.0, 1.0 / 60.0},                                // 7-point
	        9: {4.0 / 5.0, -1.0 / 5.0, 4.0 / 105.0, -1.0 / 280.0},                  // 9-point
        },
        secondDervativeCoeffs : map[int][]float64{
	        3: {-2.0, 1.0},                                                         // 3-point
	        5: {-5.0 / 2.0, 4.0 / 3.0, -1.0 / 12.0},                                // 5-point
	        7: {-49.0 / 18.0, 3.0 / 2.0, -3.0 / 20.0, 1.0 / 90.0},                  // 7-point
	        9: {-205.0 / 72.0, 8.0 / 5.0, -1.0 / 5.0, 8.0 / 315.0, -1.0 / 560.0},   // 9-point
        },
    }
}

var firstDerivativeStencils = map[string]map[int][]float64{
    "central": {
        3: {1.0 / 2.0},                                                         // 3-point (only positive side)
	    5: {2.0 / 3.0, -1.0 / 12.0},                                            // 5-point
	    7: {3.0 / 4.0, -3.0 / 20.0, 1.0 / 60.0},                                // 7-point
	    9: {4.0 / 5.0, -1.0 / 5.0, 4.0 / 105.0, -1.0 / 280.0}, 
    },
    "forward": {
        1: {-1.0, 1.0},                                               // 2-point forward difference
        3: {-3.0 / 2.0, 2.0, -1.0 / 2.0},
        5: {-11 / 6.0, 3.0, -3.0 / 2.0, 1.0 / 3.0},
        7: {-25 / 12.0, 4.0, -3.0, 4.0 / 3.0, -1.0 / 4.0},
        9: {-137 / 60.0, 5.0, -5.0, 10.0 / 3.0, -5.0 / 4.0, 1.0 / 5.0},
    },
    "backward": {
        1: {-1.0, 1.0},                                               // 2-point forward difference
        3: {1.0 / 2.0, -2.0, 3.0 / 2.0},
        5: {-1 / 3.0, 3.0 / 2.0, -3.0, 11.0 / 6.0},
        7: {1.0 / 4.0, -4.0 / 3.0, 3.0, -4.0, 25.0 / 12.0},
        9: {-1.0 / 5.0, 5.0 / 4.0, -10.0 / 3.0, 5.0, -5.0, 137.0 / 60.0},
    },
}
var secondDerivativeStencils = map[string]map[int][]float64{
    // TODO
}


func (fd *FiniteDifference) Name() string {
    return fmt.Sprintf("Finite Difference (order %d)", fd.pointFormula)
}

func (fd *FiniteDifference) differetiateAtX(f gridData.RadGrid, x float64, delta float64, order int) float64 {
    var derivative float64
    if order  > 2 {
        panic("The higher-order (>2) derivatives are not available.")
    }

    var stencilCoeffs map[string]map[int][]float64
    switch order {
    case 1: stencilCoeffs = firstDerivativeStencils
    case 2: stencilCoeffs = secondDerivativeStencils
    }

    fn := f.FunctionAtx()
    sign := 1.0
    for i, coeff := range stencilCoeffs["central"][order]{
        derivative += coeff * f[]
        derivative += sign * coeff * f(x - (i + 1) * delta)
    }
    return derivative
}


// the central difference formula:
//
//     f'(x) ≈ (f(x+δ) - f(x-δ)) / (2δ)


type GradientScalarOps interface {
	calculateGradientScalarAtX(x float64) float64 // univariate first derivative
}

type GradientScalar struct{
	fn  func(float64) float64;
	delta float64
}


func (g *GradientScalar) NewDefine(f func(float64) float64, delta float64) *GradientScalar {
	return &GradientScalar{fn: f, delta: delta}
}

func (g GradientScalar) centralDifference(f func(float64) float64, x, delta float64) float64 {
	return (f(x+delta) - f(x-delta)) / (2 * delta)
}

func (g *GradientScalar) calculateGradientScalarAtX(x float64) float64 {
    if g.delta == 0 {
        g.delta = DefaultDelta
    }
    nd := (&NumericalDifferentiation[float64]{}).NewDefine(g.fn, g.delta)
    return nd.centralDifference(nd.fn, x, nd.delta)
}

// Gradient computes the gradient vector of a scalar‑valued multivariate
// function f at the point x.

type GradientOps interface {
	GradientScalarOps
	calculateGradientAtX(x []float64) []float64  // multivariate gradient vector
}

type Gradient struct{
	fn func([]float64) float64
	delta float64
}

func (g *Gradient) NewDefine(f func([]float64) float64, delta float64) *Gradient {
	return &Gradient{fn: f, delta: delta}
}

func (g Gradient) centralDifference(f func([]float64) float64, x []float64, delta float64) []float64{
	n := len(x)
	grad := make([]float64, n)

	// make a copy of the original point so we can restore it
    orig := make([]float64, n)
    copy(orig, x)

    for j := range n {
        x[j] = orig[j] + delta
        fPlus := f(x)
        x[j] = orig[j] - delta
        fMinus := f(x)
        grad[j] = (fPlus - fMinus) / (2 * delta)
        x[j] = orig[j]
    }
	return grad
}

func (g *Gradient) calculateGradientAtX(x []float64) []float64 {
    if g.delta == 0 {
        g.delta = DefaultDelta
    }
    n := len(x)
	grad := make([]float64, n)

	// make a copy of the original point so we can restore it
    orig := make([]float64, n)
    copy(orig, x)

    nd := (&NumericalDifferentiation[[]float64]{}).NewDefine(g.fn, g.delta)

    for j := range n {
        x[j] = orig[j] + g.delta
        fPlus := g.fn(x)
        x[j] = orig[j] - g.delta
        fMinus := g.fn(x)
        grad[j] = (fPlus - fMinus) / (2 * delta)
        x[j] = orig[j]
    }
	return grad
}

func (g *Gradient) calculateGradientScalarAtX(x float64) float64 {
	if g.delta == 0 {
        g.delta = DefaultDelta
    }
	var floatchecker float64
	if reflect.TypeOf(x) != reflect.TypeOf(floatchecker) {
		panic("GradientScalarAtX() can only be used for univariate funcions.\n Use GradientAtX for multivariate functions.")
	}
	y := []float64{x}  
    result := g.centralDifference(g.fn, y, g.delta)
	return result[0]
}


type JacobianOps interface {
	GradientOps
	calculateJacobianAtX(x []float64) [][]float64  // multivariate gradient vector
}

type Jacobian struct {
	fn func([]float64) []float64
	delta float64
}

func (g *Jacobian) NewDefine(f func([]float64) []float64, delta float64) *Jacobian {
	return &Jacobian{fn: f, delta: delta}
}

// Jacobian computes the Jacobian matrix of a vector‑valued function
func (g *Jacobian) calculateJacobianAtX(x []float64) [][]float64 {
    n := len(x)
    if g.delta == 0 {
        g.delta = DefaultDelta
    }

    f0 := g.fn(x)
    m := len(f0)
    J := make([][]float64, m)
    for i := 0; i < m; i++ {
        J[i] = make([]float64, n)
    }

    orig := make([]float64, n)
    copy(orig, x)

    for j := 0; j < n; j++ {
        x[j] = orig[j] + g.delta
        fPlus := g.fn(x)
        x[j] = orig[j] - g.delta
        fMinus := g.fn(x)
        for i := 0; i < m; i++ {
            J[i][j] = (fPlus[i] - fMinus[i]) / (2 * g.delta)
        }
        x[j] = orig[j]
    }
    return J
}

