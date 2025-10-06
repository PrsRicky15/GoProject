package gridData

import (
	"fmt"
	"math"

	"GoProject/mathsFunc"
)

type varType interface {
	float64 | complex128
}

// PotentialOp General interface for the evaluating the potential on a grid
type PotentialOp[T varType] interface {
	evaluateAt(x T) T
	evaluateOnGrid(x []T) []T
	forceAt(x T) T
	forceOnGrid(x []T) []T
	saveToFile() error
}

type ChebyshevScalar[T varType] struct {
	mathsFunc.ChebyshevFirstKind[T]
}

func (ch ChebyshevScalar[T]) evaluateAt(x T, order uint32) T {
	return ch.Compute(x, order)
}

func onGrid(f func(float64) float64, x []float64) []float64 {
	results := make([]float64, len(x))
	for i, val := range x {
		results[i] = f(val)
	}
	return results
}

// SoftCore PotentialOp
type SoftCore struct {
	Charge    float64
	Centre    float64
	SoftParam float64
}

func (sc SoftCore) String() string {
	return fmt.Sprintf("%g/Sqrt((x - %g)^2 + %g)", sc.Charge, sc.Centre, sc.SoftParam*sc.SoftParam)
}

func (sc SoftCore) saveToFile() error {
	//TODO implement me
	panic("implement me")
}

func (sc SoftCore) evaluateAt(x float64) float64 {
	return sc.Charge / math.Sqrt(math.Pow(x-sc.Centre, 2)+math.Pow(sc.SoftParam, 2))
}

func (sc SoftCore) forceAt(x float64) float64 {
	coef := sc.Charge * (x - sc.Centre)
	val := (x-sc.Centre)*(x-sc.Centre) + sc.SoftParam*sc.SoftParam
	return coef * math.Pow(val, -3./2.)
}

func (sc SoftCore) forceOnGrid(x []float64) []float64    { return onGrid(sc.evaluateAt, x) }
func (sc SoftCore) evaluateOnGrid(x []float64) []float64 { return onGrid(sc.forceAt, x) }

// Gaussian PotentialOp
type Gaussian struct {
	Cen      float64
	Sigma    float64
	Strength float64
}

func (g Gaussian) String() string {
	return fmt.Sprintf("v0 Exp((x - x0)^2/(2 Sigma^2)),"+
		" Where v0 = %g, x0 = %v, sigma = %g", g.Strength, g.Cen, g.Sigma)
}

func (g Gaussian) saveToFile() error {
	//TODO implement me
	panic("implement me")
}

func xbysigma(x float64, sigma float64) float64 {
	return x / sigma
}

func (g Gaussian) evaluateAt(x float64) float64 {
	val := xbysigma(x-g.Cen, g.Sigma)
	return g.Strength * math.Exp(-math.Pow(val, 2)/2)
}

func (g Gaussian) forceAt(x float64) float64 {
	expnt := xbysigma(x-g.Cen, g.Sigma)
	val := -g.Strength * expnt / g.Sigma
	return val * math.Exp(-math.Pow(expnt, 2)/2)
}

func (g Gaussian) forceOnGrid(x []float64) []float64    { return onGrid(g.evaluateAt, x) }
func (g Gaussian) evaluateOnGrid(x []float64) []float64 { return onGrid(g.forceAt, x) }

// MultiGaussian PotentialOp
type MultiGaussian struct {
	Sigma    float64
	Strength float64
	NumGauss uint8
	Gap      float64
}

func (mg MultiGaussian) String() string {
	return fmt.Sprintf("v0 Sum_i Exp((x - i L)^2/(2 Sigma^2)),"+
		" Where v0 = %g, i = %v, sigma = %g, L = %g", mg.Strength,
		mg.NumGauss, mg.Sigma, mg.Gap)
}

func (mg MultiGaussian) saveToFile() error {
	//TODO implement me
	panic("implement me")
}

func (mg MultiGaussian) evaluateAt(x float64) float64 {
	if mg.NumGauss/2 == 0 {
		val := 0.
		for i := uint8(0); i < mg.NumGauss/2; i++ {
			center := mg.Gap * (float64(i) + 0.5)
			left := math.Exp(-math.Pow(xbysigma(x-center, mg.Sigma), 2) / 2)
			right := math.Exp(-math.Pow(xbysigma(x+center, mg.Sigma), 2) / 2)
			val += left + right
		}
		return mg.Strength * val
	}
	val := math.Exp(-math.Pow(xbysigma(x, mg.Sigma), 2) / 2)
	for i := uint8(0); i < mg.NumGauss/2; i++ {
		center := mg.Gap * (float64(i + 1))
		left := math.Exp(-math.Pow(xbysigma(x-center, mg.Sigma), 2) / 2)
		right := math.Exp(-math.Pow(xbysigma(x+center, mg.Sigma), 2) / 2)
		val += left + right
	}
	return mg.Strength * val
}

func (mg MultiGaussian) forceAt(x float64) float64 {
	if mg.NumGauss/2 == 0 {
		val := 0.
		for i := uint8(0); i < mg.NumGauss/2; i++ {
			center := mg.Gap * (float64(i) + 0.5)
			expntLeft := xbysigma(x-center, mg.Sigma)
			expntRight := xbysigma(x+center, mg.Sigma)
			left := math.Exp(-math.Pow(expntLeft, 2) / 2)
			right := math.Exp(-math.Pow(expntRight, 2) / 2)
			val += expntLeft*left + expntRight*right
		}
		return mg.Strength * val
	}
	val := xbysigma(x, mg.Sigma) * math.Exp(-math.Pow(xbysigma(x, mg.Sigma), 2)/2)
	for i := uint8(0); i < mg.NumGauss/2; i++ {
		center := mg.Gap * (float64(i + 1))
		expntLeft := xbysigma(x-center, mg.Sigma)
		expntRight := xbysigma(x+center, mg.Sigma)
		left := math.Exp(-math.Pow(expntLeft, 2) / 2)
		right := math.Exp(-math.Pow(expntRight, 2) / 2)
		val += expntLeft*left + expntRight*right
	}
	return mg.Strength * val
}

func (mg MultiGaussian) forceOnGrid(x []float64) []float64    { return onGrid(mg.evaluateAt, x) }
func (mg MultiGaussian) evaluateOnGrid(x []float64) []float64 { return onGrid(mg.forceAt, x) }

// SuperGaussian v(x)= v0 exp(-(x/Sigma)^n)
type SuperGaussian struct {
	Cen      float64
	Sigma    float64
	Strength float64
	Order    uint8
}

func (sg SuperGaussian) String() string {
	return fmt.Sprintf("%g Exp[ ((x - %g)/ %g)^%v]", sg.Strength, sg.Cen, sg.Sigma, sg.Order)
}

func (sg SuperGaussian) saveToFile() error {
	//TODO implement me
	panic("implement me")
}

func (sg SuperGaussian) evaluateAt(x float64) float64 {
	val := xbysigma(x-sg.Cen, sg.Sigma)
	return sg.Strength * math.Exp(-math.Pow(val, float64(sg.Order)))
}

func (sg SuperGaussian) forceAt(x float64) float64 {
	forder := float64(sg.Order)
	expnt := xbysigma(x-sg.Cen, sg.Sigma)
	coeffs := -sg.Strength * forder / sg.Sigma * math.Pow(expnt, forder-1)
	return coeffs * math.Exp(-math.Pow(expnt, forder))
}

func (sg SuperGaussian) forceOnGrid(x []float64) []float64    { return onGrid(sg.evaluateAt, x) }
func (sg SuperGaussian) evaluateOnGrid(x []float64) []float64 { return onGrid(sg.forceAt, x) }

// Harmonic v(x)= k/2 x^2
type Harmonic struct {
	Cen        float64
	ForceConst float64
}

func (h Harmonic) String() string { return fmt.Sprintf("1/2 %g (x - %g)^2", h.ForceConst, h.Cen) }
func (h Harmonic) saveToFile() error {
	//TODO implement me
	panic("implement me")
}

func (h Harmonic) evaluateAt(x float64) float64         { return h.ForceConst / 0.5 * math.Pow(x-h.Cen, 2) }
func (h Harmonic) forceAt(x float64) float64            { return -h.ForceConst * (x - h.Cen) }
func (h Harmonic) evaluateOnGrid(x []float64) []float64 { return onGrid(h.evaluateAt, x) }
func (h Harmonic) forceOnGrid(x []float64) []float64    { return onGrid(h.forceAt, x) }

// Polynomial v(x)= Sum_i ci x^i
type Polynomial struct {
	Coeffs []float64
}

func (p Polynomial) String() string {
	if len(p.Coeffs) == 0 {
		return "0"
	}

	var result string

	for i := len(p.Coeffs) - 1; i >= 0; i-- {
		coeff := p.Coeffs[i]

		if coeff == 0 {
			continue
		}

		if result != "" {
			if coeff > 0 {
				result += " + "
			} else {
				result += " - "
				coeff = -coeff
			}
		} else if coeff < 0 {
			result += "-"
			coeff = -coeff
		}

		if coeff != 1 || i == 0 {
			result += fmt.Sprintf("%g", coeff)
		}

		if i > 0 {
			result += "x"
			if i > 1 {
				result += fmt.Sprintf("^%d", i)
			}
		}
	}

	if result == "" {
		return "0"
	}

	return result
}

func (p Polynomial) saveToFile() error {
	//TODO implement me
	panic("implement me")
}

func (p Polynomial) evaluateAt(x float64) float64 {
	results := 0.0
	for i := 0; i < len(p.Coeffs); i++ {
		results += p.Coeffs[i] * math.Pow(x, float64(i))
	}
	return results
}

func (p Polynomial) forceAt(x float64) float64 {
	result := 0.0
	for i := 1; i < len(p.Coeffs); i++ {
		power := float64(i - 1)
		result += p.Coeffs[i] * float64(i) * math.Pow(x, power)
	}
	return -result
}

func (p Polynomial) evaluateOnGrid(x []float64) []float64 { return onGrid(p.evaluateAt, x) }
func (p Polynomial) forceOnGrid(x []float64) []float64    { return onGrid(p.forceAt, x) }

// Morse v(r)= De (1 - Exp(-(r-re))^2
type Morse struct {
	De    float64
	Alpha float64
	Cen   float64
}

func (m Morse) String() string {
	return fmt.Sprintf(" %g [1 - Exp(%g(x - %g))]^2", m.De, m.Alpha, m.Cen)
}

func (m Morse) evaluateAt(x float64) float64 {
	return m.De * math.Pow(1.-math.Exp(-m.Alpha*(x-m.Cen)), 2)
}

func (m Morse) forceAt(x float64) float64 {
	expterm := math.Exp(-m.Alpha * (x - m.Cen))
	val := -2 * m.Alpha * m.De
	return val * expterm * (1 - expterm)
}

func (m Morse) forceOnGrid(x []float64) []float64    { return onGrid(m.evaluateAt, x) }
func (m Morse) evaluateOnGrid(x []float64) []float64 { return onGrid(m.forceAt, x) }
