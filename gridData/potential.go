package gridData

import (
	"fmt"
	"math"
	"math/cmplx"
)

type VarType interface {
	float64 | complex128
}

// PotentialOp General interface for the evaluating the potential on a grid
type PotentialOp[T VarType] interface {
	EvaluateAt(x T) T
	EvaluateOnGrid(x []T) []T
	EvaluateOnGridInPlace(fn, x []T)
	ForceAt(x T) T
	ForceOnGrid(x []T) []T
	ForceOnGridInPlace(fn, x []T)
}

// Made generic to work with VarType
func onGrid[T VarType](f func(T) T, x []T) []T {
	results := make([]T, len(x))
	for i, val := range x {
		results[i] = f(val)
	}
	return results
}

// Made generic to work with VarType
func onGridInPlace[T VarType](f func(T) T, fn, x []T) {
	for i, val := range x {
		fn[i] = f(val)
	}
}

// Morse v(r)= De (1 - Exp(-(r-re))^2
type Morse[T VarType] struct {
	De    float64
	Alpha float64
	Cen   float64
}

func (m Morse[T]) String() string {
	return fmt.Sprintf(" De [1 - Exp(-a(x - x0))]^2, where, De: %g, a: %g, x0: %g", m.De, m.Alpha, m.Cen)
}

func (m Morse[T]) EvaluateAt(x T) T {
	var result any

	switch any(x).(type) {
	case float64:
		mf64 := MorseF64(m)
		result = mf64.evaluateAt(any(x).(float64))
	case complex128:
		mz64 := MorseZ64(m)
		result = mz64.evaluateAt(any(x).(complex128))
	default:
		panic("unsupported type")
	}

	return result.(T)
}

func (m Morse[T]) ForceAt(x T) T {
	var result any

	switch any(x).(type) {
	case float64:
		mf64 := MorseF64(m)
		result = mf64.forceAt(any(x).(float64))
	case complex128:
		mz64 := MorseZ64(m)
		result = mz64.forceAt(any(x).(complex128))
	default:
		panic("unsupported type")
	}

	return result.(T)
}

func (m Morse[T]) EvaluateOnGrid(x []T) []T        { return onGrid(m.EvaluateAt, x) }
func (m Morse[T]) EvaluateOnGridInPlace(fn, x []T) { onGridInPlace(m.EvaluateAt, fn, x) }
func (m Morse[T]) ForceOnGrid(x []T) []T {
	return onGrid(m.ForceAt, x)
}
func (m Morse[T]) ForceOnGridInPlace(fn, x []T) { onGridInPlace(m.ForceAt, fn, x) }

// MorseF64 For float64 specialization
type MorseF64 Morse[float64]

func (m MorseF64) evaluateAt(x float64) float64 {
	oneMinusExpnt := 1. - math.Exp(-m.Alpha*(x-m.Cen))
	return m.De * oneMinusExpnt * oneMinusExpnt
}

func (m MorseF64) forceAt(x float64) float64 {
	expterm := math.Exp(-m.Alpha * (x - m.Cen))
	val := -2 * m.Alpha * m.De
	return val * expterm * (1 - expterm)
}

// MorseZ64 For complex64 specialization
type MorseZ64 Morse[complex128]

func (m MorseZ64) evaluateAt(x complex128) complex128 {
	axmx0 := complex(m.Alpha, 0.) * (x - complex(m.Cen, 0.))
	oneMinusExpnt := 1. - cmplx.Exp(-axmx0)
	return complex(m.De, 0.) * oneMinusExpnt * oneMinusExpnt
}

func (m MorseZ64) forceAt(x complex128) complex128 {
	expterm := cmplx.Exp(-complex(m.Alpha, 0.) * (x - complex(m.Cen, 0.)))
	val := complex(-2*m.Alpha*m.De, 0.)
	return val * expterm * (1 - expterm)
}

// SoftCore with generic type parameter
type SoftCore[T VarType] struct {
	Charge    float64
	Centre    float64
	SoftParam float64
}

func (sc SoftCore[T]) String() string {
	return fmt.Sprintf("Za/Sqrt((x - x0)^2 + a^2), where Za=%g, x0=%g, a^2=%g",
		sc.Charge, sc.Centre, sc.SoftParam*sc.SoftParam)
}

func (sc SoftCore[T]) EvaluateAt(x T) T {
	var result any

	switch any(x).(type) {
	case float64:
		s := SoftCoreF64(sc)
		result = s.evaluateAt(any(x).(float64))
	case complex128:
		s := SoftCoreZ64(sc)
		result = s.evaluateAt(any(x).(complex128))
	default:
		panic("unsupported type")
	}

	return result.(T)
}

func (sc SoftCore[T]) ForceAt(x T) T {
	var result any

	switch any(x).(type) {
	case float64:
		s := SoftCoreF64(sc)
		result = s.forceAt(any(x).(float64))
	case complex128:
		s := SoftCoreZ64(sc)
		result = s.forceAt(any(x).(complex128))
	default:
		panic("unsupported type")
	}

	return result.(T)
}

func (sc SoftCore[T]) EvaluateOnGrid(x []T) []T {
	return onGrid(sc.EvaluateAt, x)
}
func (sc SoftCore[T]) ForceOnGrid(x []T) []T {
	return onGrid(sc.ForceAt, x)
}
func (sc SoftCore[T]) EvaluateOnGridInPlace(fn, x []T) { onGridInPlace(sc.EvaluateAt, fn, x) }
func (sc SoftCore[T]) ForceOnGridInPlace(fn, x []T)    { onGridInPlace(sc.ForceAt, fn, x) }

// SoftCoreF64 For float64 specialization
type SoftCoreF64 SoftCore[float64]

func (sc SoftCoreF64) evaluateAt(x float64) float64 {
	return sc.Charge / math.Sqrt(math.Pow(x-sc.Centre, 2)+math.Pow(sc.SoftParam, 2))
}

func (sc SoftCoreF64) forceAt(x float64) float64 {
	coef := sc.Charge * (x - sc.Centre)
	val := (x-sc.Centre)*(x-sc.Centre) + sc.SoftParam*sc.SoftParam
	return coef * math.Pow(val, -3./2.)
}

// SoftCoreZ64 For complex128 specialization
type SoftCoreZ64 SoftCore[complex128]

func (sc SoftCoreZ64) evaluateAt(x complex128) complex128 {
	return complex(sc.Charge, 0.) / cmplx.Sqrt(cmplx.Pow(x-complex(sc.Centre, 0.), 2)+
		complex(math.Pow(sc.SoftParam, 2), 0.))
}

func (sc SoftCoreZ64) forceAt(x complex128) complex128 {
	xmx0 := x - complex(sc.Centre, 0)
	coef := complex(sc.Charge, 0.) * xmx0
	val := xmx0*xmx0 + complex(sc.SoftParam*sc.SoftParam, 0.)
	return coef * cmplx.Pow(val, -3./2.)
}

// Gaussian PotentialOp
type Gaussian[T VarType] struct {
	Cen      float64
	Sigma    float64
	Strength float64
}

func (g Gaussian[T]) String() string {
	return fmt.Sprintf("v0 Exp((x - x0)^2/(2 Sigma^2)),"+
		" Where v0 = %g, x0 = %v, sigma = %g", g.Strength, g.Cen, g.Sigma)
}

func (g Gaussian[T]) EvaluateAt(x T) T {
	var result any

	switch any(x).(type) {
	case float64:
		gf64 := Gaussianf64(g)
		result = gf64.evaluateAt(any(x).(float64))
	case complex128:
		gz64 := GaussianZ64(g)
		result = gz64.evaluateAt(any(x).(complex128))
	default:
		panic("unsupported type")
	}

	return result.(T)
}

func (g Gaussian[T]) ForceAt(x T) T {
	var result any

	switch any(x).(type) {
	case float64:
		gf64 := Gaussianf64(g)
		result = gf64.forceAt(any(x).(float64))
	case complex128:
		gz64 := GaussianZ64(g)
		result = gz64.forceAt(any(x).(complex128))
	default:
		panic("unsupported type")
	}

	return result.(T)
}

func (g Gaussian[T]) EvaluateOnGrid(x []T) []T {
	return onGrid(g.EvaluateAt, x)
}
func (g Gaussian[T]) ForceOnGrid(x []T) []T {
	return onGrid(g.ForceAt, x)
}
func (g Gaussian[T]) EvaluateOnGridInPlace(fn, x []T) { onGridInPlace(g.EvaluateAt, fn, x) }
func (g Gaussian[T]) ForceOnGridInPlace(fn, x []T)    { onGridInPlace(g.ForceAt, fn, x) }

func xBySigma(x float64, sigma float64) float64 {
	return x / sigma
}
func xBySigmaZ64(x complex128, sigma float64) complex128 { return x / complex(sigma, 0) }

type Gaussianf64 Gaussian[float64]

func (g Gaussianf64) evaluateAt(x float64) float64 {
	expnt := xBySigma(x-g.Cen, g.Sigma)
	return g.Strength * math.Exp(-expnt*expnt/2.)
}

func (g Gaussianf64) forceAt(x float64) float64 {
	expnt := xBySigma(x-g.Cen, g.Sigma)
	val := -g.Strength * expnt / g.Sigma
	return val * math.Exp(-expnt*expnt/2)
}

type GaussianZ64 Gaussian[complex128]

func (g GaussianZ64) evaluateAt(x complex128) complex128 {
	val := xBySigmaZ64(x-complex(g.Cen, 0.), g.Sigma)
	return complex(g.Strength, 0.) * cmplx.Exp(-cmplx.Pow(val, 2)/2)
}

func (g GaussianZ64) forceAt(x complex128) complex128 {
	expnt := xBySigmaZ64(x-complex(g.Cen, 0.), g.Sigma)
	val := -complex(g.Strength/g.Sigma, 0.) * expnt
	return val * cmplx.Exp(-cmplx.Pow(expnt, 2)/complex(2, 0.))
}

// MultiGaussian PotentialOp
type MultiGaussian[T VarType] struct {
	Sigma    float64
	Strength float64
	NumGauss uint8
	Gap      float64
}

func (mg MultiGaussian[T]) String() string {
	return fmt.Sprintf("v0 Sum_i Exp((x - i L)^2/(2 Sigma^2)),"+
		" Where v0 = %g, i = %v, sigma = %g, L = %g", mg.Strength,
		mg.NumGauss, mg.Sigma, mg.Gap)
}

func (mg MultiGaussian[T]) EvaluateAt(x T) T {
	var result any

	switch any(x).(type) {
	case float64:
		gf64 := MultiGaussF64(mg)
		result = gf64.evaluateAt(any(x).(float64))
	case complex128:
		gz64 := MultiGaussZ64(mg)
		result = gz64.evaluateAt(any(x).(complex128))
	default:
		panic("unsupported type")
	}

	return result.(T)
}

func (mg MultiGaussian[T]) ForceAt(x T) T {
	var result any

	switch any(x).(type) {
	case float64:
		gf64 := MultiGaussF64(mg)
		result = gf64.forceAt(any(x).(float64))
	case complex128:
		gz64 := MultiGaussZ64(mg)
		result = gz64.forceAt(any(x).(complex128))
	default:
		panic("unsupported type")
	}

	return result.(T)
}

func (mg MultiGaussian[T]) EvaluateOnGrid(x []T) []T {
	return onGrid(mg.EvaluateAt, x)
}
func (mg MultiGaussian[T]) ForceOnGrid(x []T) []T {
	return onGrid(mg.ForceAt, x)
}
func (mg MultiGaussian[T]) EvaluateOnGridInPlace(fn, x []T) { onGridInPlace(mg.EvaluateAt, fn, x) }
func (mg MultiGaussian[T]) ForceOnGridInPlace(fn, x []T)    { onGridInPlace(mg.ForceAt, fn, x) }

type MultiGaussF64 MultiGaussian[float64]

func (mg MultiGaussF64) evaluateAt(x float64) float64 {
	if mg.NumGauss%2 == 0 {
		val := 0.
		for i := uint8(0); i < mg.NumGauss/2; i++ {
			center := mg.Gap * (float64(i) + 0.5)
			left := math.Exp(-math.Pow(xBySigma(x-center, mg.Sigma), 2) / 2)
			right := math.Exp(-math.Pow(xBySigma(x+center, mg.Sigma), 2) / 2)
			val += left + right
		}
		return mg.Strength * val
	}
	val := math.Exp(-math.Pow(xBySigma(x, mg.Sigma), 2) / 2)
	for i := uint8(0); i < mg.NumGauss/2; i++ {
		center := mg.Gap * (float64(i + 1))
		left := math.Exp(-math.Pow(xBySigma(x-center, mg.Sigma), 2) / 2)
		right := math.Exp(-math.Pow(xBySigma(x+center, mg.Sigma), 2) / 2)
		val += left + right
	}
	return mg.Strength * val
}

func (mg MultiGaussF64) forceAt(x float64) float64 {
	if mg.NumGauss%2 == 0 {
		val := 0.
		for i := uint8(0); i < mg.NumGauss/2; i++ {
			center := mg.Gap * (float64(i) + 0.5)
			expntLeft := xBySigma(x-center, mg.Sigma)
			expntRight := xBySigma(x+center, mg.Sigma)
			left := math.Exp(-math.Pow(expntLeft, 2) / 2)
			right := math.Exp(-math.Pow(expntRight, 2) / 2)
			val += expntLeft*left + expntRight*right
		}
		return mg.Strength * val
	}
	val := xBySigma(x, mg.Sigma) * math.Exp(-math.Pow(xBySigma(x, mg.Sigma), 2)/2)
	for i := uint8(0); i < mg.NumGauss/2; i++ {
		center := mg.Gap * (float64(i + 1))
		expntLeft := xBySigma(x-center, mg.Sigma)
		expntRight := xBySigma(x+center, mg.Sigma)
		left := math.Exp(-math.Pow(expntLeft, 2) / 2)
		right := math.Exp(-math.Pow(expntRight, 2) / 2)
		val += expntLeft*left + expntRight*right
	}
	return mg.Strength * val
}

type MultiGaussZ64 MultiGaussian[complex128]

func (mg MultiGaussZ64) evaluateAt(x complex128) complex128 {
	if mg.NumGauss%2 == 0 {
		val := complex(0., 0.)
		for i := uint8(0); i < mg.NumGauss/2; i++ {
			center := complex(mg.Gap*(float64(i)+0.5), 0.)
			left := cmplx.Exp(-cmplx.Pow(xBySigmaZ64(x-center, mg.Sigma), 2) / 2)
			right := cmplx.Exp(-cmplx.Pow(xBySigmaZ64(x+center, mg.Sigma), 2) / 2)
			val += left + right
		}
		return complex(mg.Strength, 0.) * val
	}

	val := cmplx.Exp(-cmplx.Pow(xBySigmaZ64(x, mg.Sigma), 2) / 2)
	for i := uint8(0); i < mg.NumGauss/2; i++ {
		center := complex(mg.Gap*(float64(i+1)), 0.)
		left := cmplx.Exp(-cmplx.Pow(xBySigmaZ64(x-center, mg.Sigma), 2) / 2)
		right := cmplx.Exp(-cmplx.Pow(xBySigmaZ64(x+center, mg.Sigma), 2) / 2)
		val += left + right
	}
	return complex(mg.Strength, 0.) * val
}

func (mg MultiGaussZ64) forceAt(x complex128) complex128 {
	if mg.NumGauss%2 == 0 {
		val := complex(0., 0.)
		for i := uint8(0); i < mg.NumGauss/2; i++ {
			center := complex(mg.Gap*(float64(i)+0.5), 0.)
			expntLeft := xBySigmaZ64(x-center, mg.Sigma)
			expntRight := xBySigmaZ64(x+center, mg.Sigma)
			left := cmplx.Exp(-cmplx.Pow(expntLeft, 2) / 2)
			right := cmplx.Exp(-cmplx.Pow(expntRight, 2) / 2)
			val += expntLeft*left + expntRight*right
		}
		return complex(mg.Strength, 0.) * val
	}
	val := xBySigmaZ64(x, mg.Sigma) * cmplx.Exp(-cmplx.Pow(xBySigmaZ64(x, mg.Sigma), 2)/2)
	for i := uint8(0); i < mg.NumGauss/2; i++ {
		center := complex(mg.Gap*(float64(i+1)), 0.)
		expntLeft := xBySigmaZ64(x-center, mg.Sigma)
		expntRight := xBySigmaZ64(x+center, mg.Sigma)
		left := cmplx.Exp(-cmplx.Pow(expntLeft, 2) / 2)
		right := cmplx.Exp(-cmplx.Pow(expntRight, 2) / 2)
		val += expntLeft*left + expntRight*right
	}
	return complex(mg.Strength, 0.) * val
}

// SuperGaussian v(x)= v0 exp(-(x/Sigma)^n)
type SuperGaussian[T VarType] struct {
	Cen      float64
	Sigma    float64
	Strength float64
	Order    uint8
}

func (sg SuperGaussian[T]) String() string {
	return fmt.Sprintf("%g Exp[ ((x - %g)/ %g)^%v]", sg.Strength, sg.Cen, sg.Sigma, sg.Order)
}

func (sg SuperGaussian[T]) EvaluateAt(x T) T {
	var result any

	switch any(x).(type) {
	case float64:
		gf64 := SupGaussF64(sg)
		result = gf64.evaluateAt(any(x).(float64))
	case complex128:
		gz64 := SupGaussZ64(sg)
		result = gz64.evaluateAt(any(x).(complex128))
	default:
		panic("unsupported type")
	}

	return result.(T)
}

func (sg SuperGaussian[T]) ForceAt(x T) T {
	var result any

	switch any(x).(type) {
	case float64:
		gf64 := SupGaussF64(sg)
		result = gf64.forceAt(any(x).(float64))
	case complex128:
		gz64 := SupGaussZ64(sg)
		result = gz64.forceAt(any(x).(complex128))
	default:
		panic("unsupported type")
	}

	return result.(T)
}

func (sg SuperGaussian[T]) EvaluateOnGrid(x []T) []T {
	return onGrid(sg.EvaluateAt, x)
}
func (sg SuperGaussian[T]) ForceOnGrid(x []T) []T {
	return onGrid(sg.ForceAt, x)
}
func (sg SuperGaussian[T]) EvaluateOnGridInPlace(fn, x []T) { onGridInPlace(sg.EvaluateAt, fn, x) }
func (sg SuperGaussian[T]) ForceOnGridInPlace(fn, x []T)    { onGridInPlace(sg.ForceAt, fn, x) }

type SupGaussF64 SuperGaussian[float64]

func (sg SupGaussF64) evaluateAt(x float64) float64 {
	val := xBySigma(x-sg.Cen, sg.Sigma)
	return sg.Strength * math.Exp(-math.Pow(val, float64(sg.Order)))
}

func (sg SupGaussF64) forceAt(x float64) float64 {
	forder := float64(sg.Order)
	expnt := xBySigma(x-sg.Cen, sg.Sigma)
	coeffs := -sg.Strength * forder / sg.Sigma * math.Pow(expnt, forder-1)
	return coeffs * math.Exp(-math.Pow(expnt, forder))
}

type SupGaussZ64 SuperGaussian[complex128]

func (sg SupGaussZ64) evaluateAt(x complex128) complex128 {
	val := xBySigmaZ64(x-complex(sg.Cen, 0.), sg.Sigma)
	return complex(sg.Strength, 0.) * cmplx.Exp(-cmplx.Pow(val, 2))
}

func (sg SupGaussZ64) forceAt(x complex128) complex128 {
	forder := float64(sg.Order)
	expnt := xBySigmaZ64(x-complex(sg.Cen, 0.), sg.Sigma)
	coeffs := -complex(sg.Strength*forder/sg.Sigma, 0.) * cmplx.Pow(expnt, complex(forder-1, 0.))
	return coeffs * cmplx.Exp(-cmplx.Pow(expnt, complex(forder, 0.)))
}

// Harmonic v(x)= k/2 x^2
type Harmonic[T VarType] struct {
	Cen        float64
	ForceConst float64
}

func (h Harmonic[T]) String() string { return fmt.Sprintf("1/2 %g (x - %g)^2", h.ForceConst, h.Cen) }

func (h Harmonic[T]) EvaluateAt(x T) T {
	var result any

	switch any(x).(type) {
	case float64:
		gf64 := HarmonicF64(h)
		result = gf64.evaluateAt(any(x).(float64))
	case complex128:
		gz64 := HarmonicZ64(h)
		result = gz64.evaluateAt(any(x).(complex128))
	default:
		panic("unsupported type")
	}

	return result.(T)
}

func (h Harmonic[T]) ForceAt(x T) T {
	var result any

	switch any(x).(type) {
	case float64:
		gf64 := HarmonicF64(h)
		result = gf64.forceAt(any(x).(float64))
	case complex128:
		gz64 := HarmonicZ64(h)
		result = gz64.forceAt(any(x).(complex128))
	default:
		panic("unsupported type")
	}

	return result.(T)
}

func (h Harmonic[T]) EvaluateOnGrid(x []T) []T {
	return onGrid(h.EvaluateAt, x)
}
func (h Harmonic[T]) ForceOnGrid(x []T) []T {
	return onGrid(h.ForceAt, x)
}
func (h Harmonic[T]) EvaluateOnGridInPlace(fn, x []T) { onGridInPlace(h.EvaluateAt, fn, x) }
func (h Harmonic[T]) ForceOnGridInPlace(fn, x []T)    { onGridInPlace(h.ForceAt, fn, x) }

type HarmonicF64 Harmonic[float64]

func (h HarmonicF64) evaluateAt(x float64) float64 { return h.ForceConst / 2 * math.Pow(x-h.Cen, 2) }
func (h HarmonicF64) forceAt(x float64) float64    { return -h.ForceConst * (x - h.Cen) }

type HarmonicZ64 Harmonic[complex128]

func (h HarmonicZ64) evaluateAt(x complex128) complex128 {
	return complex(h.ForceConst/2, 0.) * cmplx.Pow(x-complex(h.Cen, 0.), 2)
}

func (h HarmonicZ64) forceAt(x complex128) complex128 {
	return -complex(h.ForceConst, 0.) * (x - complex(h.Cen, 0.))
}

// Polynomial v(x)= Sum_i ci x^i
type Polynomial[T VarType] struct {
	Coeffs []float64
}

func (p Polynomial[T]) String() string {
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

func (p Polynomial[T]) EvaluateAt(x T) T {
	var result any

	switch any(x).(type) {
	case float64:
		f64 := PolynomialF64(p)
		result = f64.evaluateAt(any(x).(float64))
	case complex128:
		z64 := PolynomialZ64(p)
		result = z64.evaluateAt(any(x).(complex128))
	default:
		panic("unsupported type")
	}

	return result.(T)
}

func (p Polynomial[T]) ForceAt(x T) T {
	var result any

	switch any(x).(type) {
	case float64:
		f64 := PolynomialF64(p)
		result = f64.forceAt(any(x).(float64))
	case complex128:
		z64 := PolynomialZ64(p)
		result = z64.forceAt(any(x).(complex128))
	default:
		panic("unsupported type")
	}

	return result.(T)
}

func (p Polynomial[T]) EvaluateOnGrid(x []T) []T {
	return onGrid(p.EvaluateAt, x)
}
func (p Polynomial[T]) ForceOnGrid(x []T) []T {
	return onGrid(p.ForceAt, x)
}
func (p Polynomial[T]) EvaluateOnGridInPlace(fn, x []T) { onGridInPlace(p.EvaluateAt, fn, x) }
func (p Polynomial[T]) ForceOnGridInPlace(fn, x []T)    { onGridInPlace(p.ForceAt, fn, x) }

type PolynomialF64 Polynomial[float64]

func (p PolynomialF64) evaluateAt(x float64) float64 {
	if len(p.Coeffs) == 0 {
		return 0
	}
	result := p.Coeffs[len(p.Coeffs)-1]
	for i := len(p.Coeffs) - 2; i >= 0; i-- {
		result = result*x + p.Coeffs[i]
	}
	return result
}

func (p PolynomialF64) forceAt(x float64) float64 {
	result := 0.0
	for i := 1; i < len(p.Coeffs); i++ {
		power := float64(i - 1)
		result += p.Coeffs[i] * float64(i) * math.Pow(x, power)
	}
	return -result
}

type PolynomialZ64 Polynomial[complex128]

func (p PolynomialZ64) evaluateAt(x complex128) complex128 {
	if len(p.Coeffs) == 0 {
		return 0
	}
	result := complex(p.Coeffs[len(p.Coeffs)-1], 0)
	for i := len(p.Coeffs) - 2; i >= 0; i-- {
		result = result*x + complex(p.Coeffs[i], 0)
	}
	return result
}

func (p PolynomialZ64) forceAt(x complex128) complex128 {
	result := complex(0, 0)
	for i := 1; i < len(p.Coeffs); i++ {
		power := float64(i - 1)
		result += complex(p.Coeffs[i]*float64(i), 0) * cmplx.Pow(x, complex(power, 0))
	}
	return -result
}
