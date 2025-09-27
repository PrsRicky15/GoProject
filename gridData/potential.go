package gridData

import "math"

// Evaluate General interface for the evaluating the potential on a grid
type Evaluate interface {
	ForceAt(x float64) float64
	EvaluateAt(x float64) float64
	EvaluateOnGrid(x []float64) []float64
	ForceOnGrid(x []float64) []float64
}

func onGrid(f func(float64) float64, x []float64) []float64 {
	results := make([]float64, len(x))
	for i, val := range x {
		results[i] = f(val)
	}
	return results
}

// SoftCore Potential
type SoftCore struct {
	charge    float64
	centre    float64
	softParam float64
}

func (sc SoftCore) EvaluateAt(x float64) float64 {
	return sc.charge / math.Sqrt(math.Pow(x-sc.centre, 2)+math.Pow(sc.softParam, 2))
}

func (sc SoftCore) ForceAt(x float64) float64 {
	coef := sc.charge * (x - sc.centre)
	val := (x-sc.centre)*(x-sc.centre) + sc.softParam*sc.softParam
	return coef * math.Pow(val, 3./2.)
}

func (sc SoftCore) ForceOnGrid(x []float64) []float64    { return onGrid(sc.EvaluateAt, x) }
func (sc SoftCore) EvaluateOnGrid(x []float64) []float64 { return onGrid(sc.ForceAt, x) }

// Gaussian Potential
type Gaussian struct {
	cen      float64
	sigma    float64
	strength float64
}

func xbysigma(x float64, sigma float64) float64 {
	return x / sigma
}

func (g Gaussian) EvaluateAt(x float64) float64 {
	val := xbysigma(x-g.cen, g.sigma)
	return g.strength * math.Exp(-math.Pow(val, 2)/2)
}

func (g Gaussian) ForceAt(x float64) float64 {
	expnt := xbysigma(x-g.cen, g.sigma)
	val := -g.strength * expnt / g.sigma
	return val * math.Exp(-math.Pow(expnt, 2)/2)
}

func (g Gaussian) ForceOnGrid(x []float64) []float64    { return onGrid(g.EvaluateAt, x) }
func (g Gaussian) EvaluateOnGrid(x []float64) []float64 { return onGrid(g.ForceAt, x) }

// MultiGaussian Potential
type MultiGaussian struct {
	sigma    float64
	strength float64
	numGauss uint8
	gap      float64
}

func (mg MultiGaussian) EvaluateAt(x float64) float64 {
	if mg.numGauss/2 == 0 {
		val := 0.
		for i := uint8(0); i < mg.numGauss/2; i++ {
			center := mg.gap * (float64(i) + 0.5)
			left := math.Exp(-math.Pow(xbysigma(x-center, mg.sigma), 2) / 2)
			right := math.Exp(-math.Pow(xbysigma(x+center, mg.sigma), 2) / 2)
			val += left + right
		}
		return mg.strength * val
	}
	val := math.Exp(-math.Pow(xbysigma(x, mg.sigma), 2) / 2)
	for i := uint8(0); i < mg.numGauss/2; i++ {
		center := mg.gap * (float64(i + 1))
		left := math.Exp(-math.Pow(xbysigma(x-center, mg.sigma), 2) / 2)
		right := math.Exp(-math.Pow(xbysigma(x+center, mg.sigma), 2) / 2)
		val += left + right
	}
	return mg.strength * val
}

func (mg MultiGaussian) ForceAt(x float64) float64 {
	if mg.numGauss/2 == 0 {
		val := 0.
		for i := uint8(0); i < mg.numGauss/2; i++ {
			center := mg.gap * (float64(i) + 0.5)
			expntLeft := xbysigma(x-center, mg.sigma)
			expntRight := xbysigma(x+center, mg.sigma)
			left := math.Exp(-math.Pow(expntLeft, 2) / 2)
			right := math.Exp(-math.Pow(expntRight, 2) / 2)
			val += expntLeft*left + expntRight*right
		}
		return mg.strength * val
	}
	val := xbysigma(x, mg.sigma) * math.Exp(-math.Pow(xbysigma(x, mg.sigma), 2)/2)
	for i := uint8(0); i < mg.numGauss/2; i++ {
		center := mg.gap * (float64(i + 1))
		expntLeft := xbysigma(x-center, mg.sigma)
		expntRight := xbysigma(x+center, mg.sigma)
		left := math.Exp(-math.Pow(expntLeft, 2) / 2)
		right := math.Exp(-math.Pow(expntRight, 2) / 2)
		val += expntLeft*left + expntRight*right
	}
	return mg.strength * val
}

func (mg MultiGaussian) ForceOnGrid(x []float64) []float64    { return onGrid(mg.EvaluateAt, x) }
func (mg MultiGaussian) EvaluateOnGrid(x []float64) []float64 { return onGrid(mg.ForceAt, x) }

// SuperGaussian v(x)= v0 exp(-(x/sigma)^n)
type SuperGaussian struct {
	cen      float64
	sigma    float64
	strength float64
	order    uint8
}

func (sg SuperGaussian) EvaluateAt(x float64) float64 {
	val := xbysigma(x-sg.cen, sg.sigma)
	return sg.strength * math.Exp(-math.Pow(val, float64(sg.order)))
}

func (sg SuperGaussian) ForceAt(x float64) float64 {
	forder := float64(sg.order)
	expnt := xbysigma(x-sg.cen, sg.sigma)
	coef := -sg.strength * forder / sg.sigma * math.Pow(expnt, forder-1)
	return coef * math.Exp(-math.Pow(expnt, forder))
}

func (sg SuperGaussian) ForceOnGrid(x []float64) []float64    { return onGrid(sg.EvaluateAt, x) }
func (sg SuperGaussian) EvaluateOnGrid(x []float64) []float64 { return onGrid(sg.ForceAt, x) }

// Harmonic v(x)= k/2 x^2
type Harmonic struct {
	cen        float64
	forceConst float64
}

func (h Harmonic) EvaluateAt(x float64) float64         { return h.forceConst / 0.5 * math.Pow(x-h.cen, 2) }
func (h Harmonic) ForceAt(x float64) float64            { return -h.forceConst * (x - h.cen) }
func (h Harmonic) EvaluateOnGrid(x []float64) []float64 { return onGrid(h.EvaluateAt, x) }
func (h Harmonic) ForceOnGrid(x []float64) []float64    { return onGrid(h.ForceAt, x) }

// Polynomial v(x)= Sum_i ci x^i
type Polynomial struct {
	coeffs []float64
}

func (p Polynomial) EvaluateAt(x float64) float64 {
	results := 0.0
	for i := 0; i < len(p.coeffs); i++ {
		results += p.coeffs[i] * math.Pow(x, float64(i))
	}
	return results
}

func (p Polynomial) ForceAt(x float64) float64 {
	result := 0.0
	for i := 1; i < len(p.coeffs); i++ {
		power := float64(i - 1)
		result += p.coeffs[i] * float64(i) * math.Pow(x, power)
	}
	return -result
}

func (p Polynomial) EvaluateOnGrid(x []float64) []float64 { return onGrid(p.EvaluateAt, x) }
func (p Polynomial) ForceOnGrid(x []float64) []float64    { return onGrid(p.ForceAt, x) }

// Morse v(r)= De (1 - Exp(-(r-re))^2
type Morse struct {
	de    float64
	alpha float64
	cen   float64
}

func (m Morse) EvaluateAt(x float64) float64 {
	return m.de * math.Pow(1.-math.Exp(-m.alpha*(x-m.cen)), 2)
}

func (m Morse) ForceAt(x float64) float64 {
	expterm := math.Exp(-m.alpha * (x - m.cen))
	val := -2 * m.alpha * m.de
	return val * expterm * (1 - expterm)
}

func (m Morse) ForceOnGrid(x []float64) []float64    { return onGrid(m.EvaluateAt, x) }
func (m Morse) EvaluateOnGrid(x []float64) []float64 { return onGrid(m.ForceAt, x) }
