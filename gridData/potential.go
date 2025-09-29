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
	Charge    float64
	Centre    float64
	SoftParam float64
}

func (sc SoftCore) EvaluateAt(x float64) float64 {
	return sc.Charge / math.Sqrt(math.Pow(x-sc.Centre, 2)+math.Pow(sc.SoftParam, 2))
}

func (sc SoftCore) ForceAt(x float64) float64 {
	coef := sc.Charge * (x - sc.Centre)
	val := (x-sc.Centre)*(x-sc.Centre) + sc.SoftParam*sc.SoftParam
	return coef * math.Pow(val, -3./2.)
}

func (sc SoftCore) ForceOnGrid(x []float64) []float64    { return onGrid(sc.EvaluateAt, x) }
func (sc SoftCore) EvaluateOnGrid(x []float64) []float64 { return onGrid(sc.ForceAt, x) }

// Gaussian Potential
type Gaussian struct {
	Cen      float64
	Sigma    float64
	Strength float64
}

func xbysigma(x float64, sigma float64) float64 {
	return x / sigma
}

func (g Gaussian) EvaluateAt(x float64) float64 {
	val := xbysigma(x-g.Cen, g.Sigma)
	return g.Strength * math.Exp(-math.Pow(val, 2)/2)
}

func (g Gaussian) ForceAt(x float64) float64 {
	expnt := xbysigma(x-g.Cen, g.Sigma)
	val := -g.Strength * expnt / g.Sigma
	return val * math.Exp(-math.Pow(expnt, 2)/2)
}

func (g Gaussian) ForceOnGrid(x []float64) []float64    { return onGrid(g.EvaluateAt, x) }
func (g Gaussian) EvaluateOnGrid(x []float64) []float64 { return onGrid(g.ForceAt, x) }

// MultiGaussian Potential
type MultiGaussian struct {
	Sigma    float64
	Strength float64
	NumGauss uint8
	Gap      float64
}

func (mg MultiGaussian) EvaluateAt(x float64) float64 {
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

func (mg MultiGaussian) ForceAt(x float64) float64 {
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

func (mg MultiGaussian) ForceOnGrid(x []float64) []float64    { return onGrid(mg.EvaluateAt, x) }
func (mg MultiGaussian) EvaluateOnGrid(x []float64) []float64 { return onGrid(mg.ForceAt, x) }

// SuperGaussian v(x)= v0 exp(-(x/Sigma)^n)
type SuperGaussian struct {
	Cen      float64
	Sigma    float64
	Strength float64
	Order    uint8
}

func (sg SuperGaussian) EvaluateAt(x float64) float64 {
	val := xbysigma(x-sg.Cen, sg.Sigma)
	return sg.Strength * math.Exp(-math.Pow(val, float64(sg.Order)))
}

func (sg SuperGaussian) ForceAt(x float64) float64 {
	forder := float64(sg.Order)
	expnt := xbysigma(x-sg.Cen, sg.Sigma)
	coeffs := -sg.Strength * forder / sg.Sigma * math.Pow(expnt, forder-1)
	return coeffs * math.Exp(-math.Pow(expnt, forder))
}

func (sg SuperGaussian) ForceOnGrid(x []float64) []float64    { return onGrid(sg.EvaluateAt, x) }
func (sg SuperGaussian) EvaluateOnGrid(x []float64) []float64 { return onGrid(sg.ForceAt, x) }

// Harmonic v(x)= k/2 x^2
type Harmonic struct {
	Cen        float64
	ForceConst float64
}

func (h Harmonic) EvaluateAt(x float64) float64         { return h.ForceConst / 0.5 * math.Pow(x-h.Cen, 2) }
func (h Harmonic) ForceAt(x float64) float64            { return -h.ForceConst * (x - h.Cen) }
func (h Harmonic) EvaluateOnGrid(x []float64) []float64 { return onGrid(h.EvaluateAt, x) }
func (h Harmonic) ForceOnGrid(x []float64) []float64    { return onGrid(h.ForceAt, x) }

// Polynomial v(x)= Sum_i ci x^i
type Polynomial struct {
	Coeffs []float64
}

func (p Polynomial) EvaluateAt(x float64) float64 {
	results := 0.0
	for i := 0; i < len(p.Coeffs); i++ {
		results += p.Coeffs[i] * math.Pow(x, float64(i))
	}
	return results
}

func (p Polynomial) ForceAt(x float64) float64 {
	result := 0.0
	for i := 1; i < len(p.Coeffs); i++ {
		power := float64(i - 1)
		result += p.Coeffs[i] * float64(i) * math.Pow(x, power)
	}
	return -result
}

func (p Polynomial) EvaluateOnGrid(x []float64) []float64 { return onGrid(p.EvaluateAt, x) }
func (p Polynomial) ForceOnGrid(x []float64) []float64    { return onGrid(p.ForceAt, x) }

// Morse v(r)= De (1 - Exp(-(r-re))^2
type Morse struct {
	De    float64
	Alpha float64
	Cen   float64
}

func (m Morse) EvaluateAt(x float64) float64 {
	return m.De * math.Pow(1.-math.Exp(-m.Alpha*(x-m.Cen)), 2)
}

func (m Morse) ForceAt(x float64) float64 {
	expterm := math.Exp(-m.Alpha * (x - m.Cen))
	val := -2 * m.Alpha * m.De
	return val * expterm * (1 - expterm)
}

func (m Morse) ForceOnGrid(x []float64) []float64    { return onGrid(m.EvaluateAt, x) }
func (m Morse) EvaluateOnGrid(x []float64) []float64 { return onGrid(m.ForceAt, x) }
