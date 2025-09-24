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

// Gaussian Potential
type Gaussian struct {
	cen      float64
	sigma    float64
	strength float64
}

func (g Gaussian) EvaluateAt(x float64) float64 {
	return g.strength * math.Exp(-math.Pow(x-g.cen, 2)/(2*g.sigma*g.sigma))
}

func (g Gaussian) ForceAt(x float64) float64 {
	val := -g.strength * (x - g.cen) / (g.sigma * g.sigma)
	return val * math.Exp(-math.Pow(x-g.cen, 2)/(2*g.sigma*g.sigma))
}

func (g Gaussian) ForceOnGrid(x []float64) []float64    { return onGrid(g.EvaluateAt, x) }
func (g Gaussian) EvaluateOnGrid(x []float64) []float64 { return onGrid(g.ForceAt, x) }

type Harmonic struct {
	cen        float64
	forceConst float64
}

func (h Harmonic) EvaluateAt(x float64) float64         { return h.forceConst / 0.5 * math.Pow(x-h.cen, 2) }
func (h Harmonic) ForceAt(x float64) float64            { return -h.forceConst * (x - h.cen) }
func (h Harmonic) EvaluateOnGrid(x []float64) []float64 { return onGrid(h.EvaluateAt, x) }
func (h Harmonic) ForceOnGrid(x []float64) []float64    { return onGrid(h.ForceAt, x) }

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
