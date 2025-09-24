package gridData

import "math"

// Evaluate General interface for the evaluating the potential on a grid
type Evaluate interface {
	evaluateAt(x float64) float64
	evaluateOnGrid(x []float64) []float64
	forceAt(x float64) float64
}

type Potential struct {
	pot Evaluate
}

// Gaussian Potential
type Gaussian struct {
	cen      float64
	sigma    float64
	strength float64
}

func (g Gaussian) evaluateAt(x float64) float64 {
	return g.strength * math.Exp(-math.Pow(x-g.cen, 2)/(2*g.sigma*g.sigma))
}

func (g Gaussian) evaluateOnGrid(x []float64) []float64 {
	val := make([]float64, len(x))
	for i := 0; i < len(x); i++ {
		val[i] = g.evaluateAt(x[i])
	}
	return val
}

type Harmonic struct {
	cen         float64
	force_const float64
}

func (h Harmonic) evaluateAt(x float64) float64 {
	return -h.force_const / 0.5 * math.Pow(x-h.cen, 2)
}

func (h Harmonic) evaluateOnGrid(x []float64) []float64 {
	val := make([]float64, len(x))
	for i := 0; i < len(x); i++ {
		val[i] = h.evaluateAt(x[i])
	}
	return val
}
