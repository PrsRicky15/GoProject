package input

import (
	"fmt"
	"math"
)

type RGrid struct {
	rMin    float64
	rMax    float64
	nPoints uint32
	deltaR  float64
	length  float64
	deltaK  float64
	kMin    float64
	kMax    float64
	cutoffE float64
}

func NewRGrid(rMin float64, rMax float64, nPoints uint32) (*RGrid, error) {
	if nPoints == 0 {
		return nil, fmt.Errorf("number of grid points must be positive")
	}
	if rMax <= rMin {
		return nil, fmt.Errorf("rMax (%g) must be greater than rMin (%g)", rMax, rMin)
	}

	deltaR := (rMax - rMin) / float64(nPoints)
	length := rMax - rMin
	deltaK := 2 * math.Pi / length
	kMin := -math.Pi / deltaR
	kMax := math.Pi / deltaR
	cutoffE := math.Pow(kMax, 2) / 2

	return &RGrid{
		rMin:    rMin,
		rMax:    rMax,
		nPoints: nPoints,
		deltaR:  deltaR,
		length:  length,
		deltaK:  deltaK,
		kMin:    kMin,
		kMax:    kMax,
		cutoffE: cutoffE,
	}, nil
}

func NewFromLength(length float64, nPoints uint32) (*RGrid, error) {
	if length <= 0 {
		return nil, fmt.Errorf("length must be positive")
	}

	halfLength := length / 2
	return NewRGrid(-halfLength, halfLength, nPoints)
}

func (g *RGrid) Redimension(nPoints uint32) (*RGrid, error) {
	return NewRGrid(g.rMin, g.rMax, nPoints)
}

func (g *RGrid) RedimensionRange(rMin float64, rMax float64, nPoints uint32) (*RGrid, error) {
	return NewRGrid(rMin, rMax, nPoints)
}

func (g *RGrid) RedimensionLength(length float64, nPoints uint32) (*RGrid, error) {
	return NewFromLength(length, nPoints)
}

func (g *RGrid) RMin() float64    { return g.rMin }
func (g *RGrid) RMax() float64    { return g.rMax }
func (g *RGrid) NPoints() uint32  { return g.nPoints }
func (g *RGrid) DeltaR() float64  { return g.deltaR }
func (g *RGrid) Length() float64  { return g.length }
func (g *RGrid) DeltaK() float64  { return g.deltaK }
func (g *RGrid) KMin() float64    { return g.kMin }
func (g *RGrid) KMax() float64    { return g.kMax }
func (g *RGrid) CutoffE() float64 { return g.cutoffE }

func (g *RGrid) String() string {
	return fmt.Sprintf("RGrid{rMin: %.6g, rMax: %.6g, nPoints: %d, deltaR: %.6g, length: %.6g}",
		g.rMin, g.rMax, g.nPoints, g.deltaR, g.length)
}

func (g *RGrid) RValues() []float64 {
	values := make([]float64, g.nPoints)
	for i := uint32(0); i < g.nPoints; i++ {
		values[i] = g.rMin + float64(i)*g.deltaR
	}
	return values
}

func (g *RGrid) KValues() []float64 {
	values := make([]float64, g.nPoints)
	values[0] = 0.
	values[g.nPoints/2] = -float64(g.nPoints/2) * g.deltaK
	for i := uint32(1); i < g.nPoints/2; i++ {
		values[i] = -float64(i) * g.deltaK
		values[i+g.nPoints/2] = float64(g.nPoints/2-i) * g.deltaK
	}
	return values
}

func (g *RGrid) DisplayInfo() {
	fmt.Printf("Real Space - Min: %.6g, Max: %.6g, Dr: %.6g\n", g.rMin, g.rMax, g.deltaR)
	fmt.Printf("K Space    - Min: %.6g, Max: %.6g, Dk: %.6g\n", g.kMin, g.kMax, g.deltaK)
	fmt.Printf("Grid       - Length: %.6g, Points: %d, Cutoff Energy: %.6g\n",
		g.length, g.nPoints, g.cutoffE)
}

func (g *RGrid) DisplayRgrid() {
	rPoints := g.RValues()
	for ri, val := range rPoints {
		fmt.Printf("r-%d %14.7e\n", ri, val)
	}
}

func (g *RGrid) DisplayKgrid() {
	kPoints := g.KValues()
	for ki, val := range kPoints {
		fmt.Printf("k-%d %14.7e\n", ki, val)
	}
}
