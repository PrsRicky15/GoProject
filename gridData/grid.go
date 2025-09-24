package gridData

import (
	"fmt"
	"math"
)

type RadGrid struct {
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

func NewRGrid(rMin, rMax float64, nPoints uint32) (*RadGrid, error) {
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

	return &RadGrid{
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

func NewFromLength(length float64, nPoints uint32) (*RadGrid, error) {
	if length <= 0 {
		return nil, fmt.Errorf("length must be positive")
	}
	halfLength := length / 2
	return NewRGrid(-halfLength, halfLength, nPoints)
}

func (g *RadGrid) Redimension(nPoints uint32) (*RadGrid, error) {
	return NewRGrid(g.rMin, g.rMax, nPoints)
}

func (g *RadGrid) RedimensionRange(rMin, rMax float64, nPoints uint32) (*RadGrid, error) {
	return NewRGrid(rMin, rMax, nPoints)
}

func (g *RadGrid) RedimensionLength(length float64, nPoints uint32) (*RadGrid, error) {
	return NewFromLength(length, nPoints)
}

func (g *RadGrid) RMin() float64    { return g.rMin }
func (g *RadGrid) RMax() float64    { return g.rMax }
func (g *RadGrid) NPoints() uint32  { return g.nPoints }
func (g *RadGrid) DeltaR() float64  { return g.deltaR }
func (g *RadGrid) Length() float64  { return g.length }
func (g *RadGrid) DeltaK() float64  { return g.deltaK }
func (g *RadGrid) KMin() float64    { return g.kMin }
func (g *RadGrid) KMax() float64    { return g.kMax }
func (g *RadGrid) CutoffE() float64 { return g.cutoffE }

func (g *RadGrid) String() string {
	return fmt.Sprintf("RadGrid{rMin: %.6g, rMax: %.6g, nPoints: %d, deltaR: %.6g, length: %.6g}",
		g.rMin, g.rMax, g.nPoints, g.deltaR, g.length)
}

func (g *RadGrid) RValues() []float64 {
	values := make([]float64, g.nPoints)
	for i := uint32(0); i < g.nPoints; i++ {
		values[i] = g.rMin + float64(i)*g.deltaR
	}
	return values
}

func (g *RadGrid) KValues() []float64 {
	values := make([]float64, g.nPoints)
	values[0] = 0.
	values[g.nPoints/2] = -float64(g.nPoints/2) * g.deltaK
	for i := uint32(1); i < g.nPoints/2; i++ {
		values[i] = -float64(i) * g.deltaK
		values[i+g.nPoints/2] = float64(g.nPoints/2-i) * g.deltaK
	}
	return values
}

func (g *RadGrid) DisplayInfo() {
	fmt.Printf("Real Space - Min: %.6g, Max: %.6g, Dr: %.6g\n", g.RMin(), g.RMax(), g.DeltaR())
	fmt.Printf("K Space    - Min: %.6g, Max: %.6g, Dk: %.6g\n", g.KMin(), g.KMax(), g.DeltaK())
	fmt.Printf("Grid       - Length: %.6g, Points: %d, Cutoff Energy: %.6g\n",
		g.Length(), g.NPoints(), g.CutoffE())
}

func (g *RadGrid) DisplayRgrid() {
	rPoints := g.RValues()
	for ri, val := range rPoints {
		fmt.Printf("r-%d %14.7e\n", ri, val)
	}
}

func (g *RadGrid) DisplayKgrid() {
	kPoints := g.KValues()
	for ki, val := range kPoints {
		fmt.Printf("k-%d %14.7e\n", ki, val)
	}
}

// TimeGrid time-grid definition for the time-dependent differential equation solver
type TimeGrid struct {
	tMin     float64
	tMax     float64
	nPoints  uint32
	deltaT   float64
	length   float64
	dOmega   float64
	omegaMin float64
	omegaMax float64
}

func NewTimeGrid(tMin, tMax float64, nPoints uint32) (*TimeGrid, error) {
	if nPoints == 0 {
		return nil, fmt.Errorf("number of grid points must be positive")
	}
	if tMax <= tMin {
		return nil, fmt.Errorf("rMax (%g) must be greater than rMin (%g)", tMax, tMin)
	}

	deltaT := (tMax - tMin) / float64(nPoints)
	length := tMax - tMin
	deltaW := 2 * math.Pi / length
	wMin := -math.Pi / deltaT
	wMax := math.Pi / deltaT

	return &TimeGrid{
		tMin:     tMin,
		tMax:     tMax,
		nPoints:  nPoints,
		deltaT:   deltaT,
		length:   length,
		dOmega:   deltaW,
		omegaMin: wMin,
		omegaMax: wMax,
	}, nil
}

func TimeGridFromLength(length float64, nPoints uint32) (*TimeGrid, error) {
	if length <= 0 {
		return nil, fmt.Errorf("length must be positive")
	}
	return NewTimeGrid(0., length, nPoints)
}

func (g *TimeGrid) Redimension(nPoints uint32) (*TimeGrid, error) {
	return NewTimeGrid(g.tMin, g.tMax, nPoints)
}

func (g *TimeGrid) RedimensionRange(tMin, tMax float64, nPoints uint32) (*TimeGrid, error) {
	return NewTimeGrid(tMin, tMax, nPoints)
}

func (g *TimeGrid) RedimensionLength(length float64, nPoints uint32) (*TimeGrid, error) {
	return TimeGridFromLength(length, nPoints)
}

func (g *TimeGrid) TMin() float64     { return g.tMin }
func (g *TimeGrid) TMax() float64     { return g.tMax }
func (g *TimeGrid) NPoints() uint32   { return g.nPoints }
func (g *TimeGrid) DeltaT() float64   { return g.deltaT }
func (g *TimeGrid) Length() float64   { return g.length }
func (g *TimeGrid) DOmega() float64   { return g.dOmega }
func (g *TimeGrid) OmegaMin() float64 { return g.omegaMin }
func (g *TimeGrid) OmegaMax() float64 { return g.omegaMax }

func (g *TimeGrid) String() string {
	return fmt.Sprintf("TimeGrid{tMin: %.6g, tMax: %.6g, nPoints: %d, deltaT: %.6g, length: %.6g}",
		g.tMin, g.tMax, g.nPoints, g.deltaT, g.length)
}

func (g *TimeGrid) TValues() []float64 {
	values := make([]float64, g.nPoints)
	for i := uint32(0); i < g.nPoints; i++ {
		values[i] = g.tMin + float64(i)*g.deltaT
	}
	return values
}

func (g *TimeGrid) WValues() []float64 {
	values := make([]float64, g.nPoints)
	values[0] = 0.
	values[g.nPoints/2] = -float64(g.nPoints/2) * g.deltaT
	for i := uint32(1); i < g.nPoints/2; i++ {
		values[i] = -float64(i) * g.deltaT
		values[i+g.nPoints/2] = float64(g.nPoints/2-i) * g.dOmega
	}
	return values
}

func (g *TimeGrid) DisplayInfo() {
	fmt.Printf("Real Time - Min: %.6g, Max: %.6g, Dr: %.6g\n", g.TMin(), g.TMax(), g.DeltaT())
	fmt.Printf("w Space    - Min: %.6g, Max: %.6g, Dw: %.6g\n", g.OmegaMin(), g.OmegaMax(), g.DOmega())
	fmt.Printf("Grid       - Length: %.6g, Points: %d\n", g.Length(), g.NPoints())
}

func (g *TimeGrid) DisplayTimeGrid() {
	rPoints := g.TValues()
	for ri, val := range rPoints {
		fmt.Printf("t-%d %14.7e\n", ri, val)
	}
}

func (g *TimeGrid) DisplayOmegaGrid() {
	kPoints := g.WValues()
	for ki, val := range kPoints {
		fmt.Printf("w-%d %14.7e\n", ki, val)
	}
}
