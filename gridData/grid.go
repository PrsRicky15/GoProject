package gridData

import (
	"fmt"
	"math"
)

type getGridData interface {
	getMin() float64
	getdR() float64
	getNgrid() uint32
}

func displayGrid(space func() []float64) {
	points := space()
	for i, val := range points {
		fmt.Printf("r-%d %14.7e\n", i, val)
	}
}

func generatePoints(g getGridData) []float64 {
	values := make([]float64, g.getNgrid())
	for i := uint32(0); i < g.getNgrid(); i++ {
		values[i] = g.getMin() + float64(i)*g.getdR()
	}
	return values
}

func createGrid(min, max float64, nPoints uint32, paramName string) (setGrid, error) {
	if nPoints == 0 {
		return setGrid{0., 0., 0., 0., 0.},
			fmt.Errorf("number of grid points must be positive")
	}
	if max <= min {
		return setGrid{0., 0., 0., 0., 0.},
			fmt.Errorf("%sMax (%g) must be greater than %sMin (%g)", paramName, max, paramName, min)
	}
	dSpace := (max - min) / float64(nPoints)
	length := max - min
	dConjugate := 2 * math.Pi / length
	cMin := -math.Pi / dSpace
	cMax := math.Pi / dSpace
	return setGrid{length, dSpace, dConjugate, cMin, cMax}, nil
}

type setGrid struct {
	length  float64
	deltaS  float64
	deltaCS float64
	cMin    float64
	cMax    float64
}

type RadGrid struct {
	rMin     float64
	rMax     float64
	nPoints  uint32
	gridData setGrid
	cutoffE  float64
}

func NewRGrid(rMin, rMax float64, nPoints uint32) (*RadGrid, error) {
	igridData, err := createGrid(rMin, rMax, nPoints, "Rgrid")
	if err != nil {
		return nil, err
	}

	cutoffE := math.Pow(igridData.cMax, 2) / 2
	return &RadGrid{
		rMin:     rMin,
		rMax:     rMax,
		nPoints:  nPoints,
		gridData: igridData,
		cutoffE:  cutoffE,
	}, nil
}

func NewFromLength(length float64, nPoints uint32) (*RadGrid, error) {
	if length <= 0 {
		return nil, fmt.Errorf("length must be positive")
	}
	halfLength := length / 2
	return NewRGrid(-halfLength, halfLength, nPoints)
}

func (g *RadGrid) getMin() float64  { return g.rMin }
func (g *RadGrid) getdR() float64   { return g.gridData.deltaS }
func (g *RadGrid) getNgrid() uint32 { return g.nPoints }

func (g *RadGrid) redefine(rMin, rMax float64, nPoints uint32) error {
	igridData, err := createGrid(rMin, rMax, nPoints, "Rgrid")
	if err != nil {
		return err
	}
	g.rMin = rMin
	g.rMax = rMax
	g.nPoints = nPoints
	g.gridData = igridData
	g.cutoffE = math.Pow(igridData.cMax, 2) / 2
	return nil
}

func (g *RadGrid) ReDefine(nPoints uint32) error { return g.redefine(g.rMin, g.rMax, nPoints) }
func (g *RadGrid) ReDefineMinMax(rMin, rMax float64, nPoints uint32) error {
	return g.redefine(rMin, rMax, nPoints)
}
func (g *RadGrid) ReDefineLength(length float64, nPoints uint32) error {
	return g.redefine(-0.5*length, 0.5*length, nPoints)
}

func (g *RadGrid) RValues() []float64 {
	return generatePoints(g)
}

func (g *RadGrid) KValues() []float64 {
	values := make([]float64, g.nPoints)
	values[0] = 0.
	values[g.nPoints/2] = -float64(g.nPoints/2) * g.gridData.deltaCS
	for i := uint32(1); i < g.nPoints/2; i++ {
		values[i] = -float64(i) * g.gridData.deltaCS
		values[i+g.nPoints/2] = float64(g.nPoints/2-i) * g.gridData.deltaCS
	}
	return values
}

func (g *RadGrid) RMin() float64    { return g.rMin }
func (g *RadGrid) RMax() float64    { return g.rMax }
func (g *RadGrid) NPoints() uint32  { return g.nPoints }
func (g *RadGrid) DeltaR() float64  { return g.gridData.deltaS }
func (g *RadGrid) Length() float64  { return g.gridData.length }
func (g *RadGrid) DeltaK() float64  { return g.gridData.deltaCS }
func (g *RadGrid) KMin() float64    { return g.gridData.cMin }
func (g *RadGrid) KMax() float64    { return g.gridData.cMax }
func (g *RadGrid) CutoffE() float64 { return g.cutoffE }

func (g *RadGrid) String() string {
	return fmt.Sprintf("RadGrid{rMin: %.6g, rMax: %.6g, nPoints: %d, deltaR: %.6g, length: %.6g}",
		g.rMin, g.rMax, g.nPoints, g.gridData.deltaS, g.gridData.length)
}

func (g *RadGrid) DisplayInfo() {
	fmt.Printf("Real Space - Min: %.6g, Max: %.6g, Dr: %.6g\n", g.RMin(), g.RMax(), g.DeltaR())
	fmt.Printf("K Space    - Min: %.6g, Max: %.6g, Dk: %.6g\n", g.KMin(), g.KMax(), g.DeltaK())
	fmt.Printf("Grid       - Length: %.6g, Points: %d, Cutoff Energy: %.6g\n",
		g.Length(), g.NPoints(), g.CutoffE())
}

func (g *RadGrid) DisplayRgrid()                               { displayGrid(g.RValues) }
func (g *RadGrid) DisplayKgrid()                               { displayGrid(g.KValues) }
func (g *RadGrid) PotentialAt(pot Evaluate, x float64) float64 { return pot.EvaluateAt(x) }
func (g *RadGrid) ForceAt(pot Evaluate, x float64) float64     { return pot.ForceAt(x) }
func (g *RadGrid) PotentialOnGrid(pot Evaluate) []float64      { return pot.EvaluateOnGrid(g.RValues()) }
func (g *RadGrid) ForceOnGrid(pot Evaluate) []float64          { return pot.ForceOnGrid(g.RValues()) }

// TimeGrid represents a time-grid definition for time-dependent differential equation solver
type TimeGrid struct {
	tMin     float64
	tMax     float64
	length   float64
	nPoints  uint32
	deltaT   float64
	dOmega   float64
	omegaMin float64
	omegaMax float64
}

func DisplayTimeGrid(space func() []float64) {
	points := space()
	for ri, val := range points {
		fmt.Printf("t-%d %14.7e\n", ri, val)
	}
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
