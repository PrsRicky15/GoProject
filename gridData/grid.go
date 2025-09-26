package gridData

import (
	"fmt"
	"math"
)

// setGrid Setup some parameters for time- and space- grid
type setGrid struct {
	length  float64
	deltaS  float64
	deltaCS float64
	cMin    float64
	cMax    float64
}

// getGridData get the basic input and setup time- and space- grid
type getGridData interface {
	getMin() float64
	getdS() float64
	getNgrid() uint32
	getdCS() float64
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
		values[i] = g.getMin() + float64(i)*g.getdS()
	}
	return values
}

func generateConjugatePoints(g getGridData) []float64 {
	nby2 := g.getNgrid() / 2
	values := make([]float64, g.getNgrid())
	values[0] = 0.
	values[g.getNgrid()/2] = -float64(nby2) * g.getdCS()
	for i := uint32(1); i < nby2; i++ {
		values[i] = -float64(i) * g.getdCS()
		values[i+nby2] = float64(nby2-i) * g.getdCS()
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

// RadGrid Represents a real-space grid
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
func (g *RadGrid) getdS() float64   { return g.gridData.deltaS }
func (g *RadGrid) getdCS() float64  { return g.gridData.deltaCS }
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
func (g *RadGrid) RValues() []float64 { return generatePoints(g) }
func (g *RadGrid) KValues() []float64 { return generateConjugatePoints(g) }
func (g *RadGrid) RMin() float64      { return g.rMin }
func (g *RadGrid) RMax() float64      { return g.rMax }
func (g *RadGrid) NPoints() uint32    { return g.nPoints }
func (g *RadGrid) DeltaR() float64    { return g.gridData.deltaS }
func (g *RadGrid) Length() float64    { return g.gridData.length }
func (g *RadGrid) DeltaK() float64    { return g.gridData.deltaCS }
func (g *RadGrid) KMin() float64      { return g.gridData.cMin }
func (g *RadGrid) KMax() float64      { return g.gridData.cMax }
func (g *RadGrid) CutoffE() float64   { return g.cutoffE }

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
	nPoints  uint32
	gridData setGrid
}

func NewTimeGrid(tMin, tMax float64, nPoints uint32) (*TimeGrid, error) {
	igridData, err := createGrid(tMin, tMax, nPoints, "Time-grid")
	if err != nil {
		return nil, err
	}
	return &TimeGrid{
		tMin:     tMin,
		tMax:     tMax,
		nPoints:  nPoints,
		gridData: igridData,
	}, nil
}

func TimeGridFromLength(length float64, nPoints uint32) (*TimeGrid, error) {
	if length <= 0 {
		return nil, fmt.Errorf("total duration must be positive")
	}
	return NewTimeGrid(0., length, nPoints)
}

func (t *TimeGrid) getMin() float64  { return t.tMin }
func (t *TimeGrid) getdS() float64   { return t.gridData.deltaS }
func (t *TimeGrid) getdCS() float64  { return t.gridData.deltaCS }
func (t *TimeGrid) getNgrid() uint32 { return t.nPoints }
func (t *TimeGrid) redefine(tMin, tMax float64, nPoints uint32) error {
	igridData, err := createGrid(tMin, tMax, nPoints, "Time-grid")
	if err != nil {
		return err
	}
	t.tMin = tMin
	t.tMax = tMax
	t.nPoints = nPoints
	t.gridData = igridData
	return nil
}

func (t *TimeGrid) ReDefine(nPoints uint32) error { return t.redefine(t.tMin, t.tMax, nPoints) }
func (t *TimeGrid) ReDefineMinMax(tMin, tMax float64, nPoints uint32) error {
	return t.redefine(tMin, tMax, nPoints)
}
func (t *TimeGrid) ReDefineLength(length float64, nPoints uint32) error {
	return t.redefine(0., 0.5*length, nPoints)
}

func (t *TimeGrid) TMin() float64      { return t.tMin }
func (t *TimeGrid) TMax() float64      { return t.tMax }
func (t *TimeGrid) NPoints() uint32    { return t.nPoints }
func (t *TimeGrid) DeltaT() float64    { return t.gridData.deltaS }
func (t *TimeGrid) Length() float64    { return t.gridData.length }
func (t *TimeGrid) DOmega() float64    { return t.gridData.deltaCS }
func (t *TimeGrid) OmegaMin() float64  { return t.gridData.cMin }
func (t *TimeGrid) OmegaMax() float64  { return t.gridData.cMax }
func (t *TimeGrid) TValues() []float64 { return generatePoints(t) }
func (t *TimeGrid) WValues() []float64 { return generateConjugatePoints(t) }
func (t *TimeGrid) DisplayTimeGrid()   { displayGrid(t.TValues) }
func (t *TimeGrid) DisplayOmegaGrid()  { displayGrid(t.WValues) }

func (t *TimeGrid) String() string {
	return fmt.Sprintf("TimeGrid{tMin: %.6g, tMax: %.6g, nPoints: %d, deltaT: %.6g, length: %.6g}",
		t.tMin, t.tMax, t.nPoints, t.gridData.deltaS, t.gridData.length)
}

func (t *TimeGrid) DisplayInfo() {
	fmt.Printf("Real Time - Min: %.6g, Max: %.6g, Dr: %.6g\n", t.TMin(), t.TMax(), t.DeltaT())
	fmt.Printf("w Space    - Min: %.6g, Max: %.6g, Dw: %.6g\n", t.OmegaMin(), t.OmegaMax(), t.DOmega())
	fmt.Printf("Grid       - Length: %.6g, Points: %d\n", t.Length(), t.NPoints())
}
