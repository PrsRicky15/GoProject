package gridData

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
)

// RadGrid Represents a real-space grid
type RadGrid struct {
	rMin     float64
	rMax     float64
	nPoints  uint32
	gridData setGrid
	cutoffE  float64
}

// TimeGrid represents a time-grid definition for time-dependent differential equation solver
type TimeGrid struct {
	tMin       float64
	tMax       float64
	gridData   setGrid
	macroDt    float64
	microSteps uint32
	macroSteps uint32
	nPoints    uint32
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

func NewRGridFromFile(dirPath string) (*RadGrid, error) {
	info, err := os.Stat(dirPath)
	if err != nil {
		return nil, fmt.Errorf("directory check failed: %w", err)
	}
	if !info.IsDir() {
		return nil, fmt.Errorf("%s is not a directory", dirPath)
	}

	filePath := dirPath + "/rgrid.inp"

	file, err := os.Open(filePath)
	if err != nil {
		return nil, fmt.Errorf("could not open %s: %w", filePath, err)
	}

	var rmin, rmax float64
	var nPoints int

	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}

		parts := strings.Split(line, ":")
		if len(parts) != 2 {
			continue
		}

		key := strings.TrimSpace(parts[0])
		value := strings.TrimSpace(parts[1])

		switch key {
		case "rMin":
			rmin, _ = strconv.ParseFloat(value, 64)
		case "rMax":
			rmax, _ = strconv.ParseFloat(value, 64)
		case "nPoints":
			nPoints, _ = strconv.Atoi(value)
		}
	}

	return NewRGrid(rmin, rmax, uint32(nPoints))
}

func NewFromLength(length float64, nPoints uint32) (*RadGrid, error) {
	if length <= 0 {
		return nil, fmt.Errorf("length must be positive")
	}
	halfLength := length / 2
	return NewRGrid(-halfLength, halfLength, nPoints)
}

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

func (g *RadGrid) getMin() float64  { return g.rMin }
func (g *RadGrid) getdS() float64   { return g.gridData.deltaS }
func (g *RadGrid) getdCS() float64  { return g.gridData.deltaCS }
func (g *RadGrid) getNgrid() uint32 { return g.nPoints }

func (g *RadGrid) RValues() []float64                             { return generatePoints(g) }
func (g *RadGrid) KValues() []float64                             { return generateConjugatePoints(g) }
func (g *RadGrid) RMin() float64                                  { return g.rMin }
func (g *RadGrid) RMax() float64                                  { return g.rMax }
func (g *RadGrid) NPoints() uint32                                { return g.nPoints }
func (g *RadGrid) DeltaR() float64                                { return g.gridData.deltaS }
func (g *RadGrid) Length() float64                                { return g.gridData.length }
func (g *RadGrid) DeltaK() float64                                { return g.gridData.deltaCS }
func (g *RadGrid) KMin() float64                                  { return g.gridData.cMin }
func (g *RadGrid) KMax() float64                                  { return g.gridData.cMax }
func (g *RadGrid) CutoffE() float64                               { return g.cutoffE }
func (g *RadGrid) DisplayRGrid()                                  { displayGrid(g.RValues) }
func (g *RadGrid) DisplayKGrid()                                  { displayGrid(g.KValues) }
func (g *RadGrid) PotentialAt(pot PotentialOp, x float64) float64 { return pot.evaluateAt(x) }
func (g *RadGrid) ForceAt(pot PotentialOp, x float64) float64     { return pot.forceAt(x) }
func (g *RadGrid) PotentialOnGrid(pot PotentialOp) []float64      { return pot.evaluateOnGrid(g.RValues()) }
func (g *RadGrid) ForceOnGrid(pot PotentialOp) []float64          { return pot.forceOnGrid(g.RValues()) }
func (g *RadGrid) DisplayPotential(Pot PotentialOp, format string) {
	displayFunc(g, Pot, format, g.PotentialAt)
}
func (g *RadGrid) DisplayForce(Pot PotentialOp, format string) {
	displayFunc(g, Pot, format, g.ForceAt)
}

func (g *RadGrid) PrintPotentToFile(Pot PotentialOp, filename string, format string) error {
	err := functionToFile(g, Pot, filename, format, g.PotentialAt)
	return err
}

func (g *RadGrid) PrintForceToFile(Pot PotentialOp, filename string, format string) error {
	err := functionToFile(g, Pot, filename, format, g.ForceAt)
	return err
}

func (g *RadGrid) PrintVectorToFile(vec []float64, filename string, format string) error {
	err := vectorToFile(g, vec, filename, format)
	return err
}

func (g *RadGrid) String() string {
	return fmt.Sprintf("RadGrid{rMin: %.6g, rMax: %.6g, nPoints: %d, deltaK: %.6g, CutOffe: %.6g}",
		g.rMin, g.rMax, g.nPoints, g.gridData.deltaCS, g.cutoffE)
}

func (g *RadGrid) DisplayInfo() {
	fmt.Println("#-------------------------------------------------------------")
	fmt.Println("# Real-Space-grid parameters:")
	fmt.Printf("Real Space - Min:    %8.4g | Max: %8.4g | Dr: %8.4g\n", g.RMin(), g.RMax(), g.DeltaR())
	fmt.Printf("K Space    - Min:    %8.4g | Max: %8.4g | Dk: %8.4g\n", g.KMin(), g.KMax(), g.DeltaK())
	fmt.Printf("Grid       - Length: %8.4g | Points: %5d | Cutoff Energy: %8.4g\n",
		g.Length(), g.NPoints(), g.CutoffE())
	fmt.Println("#-------------------------------------------------------------")
}

func NewTimeGrid(MacroDT float64, MacroSteps, MicroSteps uint32) (*TimeGrid, error) {
	nPoints := MacroSteps * MicroSteps
	tMin := 0.
	tMax := MacroDT * float64(MacroSteps)
	igridData, err := createGrid(tMin, tMax, nPoints, "Time-grid")
	if err != nil {
		return nil, err
	}
	return &TimeGrid{
		tMin:       tMin,
		tMax:       tMax,
		gridData:   igridData,
		macroDt:    MacroDT,
		microSteps: MicroSteps,
		macroSteps: MacroSteps,
		nPoints:    nPoints,
	}, nil
}

func NewTGridFromFile(dirPath string) (*TimeGrid, error) {
	info, err := os.Stat(dirPath)
	if err != nil {
		return nil, fmt.Errorf("directory check failed: %w", err)
	}
	if !info.IsDir() {
		return nil, fmt.Errorf("%s is not a directory", dirPath)
	}

	filePath := dirPath + "/tgrid.inp"

	file, err := os.Open(filePath)
	if err != nil {
		return nil, fmt.Errorf("could not open %s: %w", filePath, err)
	}

	var MacroDt float64
	var MacroSteps int
	var MicroSteps int

	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}

		parts := strings.Split(line, ":")
		if len(parts) != 2 {
			continue
		}

		key := strings.TrimSpace(parts[0])
		value := strings.TrimSpace(parts[1])

		switch key {
		case "MacroDT":
			MacroDt, _ = strconv.ParseFloat(value, 64)
		case "MacroSteps":

			MacroSteps, _ = strconv.Atoi(value)
		case "MicroSteps":
			MicroSteps, _ = strconv.Atoi(value)
		}
	}

	return NewTimeGrid(MacroDt, uint32(MacroSteps), uint32(MicroSteps))
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

func (t *TimeGrid) PotentialAt(pot PotentialOp, x float64) float64 { return pot.evaluateAt(x) }

func (t *TimeGrid) PrintPotentToFile(Pot PotentialOp, filename string, format string) error {
	err := functionToFile(t, Pot, filename, format, t.PotentialAt)
	return err
}

func (t *TimeGrid) String() string {
	return fmt.Sprintf("TGrid{Duration: %.6g, macroDT: %.6g, macroSteps: %d, microSteps: %.6v}",
		t.tMax, t.macroDt, t.macroSteps, t.microSteps)
}

func (t *TimeGrid) DisplayInfo() {
	fmt.Println("#-------------------------------------------------------------")
	fmt.Println("# Time-grid parameters:")
	fmt.Printf("Real Time  - Min:    %8.4g | Max: %8.4g | Dr: %8.4g\n", t.TMin(), t.TMax(), t.DeltaT())
	fmt.Printf("w Space    - Min:    %8.4g | Max: %8.4g | Dw: %8.4g\n", t.OmegaMin(), t.OmegaMax(), t.DOmega())
	fmt.Printf("Grid       - Length: %8.4g | Points: %d\n", t.Length(), t.NPoints())
	fmt.Println("#-------------------------------------------------------------")
}
