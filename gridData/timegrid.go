package gridData

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

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

func (t *TimeGrid) PotentialAt(pot PotentialOp[float64], x float64) float64 { return pot.evaluateAt(x) }

func (t *TimeGrid) PrintPotentToFile(Pot PotentialOp[float64], filename string, format string) error {
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
