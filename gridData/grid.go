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

func (g *RadGrid) RMin() float64    { return g.rMin }
func (g *RadGrid) RMax() float64    { return g.rMax }
func (g *RadGrid) NPoints() uint32  { return g.nPoints }
func (g *RadGrid) DeltaR() float64  { return g.gridData.deltaS }
func (g *RadGrid) Length() float64  { return g.gridData.length }
func (g *RadGrid) DeltaK() float64  { return g.gridData.deltaCS }
func (g *RadGrid) KMin() float64    { return g.gridData.cMin }
func (g *RadGrid) KMax() float64    { return g.gridData.cMax }
func (g *RadGrid) CutoffE() float64 { return g.cutoffE }

func (g *RadGrid) RValues() []float64 { return generatePoints(g) }
func (g *RadGrid) KValues() []float64 { return generateConjugatePoints(g) }
func (g *RadGrid) DisplayRGrid()      { displayGrid(g.RValues) }
func (g *RadGrid) DisplayKGrid()      { displayGrid(g.KValues) }

func (g *RadGrid) FunctionAtx(pot Rfunc, x float64) float64 {
	return pot.EvaluateAt(x)
}

func (g *RadGrid) FunctionOnGridInPlace(pot Rfunc, f []float64) {
	for i := range f {
		x := g.rMin + g.gridData.deltaS*float64(i)
		f[i] = pot.EvaluateAt(x)
	}
}

func (g *RadGrid) PotentialAtR(pot PotentialOp[float64], x float64) float64 {
	return pot.EvaluateAt(x)
}
func (g *RadGrid) ForceAtR(pot PotentialOp[float64], x float64) float64 { return pot.ForceAt(x) }
func (g *RadGrid) PotentialAtZ(pot PotentialOp[complex128], x complex128) complex128 {
	return pot.EvaluateAt(x)
}
func (g *RadGrid) ForceAtZ(pot PotentialOp[complex128], x complex128) complex128 {
	return pot.ForceAt(x)
}

func (g *RadGrid) PotentialOnGrid(pot PotentialOp[float64]) []float64 {
	return pot.EvaluateOnGrid(g.RValues())
}

func (g *RadGrid) ForceOnGrid(pot PotentialOp[float64]) []float64 {
	return pot.ForceOnGrid(g.RValues())
}

func (g *RadGrid) DisplayPotentialRe(Pot PotentialOp[float64], format string) {
	displayFunc(g, Pot, format, g.PotentialAtR)
}
func (g *RadGrid) DisplayPotentialC(Pot PotentialOp[complex128], format string, theta float64) {
	displayFunc(g, Pot, format, g.PotentialAtZ, theta)
}
func (g *RadGrid) DisplayForceRe(Pot PotentialOp[float64], format string) {
	displayFunc(g, Pot, format, g.ForceAtR)
}
func (g *RadGrid) DisplayForceC(Pot PotentialOp[complex128], format string, theta float64) {
	displayFunc(g, Pot, format, g.ForceAtZ, theta)
}

func (g *RadGrid) PrintPotentToFileRe(Pot PotentialOp[float64], filename string, format string) error {
	err := functionToFile(g, Pot, filename, format, g.PotentialAtR)
	return err
}

func (g *RadGrid) PrintPotentToFileZ(Pot PotentialOp[complex128], filename string, format string, theta float64) error {
	err := functionToFile(g, Pot, filename, format, g.PotentialAtZ, theta)
	return err
}

func (g *RadGrid) PrintForceToFileRe(Pot PotentialOp[float64], filename string, format string) error {
	err := functionToFile(g, Pot, filename, format, g.ForceAtR)
	return err
}

func (g *RadGrid) PrintForceToFileZ(Pot PotentialOp[complex128], filename string, format string, theta float64) error {
	err := functionToFile(g, Pot, filename, format, g.ForceAtZ, theta)
	return err
}

func (g *RadGrid) PrintVectorToFileRe(vec []float64, filename string, format string) error {
	err := vectorToFile(g, vec, filename, format)
	return err
}

func (g *RadGrid) PrintVectorToFileZ(vec []complex128, filename string, format string) error {
	err := vectorToFile(g, vec, filename, format)
	return err
}

func RVecToComplexVec(rVec []float64) []complex128 {
	result := make([]complex128, len(rVec))
	for i := 0; i < len(rVec); i++ {
		result[i] = complex(rVec[i], 0)
	}
	return result
}
