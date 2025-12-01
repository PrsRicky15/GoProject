package gridData

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"os"
	"strconv"
	"strings"
)

// NDimGrid represents an n-dimensional grid
type NDimGrid struct {
	nDim     uint32    // Number of dimensions
	mins     []float64 // Minimum values for each dimension
	maxs     []float64 // Maximum values for each dimension
	nPoints  []uint32  // Number of points in each dimension
	gridData []setGrid // Grid data for each dimension
	cutoffE  []float64 // Cutoff energy for each dimension
}

func (g *NDimGrid) String() string {
	return fmt.Sprintf("NDimGrid{nDim: %d, nPoints: %v, mins: %v, maxs: %v}",
		g.nDim, g.nPoints, g.mins, g.maxs)
}

func (g *NDimGrid) DisplayInfo() {
	fmt.Println("#-------------------------------------------------------------")
	fmt.Printf("# %d-Dimensional Grid Parameters:\n", g.nDim)
	for dim := uint32(0); dim < g.nDim; dim++ {
		fmt.Printf("# Dimension %d:\n", dim)
		fmt.Printf("  Real Space - Min: %8.4g | Max: %8.4g | Dr: %8.4g\n",
			g.mins[dim], g.maxs[dim], g.gridData[dim].deltaS)
		fmt.Printf("  K Space    - Min: %8.4g | Max: %8.4g | Dk: %8.4g\n",
			g.gridData[dim].cMin, g.gridData[dim].cMax, g.gridData[dim].deltaCS)
		fmt.Printf("  Grid       - Length: %8.4g | Points: %5d | Cutoff Energy: %8.4g\n",
			g.gridData[dim].length, g.nPoints[dim], g.cutoffE[dim])
	}
	fmt.Println("#-------------------------------------------------------------")
}

// NewNDimGrid creates a new n-dimensional grid with specified parameters for each dimension
func NewNDimGrid(mins, maxs []float64, nPoints []uint32) (*NDimGrid, error) {
	nDim := len(mins)
	if len(maxs) != nDim || len(nPoints) != nDim {
		return nil, fmt.Errorf("inconsistent dimensions: mins=%d, maxs=%d, nPoints=%d",
			len(mins), len(maxs), len(nPoints))
	}

	gridData := make([]setGrid, nDim)
	cutoffE := make([]float64, nDim)

	for i := 0; i < nDim; i++ {
		gd, err := createGrid(mins[i], maxs[i], nPoints[i], fmt.Sprintf("Grid_dim%d", i))
		if err != nil {
			return nil, fmt.Errorf("error creating grid for dimension %d: %w", i, err)
		}
		gridData[i] = gd
		cutoffE[i] = math.Pow(gd.cMax, 2) / 2
	}

	return &NDimGrid{
		nDim:     uint32(nDim),
		mins:     mins,
		maxs:     maxs,
		nPoints:  nPoints,
		gridData: gridData,
		cutoffE:  cutoffE,
	}, nil
}

// NewNDimGridUniform creates a uniform n-dimensional grid (same parameters for all dimensions)
func NewNDimGridUniform(nDim uint32, min, max float64, nPoints uint32) (*NDimGrid, error) {
	mins := make([]float64, nDim)
	maxs := make([]float64, nDim)
	nPts := make([]uint32, nDim)

	for i := uint32(0); i < nDim; i++ {
		mins[i] = min
		maxs[i] = max
		nPts[i] = nPoints
	}

	return NewNDimGrid(mins, maxs, nPts)
}

// NewNDimGridFromFile reads grid parameters from a file
func NewNDimGridFromFile(dirPath string) (*NDimGrid, error) {
	info, err := os.Stat(dirPath)
	if err != nil {
		return nil, fmt.Errorf("directory check failed: %w", err)
	}
	if !info.IsDir() {
		return nil, fmt.Errorf("%s is not a directory", dirPath)
	}

	filePath := dirPath + "/ndgrid.inp"
	file, err := os.Open(filePath)
	if err != nil {
		return nil, fmt.Errorf("could not open %s: %w", filePath, err)
	}
	defer func(file *os.File) {
		err := file.Close()
		if err != nil {
			log.Fatal(err)
		}
	}(file)

	var nDim int
	var mins, maxs []float64
	var nPoints []uint32

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
		case "nDim":
			nDim, _ = strconv.Atoi(value)
			mins = make([]float64, nDim)
			maxs = make([]float64, nDim)
			nPoints = make([]uint32, nDim)
		case "mins":
			parseFloatArray(value, mins)
		case "maxs":
			parseFloatArray(value, maxs)
		case "nPoints":
			parseUintArray(value, nPoints)
		}
	}

	return NewNDimGrid(mins, maxs, nPoints)
}

// NewNDimGridFromLengths creates a grid centered at origin with specified lengths
func NewNDimGridFromLengths(lengths []float64, nPoints []uint32) (*NDimGrid, error) {
	nDim := len(lengths)
	mins := make([]float64, nDim)
	maxs := make([]float64, nDim)

	for i := 0; i < nDim; i++ {
		if lengths[i] <= 0 {
			return nil, fmt.Errorf("length must be positive for dimension %d", i)
		}
		halfLength := lengths[i] / 2
		mins[i] = -halfLength
		maxs[i] = halfLength
	}

	return NewNDimGrid(mins, maxs, nPoints)
}

func (g *NDimGrid) ReDefine(nPoints []uint32) error {
	return g.redefine(g.mins, g.maxs, nPoints)
}

func (g *NDimGrid) ReDefineMinMax(mins, maxs []float64, nPoints []uint32) error {
	return g.redefine(mins, maxs, nPoints)
}

func (g *NDimGrid) ReDefineLengths(lengths []float64, nPoints []uint32) error {
	nDim := len(lengths)
	mins := make([]float64, nDim)
	maxs := make([]float64, nDim)

	for i := 0; i < nDim; i++ {
		mins[i] = -0.5 * lengths[i]
		maxs[i] = 0.5 * lengths[i]
	}

	return g.redefine(mins, maxs, nPoints)
}

func (g *NDimGrid) redefine(mins, maxs []float64, nPoints []uint32) error {
	nDim := len(mins)
	gridData := make([]setGrid, nDim)
	cutoffE := make([]float64, nDim)

	for i := 0; i < nDim; i++ {
		gd, err := createGrid(mins[i], maxs[i], nPoints[i], fmt.Sprintf("Grid_dim%d", i))
		if err != nil {
			return err
		}
		gridData[i] = gd
		cutoffE[i] = math.Pow(gd.cMax, 2) / 2
	}

	g.nDim = uint32(nDim)
	g.mins = mins
	g.maxs = maxs
	g.nPoints = nPoints
	g.gridData = gridData
	g.cutoffE = cutoffE
	return nil
}

func (g *NDimGrid) NDim() uint32                { return g.nDim }
func (g *NDimGrid) Mins() []float64             { return g.mins }
func (g *NDimGrid) Maxs() []float64             { return g.maxs }
func (g *NDimGrid) NPoints() []uint32           { return g.nPoints }
func (g *NDimGrid) Min(dim uint32) float64      { return g.mins[dim] }
func (g *NDimGrid) Max(dim uint32) float64      { return g.maxs[dim] }
func (g *NDimGrid) NumPoints(dim uint32) uint32 { return g.nPoints[dim] }
func (g *NDimGrid) DeltaR(dim uint32) float64   { return g.gridData[dim].deltaS }
func (g *NDimGrid) Length(dim uint32) float64   { return g.gridData[dim].length }
func (g *NDimGrid) DeltaK(dim uint32) float64   { return g.gridData[dim].deltaCS }
func (g *NDimGrid) KMin(dim uint32) float64     { return g.gridData[dim].cMin }
func (g *NDimGrid) KMax(dim uint32) float64     { return g.gridData[dim].cMax }
func (g *NDimGrid) CutoffE(dim uint32) float64  { return g.cutoffE[dim] }

// TotalPoints returns the total number of grid points across all dimensions
func (g *NDimGrid) TotalPoints() uint64 {
	total := uint64(1)
	for i := uint32(0); i < g.nDim; i++ {
		total *= uint64(g.nPoints[i])
	}
	return total
}

// RValues returns grid points for a specific dimension
func (g *NDimGrid) RValues(dim uint32) []float64 {
	if dim >= g.nDim {
		return nil
	}
	result := make([]float64, g.nPoints[dim])
	for i := uint32(0); i < g.nPoints[dim]; i++ {
		result[i] = g.mins[dim] + float64(i)*g.gridData[dim].deltaS
	}
	return result
}

// KValues returns k-space points for a specific dimension
func (g *NDimGrid) KValues(dim uint32) []float64 {
	if dim >= g.nDim {
		return nil
	}
	// Implementation depends on generateConjugatePoints from your original code
	return generateConjugatePointsND(g, dim)
}

// GetPoint returns coordinates for a flat index
func (g *NDimGrid) GetPoint(flatIndex uint64) []float64 {
	coords := make([]float64, g.nDim)
	idx := flatIndex

	for dim := int(g.nDim) - 1; dim >= 0; dim-- {
		dimSize := uint64(g.nPoints[dim])
		localIdx := idx % dimSize
		coords[dim] = g.mins[dim] + float64(localIdx)*g.gridData[dim].deltaS
		idx /= dimSize
	}

	return coords
}

// GetFlatIndex converts multi-dimensional indices to flat index
func (g *NDimGrid) GetFlatIndex(indices []uint32) uint64 {
	flatIdx := uint64(0)
	stride := uint64(1)

	for dim := int(g.nDim) - 1; dim >= 0; dim-- {
		flatIdx += uint64(indices[dim]) * stride
		stride *= uint64(g.nPoints[dim])
	}

	return flatIdx
}

// FunctionOnGridInPlace evaluates a function on the entire grid
func (g *NDimGrid) FunctionOnGridInPlace(fn func([]float64) float64, f []float64) {
	if uint64(len(f)) != g.TotalPoints() {
		return
	}

	for i := uint64(0); i < g.TotalPoints(); i++ {
		coords := g.GetPoint(i)
		f[i] = fn(coords)
	}
}

// PotentialOnGrid evaluates a potential on the grid
func (g *NDimGrid) PotentialOnGrid(pot NDPotentialOp) []float64 {
	totalPts := g.TotalPoints()
	result := make([]float64, totalPts)

	for i := uint64(0); i < totalPts; i++ {
		coords := g.GetPoint(i)
		result[i] = pot.EvaluateAt(coords)
	}

	return result
}

// PrintVectorToFile writes a vector to file with coordinates
func (g *NDimGrid) PrintVectorToFile(vec []float64, filename string, format string) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer func(file *os.File) {
		err := file.Close()
		if err != nil {
			log.Fatal(err)
		}
	}(file)

	writer := bufio.NewWriter(file)
	defer func(writer *bufio.Writer) {
		err := writer.Flush()
		if err != nil {
			log.Fatal(err)
		}
	}(writer)

	for i := uint64(0); i < g.TotalPoints() && i < uint64(len(vec)); i++ {
		coords := g.GetPoint(i)

		for _, coord := range coords {
			_, err2 := fmt.Fprintf(writer, format+" ", coord)
			if err2 != nil {
				return err2
			}
		}

		_, err := fmt.Fprintf(writer, format+"\n", vec[i])
		if err != nil {
			return err
		}
	}

	return nil
}

// Helper functions
func parseFloatArray(s string, arr []float64) {
	parts := strings.Split(s, ",")
	for i, p := range parts {
		if i < len(arr) {
			arr[i], _ = strconv.ParseFloat(strings.TrimSpace(p), 64)
		}
	}
}

func parseUintArray(s string, arr []uint32) {
	parts := strings.Split(s, ",")
	for i, p := range parts {
		if i < len(arr) {
			val, _ := strconv.Atoi(strings.TrimSpace(p))
			arr[i] = uint32(val)
		}
	}
}

func generateConjugatePointsND(g *NDimGrid, dim uint32) []float64 {
	result := make([]float64, g.nPoints[dim])
	for i := uint32(0); i < g.nPoints[dim]; i++ {
		result[i] = g.gridData[dim].cMin + float64(i)*g.gridData[dim].deltaCS
	}
	return result
}

// NDPotentialOp interface for n-dimensional potentials
type NDPotentialOp interface {
	EvaluateAt(coords []float64) float64
	ForceAt(coords []float64) []float64
}
