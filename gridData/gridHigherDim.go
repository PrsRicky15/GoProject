package gridData

// Restores the NDimGrid API from the project's PrashantCommit/webdev history
// (gridHigherDim.go at 8cfc52f), with checked indexing and input ownership.
import (
	"bufio"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

// NDimGrid is a Cartesian tensor grid. Each axis samples [min,max); the last
// coordinate varies fastest, preserving the original multidimensional layout.
// It stores axis metadata, not a materialized array of all coordinate tuples.
type NDimGrid struct {
	axes    []RadGrid
	strides []uint64
	total   uint64
}

func NewNDimGrid(mins, maxs []float64, nPoints []uint32) (*NDimGrid, error) {
	if len(mins) == 0 || len(mins) != len(maxs) || len(mins) != len(nPoints) {
		return nil, fmt.Errorf("NDimGrid: nonempty, equal-length bounds and point counts required")
	}
	g := &NDimGrid{axes: make([]RadGrid, len(mins)), strides: make([]uint64, len(mins)), total: 1}
	for d := len(mins) - 1; d >= 0; d-- {
		length := maxs[d] - mins[d]
		if !ndFinite(mins[d]) || !ndFinite(maxs[d]) || !ndFinite(length) || length <= 0 || nPoints[d] == 0 {
			return nil, fmt.Errorf("NDimGrid: invalid bounds or point count on axis %d", d)
		}
		if g.total > uint64(^uint(0)>>1)/uint64(nPoints[d]) {
			return nil, fmt.Errorf("NDimGrid: total point count exceeds addressable vector length")
		}
		axis, err := NewRGrid(mins[d], maxs[d], nPoints[d])
		if err != nil {
			return nil, err
		}
		if !ndFinite(axis.DeltaR()) || axis.DeltaR() <= 0 || !ndFinite(axis.DeltaK()) || !ndFinite(axis.CutoffE()) {
			return nil, fmt.Errorf("NDimGrid: unrepresentable spacing or frequency on axis %d", d)
		}
		g.axes[d] = *axis
		g.strides[d] = g.total
		g.total *= uint64(nPoints[d])
	}
	return g, nil
}

func ndFinite(x float64) bool { return !math.IsNaN(x) && !math.IsInf(x, 0) }

func NewNDimGridUniform(nDim uint32, min, max float64, nPoints uint32) (*NDimGrid, error) {
	mins, maxs, ns := make([]float64, nDim), make([]float64, nDim), make([]uint32, nDim)
	for d := range mins {
		mins[d], maxs[d], ns[d] = min, max, nPoints
	}
	return NewNDimGrid(mins, maxs, ns)
}

func NewNDimGridFromLengths(lengths []float64, nPoints []uint32) (*NDimGrid, error) {
	mins, maxs := make([]float64, len(lengths)), make([]float64, len(lengths))
	for d, l := range lengths {
		mins[d], maxs[d] = -l/2, l/2
	}
	return NewNDimGrid(mins, maxs, nPoints)
}

// Clone returns an independent grid definition.
func (g *NDimGrid) Clone() *NDimGrid {
	return &NDimGrid{axes: append([]RadGrid(nil), g.axes...), strides: append([]uint64(nil), g.strides...), total: g.total}
}

func (g *NDimGrid) ReDefine(nPoints []uint32) error {
	return g.ReDefineMinMax(g.Mins(), g.Maxs(), nPoints)
}
func (g *NDimGrid) ReDefineMinMax(mins, maxs []float64, nPoints []uint32) error {
	next, err := NewNDimGrid(mins, maxs, nPoints)
	if err == nil {
		*g = *next
	}
	return err
}
func (g *NDimGrid) ReDefineLengths(lengths []float64, nPoints []uint32) error {
	next, err := NewNDimGridFromLengths(lengths, nPoints)
	if err == nil {
		*g = *next
	}
	return err
}

func (g *NDimGrid) NDim() uint32        { return uint32(len(g.axes)) }
func (g *NDimGrid) TotalPoints() uint64 { return g.total }
func (g *NDimGrid) Mins() []float64 {
	x := make([]float64, len(g.axes))
	for d := range x {
		x[d] = g.Min(uint32(d))
	}
	return x
}
func (g *NDimGrid) Maxs() []float64 {
	x := make([]float64, len(g.axes))
	for d := range x {
		x[d] = g.Max(uint32(d))
	}
	return x
}
func (g *NDimGrid) NPoints() []uint32 {
	x := make([]uint32, len(g.axes))
	for d := range x {
		x[d] = g.NumPoints(uint32(d))
	}
	return x
}
func (g *NDimGrid) Min(d uint32) float64      { return g.axes[d].RMin() }
func (g *NDimGrid) Max(d uint32) float64      { return g.axes[d].RMax() }
func (g *NDimGrid) NumPoints(d uint32) uint32 { return g.axes[d].NPoints() }
func (g *NDimGrid) DeltaR(d uint32) float64   { return g.axes[d].DeltaR() }
func (g *NDimGrid) Length(d uint32) float64   { return g.axes[d].Length() }
func (g *NDimGrid) DeltaK(d uint32) float64   { return g.axes[d].DeltaK() }
func (g *NDimGrid) KMin(d uint32) float64     { return g.axes[d].KMin() }
func (g *NDimGrid) KMax(d uint32) float64     { return g.axes[d].KMax() }
func (g *NDimGrid) CutoffE(d uint32) float64  { return g.axes[d].CutoffE() }
func (g *NDimGrid) Strides() []uint64         { return append([]uint64(nil), g.strides...) }
func (g *NDimGrid) VolumeElement() float64 {
	v := 1.0
	for d := range g.axes {
		v *= g.axes[d].DeltaR()
	}
	return v
}
func (g *NDimGrid) RValues(d uint32) ([]float64, error) {
	if d >= g.NDim() {
		return nil, fmt.Errorf("NDimGrid: invalid axis %d", d)
	}
	return g.axes[d].RValues(), nil
}

// KValues uses the same FFT ordering as RadGrid, including a negative Nyquist bin.
func (g *NDimGrid) KValues(d uint32) []float64 {
	if d >= g.NDim() {
		return nil
	}
	return g.axes[d].KValues()
}

// PointTo writes coordinates without allocating. Invalid indices/buffer sizes panic.
func (g *NDimGrid) PointTo(flat uint64, out []float64) {
	if flat >= g.total || len(out) != len(g.axes) {
		panic("NDimGrid: invalid point index or coordinate buffer")
	}
	for d, a := range g.axes {
		out[d] = a.RMin() + float64(flat/g.strides[d]%uint64(a.NPoints()))*a.DeltaR()
	}
}
func (g *NDimGrid) GetPoint(flat uint64) []float64 {
	x := make([]float64, len(g.axes))
	g.PointTo(flat, x)
	return x
}
func (g *NDimGrid) GetFlatIndex(indices []uint32) uint64 {
	if len(indices) != len(g.axes) || g.total == 0 {
		panic("NDimGrid: invalid index dimension")
	}
	var flat uint64
	for d, i := range indices {
		if i >= g.axes[d].NPoints() {
			panic("NDimGrid: index outside axis")
		}
		flat += uint64(i) * g.strides[d]
	}
	return flat
}

// FunctionOnGridInPlace reuses the coordinate slice: fn must not retain or alter it.
func (g *NDimGrid) FunctionOnGridInPlace(fn func([]float64) float64, out []float64) {
	if fn == nil || uint64(len(out)) != g.total {
		panic("NDimGrid: invalid function or output length")
	}
	x := make([]float64, len(g.axes))
	for i := range out {
		g.PointTo(uint64(i), x)
		out[i] = fn(x)
	}
}

// NDimPotentialOp preserves the original NDimEvaluateAt method without forcing
// one-coordinate PotentialOp implementations to invent multidimensional behavior.
type NDimPotentialOp interface{ NDimEvaluateAt([]float64) float64 }
type NDimPotentialFunc func([]float64) float64

func (f NDimPotentialFunc) NDimEvaluateAt(x []float64) float64 { return f(x) }

func (g *NDimGrid) PotentialOnGrid(pot NDimPotentialOp) []float64 {
	out := make([]float64, int(g.total))
	g.PotentialOnGridInPlace(pot, out)
	return out
}
func (g *NDimGrid) PotentialOnGridInPlace(pot NDimPotentialOp, out []float64) {
	if pot == nil {
		panic("NDimGrid: nil potential")
	}
	g.FunctionOnGridInPlace(pot.NDimEvaluateAt, out)
}

// NDimEvaluateAt extends the harmonic well by summing over Cartesian coordinates.
func (h Harmonic[T]) NDimEvaluateAt(x []T) T {
	var v T
	for _, xi := range x {
		v += h.EvaluateAt(xi)
	}
	return v
}

// SeparablePotential explicitly constructs V(x)=sum_d V_d(x_d).
// Use NDimPotentialFunc instead when coordinates are coupled.
type SeparablePotential []PotentialOp[float64]

func (p SeparablePotential) NDimEvaluateAt(x []float64) float64 {
	if len(p) != len(x) {
		panic("SeparablePotential: dimension mismatch")
	}
	v := 0.0
	for d, pot := range p {
		if pot != nil {
			v += pot.EvaluateAt(x[d])
		}
	}
	return v
}

func (g *NDimGrid) String() string {
	return fmt.Sprintf("NDimGrid{nDim: %d, nPoints: %v, mins: %v, maxs: %v}", g.NDim(), g.NPoints(), g.Mins(), g.Maxs())
}
func (g *NDimGrid) DisplayInfo() {
	fmt.Println(g)
	for d := range g.axes {
		fmt.Printf("Axis %d: ", d)
		g.axes[d].DisplayInfo()
	}
}

func (g *NDimGrid) PrintVectorToFile(vec []float64, filename, format string) (err error) {
	if uint64(len(vec)) != g.total {
		return fmt.Errorf("NDimGrid: vector length mismatch")
	}
	f, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer func() {
		if e := f.Close(); err == nil {
			err = e
		}
	}()
	w := bufio.NewWriter(f)
	x := make([]float64, g.NDim())
	for i, v := range vec {
		g.PointTo(uint64(i), x)
		for _, xi := range x {
			if _, err = fmt.Fprintf(w, format+" ", xi); err != nil {
				return err
			}
		}
		if _, err = fmt.Fprintf(w, format+"\n", v); err != nil {
			return err
		}
	}
	return w.Flush()
}

// NewNDimGridFromFile reads nDim, mins, maxs and nPoints from ndgrid.inp.
// Lists are comma-separated; key order is arbitrary and malformed input is an error.
func NewNDimGridFromFile(dirPath string) (*NDimGrid, error) {
	f, err := os.Open(filepath.Join(dirPath, "ndgrid.inp"))
	if err != nil {
		return nil, err
	}
	defer f.Close()
	fields := make(map[string]string)
	s := bufio.NewScanner(f)
	for s.Scan() {
		line := strings.TrimSpace(strings.SplitN(s.Text(), "#", 2)[0])
		if line == "" {
			continue
		}
		key, value, ok := strings.Cut(line, ":")
		key = strings.TrimSpace(key)
		if !ok || (key != "nDim" && key != "mins" && key != "maxs" && key != "nPoints") {
			return nil, fmt.Errorf("ndgrid.inp: invalid field %q", line)
		}
		if _, duplicate := fields[key]; duplicate {
			return nil, fmt.Errorf("ndgrid.inp: duplicate %s", key)
		}
		fields[key] = strings.TrimSpace(value)
	}
	if err := s.Err(); err != nil {
		return nil, err
	}
	dims, err := strconv.ParseUint(fields["nDim"], 10, 32)
	if err != nil || dims == 0 {
		return nil, fmt.Errorf("ndgrid.inp: invalid nDim")
	}
	lo, hi, ns := strings.Split(fields["mins"], ","), strings.Split(fields["maxs"], ","), strings.Split(fields["nPoints"], ",")
	if uint64(len(lo)) != dims || len(hi) != len(lo) || len(ns) != len(lo) {
		return nil, fmt.Errorf("ndgrid.inp: axis list length mismatch")
	}
	mins, maxs, points := make([]float64, len(lo)), make([]float64, len(lo)), make([]uint32, len(lo))
	for d := range lo {
		mins[d], err = strconv.ParseFloat(strings.TrimSpace(lo[d]), 64)
		if err != nil {
			return nil, err
		}
		maxs[d], err = strconv.ParseFloat(strings.TrimSpace(hi[d]), 64)
		if err != nil {
			return nil, err
		}
		n, e := strconv.ParseUint(strings.TrimSpace(ns[d]), 10, 32)
		if e != nil {
			return nil, e
		}
		points[d] = uint32(n)
	}
	return NewNDimGrid(mins, maxs, points)
}
