package gridData

import (
	"math"
	"os"
	"path/filepath"
	"reflect"
	"testing"
)

func TestNDimGridOrderingAndOwnership(t *testing.T) {
	lo, hi, ns := []float64{-1, 2, 10}, []float64{1, 8, 18}, []uint32{2, 3, 4}
	g, err := NewNDimGrid(lo, hi, ns)
	if err != nil {
		t.Fatal(err)
	}
	lo[0], hi[0], ns[0] = 99, 100, 200
	if g.TotalPoints() != 24 || !reflect.DeepEqual(g.Strides(), []uint64{12, 4, 1}) {
		t.Fatal("wrong tensor shape/strides")
	}
	flat := uint64(0)
	for i := uint32(0); i < 2; i++ {
		for j := uint32(0); j < 3; j++ {
			for k := uint32(0); k < 4; k++ {
				if g.GetFlatIndex([]uint32{i, j, k}) != flat {
					t.Fatal("changed original last-axis-fastest ordering")
				}
				want := []float64{-1 + float64(i), 2 + 2*float64(j), 10 + 2*float64(k)}
				if !reflect.DeepEqual(g.GetPoint(flat), want) {
					t.Fatalf("point %d: %v want %v", flat, g.GetPoint(flat), want)
				}
				flat++
			}
		}
	}
	g.Mins()[0], g.Maxs()[0], g.NPoints()[0] = 100, 200, 400
	copyGrid := g.Clone()
	if err := copyGrid.ReDefine([]uint32{1, 1, 1}); err != nil {
		t.Fatal(err)
	}
	if g.Min(0) != -1 || g.TotalPoints() != 24 || g.VolumeElement() != 4 {
		t.Fatal("metadata aliases caller or incorrect volume")
	}
	k := g.KValues(2)
	for i, m := range []float64{0, 1, -2, -1} {
		if math.Abs(k[i]-m*math.Pi/4) > 1e-14 {
			t.Fatal("wrong FFT ordering")
		}
	}
	out := make([]float64, 24)
	g.PotentialOnGridInPlace(NDimPotentialFunc(func(x []float64) float64 { return x[0]*x[1] + x[2] }), out)
	if out[0] != 8 || out[23] != 16 {
		t.Fatal("coupled potential sampled incorrectly")
	}
}

func TestNDimGridValidation(t *testing.T) {
	for _, tc := range []struct {
		lo, hi []float64
		ns     []uint32
	}{
		{nil, nil, nil}, {[]float64{0}, nil, []uint32{1}}, {[]float64{0}, []float64{1}, []uint32{0}},
		{[]float64{math.NaN()}, []float64{1}, []uint32{1}}, {[]float64{0}, []float64{math.Inf(1)}, []uint32{1}},
		{[]float64{1}, []float64{0}, []uint32{1}}, {[]float64{0}, []float64{1e-310}, []uint32{100}},
		{[]float64{0, 0, 0}, []float64{1, 1, 1}, []uint32{math.MaxUint32, math.MaxUint32, 2}},
	} {
		if _, err := NewNDimGrid(tc.lo, tc.hi, tc.ns); err == nil {
			t.Fatalf("accepted invalid grid: %+v", tc)
		}
	}
	g, _ := NewNDimGridUniform(2, -1, 1, 4)
	if err := g.ReDefineMinMax([]float64{0, 0}, []float64{1}, []uint32{1, 2}); err == nil || g.TotalPoints() != 16 {
		t.Fatal("failed redefine must preserve old grid")
	}
	for _, fn := range []func(){func() { g.GetPoint(16) }, func() { g.GetFlatIndex([]uint32{0}) }, func() { g.GetFlatIndex([]uint32{4, 0}) }, func() { g.PotentialOnGridInPlace(Harmonic[float64]{ForceConst: 1}, make([]float64, 2)) }} {
		func() {
			defer func() {
				if recover() == nil {
					t.Error("expected invalid indexing/buffer to panic")
				}
			}()
			fn()
		}()
	}
	if _, err := g.RValues(2); err == nil || g.KValues(2) != nil {
		t.Fatal("accepted invalid dimension")
	}
}

func TestNDimGridFileAndPotentialAdapters(t *testing.T) {
	dir := t.TempDir()
	file := filepath.Join(dir, "ndgrid.inp")
	valid := "mins: -2,-3\nnPoints: 4,6\nmaxs: 2,3 # upper bounds\nnDim: 2\n"
	for _, bad := range []string{"nDim: -1", "nDim: 2\nmins: 0,0\nmaxs: 1,1\nnPoints: -1,3", valid + "nDim: 3\n", "nDim: 1\nmins: bad\nmaxs: 1\nnPoints: 3"} {
		if err := os.WriteFile(file, []byte(bad), 0600); err != nil {
			t.Fatal(err)
		}
		if _, err := NewNDimGridFromFile(dir); err == nil {
			t.Fatalf("accepted malformed file %q", bad)
		}
	}
	if err := os.WriteFile(file, []byte(valid), 0600); err != nil {
		t.Fatal(err)
	}
	g, err := NewNDimGridFromFile(dir)
	if err != nil {
		t.Fatal(err)
	}
	if g.TotalPoints() != 24 {
		t.Fatal("wrong file grid")
	}
	h := Harmonic[float64]{ForceConst: 2, Cen: 1}
	if h.NDimEvaluateAt([]float64{2, 3}) != 5 {
		t.Fatal("harmonic extension")
	}
	p := SeparablePotential{h, Polynomial[float64]{Coeffs: []float64{3, 2}}}
	if p.NDimEvaluateAt([]float64{2, 3}) != 10 {
		t.Fatal("separable adapter")
	}
	if err := g.PrintVectorToFile([]float64{1}, filepath.Join(dir, "bad.dat"), "%g"); err == nil {
		t.Fatal("accepted wrong file vector size")
	}
	if err := g.PrintVectorToFile(g.PotentialOnGrid(h), filepath.Join(dir, "v.dat"), "%.8g"); err != nil {
		t.Fatal(err)
	}
}
