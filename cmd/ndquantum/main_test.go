package main

import (
	"encoding/csv"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"testing"

	"GoProject/Quantum"
	"GoProject/gridData"
	"gonum.org/v1/gonum/mat"
)

func TestExportedDensityNormalizationAndOrdering(t *testing.T) {
	g, err := gridData.NewNDimGrid([]float64{-1, -3}, []float64{1, 3}, []uint32{2, 3})
	if err != nil {
		t.Fatal(err)
	}
	vectors := mat.NewDense(6, 1, nil)
	for i := 0; i < 6; i++ {
		vectors.Set(i, 0, 1/math.Sqrt(6))
	}
	path := filepath.Join(t.TempDir(), "density.csv")
	if err := writeDensities(path, g, &Quantum.EigenResult{Values: []float64{0}, Vectors: vectors}); err != nil {
		t.Fatal(err)
	}
	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	rows, err := csv.NewReader(f).ReadAll()
	if err != nil {
		t.Fatal(err)
	}
	if len(rows) != 7 || len(rows[0]) != 3 {
		t.Fatal("missing density rows/columns")
	}
	var integral float64
	for i, row := range rows[1:] {
		x, _ := strconv.ParseFloat(row[0], 64)
		y, _ := strconv.ParseFloat(row[1], 64)
		rho, err := strconv.ParseFloat(row[2], 64)
		if err != nil {
			t.Fatal(err)
		}
		if x != float64(i/3)-1 || y != float64(i%3)*2-3 {
			t.Fatal("export changed original tensor ordering")
		}
		integral += rho * 2 // dx=1,dy=2
	}
	if math.Abs(integral-1) > 1e-14 {
		t.Fatalf("density integral=%g; want 1", integral)
	}
}
