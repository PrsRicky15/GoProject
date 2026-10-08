package OperatorAlgebra

import (
	"GoProject/gridData"
	"testing"
)

func TestKeDVR_Evaluate(t *testing.T) {
	rgrid, err := gridData.NewFromLength(10., 30)
	if err != nil {
		t.Fatal(err)
	}
	kinE := NewKeDVR(rgrid, 1.)
	kinetic := kinE.GetMat()
	rows, cols := kinetic.Dims()
	if rows != 30 || cols != 30 {
		t.Fatalf("kinetic matrix dimensions = %dx%d, want 30x30", rows, cols)
	}
}
