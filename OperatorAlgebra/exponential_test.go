package OperatorAlgebra

import (
	"GoProject/gridData"
	"gonum.org/v1/gonum/mat"
	"math"
	"math/cmplx"
	"testing"
)

func TestDVRCompletedExponentials(t *testing.T) {
	g, _ := gridData.NewRGrid(-5, 5, 32)
	k := NewKeDVR(g, 1)
	before := mat.DenseCopyOf(k.GetMat())
	values, vectors, err := k.RealDiagonalize()
	if err != nil {
		t.Fatal(err)
	}
	if !mat.Equal(before, k.GetMat()) {
		t.Fatal("diagonalization changed kinetic matrix")
	}
	x := make([]complex128, 32)
	for i := range x {
		x[i] = complex(vectors.At(i, 0), 0)
	}
	out, err := k.ExpIdt(.1, x)
	if err != nil {
		t.Fatal(err)
	}
	phase := cmplx.Exp(complex(0, .1*values[0]))
	for i := range out {
		if cmplx.Abs(out[i]-phase*x[i]) > 1e-12 {
			t.Fatal("incorrect exponential phase")
		}
	}
	k.ExpIdtTo(-.1, out, out)
	for i := range out {
		if cmplx.Abs(out[i]-x[i]) > 1e-12 {
			t.Fatal("aliasing or reversibility failure")
		}
	}
	realState := make([]float64, 32)
	for i := range realState {
		realState[i] = vectors.At(i, 0)
	}
	k.ExpDtTo(-.1, realState, realState)
	for i := range realState {
		if math.Abs(realState[i]-math.Exp(-.1*values[0])*vectors.At(i, 0)) > 1e-12 {
			t.Fatal("real exponential failure")
		}
	}
}
