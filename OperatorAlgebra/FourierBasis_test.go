package OperatorAlgebra

import (
	"GoProject/gridData"
	"math"
	"math/cmplx"
	"testing"
)

func TestFourierOperatorsOnPlaneWaves(t *testing.T) {
	for _, n := range []uint32{7, 8} {
		grid, err := gridData.NewFromLength(2*math.Pi, n)
		if err != nil {
			t.Fatal(err)
		}
		basis := FFTInit(grid, 2)
		t.Cleanup(basis.Clean)
		for _, mode := range []int{-3, -1, 0, 1, 3} {
			psi, out := make([]complex128, n), make([]complex128, n)
			for i, x := range grid.RValues() {
				psi[i] = cmplx.Exp(complex(0, float64(mode)*x))
			}
			basis.MomentumOp(psi, out)
			for i := range psi {
				if cmplx.Abs(out[i]-complex(float64(mode), 0)*psi[i]) > 1e-11 {
					t.Fatalf("n=%d mode=%d: incorrect momentum", n, mode)
				}
			}
			basis.LaplacianOp(psi, out)
			for i := range psi {
				if cmplx.Abs(out[i]-complex(float64(mode*mode)/4, 0)*psi[i]) > 1e-11 {
					t.Fatalf("n=%d mode=%d: incorrect kinetic energy", n, mode)
				}
			}
		}
	}
}
