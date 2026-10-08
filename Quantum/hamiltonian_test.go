package Quantum

import (
	"GoProject/gridData"
	"fmt"
	"math"
	"math/cmplx"
	"os"
	"slices"
	"testing"

	"gonum.org/v1/gonum/mat"
)

func TestHamiltonianOp(t *testing.T) {
	t.Chdir(t.TempDir())
	grid, err := gridData.NewFromLength(16, 48)
	if err != nil {
		t.Fatal(err)
	}
	h := NewHamil(grid, 1, gridData.Harmonic[float64]{ForceConst: 1})
	defer h.Close()
	kinetic := mat.DenseCopyOf(h.kinE.(*dvrKinetic).basis.GetMat())
	first := h.EvaluateOp()
	h.Mat()
	second := h.EvaluateOp()
	if !mat.EqualApprox(first, second, 1e-14) {
		t.Fatal("rebuilding H accumulated the potential")
	}
	if !mat.Equal(kinetic, h.kinE.(*dvrKinetic).basis.GetMat()) {
		t.Fatal("building H modified the cached kinetic matrix")
	}
	for i, x := range grid.RValues() {
		if math.Abs(first.At(i, i)-kinetic.At(i, i)-.5*x*x) > 1e-13 {
			t.Fatalf("potential was not added to diagonal %d", i)
		}
	}
	if _, err := os.Stat("potent.dat"); !os.IsNotExist(err) {
		t.Fatal("constructing H must not write potent.dat")
	}
}

type hamiltonianConstructor func(*gridData.RadGrid, float64, gridData.PotentialOp[float64]) (*HamiltonianOp, error)

var hamiltonianConstructors = map[string]hamiltonianConstructor{
	"DVR": NewHamiltonian, "Fourier": NewFourierHamiltonian,
}

func testHamiltonian(t *testing.T, build hamiltonianConstructor, n uint32, length, mass float64, pot gridData.PotentialOp[float64]) *HamiltonianOp {
	t.Helper()
	g, err := gridData.NewFromLength(length, n)
	if err != nil {
		t.Fatal(err)
	}
	h, err := build(g, mass, pot)
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(h.Close)
	return h
}

func TestHamiltonianApplicationAndAliasing(t *testing.T) {
	for name, build := range hamiltonianConstructors {
		for _, n := range []uint32{7, 8} {
			t.Run(fmt.Sprintf("%s/%d", name, n), func(t *testing.T) {
				h := testHamiltonian(t, build, n, 8, 2, gridData.Harmonic[float64]{ForceConst: 3})
				matrix := h.EvaluateOp()
				psiR, psiC := make([]float64, n), make([]complex128, n)
				wantR, wantC := make([]float64, n), make([]complex128, n)
				for i := range psiR {
					psiR[i] = math.Sin(float64(i) + .3)
					psiC[i] = complex(psiR[i], math.Cos(.7*float64(i)))
				}
				for i := range psiR {
					for j := range psiR {
						if math.Abs(matrix.At(i, j)-matrix.At(j, i)) > 1e-12 {
							t.Fatal("H is not symmetric")
						}
						wantR[i] += matrix.At(i, j) * psiR[j]
						wantC[i] += complex(matrix.At(i, j), 0) * psiC[j]
					}
				}
				outR, outC := make([]float64, n), make([]complex128, n)
				savedR, savedC := slices.Clone(psiR), slices.Clone(psiC)
				h.ApplyHReal(psiR, outR)
				h.ApplyHComplex(psiC, outC)
				if !slices.Equal(psiR, savedR) || !slices.Equal(psiC, savedC) {
					t.Fatal("out-of-place application changed the input")
				}
				assertStates(t, outR, outC, wantR, wantC)
				h.ApplyHReal(psiR, psiR)
				h.ApplyHComplex(psiC, psiC)
				assertStates(t, psiR, psiC, wantR, wantC)
				overlapR := append(slices.Clone(savedR), 0)
				overlapC := append(slices.Clone(savedC), 0)
				h.ApplyHReal(overlapR[:n], overlapR[1:])
				h.ApplyHComplex(overlapC[:n], overlapC[1:])
				assertStates(t, overlapR[1:], overlapC[1:], wantR, wantC)
				matrix.Set(0, 0, 1e9)
				h.ApplyHReal(savedR, outR)
				h.ApplyHComplex(savedC, outC)
				assertStates(t, outR, outC, wantR, wantC)
			})
		}
	}
}

func assertStates(t *testing.T, gotR []float64, gotC []complex128, wantR []float64, wantC []complex128) {
	t.Helper()
	for i := range wantR {
		if !(math.Abs(gotR[i]-wantR[i]) <= 1e-10) || !(cmplx.Abs(gotC[i]-wantC[i]) <= 1e-10) {
			t.Fatalf("H psi at %d: real=%g want %g, complex=%v want %v", i, gotR[i], wantR[i], gotC[i], wantC[i])
		}
	}
}

func TestHamiltonianHarmonicSpectrumAndEnergy(t *testing.T) {
	// m=2 and force constant=8 give omega=2, so E_n=(n+1/2)*2.
	for name, build := range hamiltonianConstructors {
		t.Run(name, func(t *testing.T) {
			h := testHamiltonian(t, build, 96, 16, 2, gridData.Harmonic[float64]{ForceConst: 8})
			before := h.EvaluateOp()
			values, vectors, err := h.RealDiagonalize()
			if err != nil {
				t.Fatal(err)
			}
			for state := 0; state < 4; state++ {
				want := 2*float64(state) + 1
				if !(math.Abs(values[state]-want) <= 1e-7) {
					t.Errorf("E_%d=%g, want %g", state, values[state], want)
				}
				psi, out := make([]float64, h.Dim()), make([]float64, h.Dim())
				mat.Col(psi, state, vectors)
				h.ApplyHReal(psi, out)
				for i := range psi {
					if math.Abs(out[i]-values[state]*psi[i]) > 1e-9 {
						t.Fatalf("eigenvector %d has a nonzero residual", state)
					}
				}
			}
			psi := make([]complex128, h.Dim())
			for i, x := range h.Grid().RValues() {
				psi[i] = complex(3*math.Exp(-2*x*x), 4*math.Exp(-2*x*x))
			}
			if energy := h.Energy(psi); !(math.Abs(energy-1) <= 1e-8) {
				t.Fatalf("ground-state energy=%g, want 1", energy)
			}
			if !mat.EqualApprox(before, h.EvaluateOp(), 1e-13) {
				t.Fatal("diagonalization changed H")
			}
		})
	}
}

func TestFourierHamiltonianPlaneWaves(t *testing.T) {
	for _, n := range []uint32{7, 8} {
		h := testHamiltonian(t, NewFourierHamiltonian, n, 2*math.Pi, 2, gridData.Polynomial[float64]{Coeffs: []float64{3}})
		for _, mode := range []int{-3, -1, 0, 1, 3} {
			psi, out := make([]complex128, n), make([]complex128, n)
			for i, x := range h.Grid().RValues() {
				psi[i] = cmplx.Exp(complex(0, float64(mode)*x))
			}
			h.ApplyHComplex(psi, out)
			energy := float64(mode*mode)/4 + 3
			for i := range psi {
				if cmplx.Abs(out[i]-complex(energy, 0)*psi[i]) > 1e-11 {
					t.Fatalf("n=%d mode=%d index=%d: incorrect plane-wave energy", n, mode, i)
				}
			}
		}
	}
}

func TestHamiltonianTimeDependentPotential(t *testing.T) {
	for name, build := range hamiltonianConstructors {
		t.Run(name, func(t *testing.T) {
			h := testHamiltonian(t, build, 8, 8, 1, nil)
			psiR, psiC := make([]float64, h.Dim()), make([]complex128, h.Dim())
			for i := range psiR {
				psiR[i] = float64(i + 1)
				psiC[i] = complex(psiR[i], -psiR[i])
			}
			kinR, kinC := make([]float64, h.Dim()), make([]complex128, h.Dim())
			h.ApplyHReal(psiR, kinR)
			h.ApplyHComplex(psiC, kinC)
			static := make([]float64, h.Dim())
			for i := range static {
				static[i] = 10
			}
			if err := h.SetStaticPotential(static); err != nil {
				t.Fatal(err)
			}
			static[0] = 999
			copyOfStatic := h.VStatic()
			copyOfStatic[0] = -999
			if h.VStatic()[0] != 10 {
				t.Fatal("potential aliases caller memory")
			}
			h.SetTimeDependentPotential(func(time float64, out []float64) {
				for i := range out {
					out[i] = 3 + time*float64(i+1)
				}
			})
			if !h.IsTimeDependent() {
				t.Fatal("time-dependent potential was not enabled")
			}
			for _, time := range []float64{0, .25, 2} {
				outR, outC := make([]float64, h.Dim()), make([]complex128, h.Dim())
				h.ApplyHRealAt(psiR, outR, time)
				h.ApplyHComplexAt(psiC, outC, time)
				potential := make([]float64, h.Dim())
				h.PotentialAt(time, potential)
				matrix := h.MatrixAt(time)
				var numerator, norm float64
				for i := range psiR {
					vi := 3 + time*float64(i+1)
					wantC := kinC[i] + complex(vi, 0)*psiC[i]
					if math.Abs(outR[i]-kinR[i]-vi*psiR[i]) > 1e-10 || cmplx.Abs(outC[i]-wantC) > 1e-10 || potential[i] != vi {
						t.Fatalf("incorrect time-dependent potential at t=%g, index=%d", time, i)
					}
					var fromMatrix complex128
					for j := range psiC {
						fromMatrix += complex(matrix.At(i, j), 0) * psiC[j]
					}
					if cmplx.Abs(outC[i]-fromMatrix) > 1e-10 {
						t.Fatal("time-dependent matrix differs from application")
					}
					numerator += real(cmplx.Conj(psiC[i]) * wantC)
					norm += real(cmplx.Conj(psiC[i]) * psiC[i])
				}
				if math.Abs(h.EnergyAt(psiC, time)-numerator/norm) > 1e-10 {
					t.Fatal("incorrect time-dependent energy")
				}
			}
			h.SetTimeDependentPotential(nil)
			if h.IsTimeDependent() || h.VStatic()[0] != 10 {
				t.Fatal("static potential was not restored")
			}
			h.SetTimeDependentPotential(func(_ float64, out []float64) { clear(out) })
			if err := h.SetStaticPotential(make([]float64, h.Dim())); err != nil {
				t.Fatal(err)
			}
			if h.IsTimeDependent() {
				t.Fatal("SetStaticPotential did not disable the callback")
			}
		})
	}
}

func TestHamiltonianValidationAndGridSnapshot(t *testing.T) {
	g, err := gridData.NewFromLength(8, 8)
	if err != nil {
		t.Fatal(err)
	}
	for name, build := range hamiltonianConstructors {
		t.Run(name, func(t *testing.T) {
			for _, mass := range []float64{0, -1, math.NaN(), math.Inf(1)} {
				if _, err := build(g, mass, nil); err == nil {
					t.Fatalf("accepted mass %g", mass)
				}
			}
			if _, err := build(nil, 1, nil); err == nil {
				t.Fatal("accepted nil grid")
			}
			if _, err := build(&gridData.RadGrid{}, 1, nil); err == nil {
				t.Fatal("accepted empty grid")
			}
			if _, err := build(g, 1, gridData.Polynomial[float64]{Coeffs: []float64{math.NaN()}}); err == nil {
				t.Fatal("accepted non-finite potential")
			}
			copyGrid := *g
			h, err := build(&copyGrid, 1, nil)
			if err != nil {
				t.Fatal(err)
			}
			t.Cleanup(h.Close)
			before := h.EvaluateOp()
			if err := copyGrid.ReDefineMinMax(-10, 10, 16); err != nil {
				t.Fatal(err)
			}
			if err := h.Grid().ReDefine(32); err != nil {
				t.Fatal(err)
			}
			if h.Dim() != 8 || !mat.EqualApprox(before, h.EvaluateOp(), 1e-13) {
				t.Fatal("caller grid mutation changed H")
			}
			if err := h.SetStaticPotential([]float64{1}); err == nil {
				t.Fatal("accepted wrong potential length")
			}
			bad := make([]float64, 8)
			bad[3] = math.Inf(1)
			if err := h.SetStaticPotential(bad); err == nil {
				t.Fatal("accepted non-finite potential")
			}
			mustPanic(t, func() { h.ApplyHReal(make([]float64, 7), make([]float64, 8)) })
			mustPanic(t, func() { h.ApplyHComplex(make([]complex128, 8), make([]complex128, 9)) })
			h.Close()
			h.Close()
			mustPanic(t, func() { h.ApplyHReal(make([]float64, 8), make([]float64, 8)) })
		})
	}
}

func mustPanic(t *testing.T, f func()) {
	t.Helper()
	defer func() {
		if recover() == nil {
			t.Error("expected a panic")
		}
	}()
	f()
}
