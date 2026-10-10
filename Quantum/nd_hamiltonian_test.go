package Quantum

import (
	"context"
	"fmt"
	"math"
	"math/cmplx"
	"slices"
	"testing"

	"GoProject/gridData"
	"gonum.org/v1/gonum/mat"
)

type ndConstructor func(*gridData.NDimGrid, []float64, gridData.NDimPotentialOp) (*NDimHamiltonian, error)

var ndConstructors = map[string]ndConstructor{"DVR": NewNDimHamiltonian, "Fourier": NewNDimFourierHamiltonian}

func TestNDimFourierPlaneWaves(t *testing.T) {
	for dims := 1; dims <= 4; dims++ {
		t.Run(fmt.Sprint(dims), func(t *testing.T) {
			lo, hi, ns, masses := make([]float64, dims), make([]float64, dims), make([]uint32, dims), make([]float64, dims)
			modes := []int{1, -2, 2, -3}
			energy := 3.0
			for d := range lo {
				lo[d] = float64(d) - .7
				hi[d] = lo[d] + 2*math.Pi
				ns[d] = uint32(d + 3)
				masses[d] = float64(d + 1)
				energy += float64(modes[d]*modes[d]) / (2 * masses[d])
			}
			g, err := gridData.NewNDimGrid(lo, hi, ns)
			if err != nil {
				t.Fatal(err)
			}
			h, err := NewNDimFourierHamiltonian(g, masses, gridData.NDimPotentialFunc(func([]float64) float64 { return 3 }))
			if err != nil {
				t.Fatal(err)
			}
			defer h.Close()
			psi, out := make([]complex128, h.Dim()), make([]complex128, h.Dim())
			for i := range psi {
				phase := 0.0
				for d, x := range g.GetPoint(uint64(i)) {
					phase += float64(modes[d]) * x
				}
				psi[i] = cmplx.Exp(complex(0, phase))
			}
			h.ApplyHComplex(psi, out)
			for i := range out {
				if cmplx.Abs(out[i]-complex(energy, 0)*psi[i]) > 2e-12 {
					t.Fatalf("axis kinetic mismatch at %d", i)
				}
			}
			p, err := NewNDimSplitPropagator(h)
			if err != nil {
				t.Fatal(err)
			}
			before := slices.Clone(psi)
			if err := p.Step(psi, .3, .07); err != nil {
				t.Fatal(err)
			}
			for i := range psi {
				if cmplx.Abs(psi[i]-before[i]*cmplx.Exp(complex(0, -energy*.07))) > 3e-13 {
					t.Fatal("incorrect Fourier kinetic phase")
				}
			}
		})
	}
}

func TestNDimDenseApplicationAndTimeDependence(t *testing.T) {
	for name, build := range ndConstructors {
		t.Run(name, func(t *testing.T) {
			g, _ := gridData.NewNDimGrid([]float64{1, -2, 0}, []float64{5, 1, 4}, []uint32{3, 4, 2})
			h, err := build(g, []float64{1, 2, 3}, gridData.NDimPotentialFunc(func(x []float64) float64 { return x[0]*x[1] + 2*x[2] }))
			if err != nil {
				t.Fatal(err)
			}
			defer h.Close()
			h.SetTimeDependentFunction(func(x []float64, time float64) float64 { return x[0]*x[1] + 2*x[2] + time*x[0]*x[2] })
			matrix, err := h.MatrixAt(.7)
			if err != nil {
				t.Fatal(err)
			}
			n := h.Dim()
			r, c := make([]float64, n), make([]complex128, n)
			wantR, wantC := make([]float64, n), make([]complex128, n)
			for i := range r {
				r[i] = math.Sin(float64(i))
				c[i] = complex(r[i], math.Cos(float64(i)))
			}
			for i := range r {
				for j := range r {
					if math.Abs(matrix.At(i, j)-matrix.At(j, i)) > 1e-12 {
						t.Fatal("not symmetric")
					}
					wantR[i] += matrix.At(i, j) * r[j]
					wantC[i] += complex(matrix.At(i, j), 0) * c[j]
				}
			}
			// Deliberately overlap, rather than merely reuse the identical slice.
			rOverlap, cOverlap := append(slices.Clone(r), 0), append(slices.Clone(c), 0)
			h.ApplyHRealAt(rOverlap[:n], rOverlap[1:], .7)
			h.ApplyHComplexAt(cOverlap[:n], cOverlap[1:], .7)
			assertStates(t, rOverlap[1:], cOverlap[1:], wantR, wantC)
			before := h.Grid().TotalPoints()
			if err := g.ReDefine([]uint32{1, 1, 1}); err != nil {
				t.Fatal(err)
			}
			if h.Grid().TotalPoints() != before {
				t.Fatal("grid was not copied")
			}
			h.SetTimeDependentFunction(nil)
			if h.IsTimeDependent() {
				t.Fatal("callback not cleared")
			}
		})
	}
}

func TestNDimCoupledOscillatorSpectrum(t *testing.T) {
	// V=(x²+y²)/2+c*x*y: normal-mode frequencies sqrt(1+c),sqrt(1-c).
	for name, build := range ndConstructors {
		t.Run(name, func(t *testing.T) {
			g, _ := gridData.NewNDimGridUniform(2, -6, 6, 20)
			const coupling = .3
			h, err := build(g, []float64{1, 1}, gridData.NDimPotentialFunc(func(x []float64) float64 { return .5*(x[0]*x[0]+x[1]*x[1]) + coupling*x[0]*x[1] }))
			if err != nil {
				t.Fatal(err)
			}
			defer h.Close()
			values, _, err := h.RealDiagonalize()
			if err != nil {
				t.Fatal(err)
			}
			omegaLow, omegaHigh := math.Sqrt(1-coupling), math.Sqrt(1+coupling)
			e0 := (omegaLow + omegaHigh) / 2
			for i, want := range []float64{e0, e0 + omegaLow, e0 + omegaHigh, e0 + 2*omegaLow} {
				if math.Abs(values[i]-want) > 2e-5 {
					t.Errorf("E%d=%g want %g", i, values[i], want)
				}
			}
			result, err := h.LowestStates(context.Background(), 4, EigenOptions{})
			if err != nil {
				t.Fatal(err)
			}
			for i := range result.Values {
				if math.Abs(result.Values[i]-values[i]) > 2e-8 {
					t.Fatalf("iterative/dense mismatch: %v / %v", result.Values, values[:4])
				}
			}
			checkNDEigenpairs(t, h, result, 0)
		})
	}
}

func TestNDimDegenerateThreeDimensionalSpectrum(t *testing.T) {
	g, _ := gridData.NewNDimGridUniform(3, -6, 6, 20)
	h, err := NewNDimFourierHamiltonian(g, []float64{1, 1, 1}, gridData.Harmonic[float64]{ForceConst: 1})
	if err != nil {
		t.Fatal(err)
	}
	defer h.Close()
	if _, err := h.MatrixAt(0); err == nil {
		t.Fatal("unbounded full tensor matrix allocation")
	}
	result, err := h.LowestStates(context.Background(), 4, EigenOptions{Tolerance: 1e-8})
	if err != nil {
		t.Fatal(err)
	}
	for i, want := range []float64{1.5, 2.5, 2.5, 2.5} {
		if math.Abs(result.Values[i]-want) > 2e-5 {
			t.Errorf("E%d=%g want %g", i, result.Values[i], want)
		}
	}
	checkNDEigenpairs(t, h, result, 0)
}

func checkNDEigenpairs(t *testing.T, h *NDimHamiltonian, result *EigenResult, time float64) {
	t.Helper()
	n := h.Dim()
	for s, e := range result.Values {
		x, ax := mat.Col(nil, s, result.Vectors), make([]float64, n)
		h.ApplyHRealAt(x, ax, time)
		for i := range ax {
			ax[i] -= e * x[i]
		}
		res := math.Sqrt(ndDot(ax, ax))
		if res > 1e-7 || math.Abs(res-result.Residuals[s]) > 1e-12 {
			t.Fatal("unverified eigen residual")
		}
		for j := 0; j <= s; j++ {
			dot := ndDot(x, mat.Col(nil, j, result.Vectors))
			want := 0.0
			if j == s {
				want = 1
			}
			if math.Abs(dot-want) > 1e-10 {
				t.Fatal("eigenvectors not orthonormal")
			}
		}
	}
}

func TestNDimOneDimensionAndTranslatedDVR(t *testing.T) {
	for name, build := range ndConstructors {
		t.Run(name, func(t *testing.T) {
			g, _ := gridData.NewNDimGridUniform(1, -6, 6, 24)
			h, err := build(g, []float64{2}, gridData.Harmonic[float64]{ForceConst: 3})
			if err != nil {
				t.Fatal(err)
			}
			defer h.Close()
			oneG, _ := gridData.NewFromLength(12, 24)
			one, err := hamiltonianConstructors[name](oneG, 2, gridData.Harmonic[float64]{ForceConst: 3})
			if err != nil {
				t.Fatal(err)
			}
			defer one.Close()
			m, err := h.MatrixAt(0)
			if err != nil || !mat.EqualApprox(m, one.MatrixAt(0), 1e-12) {
				t.Fatal("1D compatibility failed")
			}
			shiftG, _ := gridData.NewNDimGridUniform(1, 4, 16, 24)
			shift, err := build(shiftG, []float64{2}, gridData.Harmonic[float64]{ForceConst: 3, Cen: 10})
			if err != nil {
				t.Fatal(err)
			}
			defer shift.Close()
			shiftM, _ := shift.MatrixAt(0)
			if !mat.EqualApprox(m, shiftM, 1e-12) {
				t.Fatal("moving Cartesian box changed kinetic boundary condition")
			}
		})
	}
}

func TestNDimSolverErrors(t *testing.T) {
	g, _ := gridData.NewNDimGridUniform(2, -3, 3, 6)
	for _, m := range [][]float64{{1}, {1, -1}, {1, math.NaN()}} {
		if _, err := NewNDimHamiltonian(g, m, nil); err == nil {
			t.Fatal("accepted masses")
		}
	}
	if _, err := NewNDimHamiltonian(nil, nil, nil); err == nil {
		t.Fatal("accepted nil grid")
	}
	if _, err := NewNDimHamiltonian(g, []float64{1, 1}, gridData.NDimPotentialFunc(func([]float64) float64 { return math.Inf(1) })); err == nil {
		t.Fatal("accepted nonfinite potential")
	}
	h, err := NewNDimHamiltonian(g, []float64{1, 1}, gridData.Harmonic[float64]{ForceConst: 1})
	if err != nil {
		t.Fatal(err)
	}
	defer h.Close()
	ctx, cancel := context.WithCancel(context.Background())
	cancel()
	if _, err := h.LowestStates(ctx, 2, EigenOptions{}); err != context.Canceled {
		t.Fatal("cancellation ignored")
	}
	if result, err := h.LowestStates(context.Background(), 2, EigenOptions{MaxIterations: 1}); err == nil || result == nil {
		t.Fatal("unconverged approximation reported as success")
	}
	for _, opts := range []EigenOptions{{Tolerance: -1}, {MaxSubspace: 3}, {Time: math.NaN()}} {
		if _, err := h.LowestStates(context.Background(), 2, opts); err == nil {
			t.Fatal("accepted invalid solver options")
		}
	}
	h.SetTimeDependentFunction(func(x []float64, time float64) float64 { return .5*(x[0]*x[0]+x[1]*x[1]) + time })
	atZero, err := h.LowestStates(context.Background(), 2, EigenOptions{})
	if err != nil {
		t.Fatal(err)
	}
	atTwo, err := h.LowestStates(context.Background(), 2, EigenOptions{Time: 2})
	if err != nil {
		t.Fatal(err)
	}
	for i := range atZero.Values {
		if math.Abs(atTwo.Values[i]-atZero.Values[i]-2) > 1e-8 {
			t.Fatal("incorrect instantaneous spectrum")
		}
	}
	h.Close()
	h.Close()
	if _, err := h.LowestStates(context.Background(), 1, EigenOptions{}); err == nil {
		t.Fatal("accepted closed Hamiltonian")
	}
}

func TestNDimSplitUnitarityAndReversibility(t *testing.T) {
	for name, build := range ndConstructors {
		t.Run(name, func(t *testing.T) {
			g, _ := gridData.NewNDimGrid([]float64{-4, -3}, []float64{4, 3}, []uint32{12, 10})
			h, err := build(g, []float64{1, 2}, nil)
			if err != nil {
				t.Fatal(err)
			}
			defer h.Close()
			h.SetTimeDependentFunction(func(x []float64, time float64) float64 { return .5*(x[0]*x[0]+x[1]*x[1]) + .1*x[0]*x[1]*math.Sin(time) })
			p, err := NewNDimSplitPropagator(h)
			if err != nil {
				t.Fatal(err)
			}
			psi := make([]complex128, h.Dim())
			for i := range psi {
				x := g.GetPoint(uint64(i))
				psi[i] = cmplx.Exp(complex(-.5*(x[0]*x[0]+x[1]*x[1]), x[0]))
			}
			before := slices.Clone(psi)
			for step := 0; step < 20; step++ {
				if err := p.Step(psi, float64(step)*.01, .01); err != nil {
					t.Fatal(err)
				}
			}
			var normBefore, normAfter float64
			for i, z := range psi {
				normAfter += real(cmplx.Conj(z) * z)
				normBefore += real(cmplx.Conj(before[i]) * before[i])
			}
			if math.Abs(normAfter/normBefore-1) > 1e-12 {
				t.Fatal("split step changed norm")
			}
			for step := 20; step > 0; step-- {
				if err := p.Step(psi, float64(step)*.01, -.01); err != nil {
					t.Fatal(err)
				}
			}
			for i, z := range psi {
				if cmplx.Abs(z-before[i]) > 1e-12 {
					t.Fatal("split step not reversible")
				}
			}
			saved := slices.Clone(psi)
			if err := p.Step(psi, 0, math.NaN()); err == nil || !slices.Equal(saved, psi) {
				t.Fatal("invalid step changed state")
			}
		})
	}
}

func TestNDimAnisotropicMassSpectrum(t *testing.T) {
	g, _ := gridData.NewNDimGridFromLengths([]float64{12, 14}, []uint32{24, 28})
	// m=(2,.5), omega=(1,1.7), force constants=(2,1.445).
	h, err := NewNDimHamiltonian(g, []float64{2, .5}, gridData.NDimPotentialFunc(func(x []float64) float64 {
		return x[0]*x[0] + .5*1.445*x[1]*x[1]
	}))
	if err != nil {
		t.Fatal(err)
	}
	defer h.Close()
	result, err := h.LowestStates(context.Background(), 3, EigenOptions{})
	if err != nil {
		t.Fatal(err)
	}
	for i, want := range []float64{1.35, 2.35, 3.05} {
		if math.Abs(result.Values[i]-want) > 1e-5 {
			t.Fatalf("anisotropic energies %v", result.Values)
		}
	}
}

func TestNDimTimeDependentPhase(t *testing.T) {
	g, _ := gridData.NewNDimGridUniform(2, -math.Pi, math.Pi, 9)
	h, err := NewNDimFourierHamiltonian(g, []float64{1, 2}, nil)
	if err != nil {
		t.Fatal(err)
	}
	defer h.Close()
	h.SetTimeDependentFunction(func(_ []float64, time float64) float64 { return time })
	p, err := NewNDimSplitPropagator(h)
	if err != nil {
		t.Fatal(err)
	}
	psi := make([]complex128, h.Dim())
	for i := range psi {
		x := g.GetPoint(uint64(i))
		psi[i] = cmplx.Exp(complex(0, x[0]+2*x[1]))
	}
	initial := slices.Clone(psi)
	for step := 0; step < 10; step++ {
		if err := p.Step(psi, float64(step)*.1, .1); err != nil {
			t.Fatal(err)
		}
	}
	// Ekin=1/2+4/4=1.5; integral_0^1 V(t)dt=0.5.
	for i := range psi {
		if cmplx.Abs(psi[i]-initial[i]*cmplx.Exp(complex(0, -2))) > 1e-12 {
			t.Fatal("time-dependent midpoint phase is incorrect")
		}
	}
}

func TestNDimSplitSecondOrder(t *testing.T) {
	for name, build := range ndConstructors {
		t.Run(name, func(t *testing.T) {
			g, _ := gridData.NewNDimGridUniform(2, -4, 4, 10)
			h, err := build(g, []float64{1, 2}, gridData.NDimPotentialFunc(func(x []float64) float64 { return .5*x[0]*x[0] + x[1]*x[1] + .3*x[0]*x[1] }))
			if err != nil {
				t.Fatal(err)
			}
			defer h.Close()
			values, vectors, err := h.RealDiagonalize()
			if err != nil {
				t.Fatal(err)
			}
			initial, exact := make([]complex128, h.Dim()), make([]complex128, h.Dim())
			for i := range initial {
				x := g.GetPoint(uint64(i))
				initial[i] = cmplx.Exp(complex(-.5*((x[0]-.5)*(x[0]-.5)+x[1]*x[1]), .7*x[0]))
			}
			for j, e := range values {
				var c complex128
				for i, z := range initial {
					c += complex(vectors.At(i, j), 0) * z
				}
				c *= cmplx.Exp(complex(0, -.4*e))
				for i := range exact {
					exact[i] += complex(vectors.At(i, j), 0) * c
				}
			}
			errors := make([]float64, 2)
			for run, steps := range []int{10, 20} {
				p, err := NewNDimSplitPropagator(h)
				if err != nil {
					t.Fatal(err)
				}
				psi := slices.Clone(initial)
				dt := .4 / float64(steps)
				for step := 0; step < steps; step++ {
					if err := p.Step(psi, float64(step)*dt, dt); err != nil {
						t.Fatal(err)
					}
				}
				for i, z := range psi {
					diff := cmplx.Abs(z - exact[i])
					errors[run] += diff * diff
				}
				errors[run] = math.Sqrt(errors[run])
			}
			if ratio := errors[0] / errors[1]; ratio < 3.8 || ratio > 4.2 {
				t.Fatalf("split method not second order: errors %v", errors)
			}
		})
	}
}
