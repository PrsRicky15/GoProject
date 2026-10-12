package Quantum

import (
	EquationSolver "GoProject/ODESolver"
	"GoProject/gridData"
	"context"
	"errors"
	"math"
	"math/cmplx"
	"slices"
	"testing"
)

func TestRKHamiltonianStaticSuperposition(t *testing.T) {
	methods := []EquationSolver.RKMethod{EquationSolver.RKEuler, EquationSolver.RKHeun2, EquationSolver.RKMidpoint2, EquationSolver.RKRalston2, EquationSolver.RKRalston3, EquationSolver.RKHeun3, EquationSolver.RKWray3, EquationSolver.RKSSP3, EquationSolver.RKKutta3, EquationSolver.RKClassic4, EquationSolver.RKThreeEighths4, EquationSolver.RKRalston4, EquationSolver.RKNystrom5, EquationSolver.RKHeunEuler21, EquationSolver.RKFehlberg21, EquationSolver.RKBogackiShampine32, EquationSolver.RKFehlberg54, EquationSolver.RKCashKarp54, EquationSolver.RKDormandPrince54}
	for name, build := range hamiltonianConstructors {
		t.Run(name, func(t *testing.T) {
			g, err := gridData.NewFromLength(12, 24)
			if err != nil {
				t.Fatal(err)
			}
			h, err := build(g, 1, gridData.Harmonic[float64]{ForceConst: 1})
			if err != nil {
				t.Fatal(err)
			}
			defer h.Close()
			energies, vectors, err := h.RealDiagonalize()
			if err != nil {
				t.Fatal(err)
			}
			initial := make([]complex128, h.Dim())
			want := make([]complex128, h.Dim())
			for i := range initial {
				initial[i] = complex(vectors.At(i, 0), vectors.At(i, 1)) / complex(math.Sqrt2, 0)
				want[i] = (complex(vectors.At(i, 0), 0)*cmplx.Exp(complex(0, -.1*energies[0])) + complex(0, vectors.At(i, 1))*cmplx.Exp(complex(0, -.1*energies[1]))) / complex(math.Sqrt2, 0)
			}
			for _, m := range methods {
				t.Run(string(m), func(t *testing.T) {
					p, err := NewRKPropagator(h, m)
					if err != nil {
						t.Fatal(err)
					}
					psi := slices.Clone(initial)
					for j := 0; j < 100; j++ {
						if err = p.Step(psi, float64(j)*.001, .001); err != nil {
							t.Fatal(err)
						}
					}
					diff := 0.
					for i, v := range psi {
						d := cmplx.Abs(v - want[i])
						diff += d * d
					}
					limit := 2e-7
					if m == EquationSolver.RKEuler {
						limit = 1e-4
					}
					if math.Sqrt(diff) > limit {
						t.Fatalf("phase/amplitude error %g", math.Sqrt(diff))
					}
				})
			}
			p, err := NewAdaptiveRKPropagator(h, EquationSolver.RKDormandPrince54, EquationSolver.AdaptiveOptions{AbsoluteTolerance: 1e-12, RelativeTolerance: 1e-10})
			if err != nil {
				t.Fatal(err)
			}
			psi := slices.Clone(initial)
			r, err := p.Integrate(psi, 0, .1)
			if err != nil || r.Time != .1 {
				t.Fatal(r, err)
			}
			for i := range psi {
				if cmplx.Abs(psi[i]-want[i]) > 1e-9 {
					t.Fatal("adaptive spectral reference mismatch")
				}
			}
			bad := []complex128{1}
			if _, err = p.Integrate(bad, 0, .1); err == nil || bad[0] != 1 {
				t.Fatal("dimension guard")
			}
			fixed, _ := NewRKPropagator(h, EquationSolver.RKClassic4)
			before := slices.Clone(psi)
			h.Close()
			if err = fixed.Step(psi, 0, .1); err == nil || !slices.Equal(psi, before) {
				t.Fatal("closed fixed adapter")
			}
			if _, err = p.Integrate(psi, 0, .1); err == nil || !slices.Equal(psi, before) {
				t.Fatal("closed adaptive adapter")
			}
		})
	}
}

func TestRKNDimensionalTimeDependentHamiltonian(t *testing.T) {
	g, err := gridData.NewNDimGridFromLengths([]float64{2 * math.Pi, 2 * math.Pi, 2 * math.Pi, 2 * math.Pi}, []uint32{4, 4, 4, 4})
	if err != nil {
		t.Fatal(err)
	}
	h, err := NewNDimFourierHamiltonian(g, []float64{1, 1, 1, 1}, nil)
	if err != nil {
		t.Fatal(err)
	}
	defer h.Close()
	h.SetTimeDependentPotential(func(t float64, out []float64) {
		for i := range out {
			out[i] = t
		}
	})
	initial := make([]complex128, h.Dim())
	x := make([]float64, 4)
	for i := range initial {
		g.PointTo(uint64(i), x)
		phase := 0.
		for _, v := range x {
			phase += v
		}
		initial[i] = cmplx.Exp(complex(0, phase)) / complex(math.Sqrt(float64(h.Dim())), 0)
	}
	check := func(psi []complex128) {
		t.Helper()
		phase := cmplx.Exp(complex(0, -(2*.2 + .2*.2/2)))
		for i := range psi {
			if cmplx.Abs(psi[i]-initial[i]*phase) > 2e-7 {
				t.Fatal("time-dependent 4D propagation mismatch")
			}
		}
	}
	p, err := NewNDimRKPropagator(h, EquationSolver.RKClassic4)
	if err != nil {
		t.Fatal(err)
	}
	psi := slices.Clone(initial)
	for j := 0; j < 20; j++ {
		if err = p.Step(psi, float64(j)*.01, .01); err != nil {
			t.Fatal(err)
		}
	}
	check(psi)
	for _, m := range []EquationSolver.RKMethod{EquationSolver.RKHeunEuler21, EquationSolver.RKFehlberg21, EquationSolver.RKBogackiShampine32, EquationSolver.RKFehlberg54, EquationSolver.RKCashKarp54, EquationSolver.RKDormandPrince54} {
		t.Run(string(m), func(t *testing.T) {
			a, err := NewNDimAdaptiveRKPropagator(h, m, EquationSolver.AdaptiveOptions{AbsoluteTolerance: 1e-12, RelativeTolerance: 1e-10, InitialStep: .2})
			if err != nil {
				t.Fatal(err)
			}
			psi := slices.Clone(initial)
			r, err := a.Integrate(psi, 0, .2)
			if err != nil || r.Time != .2 {
				t.Fatal(r, err)
			}
			check(psi)
			ctx, cancel := context.WithCancel(context.Background())
			cancel()
			before := slices.Clone(psi)
			if _, err = a.IntegrateContext(ctx, psi, .2, .3); !errors.Is(err, context.Canceled) || !slices.Equal(psi, before) {
				t.Fatal("canceled adapter")
			}
		})
	}
}
