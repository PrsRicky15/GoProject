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

func TestImplicitHamiltonianStaticPhaseAndNorm(t *testing.T) {
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
			p, err := NewImplicitPropagator(h, EquationSolver.KrylovOptions{})
			if err != nil {
				t.Fatal(err)
			}
			psi := make([]complex128, h.Dim())
			for i := range psi {
				psi[i] = complex(vectors.At(i, 0), vectors.At(i, 1)) / complex(math.Sqrt(2), 0)
			}
			for i := 0; i < 20; i++ {
				if err = p.Step(psi, float64(i)*.01, .01); err != nil {
					t.Fatal(err)
				}
			}
			errorNorm, norm := 0., 0.
			for i, v := range psi {
				want := (complex(vectors.At(i, 0), 0)*cmplx.Exp(complex(0, -.2*energies[0])) + complex(0, vectors.At(i, 1))*cmplx.Exp(complex(0, -.2*energies[1]))) / complex(math.Sqrt(2), 0)
				errorNorm += cmplx.Abs(v-want) * cmplx.Abs(v-want)
				norm += cmplx.Abs(v) * cmplx.Abs(v)
			}
			if math.Sqrt(errorNorm) > 5e-6 || math.Abs(norm-1) > 1e-9 {
				t.Fatalf("phase error %g norm %g", math.Sqrt(errorNorm), norm)
			}
			before := slices.Clone(psi)
			h.Close()
			if err = p.Step(psi, .2, .01); err == nil || !slices.Equal(psi, before) {
				t.Fatal("closed Hamiltonian should fail atomically")
			}
		})
	}
}

func TestImplicitNDimensionalTimeDependentPropagation(t *testing.T) {
	// 2048 entries exceeds the nonlinear dense real-block solver limit.
	g, err := gridData.NewNDimGridFromLengths([]float64{2 * math.Pi, 2 * math.Pi, 2 * math.Pi, 2 * math.Pi}, []uint32{8, 8, 8, 4})
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
	p, err := NewNDimImplicitPropagator(h, EquationSolver.KrylovOptions{Restart: 8})
	if err != nil {
		t.Fatal(err)
	}
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
	var errs []float64
	for _, steps := range []int{10, 20} {
		psi := slices.Clone(initial)
		dt := .2 / float64(steps)
		for i := 0; i < steps; i++ {
			if err = p.Step(psi, float64(i)*dt, dt); err != nil {
				t.Fatal(err)
			}
		}
		phase := cmplx.Exp(complex(0, -(2*.2 + .2*.2/2)))
		norm, diff := 0., 0.
		for i, v := range psi {
			norm += cmplx.Abs(v) * cmplx.Abs(v)
			e := cmplx.Abs(v - phase*initial[i])
			diff += e * e
		}
		if math.Abs(norm-1) > 1e-9 {
			t.Fatalf("norm drift %g", norm-1)
		}
		errs = append(errs, math.Sqrt(diff))
	}
	if errs[0]/errs[1] < 3.9 {
		t.Fatalf("time-dependent order errors %v", errs)
	}
	ctx, cancel := context.WithCancel(context.Background())
	cancel()
	psi := slices.Clone(initial)
	if err = p.StepContext(ctx, psi, 0, .1); !errors.Is(err, context.Canceled) || !slices.Equal(psi, initial) {
		t.Fatal("cancellation contract")
	}
}
