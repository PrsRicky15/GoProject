package EquationSolver

import (
	"context"
	"errors"
	"math"
	"math/cmplx"
	"slices"
	"testing"
)

func TestVectorOnlyCoupledRHS(t *testing.T) {
	rhs := func(t float64, y, dy []float64) error { dy[0] = y[1]; dy[1] = -y[0]; return nil }
	s, err := NewVectorImplicit(rhs, ImplicitMidpointMethod, ImplicitOptions{})
	if err != nil {
		t.Fatal(err)
	}
	y := []float64{1, 0}
	for i := 0; i < 100; i++ {
		if err = s.Step(y, float64(i)*.01, .01); err != nil {
			t.Fatal(err)
		}
	}
	if math.Hypot(y[0]-math.Cos(1), y[1]+math.Sin(1)) > 1e-5 {
		t.Fatal("whole-vector rotation", y)
	}
	sentinel := errors.New("model failure")
	s, _ = NewVectorImplicit(func(t float64, y, dy []float64) error { return sentinel }, ImplicitEulerMethod, ImplicitOptions{})
	before := slices.Clone(y)
	if err = s.Step(y, 0, .1); !errors.Is(err, sentinel) || !slices.Equal(y, before) {
		t.Fatal("callback failure contract", err)
	}
}

func TestComplexNonlinearRealBlock(t *testing.T) {
	// Conjugation is not complex analytic. Re and Im evolve with opposite signs.
	s, err := NewComplexImplicit(func(t float64, y, dy []complex128) error {
		for i, v := range y {
			dy[i] = cmplx.Conj(v)
		}
		return nil
	}, ImplicitEulerMethod, ImplicitOptions{})
	if err != nil {
		t.Fatal(err)
	}
	y := []complex128{1 + 2i, 3 - 4i}
	if err = s.Step(y, 0, .1); err != nil {
		t.Fatal(err)
	}
	if cmplx.Abs(y[0]-complex(1/.9, 2/1.1)) > 1e-8 || cmplx.Abs(y[1]-complex(3/.9, -4/1.1)) > 1e-8 {
		t.Fatal("real block coupling", y)
	}
	// A genuinely nonlinear, globally coupled RHS: density of the whole vector.
	s, _ = NewComplexImplicit(func(t float64, y, dy []complex128) error {
		rho := 0.
		for _, v := range y {
			rho += cmplx.Abs(v) * cmplx.Abs(v)
		}
		for i, v := range y {
			dy[i] = -complex(rho, 0) * v
		}
		return nil
	}, ImplicitEulerMethod, ImplicitOptions{})
	y = []complex128{1, 1i}
	before := slices.Clone(y)
	if err = s.Step(y, 0, .2); err != nil {
		t.Fatal(err)
	}
	rho := cmplx.Abs(y[0])*cmplx.Abs(y[0]) + cmplx.Abs(y[1])*cmplx.Abs(y[1])
	for i := range y {
		if cmplx.Abs(y[i]-before[i]+complex(.2*rho, 0)*y[i]) > 1e-8 {
			t.Fatal("nonlinear vector residual")
		}
	}
}

func complexMethodFactor(method ImplicitMethod, z complex128) complex128 {
	switch method {
	case ImplicitEulerMethod:
		return 1 / (1 - z)
	case ImplicitSDIRK2Method:
		g := complex(1-1/math.Sqrt(2), 0)
		return (1 + (1-2*g)*z) / ((1 - g*z) * (1 - g*z))
	default:
		return (1 + z/2) / (1 - z/2)
	}
}

func TestComplexLinearCouplingAndMethods(t *testing.T) {
	// H = sigma_x; its eigenvectors (1,1) and (1,-1) have energies +/-1.
	rhs := func(t float64, y, dy []complex128) error { dy[0] = -1i * y[1]; dy[1] = -1i * y[0]; return nil }
	for _, method := range []ImplicitMethod{ImplicitEulerMethod, ImplicitTrapezoidalMethod, ImplicitMidpointMethod, ImplicitSDIRK2Method} {
		s, err := NewComplexLinearImplicit(2, rhs, method, KrylovOptions{})
		if err != nil {
			t.Fatal(err)
		}
		y := []complex128{1, 0}
		if err = s.Step(y, 0, .2); err != nil {
			t.Fatal(err)
		}
		a, b := complexMethodFactor(method, -.2i), complexMethodFactor(method, .2i)
		want := []complex128{(a + b) / 2, (a - b) / 2}
		for i := range y {
			if cmplx.Abs(y[i]-want[i]) > 1e-11 {
				t.Fatalf("method %d: got %v want %v", method, y, want)
			}
		}
		if method == ImplicitMidpointMethod && math.Abs(complexNorm(y)-1) > 1e-12 {
			t.Fatal("midpoint norm")
		}
	}
}

func TestComplexLinearTimeDependentOrder(t *testing.T) {
	// Noncommuting H(t)=cos(t)*sigma_x+sin(t)*sigma_y. Exact evolution from
	// rotating-frame H_eff=sigma_x-sigma_z/2, with U=exp(-i*t*sigma_z/2).
	rhs := func(t float64, y, dy []complex128) error {
		dy[0] = -1i * cmplx.Exp(complex(0, -t)) * y[1]
		dy[1] = -1i * cmplx.Exp(complex(0, t)) * y[0]
		return nil
	}
	a := math.Sqrt(1.25)
	q := math.Sin(a) / a
	want := []complex128{cmplx.Exp(-.5i) * complex(math.Cos(a), .5*q), cmplx.Exp(.5i) * complex(0, -q)}
	for _, method := range []ImplicitMethod{ImplicitEulerMethod, ImplicitTrapezoidalMethod, ImplicitMidpointMethod, ImplicitSDIRK2Method} {
		var errs []float64
		for _, steps := range []int{20, 40} {
			s, _ := NewComplexLinearImplicit(2, rhs, method, KrylovOptions{})
			y := []complex128{1, 0}
			dt := 1 / float64(steps)
			for i := 0; i < steps; i++ {
				if err := s.Step(y, float64(i)*dt, dt); err != nil {
					t.Fatal(err)
				}
			}
			errs = append(errs, math.Hypot(cmplx.Abs(y[0]-want[0]), cmplx.Abs(y[1]-want[1])))
			if method == ImplicitMidpointMethod && math.Abs(complexNorm(y)-1) > 1e-10 {
				t.Fatal("time-dependent norm")
			}
		}
		threshold := 3.7
		if method == ImplicitEulerMethod {
			threshold = 1.8
		}
		if errs[0]/errs[1] < threshold {
			t.Fatalf("method %d errors %v", method, errs)
		}
	}
}

func TestGMRESRestartAndLargeVector(t *testing.T) {
	const n = 2048
	// A nonnormal, complex, upper-bidiagonal equation; every row reads a neighbor.
	rhs := func(t float64, y, dy []complex128) error {
		for i := range y {
			dy[i] = complex(-2, float64(i%7)/10) * y[i]
			if i+1 < len(y) {
				dy[i] += (.4 + .2i) * y[i+1]
			}
		}
		return nil
	}
	s, err := NewComplexLinearImplicit(n, rhs, ImplicitEulerMethod, KrylovOptions{Restart: 3, MaxIterations: 100})
	if err != nil {
		t.Fatal(err)
	}
	y := make([]complex128, n)
	for i := range y {
		y[i] = complex(1, float64(i%3))
	}
	before := slices.Clone(y)
	if err = s.Step(y, 0, .5); err != nil {
		t.Fatal(err)
	}
	out := make([]complex128, n)
	rhs(.5, y, out)
	for i := range y {
		out[i] = y[i] - .5*out[i] - before[i]
	}
	if complexNorm(out) > 1e-8 || s.Stats().Restarts == 0 {
		t.Fatalf("restart residual %g stats %+v", complexNorm(out), s.Stats())
	}
	if allocs := testing.AllocsPerRun(10, func() {
		copy(y, before)
		if err := s.Step(y, 0, .5); err != nil {
			panic(err)
		}
	}); allocs != 0 {
		t.Fatalf("warm propagation allocates %g", allocs)
	}
}

func TestComplexVectorFailures(t *testing.T) {
	sentinel := errors.New("operator unavailable")
	s, _ := NewComplexLinearImplicit(2, func(t float64, y, dy []complex128) error { return sentinel }, ImplicitMidpointMethod, KrylovOptions{})
	y := []complex128{1, 1i}
	before := slices.Clone(y)
	if err := s.Step(y, 0, .1); !errors.Is(err, sentinel) || !slices.Equal(y, before) {
		t.Fatal("operator failure", err)
	}
	ctx, cancel := context.WithCancel(context.Background())
	cancel()
	if err := s.StepContext(ctx, y, 0, .1); !errors.Is(err, context.Canceled) {
		t.Fatal("cancellation", err)
	}
	rhs := func(t float64, y, dy []complex128) error { dy[0] = -1i * y[1]; dy[1] = -1i * y[0]; return nil }
	s, _ = NewComplexLinearImplicit(2, rhs, ImplicitEulerMethod, KrylovOptions{Restart: 1, MaxIterations: 1})
	y = []complex128{1, 0}
	before = slices.Clone(y)
	if err := s.Step(y, 0, 1); err == nil || !slices.Equal(y, before) {
		t.Fatal("unconverged result committed")
	}
	s, _ = NewComplexLinearImplicit(2, func(t float64, y, dy []complex128) error {
		if t > .5 {
			return sentinel
		}
		copy(dy, y)
		return nil
	}, ImplicitSDIRK2Method, KrylovOptions{})
	if err := s.Step(y, 0, 1); !errors.Is(err, sentinel) || !slices.Equal(y, before) {
		t.Fatal("second stage failure", err)
	}
	if _, err := NewComplexLinearImplicit(1<<20, rhs, ImplicitMidpointMethod, KrylovOptions{Restart: 256}); err == nil {
		t.Fatal("workspace limit")
	}
	if err := s.Step(y, 0, math.NaN()); err == nil {
		t.Fatal("non-finite dt")
	}
}
