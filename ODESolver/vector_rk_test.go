package EquationSolver

import (
	"context"
	"errors"
	"math"
	"math/cmplx"
	"slices"
	"testing"
)

var allRKMethods = []RKMethod{RKEuler, RKHeun2, RKMidpoint2, RKRalston2, RKRalston3, RKHeun3, RKWray3, RKSSP3, RKKutta3, RKClassic4, RKThreeEighths4, RKRalston4, RKNystrom5, RKHeunEuler21, RKFehlberg21, RKBogackiShampine32, RKFehlberg54, RKCashKarp54, RKDormandPrince54}
var adaptiveRKMethods = []RKMethod{RKHeunEuler21, RKFehlberg21, RKBogackiShampine32, RKFehlberg54, RKCashKarp54, RKDormandPrince54}

// Nonlinear, coupled and non-autonomous; exact solution (exp(t),exp(-t)).
func manufacturedRK(t float64, y, dy []float64) error {
	dy[0] = y[0]*y[1] + math.Exp(t) - 1
	dy[1] = y[0]*y[0] + y[1] - math.Exp(2*t) - 2*math.Exp(-t)
	return nil
}

// A nonlinear complex equation involving the entire vector and conjugates.
// On the unit sphere it has exact solution (cos(theta),i*sin(theta)).
func complexCoupledRK(t float64, y, dy []complex128) error {
	w := 1 + t + .1*(real(y[0]*cmplx.Conj(y[0])+y[1]*cmplx.Conj(y[1]))-1)
	dy[0] = complex(0, w) * y[1]
	dy[1] = complex(0, w) * y[0]
	return nil
}
func TestRKAllTableauxAndObservedOrders(t *testing.T) {
	for _, m := range allRKMethods {
		t.Run(string(m), func(t *testing.T) {
			q, _ := rkTable(m)
			for i, row := range q.a {
				sum := 0.
				for _, v := range row {
					sum += v
				}
				if len(row) > i || math.Abs(sum-q.c[i]) > 2e-14 {
					t.Fatalf("inconsistent stage %d: %g vs %g", i, sum, q.c[i])
				}
			}
			for _, weights := range [][]float64{q.b, q.low} {
				if weights == nil {
					continue
				}
				sum := 0.
				for _, v := range weights {
					sum += v
				}
				if math.Abs(sum-1) > 2e-14 {
					t.Fatalf("weights sum %g", sum)
				}
			}
			for _, complexState := range []bool{false, true} {
				var errs [2]float64
				for level, n := range []int{40, 80} {
					h := 1. / float64(n)
					if complexState {
						s, err := NewComplexRK(complexCoupledRK, m)
						if err != nil {
							t.Fatal(err)
						}
						y := []complex128{1, 0}
						for j := 0; j < n; j++ {
							if err = s.Step(y, float64(j)*h, h); err != nil {
								t.Fatal(err)
							}
						}
						errs[level] = math.Hypot(cmplx.Abs(y[0]-complex(math.Cos(1.5), 0)), cmplx.Abs(y[1]-complex(0, math.Sin(1.5))))
					} else {
						s, err := NewVectorRK(manufacturedRK, m)
						if err != nil {
							t.Fatal(err)
						}
						y := []float64{1, 1}
						for j := 0; j < n; j++ {
							if err = s.Step(y, float64(j)*h, h); err != nil {
								t.Fatal(err)
							}
						}
						errs[level] = math.Hypot(y[0]-math.E, y[1]-1/math.E)
					}
				}
				ratio := errs[0] / errs[1]
				want := math.Pow(2, float64(q.order))
				if ratio < .65*want || ratio > 1.6*want {
					t.Fatalf("complex=%v order=%d errors=%v ratio=%g", complexState, q.order, errs, ratio)
				}
			}
		})
	}
}

func TestAllAdaptiveRKRealComplexAndBackward(t *testing.T) {
	for _, m := range adaptiveRKMethods {
		t.Run(string(m), func(t *testing.T) {
			var errs [2]float64
			var calls [2]int
			for level, tol := range []float64{1e-7, 1e-10} {
				s, err := NewVectorAdaptiveRK(manufacturedRK, m, AdaptiveOptions{AbsoluteTolerance: tol, RelativeTolerance: tol, InitialStep: 1})
				if err != nil {
					t.Fatal(err)
				}
				y := []float64{1, 1}
				r, err := s.Integrate(y, 0, 1)
				if err != nil {
					t.Fatal(err)
				}
				if r.Time != 1 || r.AcceptedSteps == 0 || r.RejectedSteps == 0 {
					t.Fatalf("missing endpoint or rejection: %+v", r)
				}
				errs[level] = math.Hypot(y[0]-math.E, y[1]-1/math.E)
				calls[level] = r.RHSEvaluations
			}
			if errs[1] > 5e-6 || errs[1] >= errs[0]/5 || calls[1] <= calls[0] {
				t.Fatalf("accuracy/control errors %v evaluations %v", errs, calls)
			}
			for _, end := range []float64{.7, -.7} {
				s, err := NewComplexAdaptiveRK(complexCoupledRK, m, AdaptiveOptions{AbsoluteTolerance: 1e-11, RelativeTolerance: 1e-10, InitialStep: 1})
				if err != nil {
					t.Fatal(err)
				}
				y := []complex128{1, 0}
				r, err := s.Integrate(y, 0, end)
				if err != nil {
					t.Fatal(err)
				}
				theta := end + end*end/2
				if r.Time != end || r.NextStep*end <= 0 || cmplx.Abs(y[0]-complex(math.Cos(theta), 0))+cmplx.Abs(y[1]-complex(0, math.Sin(theta))) > 2e-6 {
					t.Fatalf("end %g state %v result %+v", end, y, r)
				}
			}
		})
	}
}

func TestRKFailureAtomicity(t *testing.T) {
	boom := errors.New("RHS failure")
	for _, m := range allRKMethods {
		t.Run(string(m), func(t *testing.T) {
			s, _ := NewComplexRK(func(t float64, y, dy []complex128) error { y[0] = 123; return boom }, m)
			y := []complex128{1, 2i}
			before := slices.Clone(y)
			if err := s.Step(y, 0, .1); !errors.Is(err, boom) || !slices.Equal(y, before) {
				t.Fatal("RHS failure mutated state", err, y)
			}
		})
	}
	for _, m := range adaptiveRKMethods {
		t.Run(string(m), func(t *testing.T) {
			s, _ := NewVectorAdaptiveRK(func(t float64, y, dy []float64) error {
				if t > .25 {
					return boom
				}
				dy[0] = 1
				return nil
			}, m, AdaptiveOptions{MaxStep: .1})
			y := []float64{2}
			r, err := s.Integrate(y, 0, 1)
			if !errors.Is(err, boom) || y[0] != 2 || r.Time != 0 || r.AcceptedSteps == 0 {
				t.Fatalf("partial interval committed: %v %+v %v", y, r, err)
			}
		})
	}
	for _, o := range []AdaptiveOptions{{MaxAttempts: 1, InitialStep: 1, AbsoluteTolerance: 1e-14, RelativeTolerance: 1e-14}, {MinStep: 1, InitialStep: 1, AbsoluteTolerance: 1e-14, RelativeTolerance: 1e-14}} {
		s, _ := NewVectorAdaptiveRK(manufacturedRK, RKDormandPrince54, o)
		y := []float64{1, 1}
		r, err := s.Integrate(y, 0, 1)
		if err == nil || !slices.Equal(y, []float64{1, 1}) || r.Time != 0 {
			t.Fatal("limit should fail atomically")
		}
	}
	ctx, cancel := context.WithCancel(context.Background())
	s, _ := NewComplexAdaptiveRK(func(t float64, y, dy []complex128) error { dy[0] = 1; cancel(); return nil }, RKHeunEuler21, AdaptiveOptions{})
	y := []complex128{2}
	_, err := s.IntegrateContext(ctx, y, 0, 1)
	if !errors.Is(err, context.Canceled) || y[0] != 2 {
		t.Fatal("cancellation must roll back", err)
	}
	bad, _ := NewVectorRK(func(t float64, y, dy []float64) error { dy[0] = math.NaN(); return nil }, RKEuler)
	z := []float64{1}
	if err = bad.Step(z, 0, 1); err == nil || z[0] != 1 {
		t.Fatal("nonfinite stage must fail")
	}
}

func TestRKValidationAndNoAllocationAfterWarmup(t *testing.T) {
	if _, err := NewVectorRK(nil, RKEuler); err == nil {
		t.Fatal("nil RHS")
	}
	if _, err := NewVectorRK(manufacturedRK, "unknown"); err == nil {
		t.Fatal("unknown method")
	}
	if _, err := NewVectorAdaptiveRK(manufacturedRK, RKClassic4, AdaptiveOptions{}); err == nil {
		t.Fatal("nonembedded method")
	}
	for _, o := range []AdaptiveOptions{{RelativeTolerance: -1}, {AbsoluteTolerance: math.NaN()}, {MinStep: 2, MaxStep: 1}, {MaxAttempts: -1}} {
		if _, err := NewVectorAdaptiveRK(manufacturedRK, RKHeunEuler21, o); err == nil {
			t.Fatal("invalid options accepted")
		}
	}
	s, _ := NewVectorRK(manufacturedRK, RKClassic4)
	for _, tt := range [][2]float64{{math.NaN(), 1}, {math.MaxFloat64, math.MaxFloat64}, {1e20, 1}} {
		if err := s.Step([]float64{1, 1}, tt[0], tt[1]); err == nil {
			t.Fatal("invalid time accepted", tt)
		}
	}
	for _, m := range adaptiveRKMethods {
		calls := 0
		rhs := func(t float64, y, dy []complex128) error {
			calls++
			for i := range dy {
				dy[i] = 1
			}
			return nil
		}
		a, _ := NewComplexAdaptiveRK(rhs, m, AdaptiveOptions{InitialStep: .125, MaxStep: .125})
		y := make([]complex128, 2048)
		r, err := a.Integrate(y, 0, 1)
		if err != nil {
			t.Fatal(err)
		}
		q, _ := rkTable(m)
		want := 8 * len(q.b)
		if q.fsal {
			want -= 7
		}
		if calls != want || r.RHSEvaluations != want || r.AcceptedSteps != 8 || cmplx.Abs(y[0]-1) > 1e-14 {
			t.Fatalf("%s calls %d want %d result %+v", m, calls, want, r)
		}
		if allocs := testing.AllocsPerRun(20, func() {
			clear(y)
			if _, err := a.Integrate(y, 0, 1); err != nil {
				panic(err)
			}
		}); allocs != 0 {
			t.Fatalf("%s adaptive allocations %g", m, allocs)
		}
		f, _ := NewComplexRK(rhs, m)
		_ = f.Step(y, 0, .1)
		if allocs := testing.AllocsPerRun(20, func() {
			if err := f.Step(y, 0, .1); err != nil {
				panic(err)
			}
		}); allocs != 0 {
			t.Fatalf("%s fixed allocations %g", m, allocs)
		}
	}
}

func TestRKEmbeddedLowerOrders(t *testing.T) {
	for _, m := range adaptiveRKMethods {
		t.Run(string(m), func(t *testing.T) {
			var errs [2]float64
			order := 0
			steps := []int{40, 80}
			// Fehlberg's low-order weights are very close to the order-two
			// weights, so resolve its small first-order term with finer steps.
			if m == RKFehlberg21 {
				steps = []int{2048, 4096}
			}
			for level, n := range steps {
				s, _ := NewVectorRK(func(t float64, y, dy []float64) error { dy[0] = y[1]; dy[1] = y[0]; return nil }, m)
				s.q.b = s.q.low
				order = s.q.lowOrder
				y := []float64{2, 0}
				h := 1. / float64(n)
				for i := 0; i < n; i++ {
					if err := s.Step(y, float64(i)*h, h); err != nil {
						t.Fatal(err)
					}
				}
				errs[level] = math.Hypot(y[0]-(math.E+1/math.E), y[1]-(math.E-1/math.E))
			}
			ratio := errs[0] / errs[1]
			want := math.Pow(2, float64(order))
			if ratio < .65*want || ratio > 1.6*want {
				t.Fatalf("lower order %d errors %v ratio %g", order, errs, ratio)
			}
		})
	}
}

func TestCompletedLegacyRKGridMethods(t *testing.T) {
	bridge := &vectorRHSBridge{rhs: manufacturedRK}
	methods := []struct {
		method RKMethod
		old    ODESolver
	}{
		{RKHeun3, &Huens3Explicit{timeFunc: bridge, delTime: .01}},
		{RKWray3, &VDHouwenExplicit{timeFunc: bridge, delTime: .01}},
		{RKSSP3, &SSPRungeKutta3{timeFunc: bridge, delTime: .01}},
		{RKKutta3, &RungeKutta3Explicit{timeFunc: bridge, delTime: .01}},
		{RKThreeEighths4, &RungeKutta38{timeFunc: bridge, deltaTime: .01}},
		{RKRalston4, &Ralston4Order{timeFunc: bridge, deltaTime: .01}},
	}
	for _, m := range methods {
		y, z := []float64{1, 1}, []float64{1, 1}
		s, _ := NewVectorRK(manufacturedRK, m.method)
		for i := 0; i < 10; i++ {
			if err := m.old.NextStepOnGrid(y, float64(i)*.01); err != nil {
				t.Fatal(err)
			}
			if err := s.Step(z, float64(i)*.01, .01); err != nil {
				t.Fatal(err)
			}
		}
		for i := range y {
			if math.Abs(y[i]-z[i]) > 1e-14 {
				t.Fatalf("%s old/new disagreement", m.method)
			}
		}
	}
}
