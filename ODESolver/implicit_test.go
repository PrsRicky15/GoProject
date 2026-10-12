package EquationSolver

import (
	"errors"
	"math"
	"slices"
	"testing"

	"GoProject/gridData"
)

type implicitRHS struct {
	scalar func(float64, float64) float64
	grid   func([]float64, []float64, float64)
}

func (f implicitRHS) EvaluateAt(x, t float64) float64 { return f.scalar(x, t) }
func (f implicitRHS) EvaluateOnRGrid(x []float64, t float64) []float64 {
	out := make([]float64, len(x))
	f.EvaluateOnRGridInPlace(x, out, t)
	return out
}
func (f implicitRHS) EvaluateOnRGridInPlace(x, out []float64, t float64) {
	if f.grid != nil {
		f.grid(x, out, t)
		return
	}
	for i, v := range x {
		out[i] = f.scalar(v, t)
	}
}

type implicitTestSolver interface {
	ODESolver
	Stats() ImplicitStats
}

func implicitForTest(name string, dt float64, f gridData.TDPotentialOp, o ImplicitOptions) implicitTestSolver {
	switch name {
	case "EulerFP":
		return &EulerImplicitFixPoint{TdFunc: f, DeltaT: dt, Options: o}
	case "HeunFP":
		return &HeunsImplicitFixPoint{TdFunc: f, DeltaT: dt, Options: o}
	case "MidpointFP":
		return &MidPointFixPoint{TdFunc: f, DeltaT: dt, Options: o}
	case "PredictorCorrector":
		return &PredictorCorrector{TdFunc: f, DeltaT: dt, Options: o}
	case "EulerNewton":
		return &EulerImplicitNewton{TdFunc: f, DeltaT: dt, Options: o}
	case "HeunNewton":
		return &HeunsImplicitNewton{TdFunc: f, DeltaT: dt, Options: o}
	case "MidpointNewton":
		return &MidPointNewton{TdFunc: f, DeltaT: dt, Options: o}
	case "SDIRK2":
		return &SDIRK2{TdFunc: f, DeltaT: dt, Options: o}
	default:
		panic(name)
	}
}

var implicitNames = []string{"EulerFP", "HeunFP", "MidpointFP", "PredictorCorrector", "EulerNewton", "HeunNewton", "MidpointNewton", "SDIRK2"}

func TestImplicitOrderAndScalarGridAgreement(t *testing.T) {
	for _, name := range implicitNames {
		t.Run(name, func(t *testing.T) {
			var errs []float64
			for _, steps := range []int{10, 20} {
				dt := 1 / float64(steps)
				s := implicitForTest(name, dt, exponentialODE{}, ImplicitOptions{AbsoluteTolerance: 1e-13, RelativeTolerance: 1e-13})
				x := 1.
				vector := []float64{1, 2}
				for i := 0; i < steps; i++ {
					var err error
					x, err = s.NextStep(x, float64(i)*dt)
					if err != nil {
						t.Fatal(err)
					}
					if err = s.NextStepOnGrid(vector, float64(i)*dt); err != nil {
						t.Fatal(err)
					}
					if math.Abs(x-vector[0]) > 1e-11 {
						t.Fatalf("scalar/grid disagree: %g %g", x, vector[0])
					}
				}
				errs = append(errs, math.Abs(x-(2*math.E-2)))
			}
			minRatio := 3.5
			if name == "EulerFP" || name == "EulerNewton" {
				minRatio = 1.8
			}
			if errs[0]/errs[1] < minRatio {
				t.Fatalf("unexpected convergence order: errors %v", errs)
			}
		})
	}
}

func TestImplicitCoupledStiffSystem(t *testing.T) {
	f := implicitRHS{grid: func(x, y []float64, _ float64) { y[0] = -1000*x[0] + 999*x[1]; y[1] = -x[1] }}
	jac := func(_ []float64, _ float64, j []float64) error { j[0] = -1000; j[1] = 999; j[3] = -1; return nil }
	for _, name := range []string{"EulerNewton", "HeunNewton", "MidpointNewton", "SDIRK2"} {
		stability := func(z float64) float64 {
			switch name {
			case "EulerNewton":
				return 1 / (1 - z)
			case "SDIRK2":
				g := 1 - 1/math.Sqrt(2)
				return (1 + (1-2*g)*z) / math.Pow(1-g*z, 2)
			default:
				return (1 + z/2) / (1 - z/2)
			}
		}
		for _, analytic := range []bool{false, true} {
			o := ImplicitOptions{}
			if analytic {
				o.Jacobian = jac
			}
			s := implicitForTest(name, .1, f, o)
			x := []float64{1, 2}
			if err := s.NextStepOnGrid(x, 0); err != nil {
				t.Fatalf("%s analytic=%v: %v", name, analytic, err)
			}
			want := []float64{2*stability(-.1) - stability(-100), 2 * stability(-.1)}
			for i := range x {
				if math.Abs(x[i]-want[i]) > 1e-8 {
					t.Fatalf("%s coupled output %v want %v", name, x, want)
				}
			}
		}
	}
}

func TestImplicitFailuresPreserveState(t *testing.T) {
	for _, name := range implicitNames {
		for _, tc := range []struct {
			f       implicitRHS
			dt      float64
			options ImplicitOptions
		}{
			{implicitRHS{scalar: func(x, t float64) float64 { return math.NaN() }}, .1, ImplicitOptions{}},
			{implicitRHS{scalar: func(x, t float64) float64 { return x*x + 1 }}, 1, ImplicitOptions{MaxIterations: 1}},
			{implicitRHS{scalar: func(x, t float64) float64 { return x }}, math.Inf(1), ImplicitOptions{}},
			{implicitRHS{scalar: func(x, t float64) float64 { return x }}, .1, ImplicitOptions{AbsoluteTolerance: -1}},
		} {
			s := implicitForTest(name, tc.dt, tc.f, tc.options)
			x := []float64{1, 2}
			if err := s.NextStepOnGrid(x, 0); err == nil {
				t.Fatalf("%s accepted invalid step", name)
			}
			if !slices.Equal(x, []float64{1, 2}) {
				t.Fatalf("%s mutated a failed grid step", name)
			}
			if got, err := s.NextStep(1, 0); err == nil || got != 1 {
				t.Fatalf("%s scalar failure: %g %v", name, got, err)
			}
		}
	}
	// Exactly singular matrix: I - h J = 0, with nonzero residual.
	s := &EulerImplicitNewton{TdFunc: implicitRHS{scalar: func(x, t float64) float64 { return x }}, DeltaT: 1, Options: ImplicitOptions{Jacobian: func(x []float64, t float64, j []float64) error {
		for i := range x {
			j[i*len(x)+i] = 1
		}
		return nil
	}}}
	x := []float64{1, 2}
	if err := s.NextStepOnGrid(x, 0); err == nil || !slices.Equal(x, []float64{1, 2}) {
		t.Fatal("singular solve must fail atomically")
	}
	sentinel := errors.New("Jacobian unavailable")
	s.Options.Jacobian = func([]float64, float64, []float64) error { return sentinel }
	if err := s.NextStepOnGrid(x, 0); !errors.Is(err, sentinel) {
		t.Fatalf("callback error not preserved: %v", err)
	}
}

func TestImplicitZeroResizeAndLegacyPredictor(t *testing.T) {
	f := implicitRHS{scalar: func(x, t float64) float64 { return -x }}
	for _, name := range implicitNames {
		s := implicitForTest(name, .01, f, ImplicitOptions{})
		for _, n := range []int{1, 5, 2, 5} {
			x := make([]float64, n)
			if err := s.NextStepOnGrid(x, 0); err != nil {
				t.Fatal(name, err)
			}
			for i := range x {
				x[i] = float64(i + 1)
			}
			if err := s.NextStepOnGrid(x, 0); err != nil {
				t.Fatal(name, err)
			}
			if x[0] >= 1 || x[0] <= 0 {
				t.Fatal("resize lost initial state")
			}
		}
		zero := implicitForTest(name, 0, f, ImplicitOptions{})
		if y, err := zero.NextStep(2, 3); err != nil || y != 2 {
			t.Fatal("zero dt")
		}
		back := implicitForTest(name, -.01, f, ImplicitOptions{})
		if y, err := back.NextStep(1, 1); err != nil || y <= 1 {
			t.Fatal("backward step")
		}
	}
	x, fv, p := []float64{1, 2}, []float64{3, 4}, make([]float64, 2)
	(&HeunsImplicitNewton{DeltaT: .5}).PredictIni(x, fv, p)
	if !slices.Equal(p, []float64{2.5, 4}) {
		t.Fatal("Newton predictor not written")
	}
	clear(p)
	(&HeunsImplicitFixPoint{DeltaT: .5}).PredictIni(x, fv, p)
	if !slices.Equal(p, []float64{2.5, 4}) {
		t.Fatal("fixed-point predictor not written")
	}
	m := &MidPointNewton{TdFunc: f, DeltaT: .1}
	a, e := m.NextStep(1, 0)
	b, e2 := m.NextStepIm(1, 0)
	if e != nil || e2 != nil || a != b {
		t.Fatal("legacy midpoint alias")
	}
}

func TestImplicitNonlinearDampingAndReuse(t *testing.T) {
	// The undamped first Newton update is 10 from x=0; backtracking is needed.
	f := implicitRHS{scalar: func(x, t float64) float64 { return 10 - x*x*x }}
	s := &EulerImplicitNewton{TdFunc: f, DeltaT: 1, Options: ImplicitOptions{Derivative: func(x, t float64) float64 { return -3 * x * x }}}
	y, err := s.NextStep(0, 0)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(y+y*y*y-10) > 1e-8 || s.Stats().Backtracks == 0 {
		t.Fatalf("damping failed: x=%g stats=%+v", y, s.Stats())
	}
	f = implicitRHS{scalar: func(x, t float64) float64 { return -x * x * x }}
	jac := func(x []float64, t float64, j []float64) error {
		for i, v := range x {
			j[i*len(x)+i] = -3 * v * v
		}
		return nil
	}
	for _, refresh := range []int{1, 4} {
		s := &EulerImplicitNewton{TdFunc: f, DeltaT: .1, Options: ImplicitOptions{Jacobian: jac, JacobianRefresh: refresh}}
		x := []float64{1, 2}
		if err := s.NextStepOnGrid(x, 0); err != nil {
			t.Fatal(err)
		}
		if refresh == 4 && s.Stats().Factorizations >= s.Stats().Iterations {
			t.Fatalf("no factor reuse: %+v", s.Stats())
		}
		for i, v := range x {
			if math.Abs(v+.1*v*v*v-float64(i+1)) > 1e-7 {
				t.Fatal("nonlinear residual")
			}
		}
	}
}

func TestSDIRKStiffDecayAndStageFailure(t *testing.T) {
	s := &SDIRK2{TdFunc: implicitRHS{scalar: func(x, t float64) float64 { return -1e6 * x }}, DeltaT: 1}
	y, err := s.NextStep(1, 0)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(y) > 1e-5 {
		t.Fatal("SDIRK must damp a very stiff decaying mode", y)
	}
	s.TdFunc = implicitRHS{scalar: func(x, t float64) float64 {
		if t > .5 {
			return math.NaN()
		}
		return -x
	}}
	x := []float64{1, 2}
	if err := s.NextStepOnGrid(x, 0); err == nil || !slices.Equal(x, []float64{1, 2}) {
		t.Fatal("second-stage failure changed state")
	}
	if !math.IsInf(s.Stats().Residual, 1) {
		t.Fatal("failed second stage reported a converged residual")
	}
}

func TestImplicitDiagonalFastPath(t *testing.T) {
	f := implicitRHS{scalar: func(x, t float64) float64 { return -100 * x }}
	diag := func(x []float64, t float64, d []float64) error {
		for i := range d {
			d[i] = -100
		}
		return nil
	}
	s := &EulerImplicitNewton{TdFunc: f, DeltaT: .1, Options: ImplicitOptions{DiagonalJacobian: diag}}
	x := make([]float64, 2048)
	for i := range x {
		x[i] = 1
	}
	if err := s.NextStepOnGrid(x, 0); err != nil {
		t.Fatal(err)
	}
	for _, v := range x {
		if math.Abs(v-1./11) > 1e-12 {
			t.Fatal("diagonal solve")
		}
	}
	if allocations := testing.AllocsPerRun(20, func() {
		for i := range x {
			x[i] = 1
		}
		if err := s.NextStepOnGrid(x, 0); err != nil {
			panic(err)
		}
	}); allocations != 0 {
		t.Fatalf("warm diagonal step allocated %g times", allocations)
	}
	s.Options = ImplicitOptions{}
	if err := s.NextStepOnGrid(x, 0); err == nil {
		t.Fatal("missing dense size guard")
	}
}

func BenchmarkImplicitNewton(b *testing.B) {
	f := implicitRHS{scalar: func(x, t float64) float64 { return -100 * x }}
	for _, mode := range []string{"finite-difference", "analytic-dense", "analytic-diagonal"} {
		b.Run(mode, func(b *testing.B) {
			o := ImplicitOptions{}
			if mode == "analytic-dense" {
				o.Jacobian = func(x []float64, t float64, j []float64) error {
					for i := range x {
						j[i*len(x)+i] = -100
					}
					return nil
				}
			}
			if mode == "analytic-diagonal" {
				o.DiagonalJacobian = func(x []float64, t float64, d []float64) error {
					for i := range d {
						d[i] = -100
					}
					return nil
				}
			}
			s := &EulerImplicitNewton{TdFunc: f, DeltaT: .1, Options: o}
			x := make([]float64, 64)
			for i := range x {
				x[i] = 1
			}
			if err := s.NextStepOnGrid(x, 0); err != nil {
				b.Fatal(err)
			}
			b.ReportAllocs()
			b.ResetTimer()
			for k := 0; k < b.N; k++ {
				for i := range x {
					x[i] = 1
				}
				if err := s.NextStepOnGrid(x, 0); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}
