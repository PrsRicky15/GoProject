package simulation

import (
	"context"
	"errors"
	"math"
	"strings"
	"testing"
)

func TestWebVectorMethods(t *testing.T) {
	methods := []string{}
	for m := range fixedRK {
		methods = append(methods, m)
	}
	for m := range adaptiveRK {
		methods = append(methods, m)
	}
	for m := range implicitODE {
		methods = append(methods, m)
	}
	for _, kind := range []string{"quantum", "classical", "molecular"} {
		for _, method := range methods {
			t.Run(kind+"/"+method, func(t *testing.T) {
				r := testRequest(kind)
				r.Method = method
				r.Points = 32
				r.Dt = .001
				r.Steps = 100
				if kind == "molecular" {
					r.Potential = "morse"
					r.X0 = .2
				}
				result, err := Run(context.Background(), r)
				if err != nil {
					t.Fatal(err)
				}
				if metric(result, "RHS evaluations") <= 0 || metric(result, "Accepted solver steps") < 100 {
					t.Fatal("missing solver diagnostics", result.Metrics)
				}
				var got, want float64
				if kind == "quantum" {
					got = metric(result, "Final mean position")
					want = math.Cos(.1)
				} else if kind == "classical" {
					got = metric(result, "Final position")
					want = math.Cos(.1)
				} else {
					return
				}
				if math.Abs(got-want) > .001 {
					t.Fatalf("incorrect coupled propagation %g want %g", got, want)
				}
			})
		}
	}
}

func TestWebDrivenAdaptiveAndImplicit(t *testing.T) {
	for _, basis := range []string{"dvr", "fourier"} {
		for _, method := range []string{"dormand-prince54", "cash-karp54", "implicit-midpoint", "trapezoidal", "sdirk2"} {
			t.Run(basis+"/"+method, func(t *testing.T) {
				r := testRequest("quantum")
				r.Basis = basis
				r.Method = method
				r.X0 = 0
				r.DriveAmplitude = .2
				r.DriveFrequency = 1
				r.Dt = .01
				r.Steps = 100
				r.AbsoluteTolerance = 1e-11
				r.RelativeTolerance = 1e-9
				result, err := Run(context.Background(), r)
				if err != nil {
					t.Fatal(err)
				}
				want := .1 * (math.Cos(1) - math.Sin(1))
				if math.Abs(metric(result, "Final mean position")-want) > 2e-5 {
					t.Fatal("driven reference mismatch", result.Metrics)
				}
				if method == "implicit-midpoint" && metric(result, "Maximum norm drift") > 1e-7 {
					t.Fatal("midpoint norm drift")
				}
				times := result.Charts[0].Y
				if times[len(times)-1] != 1 {
					t.Fatal("wrong endpoint")
				}
			})
		}
	}
}

func TestWebAdaptiveRejectionsAndValidation(t *testing.T) {
	// The most work-intensive default adaptive selection must fit the web
	// budget; low-order error control needs far more RHS calls than RK5.
	defaultQuantum := testRequest("quantum")
	defaultQuantum.Method = "heun-euler21"
	defaultQuantum.Sigma = .7
	defaultQuantum.Steps = 400
	if _, err := Run(context.Background(), defaultQuantum); err != nil {
		t.Fatalf("default adaptive quantum job: %v", err)
	}
	r := testRequest("classical")
	r.Method = "heun-euler21"
	r.Dt = .05
	r.Steps = 2
	r.AbsoluteTolerance = 1e-10
	r.RelativeTolerance = 1e-10
	result, err := Run(context.Background(), r)
	if err != nil {
		t.Fatal(err)
	}
	if metric(result, "Rejected solver steps") < 1 || metric(result, "Accepted solver steps") <= 2 {
		t.Fatal("expected adaptive subdivision", result.Metrics)
	}
	times := result.Charts[0].Series[0].X
	if times[len(times)-1] != .1 {
		t.Fatal("rejection shifted endpoint")
	}
	for _, value := range []float64{-1, 1e-15, .1, math.NaN(), math.Inf(1)} {
		bad := r
		bad.AbsoluteTolerance = value
		if Validate(bad) == nil {
			t.Fatal("bad tolerance accepted", value)
		}
	}
	r.Method = "unrecognized"
	if Validate(r) == nil {
		t.Fatal("unknown method accepted")
	}
	r = testRequest("md")
	r.Method = "dormand-prince54"
	if Validate(r) == nil {
		t.Fatal("unsupported fluid method accepted")
	}
}

func TestWebSolverBudgetAndCancellation(t *testing.T) {
	for _, method := range []string{"rk4", "dormand-prince54", "implicit-midpoint"} {
		r := testRequest("classical")
		r.Method = method
		r.AbsoluteTolerance = 1e-9
		r.RelativeTolerance = 1e-7
		state := []float64{1, 0}
		w := &solverWork{limit: 1}
		step, err := realVectorStepper(context.Background(), r, func(t float64, y, dy []float64) error { dy[0] = y[1]; dy[1] = -y[0]; return nil }, state, w)
		if err != nil {
			t.Fatal(err)
		}
		if err = step(0, .01); err == nil || !strings.Contains(err.Error(), "work limit") || state[0] != 1 || state[1] != 0 {
			t.Fatalf("%s budget not atomic: %v %v", method, state, err)
		}
		ctx, cancel := context.WithCancel(context.Background())
		cancel()
		w = &solverWork{limit: 100}
		step, err = realVectorStepper(ctx, r, func(t float64, y, dy []float64) error { return nil }, state, w)
		if err != nil {
			t.Fatal(err)
		}
		if err = step(0, .01); !errors.Is(err, context.Canceled) {
			t.Fatalf("%s cancellation: %v", method, err)
		}
	}
}
