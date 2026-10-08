package simulation

import (
	"context"
	"math"
	"testing"
)

func testRequest(kind string) Request {
	return Request{Kind: kind, Basis: "dvr", Method: "split", Potential: "harmonic", Mass: 1, Strength: 1, Alpha: .5, HalfWidth: 8, Points: 64, States: 4, Dt: .01, Steps: 100, X0: 1, Sigma: math.Sqrt(.5), Particles: 16, Box: 6, Temperature: .2, Seed: 42}
}
func metric(r *Result, label string) float64 {
	for _, m := range r.Metrics {
		if m.Label == label {
			return m.Value
		}
	}
	return math.NaN()
}

func TestHarmonicSpectrum(t *testing.T) {
	for _, basis := range []string{"dvr", "fourier"} {
		r := testRequest("spectrum")
		r.Basis = basis
		result, err := Run(context.Background(), r)
		if err != nil {
			t.Fatal(err)
		}
		for i := 0; i < 4; i++ {
			if math.Abs(result.Metrics[i].Value-(float64(i)+.5)) > 1e-6 {
				t.Fatalf("%s E%d = %g", basis, i, result.Metrics[i].Value)
			}
		}
		if metric(result, "Largest eigenpair residual") > 1e-10 {
			t.Fatal("large eigenpair residual")
		}
	}
}

func TestQuantumMotion(t *testing.T) {
	for _, method := range []string{"split", "rk4", "ralston3", "nystrom5"} {
		r := testRequest("quantum")
		r.Method = method
		r.Dt = .002
		r.Steps = 500
		result, err := Run(context.Background(), r)
		if err != nil {
			t.Fatalf("%s: %v", method, err)
		}
		if math.Abs(metric(result, "Final mean position")-math.Cos(1)) > 3e-5 {
			t.Fatalf("%s mean position: %g", method, metric(result, "Final mean position"))
		}
		if metric(result, "Maximum norm drift") > 1e-6 {
			t.Fatalf("%s norm drift %g", method, metric(result, "Maximum norm drift"))
		}
	}
}

func TestDrivenHamiltonian(t *testing.T) {
	r := testRequest("quantum")
	r.DriveAmplitude = .2
	r.DriveFrequency = 1
	r.X0 = 0
	r.Dt = .005
	r.Steps = 200
	result, err := Run(context.Background(), r)
	if err != nil {
		t.Fatal(err)
	}
	// x''+x=-A sin(t), x(0)=x'(0)=0.
	want := .1 * (math.Cos(1) - math.Sin(1))
	if math.Abs(metric(result, "Final mean position")-want) > 1e-5 {
		t.Fatalf("driven <x>: got %g want %g", metric(result, "Final mean position"), want)
	}
}

func TestClassicalMassAndMethods(t *testing.T) {
	for _, method := range []string{"verlet", "leapfrog", "yoshida"} {
		r := testRequest("classical")
		r.Method = method
		r.Mass = 2
		r.Dt = .005
		r.Steps = 200
		result, err := Run(context.Background(), r)
		if err != nil {
			t.Fatal(err)
		}
		if math.Abs(metric(result, "Final position")-math.Cos(1/math.Sqrt(2))) > 2e-6 {
			t.Fatalf("%s incorrect mass handling", method)
		}
		if metric(result, "Maximum absolute energy drift") > 1e-5 {
			t.Fatal("excessive classical energy drift")
		}
	}
	r := testRequest("molecular")
	r.Potential = "morse"
	r.Method = "verlet"
	r.X0 = .2
	_, err := Run(context.Background(), r)
	if err != nil {
		t.Fatal(err)
	}
}

func TestMDAndValidation(t *testing.T) {
	r := testRequest("md")
	r.Method = "verlet"
	r.Dt = .001
	r.Steps = 1000
	result, err := Run(context.Background(), r)
	if err != nil {
		t.Fatal(err)
	}
	if len(result.Frames) < 2 || metric(result, "Total momentum magnitude") > 1e-10 || metric(result, "Maximum absolute energy drift") > 1e-3 {
		t.Fatalf("MD diagnostics: %+v", result.Metrics)
	}
	bad := testRequest("quantum")
	bad.Points = 10000
	if _, err := Run(context.Background(), bad); err == nil {
		t.Fatal("accepted unbounded grid")
	}
	bad = testRequest("quantum")
	bad.Method = "nystrom5"
	bad.Dt = .05
	if _, err := Run(context.Background(), bad); err == nil {
		t.Fatal("accepted unstable explicit step")
	}
	ctx, cancel := context.WithCancel(context.Background())
	cancel()
	if _, err := Run(ctx, testRequest("spectrum")); err != context.Canceled {
		t.Fatal("cancellation ignored")
	}
}
