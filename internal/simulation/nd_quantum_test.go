package simulation

import (
	"context"
	"fmt"
	"math"
	"testing"
)

func testNDRequest(dims int) Request {
	r := testRequest("nd-spectrum")
	r.Basis = "fourier"
	r.Sigma = 1.2
	r.Alpha = .25
	r.ND = &NDRequest{Axes: make([]NDAxis, dims), Coupling: .2, Softening: .7, Separation: 3, Wavenumber: 1, Order: 4, AxisX: 0, AxisY: 1}
	for d := range r.ND.Axes {
		r.ND.Axes[d] = NDAxis{Points: 24, HalfWidth: 6, Mass: 1, Scale: 1}
	}
	return r
}

func TestNDWebPotentialPresets(t *testing.T) {
	for _, dims := range []int{2, 3, 4} {
		for _, basis := range []string{"fourier", "dvr"} {
			for _, tc := range []struct {
				name     string
				strength float64
			}{
				{"harmonic", 1}, {"coupled-harmonic", 1}, {"morse", 2}, {"double-well", .2}, {"gaussian", 5}, {"super-gaussian", 5}, {"soft-core", 3}, {"double-gaussian", 4}, {"lattice", 2}, {"free", 1},
			} {
				t.Run(fmt.Sprintf("%dD/%s/%s", dims, basis, tc.name), func(t *testing.T) {
					r := testNDRequest(dims)
					for d := range r.ND.Axes {
						if dims == 3 {
							r.ND.Axes[d].Points = 16
						}
						if dims == 4 {
							r.ND.Axes[d].Points = 8
						}
					}
					r.Basis = basis
					r.Potential = tc.name
					r.Strength = tc.strength
					if tc.name == "super-gaussian" {
						r.Sigma = 2
					}
					result, err := Run(context.Background(), r)
					if err != nil {
						t.Fatal(err)
					}
					if len(result.Charts) != 3 || len(result.Metrics) != 7 {
						t.Fatal("incomplete nD response")
					}
					if dims == 2 && tc.name == "coupled-harmonic" {
						want := .5 * (math.Sqrt(1+.2) + math.Sqrt(1-.2))
						if math.Abs(result.Metrics[0].Value-want) > 1e-5 {
							t.Fatalf("coupled ground energy %g want %g", result.Metrics[0].Value, want)
						}
					}
					checkNDProjectionNorms(t, r, result)
				})
			}
		}
	}
}

func TestNDHiddenCoordinateIntegration(t *testing.T) {
	for _, dims := range []int{3, 4} {
		r := testNDRequest(dims)
		r.States = 2
		r.ND.State = 1
		r.ND.AxisX = dims - 1
		r.ND.AxisY = 0
		for d := range r.ND.Axes {
			r.ND.Axes[d].Points = 4 + d
			r.ND.Axes[d].HalfWidth = 3 + .5*float64(d)
		}
		result, err := Run(context.Background(), r)
		if err != nil {
			t.Fatal(err)
		}
		checkNDProjectionNorms(t, r, result)
		v := result.Volume
		if v == nil {
			t.Fatal("missing 3D view")
		}
		if len(v.Values) != 4*5*6 {
			t.Fatal("volume contains the wrong axes")
		}
		dv := 1.0
		for _, a := range r.ND.Axes[:3] {
			dv *= 2 * a.HalfWidth / float64(a.Points)
		}
		norm := 0.0
		for _, rho := range v.Values {
			norm += rho * dv
		}
		if math.Abs(norm-1) > 1e-10 {
			t.Fatalf("volume integral = %g", norm)
		}
		if v.X[0] != -3 || v.Y[0] != -3.5 || v.Z[0] != -4 || v.X[1] != v.X[0] || v.Z[1] <= v.Z[0] {
			t.Fatal("volume coordinates violate original grid ordering")
		}
	}
}

func checkNDProjectionNorms(t *testing.T, r Request, result *Result) {
	t.Helper()
	ax, ay := r.ND.Axes[r.ND.AxisX], r.ND.Axes[r.ND.AxisY]
	dx, dy := 2*ax.HalfWidth/float64(ax.Points), 2*ay.HalfWidth/float64(ay.Points)
	for _, s := range result.Charts[0].Series {
		norm := 0.0
		for _, rho := range s.Y {
			norm += rho * dx
		}
		if math.Abs(norm-1) > 1e-10 {
			t.Fatalf("1D projection integral=%g", norm)
		}
	}
	heat := result.Charts[1]
	if len(heat.X) != ax.Points || len(heat.Y) != ay.Points {
		t.Fatal("swapped display axes")
	}
	norm := 0.0
	for _, row := range heat.Z {
		for _, rho := range row {
			norm += rho * dx * dy
		}
	}
	if math.Abs(norm-1) > 1e-10 {
		t.Fatalf("2D projection integral=%g", norm)
	}
}

func TestNDValidationAndPotentialSlice(t *testing.T) {
	for _, change := range []func(*Request){
		func(r *Request) { r.ND = nil }, func(r *Request) { r.ND.Axes = r.ND.Axes[:1] }, func(r *Request) { r.ND.Axes[0].Points = 65 }, func(r *Request) { r.ND.AxisY = 0 }, func(r *Request) { r.ND.State = 4 }, func(r *Request) { r.ND.Coupling = .9 }, func(r *Request) { r.ND.Axes[0].Scale = math.NaN() }, func(r *Request) { r.Potential = "unknown" }, func(r *Request) { r.ND.Order = 3 },
	} {
		r := testNDRequest(2)
		change(&r)
		if err := Validate(r); err == nil {
			t.Fatal("accepted invalid nD request")
		}
	}
	if err := Validate(testNDRequest(4)); err == nil {
		t.Fatal("accepted oversized tensor grid")
	}
	r := testNDRequest(3)
	r.States = 1
	r.Strength = 2
	r.ND.AxisX = 2
	r.ND.AxisY = 0
	for d := range r.ND.Axes {
		r.ND.Axes[d].Points = 8
		r.ND.Axes[d].Center = .25 * float64(d)
		r.ND.Axes[d].Scale = .5 + float64(d)*.25
	}
	result, err := Run(context.Background(), r)
	if err != nil {
		t.Fatal(err)
	}
	cut := result.Charts[2]
	a, b := r.ND.Axes[2], r.ND.Axes[0]
	for iy, y := range cut.Y {
		for ix, x := range cut.X {
			want := math.Pow(a.Scale*(x-a.Center), 2) + math.Pow(b.Scale*(y-b.Center), 2)
			if math.Abs(cut.Z[iy][ix]-want) > 1e-12 {
				t.Fatal("potential cut does not hold hidden coordinates at centers")
			}
		}
	}
	ctx, cancel := context.WithCancel(context.Background())
	cancel()
	if _, err := Run(ctx, r); err != context.Canceled {
		t.Fatal("cancellation ignored")
	}
}
