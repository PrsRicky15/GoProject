package plotmath

import (
	"context"
	"encoding/json"
	"math"
	"os"
	"strings"
	"testing"
)

func TestExpressions(t *testing.T) {
	// The TypeScript calculator runs these same fixtures to prevent grammar drift.
	data, err := os.ReadFile("../../testdata/math-expressions.json")
	if err != nil {
		t.Fatal(err)
	}
	var cases []struct {
		Expression string  `json:"expression"`
		Expected   float64 `json:"expected"`
	}
	if err := json.Unmarshal(data, &cases); err != nil {
		t.Fatal(err)
	}
	for _, tc := range cases {
		t.Run(tc.Expression, func(t *testing.T) {
			fn, err := Compile(tc.Expression, 0)
			if err != nil {
				t.Fatal(err)
			}
			if got := fn(0, 0, 0); math.Abs(got-tc.Expected) > 1e-12 {
				t.Fatalf("got %g, want %g", got, tc.Expected)
			}
		})
	}
}

func TestRejectedExpressions(t *testing.T) {
	for _, source := range []string{"globalThis", "constructor(1)", "Math.sin(1)", "alert(1)", "x=1", "2x", "sin()", "1+", "pow(2)", "min(1,2,3)", "(2", "2)", "1;2", "2[0]", "", "1e999", strings.Repeat("1+", 200) + "1", strings.Repeat(" ", 1025)} {
		if _, err := Compile(source, 1); err == nil {
			t.Errorf("accepted %q", source)
		}
	}
	if _, err := Compile("y", 1); err == nil {
		t.Error("accepted y for 1D")
	}
}

func request(mode string, n int, expression string) Request {
	return Request{Mode: mode, Resolution: n, Expression: expression, Bounds: map[string][]float64{"x": {-2, 2}, "y": {-2, 2}, "z": {-2, 2}}, Level: 1, Style: "heatmap"}
}

func TestSampling(t *testing.T) {
	ctx := context.Background()
	line, err := Sample(ctx, request("line", 801, "exp(-x^2/2)"))
	if err != nil {
		t.Fatal(err)
	}
	if line.Count != 801 || *line.Values[400] != 1 || line.Coords["x"][0] != -2 || line.Coords["x"][800] != 2 {
		t.Fatal("incorrect Gaussian or endpoints")
	}
	for _, mode := range []string{"heatmap", "surface", "iso"} {
		r := request(mode, 10, "x+10*y")
		if mode == "iso" {
			r.Expression += "+100*z"
		}
		result, err := Sample(ctx, r)
		if err != nil {
			t.Fatal(err)
		}
		for i, v := range result.Values {
			want := result.Coords["x"][i%10] + 10*result.Coords["y"][(i/10)%10]
			if mode == "iso" {
				want += 100 * result.Coords["z"][i/100]
			}
			if v == nil || math.Abs(*v-want) > 1e-10 {
				t.Fatalf("%s sample %d transposed", mode, i)
			}
		}
	}
	gap, err := Sample(ctx, request("line", 21, "1/x"))
	if err != nil {
		t.Fatal(err)
	}
	if gap.Invalid != 1 || gap.Values[10] != nil {
		t.Fatal("missing null at singularity")
	}
	if _, err := json.Marshal(gap); err != nil {
		t.Fatal("non-finite JSON", err)
	}
}

func TestSamplingLimits(t *testing.T) {
	for _, mode := range []string{"line", "heatmap", "surface", "iso"} {
		for _, n := range []int{0, -1, 1000000} {
			if _, err := Sample(context.Background(), request(mode, n, "x")); err == nil {
				t.Errorf("accepted %s grid %d", mode, n)
			}
		}
	}
	bad := []Request{request("iso", 31, "1/x"), request("line", 21, "sqrt(-1)"), request("other", 21, "x")}
	r := request("line", 21, "x")
	r.Bounds["x"] = []float64{1, 1}
	bad = append(bad, r)
	r = request("line", 21, "x")
	r.Bounds["x"] = []float64{1}
	bad = append(bad, r)
	r = request("iso", 31, "x^2+y^2+z^2")
	r.Level = 1000
	bad = append(bad, r)
	for _, r := range bad {
		if _, err := Sample(context.Background(), r); err == nil {
			t.Errorf("accepted invalid request %+v", r)
		}
	}
	ctx, cancel := context.WithCancel(context.Background())
	cancel()
	if _, err := Sample(ctx, request("iso", 60, "x^2+y^2+z^2")); err != context.Canceled {
		t.Fatalf("cancellation: %v", err)
	}
}
