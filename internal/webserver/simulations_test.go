package webserver

import (
	"GoProject/internal/simulation"
	"encoding/json"
	"net/http/httptest"
	"strings"
	"testing"
)

func TestSimulationAPI(t *testing.T) {
	handler := Handler(t.TempDir())
	valid := `{"kind":"spectrum","basis":"dvr","potential":"harmonic","mass":1,"strength":1,"alpha":0.5,"halfWidth":8,"points":64,"states":4}`
	cases := []struct {
		body, contentType, origin string
		status                    int
	}{
		{valid, "application/json", "", 200},
		{valid, "text/plain", "", 415},
		{valid, "application/json", "https://untrusted.example", 403},
		{valid + ` {}`, "application/json", "", 400},
		{`{"kind":"spectrum","typo":1}`, "application/json", "", 400},
		{strings.Replace(valid, `"points":64`, `"points":1000000`, 1), "application/json", "", 400},
		{`{"kind":"` + strings.Repeat("x", 9000) + `"}`, "application/json", "", 413},
	}
	for _, tc := range cases {
		r := httptest.NewRequest("POST", "http://localhost/api/simulations", strings.NewReader(tc.body))
		r.Header.Set("Content-Type", tc.contentType)
		if tc.origin != "" {
			r.Header.Set("Origin", tc.origin)
		}
		w := httptest.NewRecorder()
		handler.ServeHTTP(w, r)
		if w.Code != tc.status {
			t.Fatalf("got %d want %d: %s", w.Code, tc.status, w.Body)
		}
		if tc.status == 200 {
			var result simulation.Result
			if err := json.Unmarshal(w.Body.Bytes(), &result); err != nil {
				t.Fatal(err)
			}
			if len(result.Metrics) != 5 || len(result.Charts) != 2 || w.Header().Get("Cache-Control") != "no-store" {
				t.Fatal("incomplete response")
			}
		}
	}
}

func TestNDimensionalSimulationAPI(t *testing.T) {
	handler := Handler(t.TempDir())
	valid := `{"kind":"nd-spectrum","basis":"fourier","potential":"gaussian","strength":5,"sigma":1.2,"alpha":0.25,"states":2,"nd":{"axes":[{"points":16,"halfWidth":6,"mass":1,"scale":1,"center":0},{"points":16,"halfWidth":6,"mass":1,"scale":1,"center":0}],"coupling":0.2,"softening":0.7,"separation":3,"wavenumber":1,"order":4,"state":0,"axisX":0,"axisY":1}}`
	for _, tc := range []struct {
		body   string
		status int
	}{
		{valid, 200}, {strings.Replace(valid, `"scale":1`, `"scale":0`, 1), 400}, {strings.Replace(valid, `"order":4`, `"typo":4`, 1), 400}, {strings.Replace(valid, `"axisY":1`, `"axisY":0`, 1), 400},
	} {
		r := httptest.NewRequest("POST", "http://localhost/api/simulations", strings.NewReader(tc.body))
		r.Header.Set("Content-Type", "application/json")
		w := httptest.NewRecorder()
		handler.ServeHTTP(w, r)
		if w.Code != tc.status {
			t.Fatalf("got %d want %d: %s", w.Code, tc.status, w.Body)
		}
		if tc.status == 200 {
			var result simulation.Result
			if err := json.Unmarshal(w.Body.Bytes(), &result); err != nil {
				t.Fatal(err)
			}
			if result.Request.ND == nil || len(result.Charts) != 3 || result.Charts[2].ValueLabel != "V" || result.Metrics[0].Value >= 0 {
				t.Fatal("invalid nD well result")
			}
		}
	}
}

func TestCustomPotentialAPI(t *testing.T) {
	handler := Handler(t.TempDir())
	valid := `{"kind":"nd-spectrum","basis":"fourier","potential":"custom","states":1,"nd":{"expression":"0.5*(x1^2+x2^2)","axes":[{"points":16,"halfWidth":6,"mass":1,"scale":1},{"points":16,"halfWidth":6,"mass":1,"scale":1}],"axisX":0,"axisY":1}}`
	for _, tc := range []struct {
		expression string
		status     int
	}{
		{"0.5*(x1^2+x2^2)", 200}, {"x3", 400}, {"", 400}, {"sqrt(-1)", 400}, {"x1;exit()", 400},
	} {
		body := strings.Replace(valid, "0.5*(x1^2+x2^2)", tc.expression, 1)
		r := httptest.NewRequest("POST", "http://localhost/api/simulations", strings.NewReader(body))
		r.Header.Set("Content-Type", "application/json")
		w := httptest.NewRecorder()
		handler.ServeHTTP(w, r)
		if w.Code != tc.status {
			t.Fatalf("%q: got %d want %d: %s", tc.expression, w.Code, tc.status, w.Body)
		}
		if tc.status == 200 {
			var result simulation.Result
			if err := json.Unmarshal(w.Body.Bytes(), &result); err != nil {
				t.Fatal(err)
			}
			if result.Request.ND.Expression != tc.expression || len(result.Charts) != 3 {
				t.Fatal("incomplete custom result")
			}
		} else if !strings.Contains(w.Body.String(), "custom potential") {
			t.Fatal("missing formula error")
		}
	}
}
