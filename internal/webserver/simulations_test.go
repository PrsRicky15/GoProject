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
