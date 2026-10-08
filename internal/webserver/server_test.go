package webserver

import (
	"encoding/json"
	"net/http"
	"net/http/httptest"
	"strings"
	"testing"

	"GoProject/internal/plotmath"
)

const validRequest = `{"mode":"line","expression":"1/x","resolution":21,"bounds":{"x":[-2,2]},"level":0,"style":"heatmap"}`

func TestPlotAPI(t *testing.T) {
	handler := Handler(t.TempDir())
	r := httptest.NewRequest("POST", "http://localhost/api/plot", strings.NewReader(validRequest))
	r.Header.Set("Content-Type", "application/json")
	w := httptest.NewRecorder()
	handler.ServeHTTP(w, r)
	if w.Code != 200 {
		t.Fatalf("%d: %s", w.Code, w.Body)
	}
	var result plotmath.Result
	if err := json.Unmarshal(w.Body.Bytes(), &result); err != nil {
		t.Fatal(err)
	}
	if result.Count != 21 || result.Invalid != 1 || result.Values[10] != nil {
		t.Fatal("invalid result", result)
	}
	if w.Header().Get("Cache-Control") != "no-store" {
		t.Fatal("plot responses must not be cached")
	}
}

func TestAPIRejectsInvalidRequests(t *testing.T) {
	handler := Handler(t.TempDir())
	cases := []struct {
		body, contentType, origin string
		status                    int
	}{
		{validRequest, "text/plain", "", 415},
		{validRequest, "application/json", "https://untrusted.example", 403},
		{`{"mode":"line","unknown":1}`, "application/json", "", 400},
		{validRequest + ` {}`, "application/json", "", 400},
		{strings.Replace(validRequest, `"resolution":21`, `"resolution":1e9`, 1), "application/json", "", 400},
		{strings.Replace(validRequest, `"resolution":21`, `"resolution":20.5`, 1), "application/json", "", 400},
		{`{"expression":"` + strings.Repeat("x", 9000) + `"}`, "application/json", "", 413},
		{strings.Replace(validRequest, `"1/x"`, `"globalThis"`, 1), "application/json", "", 400},
		{validRequest, "application/json", "http://localhost", 200},
	}
	for _, tc := range cases {
		r := httptest.NewRequest("POST", "http://localhost/api/plot", strings.NewReader(tc.body))
		r.Header.Set("Content-Type", tc.contentType)
		if tc.origin != "" {
			r.Header.Set("Origin", tc.origin)
		}
		w := httptest.NewRecorder()
		handler.ServeHTTP(w, r)
		if w.Code != tc.status {
			t.Errorf("want %d, got %d: %s", tc.status, w.Code, w.Body)
		}
	}
	for _, method := range []string{http.MethodGet, http.MethodPut, http.MethodDelete} {
		w := httptest.NewRecorder()
		handler.ServeHTTP(w, httptest.NewRequest(method, "/api/plot", nil))
		if w.Code < 400 {
			t.Errorf("accepted %s", method)
		}
	}
}
