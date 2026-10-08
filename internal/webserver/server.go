package webserver

import (
	"context"
	"encoding/json"
	"errors"
	"io"
	"mime"
	"net/http"
	"os"
	"path/filepath"
	"time"

	"GoProject/internal/plotmath"
)

// Handler serves the built frontend and a bounded, same-origin plotting API.
func Handler(dist string) http.Handler {
	mux := http.NewServeMux()
	registerSimulations(mux)
	slots := make(chan struct{}, 2)
	mux.HandleFunc("POST /api/plot", func(w http.ResponseWriter, r *http.Request) {
		mediaType, _, err := mime.ParseMediaType(r.Header.Get("Content-Type"))
		if err != nil || mediaType != "application/json" {
			fail(w, 415, "Send application/json.")
			return
		}
		select {
		case slots <- struct{}{}:
			defer func() { <-slots }()
		default:
			fail(w, 429, "Plot engine is busy. Try again shortly.")
			return
		}
		r.Body = http.MaxBytesReader(w, r.Body, 8192)
		defer r.Body.Close()
		decoder := json.NewDecoder(r.Body)
		decoder.DisallowUnknownFields()
		var request plotmath.Request
		if err := decoder.Decode(&request); err != nil {
			var oversized *http.MaxBytesError
			if errors.As(err, &oversized) {
				fail(w, 413, "Plot request is too large.")
			} else {
				fail(w, 400, "Invalid plot request. Check the fields and numeric values.")
			}
			return
		}
		if err := decoder.Decode(new(any)); err != io.EOF {
			fail(w, 400, "Send exactly one plot request.")
			return
		}
		ctx, cancel := context.WithTimeout(r.Context(), 5*time.Second)
		defer cancel()
		result, err := plotmath.Sample(ctx, request)
		if err != nil {
			if errors.Is(err, context.Canceled) {
				return
			}
			if errors.Is(err, context.DeadlineExceeded) {
				fail(w, 408, "Plot exceeded the time limit. Reduce samples or simplify the expression.")
				return
			}
			fail(w, 400, err.Error())
			return
		}
		w.Header().Set("Content-Type", "application/json")
		w.Header().Set("Cache-Control", "no-store")
		_ = json.NewEncoder(w).Encode(result)
	})
	mux.HandleFunc("GET /api/health", func(w http.ResponseWriter, r *http.Request) {
		w.Header().Set("Content-Type", "application/json")
		_, _ = io.WriteString(w, `{"engine":"go","status":"ok"}`)
	})
	mux.HandleFunc("/api/", func(w http.ResponseWriter, r *http.Request) {
		fail(w, 404, "Unknown API endpoint or unsupported method.")
	})
	mux.HandleFunc("GET /calculator/", func(w http.ResponseWriter, r *http.Request) {
		http.Redirect(w, r, "/#tools", http.StatusTemporaryRedirect)
	})
	files := http.FileServer(http.Dir(dist))
	mux.HandleFunc("/", func(w http.ResponseWriter, r *http.Request) {
		if r.Method != http.MethodGet && r.Method != http.MethodHead {
			w.Header().Set("Allow", "GET, HEAD")
			fail(w, http.StatusMethodNotAllowed, "Use GET to load the frontend.")
			return
		}
		if _, err := os.Stat(filepath.Join(dist, "index.html")); err != nil {
			http.Error(w, "Build the frontend first: cd webpage/web-app && bun run build", 503)
			return
		}
		files.ServeHTTP(w, r)
	})
	protected := http.NewCrossOriginProtection().Handler(mux)
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		w.Header().Set("X-Content-Type-Options", "nosniff")
		w.Header().Set("Referrer-Policy", "same-origin")
		protected.ServeHTTP(w, r)
	})
}

func fail(w http.ResponseWriter, status int, message string) {
	w.Header().Set("Content-Type", "application/json")
	w.Header().Set("Cache-Control", "no-store")
	w.WriteHeader(status)
	_ = json.NewEncoder(w).Encode(map[string]string{"error": message})
}
