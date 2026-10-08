package webserver

import (
	"GoProject/internal/simulation"
	"context"
	"encoding/json"
	"errors"
	"io"
	"mime"
	"net/http"
	"time"
)

func registerSimulations(mux *http.ServeMux) {
	slots := make(chan struct{}, 1)
	mux.HandleFunc("POST /api/simulations", func(w http.ResponseWriter, r *http.Request) {
		media, _, err := mime.ParseMediaType(r.Header.Get("Content-Type"))
		if err != nil || media != "application/json" {
			fail(w, 415, "Send application/json.")
			return
		}
		select {
		case slots <- struct{}{}:
			defer func() { <-slots }()
		default:
			fail(w, 429, "A calculation is already running. Try again when it finishes.")
			return
		}
		r.Body = http.MaxBytesReader(w, r.Body, 8192)
		defer r.Body.Close()
		decoder := json.NewDecoder(r.Body)
		decoder.DisallowUnknownFields()
		var request simulation.Request
		if err := decoder.Decode(&request); err != nil {
			var oversized *http.MaxBytesError
			if errors.As(err, &oversized) {
				fail(w, 413, "Calculation request is too large.")
			} else {
				fail(w, 400, "Invalid calculation request. Check numeric values and field names.")
			}
			return
		}
		if err := decoder.Decode(new(any)); err != io.EOF {
			fail(w, 400, "Send exactly one calculation.")
			return
		}
		ctx, cancel := context.WithTimeout(r.Context(), 10*time.Second)
		defer cancel()
		result, err := simulation.Run(ctx, request)
		if err != nil {
			if errors.Is(err, context.Canceled) {
				return
			}
			if errors.Is(err, context.DeadlineExceeded) {
				fail(w, 408, "Calculation exceeded the time limit. Reduce the grid or number of steps.")
				return
			}
			fail(w, 400, err.Error())
			return
		}
		w.Header().Set("Content-Type", "application/json")
		w.Header().Set("Cache-Control", "no-store")
		_ = json.NewEncoder(w).Encode(result)
	})
}
