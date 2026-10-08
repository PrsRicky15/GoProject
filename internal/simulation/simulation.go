// Package simulation adapts the repository's numerical solvers to bounded web jobs.
package simulation

import (
	"GoProject/gridData"
	"context"
	"fmt"
	"math"
)

type Request struct {
	Kind           string  `json:"kind"`
	Basis          string  `json:"basis"`
	Method         string  `json:"method"`
	Potential      string  `json:"potential"`
	Mass           float64 `json:"mass"`
	Strength       float64 `json:"strength"`
	Alpha          float64 `json:"alpha"`
	Center         float64 `json:"center"`
	HalfWidth      float64 `json:"halfWidth"`
	Points         int     `json:"points"`
	States         int     `json:"states"`
	Dt             float64 `json:"dt"`
	Steps          int     `json:"steps"`
	X0             float64 `json:"x0"`
	P0             float64 `json:"p0"`
	Sigma          float64 `json:"sigma"`
	DriveAmplitude float64 `json:"driveAmplitude"`
	DriveFrequency float64 `json:"driveFrequency"`
	Particles      int     `json:"particles"`
	Box            float64 `json:"box"`
	Temperature    float64 `json:"temperature"`
	Seed           int     `json:"seed"`
}
type Series struct {
	Name string    `json:"name"`
	X    []float64 `json:"x"`
	Y    []float64 `json:"y"`
}
type Chart struct {
	Title  string      `json:"title"`
	XLabel string      `json:"xLabel"`
	YLabel string      `json:"yLabel"`
	Series []Series    `json:"series,omitempty"`
	X      []float64   `json:"x,omitempty"`
	Y      []float64   `json:"y,omitempty"`
	Z      [][]float64 `json:"z,omitempty"`
}
type Metric struct {
	Label string  `json:"label"`
	Value float64 `json:"value"`
}
type Frame struct {
	Time float64   `json:"time"`
	X    []float64 `json:"x"`
	Y    []float64 `json:"y"`
}
type Result struct {
	Request Request  `json:"request"`
	Title   string   `json:"title"`
	Units   string   `json:"units"`
	Metrics []Metric `json:"metrics"`
	Charts  []Chart  `json:"charts"`
	Notes   []string `json:"notes"`
	Frames  []Frame  `json:"frames,omitempty"`
}

func bounded(v, min, max float64) bool {
	return !math.IsNaN(v) && !math.IsInf(v, 0) && v >= min && v <= max
}
func Validate(r Request) error {
	if r.Kind == "md" {
		if r.Method != "verlet" || r.Particles < 4 || r.Particles > 64 || !bounded(r.Box, 5, 30) || !bounded(r.Temperature, 0, 2) || r.Seed < 0 || r.Seed > 1000000 || !bounded(r.Dt, .00001, .01) || r.Steps < 1 || r.Steps > 2000 {
			return fmt.Errorf("MD requires velocity-Verlet, 4–64 particles, box 5–30, temperature 0–2, seed 0–1,000,000, dt 0.00001–0.01, and 1–2,000 steps")
		}
		return nil
	}
	if r.Kind != "spectrum" && r.Kind != "quantum" && r.Kind != "classical" && r.Kind != "molecular" {
		return fmt.Errorf("choose a supported calculation")
	}
	if !bounded(r.Mass, .1, 100) || !bounded(r.Strength, .01, 100) || !bounded(r.Center, -10, 10) || !bounded(r.Alpha, .05, 2) {
		return fmt.Errorf("mass must be 0.1–100, strength 0.01–100, center -10–10, and Morse alpha 0.05–2")
	}
	if r.Potential != "harmonic" && r.Potential != "morse" && r.Potential != "double-well" && r.Potential != "free" {
		return fmt.Errorf("unknown potential")
	}
	if r.Kind == "molecular" && r.Potential != "morse" {
		return fmt.Errorf("molecular vibration uses the Morse bond potential")
	}
	if r.Kind == "spectrum" || r.Kind == "quantum" {
		if r.Basis != "dvr" && r.Basis != "fourier" {
			return fmt.Errorf("select DVR or Fourier basis")
		}
		if r.Points < 32 || r.Points > 192 || !bounded(r.HalfWidth, 2, 20) {
			return fmt.Errorf("use 32–192 grid points and a half-width of 2–20")
		}
		if r.States < 1 || r.States > 8 {
			return fmt.Errorf("request 1–8 eigenstates")
		}
	}
	if r.Kind != "spectrum" {
		if !bounded(r.Dt, 0.00001, .05) || r.Steps < 1 || r.Steps > 2000 {
			return fmt.Errorf("use a time step of 0.00001–0.05 and 1–2,000 steps")
		}
		if !bounded(r.X0, -10, 10) || !bounded(r.P0, -20, 20) {
			return fmt.Errorf("initial position must be -10–10 and momentum -20–20")
		}
	}
	if r.Kind == "quantum" {
		if r.Method != "split" && r.Method != "rk4" && r.Method != "ralston3" && r.Method != "nystrom5" {
			return fmt.Errorf("select a supported quantum propagator")
		}
		if !bounded(r.Sigma, .1, 4) || math.Abs(r.X0)+3*r.Sigma >= r.HalfWidth {
			return fmt.Errorf("the initial wave packet needs at least three sigma of space to each boundary")
		}
		if !bounded(r.DriveAmplitude, -5, 5) || !bounded(r.DriveFrequency, 0, 10) {
			return fmt.Errorf("drive amplitude must be -5–5 and frequency 0–10")
		}
		if r.Points*r.Points*r.Steps > 40000000 {
			return fmt.Errorf("this job is too large; reduce grid points or time steps")
		}
	} else if r.Kind != "spectrum" && r.Method != "verlet" && r.Method != "leapfrog" && r.Method != "yoshida" {
		return fmt.Errorf("select velocity-Verlet, leapfrog or Yoshida")
	}
	return nil
}

func potential(r Request) gridData.PotentialOp[float64] {
	switch r.Potential {
	case "harmonic":
		return gridData.Harmonic[float64]{Cen: r.Center, ForceConst: r.Strength}
	case "morse":
		return gridData.Morse[float64]{Cen: r.Center, De: r.Strength, Alpha: r.Alpha}
	case "double-well":
		return gridData.Polynomial[float64]{Coeffs: []float64{0, 0, -r.Strength, 0, .5 * r.Strength}}
	default:
		return gridData.Polynomial[float64]{Coeffs: []float64{0}}
	}
}

func Run(ctx context.Context, r Request) (*Result, error) {
	if err := Validate(r); err != nil {
		return nil, err
	}
	if err := ctx.Err(); err != nil {
		return nil, err
	}
	var result *Result
	var err error
	if r.Kind == "spectrum" || r.Kind == "quantum" {
		result, err = quantum(ctx, r)
	} else if r.Kind == "md" {
		result, err = molecularDynamics(ctx, r)
	} else {
		result, err = trajectory(ctx, r)
	}
	if err != nil {
		return nil, err
	}
	if err = ctx.Err(); err != nil {
		return nil, err
	}
	for _, m := range result.Metrics {
		if !bounded(m.Value, -math.MaxFloat64, math.MaxFloat64) {
			return nil, fmt.Errorf("calculation became non-finite; reduce the step size or change the initial state")
		}
	}
	for _, f := range result.Frames {
		for _, a := range [][]float64{f.X, f.Y} {
			for _, v := range a {
				if !bounded(v, -math.MaxFloat64, math.MaxFloat64) {
					return nil, fmt.Errorf("non-finite particle position")
				}
			}
		}
	}
	for _, c := range result.Charts {
		arrays := [][]float64{c.X, c.Y}
		arrays = append(arrays, c.Z...)
		for _, s := range c.Series {
			arrays = append(arrays, s.X, s.Y)
		}
		for _, a := range arrays {
			for _, v := range a {
				if !bounded(v, -math.MaxFloat64, math.MaxFloat64) {
					return nil, fmt.Errorf("calculation diverged; reduce the step size")
				}
			}
		}
	}
	return result, nil
}
