package simulation

import (
	EquationSolver "GoProject/ODESolver"
	"GoProject/gridData"
	"context"
	"fmt"
	"math"
)

// Legacy symplectic integrators expect acceleration in ForceAt; adapt F/m here.
type acceleration struct {
	gridData.PotentialOp[float64]
	mass float64
}

func (a acceleration) ForceAt(x float64) float64 { return a.PotentialOp.ForceAt(x) / a.mass }
func (a acceleration) ForceOnGrid(x []float64) []float64 {
	out := make([]float64, len(x))
	for i, v := range x {
		out[i] = a.ForceAt(v)
	}
	return out
}

func trajectory(ctx context.Context, r Request) (*Result, error) {
	pot := potential(r)
	acc := acceleration{pot, r.Mass}
	var step func(float64, float64) (float64, float64)
	switch r.Method {
	case "verlet":
		step = (&EquationSolver.VelocityVerlet{}).NewDef(acc, r.Dt).NextStep
	case "leapfrog":
		step = (&EquationSolver.LeapFrog{}).NewDef(acc, r.Dt).NextStep
	case "yoshida":
		step = (&EquationSolver.Yoshida{}).NewDef(acc, r.Dt).NextStep
	}
	x, v := r.X0, r.P0/r.Mass
	initialEnergy := .5*r.Mass*v*v + pot.EvaluateAt(x)
	if !bounded(initialEnergy, -1e8, 1e8) {
		return nil, fmt.Errorf("initial energy is too large for this interactive calculation")
	}
	stride := (r.Steps + 399) / 400
	times, positions, momenta, energies := []float64{}, []float64{}, []float64{}, []float64{}
	maxDrift := 0.
	for i := 0; i <= r.Steps; i++ {
		if err := ctx.Err(); err != nil {
			return nil, err
		}
		energy := .5*r.Mass*v*v + pot.EvaluateAt(x)
		if !bounded(energy, -1e12, 1e12) || !bounded(x, -1e6, 1e6) {
			return nil, fmt.Errorf("trajectory diverged; reduce the step size or initial energy")
		}
		maxDrift = math.Max(maxDrift, math.Abs(energy-initialEnergy))
		if i%stride == 0 || i == r.Steps {
			times = append(times, float64(i)*r.Dt)
			positions = append(positions, x)
			momenta = append(momenta, r.Mass*v)
			energies = append(energies, energy)
		}
		if i < r.Steps {
			x, v = step(x, v)
		}
	}
	result := &Result{Request: r, Title: "Classical trajectory", Units: "Atomic units (hbar = 1)", Metrics: []Metric{{"Initial energy", initialEnergy}, {"Final energy", energies[len(energies)-1]}, {"Maximum absolute energy drift", maxDrift}, {"Final position", x}}, Charts: []Chart{
		{Title: "Position over time", XLabel: "Time", YLabel: "Position", Series: []Series{{Name: "x(t)", X: times, Y: positions}}},
		{Title: "Phase space", XLabel: "Position", YLabel: "Momentum", Series: []Series{{Name: "Trajectory", X: positions, Y: momenta}}},
		{Title: "Total energy", XLabel: "Time", YLabel: "Energy", Series: []Series{{Name: "E(t)", X: times, Y: energies}}},
	}, Notes: []string{"Uses the repository's symplectic integrator with acceleration F/m. Output is sampled to at most 401 frames; all requested integration steps are executed.", "Decrease dt and compare energy drift to check convergence. No thermostat or artificial energy correction is applied."}}
	if r.Kind == "molecular" {
		result.Title = "Morse bond vibration"
		result.Notes = append(result.Notes, "One-dimensional diatomic bond coordinate, with mass interpreted as reduced mass. This is not a many-particle molecular dynamics simulation.")
	}
	return result, nil
}
