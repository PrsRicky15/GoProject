package simulation

import (
	"GoProject/classical"
	"context"
	"math"
)

func molecularDynamics(ctx context.Context, r Request) (*Result, error) {
	s, err := classical.NewLJSystem(r.Particles, r.Box, r.Temperature, uint64(r.Seed))
	if err != nil {
		return nil, err
	}
	initial := s.Kinetic() + s.Potential
	maxDrift := 0.
	times, energies, temperatures := []float64{}, []float64{}, []float64{}
	frames := []Frame{}
	stride := (r.Steps + 149) / 150
	for i := 0; i <= r.Steps; i++ {
		if err := ctx.Err(); err != nil {
			return nil, err
		}
		energy := s.Kinetic() + s.Potential
		maxDrift = math.Max(maxDrift, math.Abs(energy-initial))
		if i%stride == 0 || i == r.Steps {
			t := float64(i) * r.Dt
			times = append(times, t)
			energies = append(energies, energy)
			temperatures = append(temperatures, s.Temperature())
			frames = append(frames, Frame{Time: t, X: append([]float64(nil), s.X...), Y: append([]float64(nil), s.Y...)})
		}
		if i < r.Steps {
			if err := s.Step(r.Dt); err != nil {
				return nil, err
			}
		}
	}
	px, py := 0., 0.
	for i := range s.VX {
		px += s.VX[i]
		py += s.VY[i]
	}
	return &Result{Request: r, Title: "Lennard-Jones molecular dynamics", Units: "Reduced LJ units: mass = sigma = epsilon = kB = 1", Metrics: []Metric{{"Initial energy", initial}, {"Final energy", energies[len(energies)-1]}, {"Maximum absolute energy drift", maxDrift}, {"Total momentum magnitude", math.Hypot(px, py)}}, Charts: []Chart{{Title: "Total energy", XLabel: "Time", YLabel: "Energy", Series: []Series{{Name: "E(t)", X: times, Y: energies}}}, {Title: "Temperature", XLabel: "Time", YLabel: "Temperature", Series: []Series{{Name: "T(t)", X: times, Y: temperatures}}}}, Frames: frames, Notes: []string{"2D monatomic model with periodic boundaries, minimum-image distances, a force-shifted cutoff min(2.5, L/2), and velocity-Verlet integration.", "Seeded lattice initialization. Initial center-of-mass momentum is removed and velocities are scaled to the requested temperature using 2N−2 degrees of freedom. The subsequent run is NVE, with no thermostat.", "This is a model fluid, not a chemically parameterized molecular force field. Reduce dt and compare energy drift before interpreting results."}}, nil
}
