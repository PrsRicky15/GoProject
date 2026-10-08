package simulation

import (
	EquationSolver "GoProject/ODESolver"
	"GoProject/Quantum"
	"GoProject/gridData"
	"context"
	"fmt"
	"math"
	"math/cmplx"
)

// Store real and imaginary parts in a coupled ODE vector for existing RK solvers.
type schrodinger struct {
	h         *Quantum.HamiltonianOp
	psi, hpsi []complex128
}

func (s *schrodinger) EvaluateAt(float64, float64) float64 {
	panic("Schrodinger evolution requires the coupled grid operation")
}
func (s *schrodinger) EvaluateOnRGrid(x []float64, t float64) []float64 {
	out := make([]float64, len(x))
	s.EvaluateOnRGridInPlace(x, out, t)
	return out
}
func (s *schrodinger) EvaluateOnRGridInPlace(x, out []float64, t float64) {
	n := s.h.Dim()
	for i := 0; i < n; i++ {
		s.psi[i] = complex(x[i], x[n+i])
	}
	s.h.ApplyHComplexAt(s.psi, s.hpsi, t)
	for i, v := range s.hpsi {
		out[i] = imag(v)
		out[n+i] = -real(v)
	}
}

func quantum(ctx context.Context, r Request) (*Result, error) {
	grid, err := gridData.NewRGrid(-r.HalfWidth, r.HalfWidth, uint32(r.Points))
	if err != nil {
		return nil, err
	}
	var h *Quantum.HamiltonianOp
	if r.Basis == "fourier" {
		h, err = Quantum.NewFourierHamiltonian(grid, r.Mass, potential(r))
	} else {
		h, err = Quantum.NewHamiltonian(grid, r.Mass, potential(r))
	}
	if err != nil {
		return nil, err
	}
	defer h.Close()
	x := grid.RValues()
	dx := grid.DeltaR()
	v0 := h.VStatic()
	for _, v := range v0 {
		if !bounded(v, -1e6, 1e6) {
			return nil, fmt.Errorf("potential is too steep; reduce the domain, Morse alpha or strength")
		}
	}
	result := &Result{Request: r, Units: "Atomic units (hbar = 1)", Metrics: []Metric{}, Charts: []Chart{}, Notes: []string{"Grid uses [-L, L), with spacing 2L/N. Fourier boundaries are periodic; DVR uses the full-line sinc kinetic matrix truncated to this symmetric domain."}}
	if r.Kind == "spectrum" {
		values, vectors, err := h.RealDiagonalize()
		if err != nil {
			return nil, err
		}
		chart := Chart{Title: "Eigenstate probability densities", XLabel: "Position", YLabel: "Probability density"}
		maxResidual := 0.
		for j := 0; j < r.States; j++ {
			if err := ctx.Err(); err != nil {
				return nil, err
			}
			psi, density, hpsi := make([]float64, r.Points), make([]float64, r.Points), make([]float64, r.Points)
			for i := range psi {
				psi[i] = vectors.At(i, j)
				density[i] = psi[i] * psi[i] / dx
			}
			h.ApplyHReal(psi, hpsi)
			sum := 0.
			for i, v := range psi {
				d := hpsi[i] - values[j]*v
				sum += d * d
			}
			maxResidual = math.Max(maxResidual, math.Sqrt(sum))
			chart.Series = append(chart.Series, Series{Name: fmt.Sprintf("n = %d", j), X: x, Y: density})
			result.Metrics = append(result.Metrics, Metric{fmt.Sprintf("Energy E%d", j), values[j]})
		}
		result.Title = "Hamiltonian energy levels"
		result.Metrics = append(result.Metrics, Metric{"Largest eigenpair residual", maxResidual})
		result.Charts = append(result.Charts, chart, Chart{Title: "Potential", XLabel: "Position", YLabel: "Energy", Series: []Series{{Name: "V(x)", X: x, Y: v0}}})
		result.Notes = append(result.Notes, "Eigenstates are normalized so sum(|psi|²) dx = 1. Finite-domain states above a dissociation threshold are box states, not necessarily bound molecular states. Increase N and L to check convergence.")
		return result, nil
	}
	if r.DriveAmplitude != 0 {
		h.SetTimeDependentPotential(func(t float64, out []float64) {
			for i, xi := range x {
				out[i] = v0[i] + r.DriveAmplitude*xi*math.Sin(r.DriveFrequency*t)
			}
		})
	}
	psi := make([]complex128, r.Points)
	norm := 0.
	for i, xi := range x {
		psi[i] = complex(math.Exp(-math.Pow((xi-r.X0)/r.Sigma, 2)/2), 0) * cmplx.Exp(complex(0, r.P0*xi))
		norm += cmplx.Abs(psi[i]) * cmplx.Abs(psi[i]) * dx
	}
	for i := range psi {
		psi[i] /= complex(math.Sqrt(norm), 0)
	}
	var step func(float64) error
	if r.Method == "split" {
		prop, err := Quantum.NewSplitPropagator(h)
		if err != nil {
			return nil, err
		}
		step = func(t float64) error { return prop.Step(psi, t, r.Dt) }
	} else {
		// Conservative spectral bound keeps explicit methods in a useful regime.
		matrix := h.MatrixAt(0)
		bound := 0.
		for i := 0; i < r.Points; i++ {
			row := 0.
			for j := 0; j < r.Points; j++ {
				row += math.Abs(matrix.At(i, j))
			}
			bound = math.Max(bound, row)
		}
		bound += math.Abs(r.DriveAmplitude) * r.HalfWidth
		limit := 1.
		if r.Method == "rk4" {
			limit = 2.
		}
		if r.Dt*bound > limit {
			return nil, fmt.Errorf("time step too large for %s; use dt <= %.6g or choose split-operator", r.Method, limit/bound)
		}
		rhs := &schrodinger{h: h, psi: make([]complex128, r.Points), hpsi: make([]complex128, r.Points)}
		var solver interface {
			NextStepOnGrid([]float64, float64) error
		}
		switch r.Method {
		case "rk4":
			solver = (&EquationSolver.RungeKutta4Explicit{}).NewDef(r.Dt, rhs)
		case "ralston3":
			solver = (&EquationSolver.Ralston3Order{}).NewDefine(r.Dt, rhs)
		case "nystrom5":
			solver = (&EquationSolver.Nystrom5Explicit{}).NewDefine(r.Dt, rhs)
		}
		state := make([]float64, 2*r.Points)
		step = func(t float64) error {
			for i, z := range psi {
				state[i] = real(z)
				state[r.Points+i] = imag(z)
			}
			if err := solver.NextStepOnGrid(state, t); err != nil {
				return err
			}
			for i := range psi {
				psi[i] = complex(state[i], state[r.Points+i])
			}
			return nil
		}
	}
	times, norms, means, energies := []float64{}, []float64{}, []float64{}, []float64{}
	densities := [][]float64{}
	stride := (r.Steps + 149) / 150
	maxNormDrift := 0.
	initialEnergy := h.EnergyAt(psi, 0)
	for i := 0; i <= r.Steps; i++ {
		if err := ctx.Err(); err != nil {
			return nil, err
		}
		t := float64(i) * r.Dt
		norm, mean := 0., 0.
		for j, z := range psi {
			prob := real(z)*real(z) + imag(z)*imag(z)
			norm += prob * dx
			mean += x[j] * prob * dx
		}
		if !bounded(norm, .5, 1.5) {
			return nil, fmt.Errorf("wavefunction norm drifted excessively; reduce dt")
		}
		maxNormDrift = math.Max(maxNormDrift, math.Abs(norm-1))
		if i%stride == 0 || i == r.Steps {
			density := make([]float64, r.Points)
			for j, z := range psi {
				density[j] = real(z)*real(z) + imag(z)*imag(z)
			}
			times = append(times, t)
			norms = append(norms, norm)
			means = append(means, mean/norm)
			energies = append(energies, h.EnergyAt(psi, t))
			densities = append(densities, density)
		}
		if i < r.Steps {
			if err := step(t); err != nil {
				return nil, err
			}
		}
	}
	result.Title = "Wave-packet propagation"
	result.Metrics = []Metric{{"Initial energy", initialEnergy}, {"Final energy", energies[len(energies)-1]}, {"Maximum norm drift", maxNormDrift}, {"Final mean position", means[len(means)-1]}}
	result.Charts = []Chart{{Title: "Probability density over time", XLabel: "Position", YLabel: "Time", X: x, Y: times, Z: densities}, {Title: "Mean position", XLabel: "Time", YLabel: "Mean position", Series: []Series{{Name: "<x>", X: times, Y: means}}}, {Title: "Wavefunction norm", XLabel: "Time", YLabel: "Norm", Series: []Series{{Name: "Norm", X: times, Y: norms}}}, {Title: "Energy expectation", XLabel: "Time", YLabel: "Energy", Series: []Series{{Name: "<H(t)>", X: times, Y: energies}}}}
	result.Notes = append(result.Notes, "Initial state: exp(-(x-x0)²/(2 sigma²)) exp(i p0 x), normalized on the grid. Drive: A x sin(omega t). A driven Hamiltonian need not conserve energy.", "No absorbing boundary is applied. Watch for reflections or periodic wrapping. No stepwise renormalization hides numerical error. Compare smaller dt and larger grids.")
	return result, nil
}
