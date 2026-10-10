package simulation

import (
	"context"
	"fmt"
	"math"

	"GoProject/Quantum"
	"GoProject/gridData"
	"GoProject/internal/plotmath"
)

type NDAxis struct {
	Points    int     `json:"points"`
	HalfWidth float64 `json:"halfWidth"`
	Mass      float64 `json:"mass"`
	Center    float64 `json:"center"`
	Scale     float64 `json:"scale"`
}
type NDRequest struct {
	Expression string   `json:"expression,omitempty"`
	Axes       []NDAxis `json:"axes"`
	Coupling   float64  `json:"coupling"`
	Softening  float64  `json:"softening"`
	Separation float64  `json:"separation"`
	Wavenumber float64  `json:"wavenumber"`
	Order      int      `json:"order"`
	State      int      `json:"state"`
	AxisX      int      `json:"axisX"`
	AxisY      int      `json:"axisY"`
}
type DensityVolume struct {
	Title  string    `json:"title"`
	Labels []string  `json:"labels"`
	X      []float64 `json:"x"`
	Y      []float64 `json:"y"`
	Z      []float64 `json:"z"`
	Values []float64 `json:"values"`
}

const maxNDWebPoints = 8192

func validateND(r Request) error {
	n := r.ND
	if n == nil || len(n.Axes) < 2 || len(n.Axes) > 4 {
		return fmt.Errorf("choose 2–4 coordinates for web calculations")
	}
	if r.Basis != "dvr" && r.Basis != "fourier" {
		return fmt.Errorf("select DVR or Fourier basis")
	}
	if r.States < 1 || r.States > 4 || n.State < 0 || n.State >= r.States {
		return fmt.Errorf("request 1–4 eigenstates and select a density state within that range")
	}
	if n.AxisX < 0 || n.AxisY < 0 || n.AxisX >= len(n.Axes) || n.AxisY >= len(n.Axes) || n.AxisX == n.AxisY {
		return fmt.Errorf("select two different projection axes")
	}
	total := 1
	for d, a := range n.Axes {
		if a.Points < 4 || a.Points > 64 || !bounded(a.HalfWidth, 2, 12) || !bounded(a.Mass, .1, 10) || !bounded(a.Center, -a.HalfWidth/2, a.HalfWidth/2) || !bounded(a.Scale, .25, 2) {
			return fmt.Errorf("axis %d requires 4–64 points, half-width 2–12, mass 0.1–10, center within half the half-width, and scale 0.25–2", d+1)
		}
		total *= a.Points
		if total > maxNDWebPoints {
			return fmt.Errorf("the full tensor grid exceeds %d points; reduce points per axis", maxNDWebPoints)
		}
	}
	if r.Potential == "custom" {
		_, err := plotmath.CompileND(n.Expression, len(n.Axes))
		if err != nil {
			return fmt.Errorf("custom potential: %w", err)
		}
		return nil
	}
	if !bounded(r.Strength, .01, 30) || !bounded(r.Alpha, .05, 1) || !bounded(r.Sigma, .3, 4) || !bounded(n.Coupling, -.45, .45) || !bounded(n.Softening, .2, 3) || !bounded(n.Separation, .2, 6) || !bounded(n.Wavenumber, .1, 3) || (n.Order != 2 && n.Order != 4 && n.Order != 6 && n.Order != 8) {
		return fmt.Errorf("invalid potential parameters; check strength, width, coupling, softening, separation or order")
	}
	switch r.Potential {
	case "harmonic", "coupled-harmonic", "morse", "double-well", "gaussian", "super-gaussian", "soft-core", "double-gaussian", "lattice", "free":
		return nil
	default:
		return fmt.Errorf("unknown multidimensional potential")
	}
}

// The closure owns reusable q[d]=scale[d]*(x[d]-center[d]) and is job-local.
func ndPotential(r Request) gridData.NDimPotentialOp {
	n := r.ND
	q := make([]float64, len(n.Axes))
	var base gridData.NDimPotentialOp
	switch r.Potential {
	case "harmonic":
		base = gridData.Harmonic[float64]{ForceConst: r.Strength}
	case "morse":
		base = gridData.Morse[float64]{De: r.Strength, Alpha: r.Alpha}
	case "double-well":
		base = gridData.Polynomial[float64]{Coeffs: []float64{r.Strength / 4, 0, -r.Strength / 2, 0, r.Strength / 4}}
	case "gaussian":
		base = gridData.Gaussian[float64]{Strength: -r.Strength, Sigma: r.Sigma}
	case "super-gaussian":
		base = gridData.SuperGaussian[float64]{Strength: -r.Strength, Sigma: r.Sigma, Order: uint8(n.Order)}
	case "soft-core":
		base = gridData.SoftCore[float64]{Charge: -r.Strength, SoftParam: n.Softening}
	}
	gauss := gridData.Gaussian[float64]{Strength: -r.Strength, Sigma: r.Sigma}
	return gridData.NDimPotentialFunc(func(x []float64) float64 {
		for d, a := range n.Axes {
			q[d] = a.Scale * (x[d] - a.Center)
		}
		if base != nil {
			return base.NDimEvaluateAt(q)
		}
		switch r.Potential {
		case "coupled-harmonic":
			v := 0.0
			for d, qi := range q {
				v += .5 * r.Strength * qi * qi
				if d > 0 {
					v += r.Strength * n.Coupling * q[d-1] * qi
				}
			}
			return v
		case "double-gaussian":
			q[0] -= n.Separation / 2
			v := gauss.NDimEvaluateAt(q)
			q[0] += n.Separation
			return v + gauss.NDimEvaluateAt(q)
		case "lattice":
			v := 0.0
			for _, qi := range q {
				s := math.Sin(n.Wavenumber * qi)
				v += r.Strength * s * s
			}
			return v
		default:
			return 0
		}
	})
}

func ndSpectrum(ctx context.Context, r Request) (*Result, error) {
	n := r.ND
	dims := len(n.Axes)
	lengths, masses, points := make([]float64, dims), make([]float64, dims), make([]uint32, dims)
	for d, a := range n.Axes {
		lengths[d] = 2 * a.HalfWidth
		masses[d] = a.Mass
		points[d] = uint32(a.Points)
	}
	g, err := gridData.NewNDimGridFromLengths(lengths, points)
	if err != nil {
		return nil, err
	}
	pot := ndPotential(r)
	if r.Potential == "custom" {
		evaluate, err := plotmath.CompileND(n.Expression, dims)
		if err != nil {
			return nil, fmt.Errorf("custom potential: %w", err)
		}
		pot = gridData.NDimPotentialFunc(evaluate)
	}
	x := make([]float64, dims)
	// Reject steep models before allocating FFT plans or starting the eigensolver.
	for i := uint64(0); i < g.TotalPoints(); i++ {
		if i%256 == 0 {
			if err := ctx.Err(); err != nil {
				return nil, err
			}
		}
		g.PointTo(i, x)
		if !bounded(pot.NDimEvaluateAt(x), -1e5, 1e5) {
			if r.Potential == "custom" {
				return nil, fmt.Errorf("custom potential is undefined or exceeds |V| = 100000 at coordinates %v; check divisions, function domains and coefficients", x)
			}
			return nil, fmt.Errorf("potential is too steep; reduce strength, coordinate scale, Morse alpha or the domain")
		}
	}
	var h *Quantum.NDimHamiltonian
	if r.Basis == "fourier" {
		h, err = Quantum.NewNDimFourierHamiltonian(g, masses, pot)
	} else {
		h, err = Quantum.NewNDimHamiltonian(g, masses, pot)
	}
	if err != nil {
		return nil, err
	}
	defer h.Close()
	// A wider restarted subspace resolves clustered lattice levels without
	// relaxing the residual tolerance. The web tensor cap bounds basis storage.
	eig, err := h.LowestStates(ctx, r.States, Quantum.EigenOptions{Tolerance: 1e-7, MaxIterations: 500, MaxSubspace: min(h.Dim(), 96)})
	if err != nil {
		return nil, err
	}
	result := &Result{Request: r, Title: fmt.Sprintf("%dD Hamiltonian energy levels", dims), Units: "Atomic units (hbar = 1)", Metrics: []Metric{}, Charts: []Chart{}, Notes: []string{
		"Energies and densities are computed in Go on the full tensor grid. Coordinate scale means q[d] = scale[d] × (x[d] − center[d]).",
		"Density plots integrate over hidden coordinates. The potential plot is a slice with hidden coordinates fixed at their centers, not an integral.",
		"Fourier boundaries are periodic. Sinc DVR is a truncated whole-line basis. Increase resolution and domain to check convergence; bounded web grids can be coarse, especially in 4D.",
		"Small eigenpair residuals verify the discrete calculation, not spatial convergence. Finite wells, free particles and lattices have box-dependent states; these are not automatically bound molecular states.",
	}}
	if r.Potential == "custom" {
		result.Notes[0] = "Custom potential V = " + n.Expression + ". Go evaluates x1, x2, … as the actual grid coordinates; preset strength, scale and center transformations are not applied."
	}
	maxResidual := 0.0
	for s, e := range eig.Values {
		result.Metrics = append(result.Metrics, Metric{fmt.Sprintf("Energy E%d", s), e})
		maxResidual = math.Max(maxResidual, eig.Residuals[s])
	}
	result.Metrics = append(result.Metrics, Metric{"Largest eigenpair residual", maxResidual}, Metric{"Tensor grid points", float64(g.TotalPoints())}, Metric{"Solver iterations", float64(eig.Iterations)})
	ax, ay := uint32(n.AxisX), uint32(n.AxisY)
	xs, _ := g.RValues(ax)
	ys, _ := g.RValues(ay)
	strides := g.Strides()
	marginal := Chart{Title: fmt.Sprintf("Eigenstate densities along x%d", ax+1), XLabel: fmt.Sprintf("x%d", ax+1), YLabel: "Marginal probability density"}
	heat := Chart{Title: fmt.Sprintf("State %d · marginal density", n.State), XLabel: fmt.Sprintf("x%d", ax+1), YLabel: fmt.Sprintf("x%d", ay+1), X: xs, Y: ys, Z: make([][]float64, len(ys))}
	cut := Chart{Title: "Potential slice · hidden coordinates at centers", XLabel: heat.XLabel, YLabel: heat.YLabel, X: xs, Y: ys, Z: make([][]float64, len(ys)), ValueLabel: "V"}
	for j := range ys {
		heat.Z[j] = make([]float64, len(xs))
		cut.Z[j] = make([]float64, len(xs))
	}
	for s := range eig.Values {
		density := make([]float64, len(xs))
		for i := 0; i < h.Dim(); i++ {
			ix := int(uint64(i) / strides[ax] % uint64(len(xs)))
			v := eig.Vectors.At(i, s)
			prob := v * v
			density[ix] += prob / g.DeltaR(ax)
			if s == n.State {
				iy := int(uint64(i) / strides[ay] % uint64(len(ys)))
				heat.Z[iy][ix] += prob / (g.DeltaR(ax) * g.DeltaR(ay))
			}
		}
		marginal.Series = append(marginal.Series, Series{Name: fmt.Sprintf("n = %d", s), X: xs, Y: density})
	}
	for d, a := range n.Axes {
		x[d] = a.Center
	}
	for iy, y := range ys {
		for ix, xval := range xs {
			x[ax], x[ay] = xval, y
			cut.Z[iy][ix] = pot.NDimEvaluateAt(x)
			if !bounded(cut.Z[iy][ix], -1e5, 1e5) {
				return nil, fmt.Errorf("potential slice is undefined or exceeds |V| = 100000 at coordinates %v; adjust the formula or hidden-coordinate slice positions", x)
			}
		}
	}
	result.Charts = append(result.Charts, marginal, heat, cut)
	if dims >= 3 {
		result.Volume = ndDensityVolume(g, eig, n.State)
		if dims == 4 {
			result.Notes = append(result.Notes, "The 3D view integrates over x4 and displays x1, x2, x3; it is not a four-dimensional geometric plot.")
		}
	}
	return result, nil
}

func ndDensityVolume(g *gridData.NDimGrid, eig *Quantum.EigenResult, state int) *DensityVolume {
	ns := g.NPoints()
	n3 := int(ns[0] * ns[1] * ns[2])
	hidden := int(g.TotalPoints()) / n3
	v := &DensityVolume{Title: fmt.Sprintf("State %d · 3D marginal density", state), Labels: []string{"x1", "x2", "x3"}, X: make([]float64, n3), Y: make([]float64, n3), Z: make([]float64, n3), Values: make([]float64, n3)}
	dv := g.DeltaR(0) * g.DeltaR(1) * g.DeltaR(2)
	x := make([]float64, g.NDim())
	for i := 0; i < n3; i++ {
		g.PointTo(uint64(i*hidden), x)
		v.X[i], v.Y[i], v.Z[i] = x[0], x[1], x[2]
		for j := 0; j < hidden; j++ {
			psi := eig.Vectors.At(i*hidden+j, state)
			v.Values[i] += psi * psi / dv
		}
	}
	return v
}
