// ndquantum solves an n-dimensional coupled harmonic model using the Go library.
// Run from the module root: go run ./cmd/ndquantum -dims 3 -points 20 -states 4
package main

import (
	"context"
	"encoding/csv"
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"os"
	"os/signal"
	"strconv"
	"strings"
	"time"

	"GoProject/Quantum"
	"GoProject/gridData"
)

func main() {
	if err := run(); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

func run() error {
	dims := flag.Int("dims", 3, "number of Cartesian coordinates")
	points := flag.String("points", "20", "points per axis: one value or a comma-separated list")
	lengths := flag.String("lengths", "12", "box lengths: one value or a comma-separated list")
	massInput := flag.String("masses", "1", "masses: one value or a comma-separated list")
	frequencies := flag.String("frequencies", "1", "uncoupled harmonic frequencies: one value or a comma-separated list")
	coupling := flag.Float64("coupling", 0, "coefficient c in c*sum_d x[d]*x[d+1] (open coordinate chain)")
	basis := flag.String("basis", "fourier", "kinetic basis: fourier or dvr")
	states := flag.Int("states", 4, "number of lowest states")
	tolerance := flag.Float64("tol", 1e-9, "eigenpair relative residual tolerance")
	subspace := flag.Int("subspace", 0, "maximum search subspace (0 selects default)")
	iterations := flag.Int("iterations", 500, "maximum eigensolver iterations")
	timeout := flag.Duration("timeout", 2*time.Minute, "calculation timeout")
	csvPath := flag.String("csv", "", "optional output file containing coordinates and eigenstate densities")
	flag.Parse()
	if *dims < 1 || *dims > 32 {
		return fmt.Errorf("example CLI supports 1..32 coordinates; the library has no fixed dimension limit")
	}
	nsFloat, err := expand(*points, *dims)
	if err != nil {
		return fmt.Errorf("points: %w", err)
	}
	ls, err := expand(*lengths, *dims)
	if err != nil {
		return fmt.Errorf("lengths: %w", err)
	}
	ms, err := expand(*massInput, *dims)
	if err != nil {
		return fmt.Errorf("masses: %w", err)
	}
	ws, err := expand(*frequencies, *dims)
	if err != nil {
		return fmt.Errorf("frequencies: %w", err)
	}
	ns := make([]uint32, *dims)
	for d, n := range nsFloat {
		if !(n >= 1 && n <= Quantum.NDimMaxPoints) || n != float64(uint32(n)) {
			return fmt.Errorf("points must be positive integers within the grid limit")
		}
		ns[d] = uint32(n)
		if !(ws[d] > 0) {
			return fmt.Errorf("frequencies must be positive")
		}
	}
	if *timeout <= 0 {
		return fmt.Errorf("timeout must be positive")
	}
	g, err := gridData.NewNDimGridFromLengths(ls, ns)
	if err != nil {
		return err
	}
	pot := gridData.NDimPotentialFunc(func(x []float64) float64 {
		v := 0.0
		for d, xi := range x {
			v += .5 * ms[d] * ws[d] * ws[d] * xi * xi
			if d > 0 {
				v += *coupling * x[d-1] * xi
			}
		}
		return v
	})
	var h *Quantum.NDimHamiltonian
	switch *basis {
	case "fourier":
		h, err = Quantum.NewNDimFourierHamiltonian(g, ms, pot)
	case "dvr":
		h, err = Quantum.NewNDimHamiltonian(g, ms, pot)
	default:
		return fmt.Errorf("basis must be fourier or dvr")
	}
	if err != nil {
		return err
	}
	defer h.Close()
	ctx, stop := signal.NotifyContext(context.Background(), os.Interrupt)
	defer stop()
	ctx, cancel := context.WithTimeout(ctx, *timeout)
	defer cancel()
	result, err := h.LowestStates(ctx, *states, Quantum.EigenOptions{Tolerance: *tolerance, MaxSubspace: *subspace, MaxIterations: *iterations})
	if err != nil {
		return err
	}
	if *csvPath != "" {
		if err := writeDensities(*csvPath, g, result); err != nil {
			return err
		}
	}
	encoder := json.NewEncoder(os.Stdout)
	encoder.SetIndent("", "  ")
	return encoder.Encode(struct {
		Dimensions int       `json:"dimensions"`
		Points     []uint32  `json:"pointsPerAxis"`
		Total      uint64    `json:"totalPoints"`
		Basis      string    `json:"basis"`
		Energies   []float64 `json:"energies"`
		Residuals  []float64 `json:"residuals"`
		Iterations int       `json:"iterations"`
	}{*dims, ns, g.TotalPoints(), *basis, result.Values, result.Residuals, result.Iterations})
}

func expand(s string, n int) ([]float64, error) {
	fields := strings.Split(s, ",")
	if len(fields) != 1 && len(fields) != n {
		return nil, fmt.Errorf("expected one value or %d values", n)
	}
	out := make([]float64, n)
	for d := range out {
		i := d
		if len(fields) == 1 {
			i = 0
		}
		value, err := strconv.ParseFloat(strings.TrimSpace(fields[i]), 64)
		if err != nil {
			return nil, err
		}
		out[d] = value
	}
	return out, nil
}

func writeDensities(path string, g *gridData.NDimGrid, result *Quantum.EigenResult) (err error) {
	dv := g.VolumeElement()
	if dv <= 0 || math.IsInf(dv, 0) || math.IsNaN(dv) {
		return fmt.Errorf("density export: unrepresentable volume element; rescale coordinate units")
	}
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer func() {
		if e := f.Close(); err == nil {
			err = e
		}
	}()
	w := csv.NewWriter(f)
	row := make([]string, int(g.NDim())+len(result.Values))
	for d := 0; d < int(g.NDim()); d++ {
		row[d] = fmt.Sprintf("x%d", d)
	}
	for s := range result.Values {
		row[int(g.NDim())+s] = fmt.Sprintf("density%d", s)
	}
	if err := w.Write(row); err != nil {
		return err
	}
	x := make([]float64, g.NDim())
	for i := 0; i < int(g.TotalPoints()); i++ {
		g.PointTo(uint64(i), x)
		for d, xi := range x {
			row[d] = strconv.FormatFloat(xi, 'g', 17, 64)
		}
		for s := range result.Values {
			v := result.Vectors.At(i, s)
			row[len(x)+s] = strconv.FormatFloat(v*v/dv, 'g', 17, 64)
		}
		if err := w.Write(row); err != nil {
			return err
		}
	}
	w.Flush()
	return w.Error()
}
