package Quantum

import (
	"context"
	"fmt"
	"math"
	"math/rand/v2"

	"gonum.org/v1/gonum/mat"
)

// EigenOptions controls a restarted block Krylov/Rayleigh-Ritz calculation.
// Zero fields select defaults. Time selects an instantaneous H(t), sampled once.
type EigenOptions struct {
	Time          float64
	Tolerance     float64
	MaxIterations int
	MaxSubspace   int
}

type EigenResult struct {
	Values     []float64
	Vectors    *mat.Dense // normalized column eigenvectors, in NDimGrid ordering
	Residuals  []float64  // ||H*v-lambda*v||_2, independently recomputed
	Iterations int
}

// LowestStates computes a few low eigenpairs without assembling the full H.
// A block of deterministic random starts handles degeneracies; double
// orthogonalization and thick restarts control loss of orthogonality/storage.
// Nonconvergence returns the latest approximations AND an error. Residual
// convergence is not a proof that every lower level has been found: verify
// difficult/clustered spectra by increasing the subspace and grid resolution.
func (h *NDimHamiltonian) LowestStates(ctx context.Context, count int, opts EigenOptions) (*EigenResult, error) {
	if ctx == nil {
		return nil, fmt.Errorf("LowestStates: nil context")
	}
	if err := ctx.Err(); err != nil {
		return nil, err
	}
	n := h.Dim()
	if h.core.closed || count < 1 || count > min(n, 64) || !finite(opts.Time) {
		return nil, fmt.Errorf("LowestStates: closed operator, invalid time or state count (maximum min(N,64))")
	}
	if opts.Tolerance == 0 {
		opts.Tolerance = 1e-9
	}
	if !finite(opts.Tolerance) || opts.Tolerance <= 0 {
		return nil, fmt.Errorf("LowestStates: tolerance must be finite and positive")
	}
	if opts.MaxIterations == 0 {
		opts.MaxIterations = 500
	}
	if opts.MaxIterations < 1 {
		return nil, fmt.Errorf("LowestStates: iteration limit must be positive")
	}
	block := min(n, count+2)
	if opts.MaxSubspace == 0 {
		opts.MaxSubspace = min(n, max(32, 8*block))
	}
	if opts.MaxSubspace < min(n, 2*block) || opts.MaxSubspace > n || opts.MaxSubspace > 512 {
		return nil, fmt.Errorf("LowestStates: subspace must be between %d and %d", min(n, 2*block), min(n, 512))
	}
	if int64(n)*int64(opts.MaxSubspace) > 16<<20 {
		return nil, fmt.Errorf("LowestStates: basis storage limit exceeded; reduce grid or subspace")
	}
	potential := make([]float64, n)
	h.PotentialAt(opts.Time, potential)
	for _, v := range potential {
		if !finite(v) {
			return nil, fmt.Errorf("LowestStates: non-finite potential")
		}
	}
	apply := func(x []float64) ([]float64, error) {
		if err := ctx.Err(); err != nil {
			return nil, err
		}
		y := make([]float64, n)
		h.kinetic.applyReal(x, y)
		for i := range y {
			y[i] += potential[i] * x[i]
			if !finite(y[i]) {
				return nil, fmt.Errorf("LowestStates: non-finite Hamiltonian application")
			}
		}
		return y, nil
	}
	var basis, images [][]float64
	add := func(v []float64) (bool, error) {
		for pass := 0; pass < 2; pass++ {
			for _, q := range basis {
				c := ndDot(q, v)
				for i := range v {
					v[i] -= c * q[i]
				}
			}
		}
		norm := math.Sqrt(ndDot(v, v))
		if norm <= 1e-13 {
			return false, nil
		}
		for i := range v {
			v[i] /= norm
		}
		av, err := apply(v)
		if err != nil {
			return false, err
		}
		basis, images = append(basis, v), append(images, av)
		return true, nil
	}
	rng := rand.New(rand.NewPCG(19, 37))
	randomVector := func() []float64 {
		v := make([]float64, n)
		for i := range v {
			v[i] = rng.NormFloat64()
		}
		return v
	}
	for len(basis) < block {
		if _, err := add(randomVector()); err != nil {
			return nil, err
		}
	}
	projected := mat.NewSymDense(opts.MaxSubspace, nil)
	previousSize := 0
	var result *EigenResult
	for iteration := 1; iteration <= opts.MaxIterations; iteration++ {
		if err := ctx.Err(); err != nil {
			return result, err
		}
		m := len(basis)
		for j := previousSize; j < m; j++ {
			for i := 0; i <= j; i++ {
				projected.SetSym(i, j, ndDot(basis[i], images[j]))
			}
		}
		previousSize = m
		var eig mat.EigenSym
		if !eig.Factorize(projected.SliceSym(0, m), true) {
			return result, fmt.Errorf("LowestStates: projected diagonalization failed")
		}
		values := eig.Values(nil)
		var rotations mat.Dense
		eig.VectorsTo(&rotations)
		kept := min(block, m)
		ritz, residuals := make([][]float64, kept), make([][]float64, kept)
		result = &EigenResult{Values: append([]float64(nil), values[:count]...), Vectors: mat.NewDense(n, count, nil), Residuals: make([]float64, count), Iterations: iteration}
		converged := true
		for s := 0; s < kept; s++ {
			x := make([]float64, n)
			for j, q := range basis {
				c := rotations.At(j, s)
				for i := range x {
					x[i] += c * q[i]
				}
			}
			ax, err := apply(x)
			if err != nil {
				return result, err
			}
			for i := range ax {
				ax[i] -= values[s] * x[i]
			}
			ritz[s], residuals[s] = x, ax
			if s < count {
				res := math.Sqrt(ndDot(ax, ax))
				result.Residuals[s] = res
				result.Vectors.SetCol(s, x)
				if res > opts.Tolerance*math.Max(1, math.Abs(values[s])) {
					converged = false
				}
			}
		}
		if converged {
			return result, nil
		}
		if m+kept > opts.MaxSubspace {
			basis, images = nil, nil
			previousSize = 0
			for _, v := range ritz {
				if _, err := add(v); err != nil {
					return result, err
				}
			}
		}
		before := len(basis)
		for _, r := range residuals {
			if len(basis) == opts.MaxSubspace {
				break
			}
			if _, err := add(r); err != nil {
				return result, err
			}
		}
		if len(basis) == before {
			ok, err := add(randomVector())
			if err != nil {
				return result, err
			}
			if !ok {
				return result, fmt.Errorf("LowestStates: search space stagnated before convergence")
			}
		}
	}
	return result, fmt.Errorf("LowestStates: did not converge within %d iterations (residuals %v)", opts.MaxIterations, result.Residuals)
}

func ndDot(x, y []float64) float64 {
	var sum float64
	for i, v := range x {
		sum += v * y[i]
	}
	return sum
}
