package EquationSolver

import (
	"GoProject/gridData"
	"fmt"
	"math"
)

// explicitGrid advances a coupled vector with a full explicit RK tableau.
// State is committed only after all stages and output have finite values.
func explicitGrid(f gridData.TDPotentialOp, x []float64, t, dt float64, a [][]float64, b, c []float64) error {
	if f == nil || math.IsNaN(dt) || math.IsInf(dt, 0) {
		return fmt.Errorf("invalid Runge-Kutta definition")
	}
	k := make([][]float64, len(b))
	work := make([]float64, len(x))
	for stage := range b {
		copy(work, x)
		for j, coef := range a[stage] {
			for i := range work {
				work[i] += dt * coef * k[j][i]
			}
		}
		k[stage] = make([]float64, len(x))
		f.EvaluateOnRGridInPlace(work, k[stage], t+dt*c[stage])
		for _, v := range k[stage] {
			if math.IsNaN(v) || math.IsInf(v, 0) {
				return fmt.Errorf("non-finite Runge-Kutta stage")
			}
		}
	}
	copy(work, x)
	for j, coef := range b {
		for i := range work {
			work[i] += dt * coef * k[j][i]
		}
	}
	for _, v := range work {
		if math.IsNaN(v) || math.IsInf(v, 0) {
			return fmt.Errorf("non-finite Runge-Kutta result")
		}
	}
	copy(x, work)
	return nil
}

func (r *Ralston2order) NextStepOnGrid(x []float64, t float64) error {
	return explicitGrid(r.timeFunc, x, t, r.delTime, [][]float64{{}, {2. / 3}}, []float64{1. / 4, 3. / 4}, []float64{0, 2. / 3})
}
func (r *Nystrom5Explicit) NextStepOnGrid(x []float64, t float64) error {
	return explicitGrid(r.timeFunc, x, t, r.deltaTime, [][]float64{{}, {1. / 3}, {4. / 25, 6. / 25}, {1. / 4, -3, 15. / 4}, {2. / 27, 10. / 9, -50. / 81, 8. / 81}, {2. / 25, 12. / 25, 2. / 15, 8. / 75}}, []float64{23. / 192, 0, 125. / 192, 0, -27. / 64, 125. / 192}, []float64{0, 1. / 3, 2. / 5, 1, 2. / 3, 4. / 5})
}
