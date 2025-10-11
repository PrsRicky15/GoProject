package mathsFunc

import (
	"GoProject/gridData"
	"fmt"
	"math"
	"slices"

	"gonum.org/v1/gonum/blas/blas64"
)

const (
	maxIter   = 20
	tolerance = 1e-7
	delta     = 1e-6
)

// ODESolverExplicit interface for different solving methods
// dx(t)/dt = f(x,t)
type ODESolverExplicit interface {
	NextStepEx(x, time float64) float64
	NextStepExOnGrid(x []float64, time float64)
}

// ODESolverImplicit interface for different solving methods
// dx(t)/dt = f(x,t)
type ODESolverImplicit interface {
	NextStepIm(x, time float64) (float64, error)
	NextStepImOnGrid(x []float64, time float64) error
}

type ODESolverNewtonIm interface {
	NextStepImNewton(x, time float64) (float64, error)
	NextStepImOnGridNewton(x []float64, time float64) error
}

// EulerSolver implements the basic Euler method
type EulerSolver struct {
	TdFunc gridData.TDPotentialOp
	DeltaT float64
	fxt    []float64
	xPre   []float64
	xNew   []float64
	diff   []float64
}

func (e *EulerSolver) Name() string {
	return "Euler Method"
}

func (e *EulerSolver) NextStepEx(xt, t float64) float64 {
	fxt := (e.TdFunc).EvaluateAt(xt, t)
	return xt + e.DeltaT*fxt
}

func (e *EulerSolver) NextStepIm(xt, t float64) (float64, error) {
	xPre := xt

	for iter := 0; iter < maxIter; iter++ {
		fxt := e.TdFunc.EvaluateAt(xPre, t+e.DeltaT)
		xNew := xt + e.DeltaT*fxt

		if math.Abs(xNew-xPre) <= tolerance {
			return xNew, nil
		}
		xPre = xNew
	}

	return xt, fmt.Errorf("implicit step did not converge after %d iterations", maxIter)
}

func (e *EulerSolver) NextStepExOnGrid(xt []float64, t float64) {
	nPoints := len(xt)
	if nPoints != len(e.fxt) {
		e.fxt = make([]float64, len(xt))
	}

	(e.TdFunc).EvaluateOnRGridInPlace(xt, e.fxt, t)
	blas64.Axpy(e.DeltaT,
		blas64.Vector{N: nPoints, Data: e.fxt, Inc: 1},
		blas64.Vector{N: nPoints, Data: xt, Inc: 1},
	)
}

func (e *EulerSolver) allocate(nPoints int) {
	if nPoints != len(e.fxt) ||
		nPoints != len(e.xNew) ||
		nPoints != len(e.xPre) ||
		nPoints != len(e.diff) {
		e.fxt = make([]float64, nPoints)
		e.xPre = make([]float64, nPoints)
		e.xNew = make([]float64, nPoints)
		e.diff = make([]float64, nPoints)
	}
}

func (e *EulerSolver) NextStepImOnGrid(xt []float64, t float64) error {
	nPoints := len(xt)
	e.allocate(nPoints)
	e.xPre = slices.Clone(xt)

	for iter := 0; iter < maxIter; iter++ {
		(e.TdFunc).EvaluateOnRGridInPlace(e.xPre, e.fxt, t+e.DeltaT)

		for i := range e.xNew {
			e.xNew[i] = xt[i] + e.DeltaT*e.fxt[i]
			e.diff[i] = e.xNew[i] - e.xPre[i]
		}

		err := blas64.Nrm2(blas64.Vector{N: nPoints, Data: e.diff, Inc: 1})

		if err <= tolerance {
			copy(xt, e.xNew)
			return nil
		}
		e.xPre = slices.Clone(e.xNew)
	}
	return fmt.Errorf("implicit step did not converge after %d iterations", maxIter)
}

func (e *EulerSolver) NextStepImNewton(xt, t float64) (float64, error) {
	xNext := xt
	for iter := 0; iter < maxIter; iter++ {
		fxt := e.TdFunc.EvaluateAt(xNext, t+e.DeltaT)
		Gxt := xNext - (xt + e.DeltaT*fxt)

		if math.Abs(Gxt) < tolerance {
			return xNext, nil
		}

		derGxt := (e.TdFunc.EvaluateAt(xNext+delta, t+e.DeltaT) - fxt) / delta
		jacobian := 1 - e.DeltaT*derGxt

		if math.Abs(jacobian) < 1e-14 {
			xNext = xt + e.DeltaT*fxt
			continue
		}

		xNext -= Gxt / jacobian
	}
	return xt, fmt.Errorf("implicit step did not converge after %d iterations", maxIter)
}

func (e *EulerSolver) NextStepImOnGridNewton(xt []float64, t float64) error {
	nPoints := len(xt)
	e.allocate(nPoints)
	e.xNew = slices.Clone(xt)
	for iter := 0; iter < maxIter; iter++ {
		e.TdFunc.EvaluateOnRGridInPlace(e.xNew, e.fxt, t+e.DeltaT)

		e.xPre = slices.Clone(xt)
		blas64.Axpy(e.DeltaT,
			blas64.Vector{N: nPoints, Data: e.fxt, Inc: 1},
			blas64.Vector{N: nPoints, Data: e.xPre, Inc: 1},
		)

		for i := range e.xNew {
			e.diff[i] = e.xNew[i] - e.xPre[i]
		}

		residualNorm := blas64.Nrm2(blas64.Vector{N: nPoints, Data: e.diff, Inc: 1})
		if residualNorm < tolerance {
			copy(xt, e.xNew)
			return nil
		}

		for i := range e.xNew {
			spaceDerivative := (e.TdFunc.EvaluateAt(e.xNew[i]+delta, t+e.DeltaT) - e.fxt[i]) / delta
			e.xPre[i] = 1 - e.DeltaT*spaceDerivative
		}

		JacobianNorm := blas64.Nrm2(blas64.Vector{N: nPoints, Data: e.xPre, Inc: 1})
		if JacobianNorm < 1e-14 {
			e.TdFunc.EvaluateOnRGridInPlace(xt, e.fxt, t+e.DeltaT)
			for i := range e.xNew {
				e.xNew[i] = xt[i] + e.DeltaT*e.fxt[i]
			}
			continue
		}

		for i := range e.xNew {
			e.xNew[i] -= e.diff[i] / e.xPre[i]
		}
	}
	return fmt.Errorf("implicit step did not converge after %d iterations", maxIter)
}

// HeunsSolver implements the basic Euler method
type HeunsSolver struct {
	TdFunc    gridData.TDPotentialOp
	DeltaT    float64
	predictor []float64
	corrector []float64
	fxt       []float64
	jacobian  []float64
	diff      []float64
}

func (h *HeunsSolver) Name() string {
	return "Heun's Method (Improved Euler method)"
}

func (h *HeunsSolver) NextStepEx(xt, t float64) float64 {
	fxt := (h.TdFunc).EvaluateAt(xt, t)
	predictor := xt + fxt*h.DeltaT
	slope2 := fxt + (h.TdFunc).EvaluateAt(predictor, t+h.DeltaT)
	return xt + 0.5*h.DeltaT*slope2
}

func (h *HeunsSolver) PredictIni(xt, fxt, xP []float64) {
	nPoints := len(xP)
	xP = slices.Clone(xt)
	vecFxt := blas64.Vector{N: nPoints, Data: fxt, Inc: 1}
	xNew := blas64.Vector{N: nPoints, Data: xP, Inc: 1}
	blas64.Axpy(h.DeltaT, vecFxt, xNew)
}

func (h *HeunsSolver) NextStepExOnGrid(xt []float64, t float64) {
	nPoints := len(xt)

	if nPoints != len(h.predictor) {
		h.predictor = make([]float64, nPoints)
		h.corrector = make([]float64, nPoints)
		h.fxt = make([]float64, nPoints)
	}

	(h.TdFunc).EvaluateOnRGridInPlace(xt, h.fxt, t)
	h.PredictIni(xt, h.fxt, h.predictor)
	(h.TdFunc).EvaluateOnRGridInPlace(h.predictor, h.corrector, t)

	for i := range xt {
		h.corrector[i] += h.fxt[i]
	}

	vecFxt := blas64.Vector{N: nPoints, Data: h.corrector, Inc: 1}
	xNew := blas64.Vector{N: nPoints, Data: xt, Inc: 1}
	blas64.Axpy(0.5*h.DeltaT, vecFxt, xNew)
}

func (h *HeunsSolver) NextStepIm(xt, t float64) (float64, error) {
	fxt := (h.TdFunc).EvaluateAt(xt, t)
	predictor := xt + h.DeltaT*fxt
	corrector := predictor

	for i := 0; i < maxIter; i++ {
		trapezoid := fxt + (h.TdFunc).EvaluateAt(corrector, t+h.DeltaT)
		correctorNew := xt + 0.5*h.DeltaT*trapezoid

		var err float64
		if math.Abs(correctorNew) > 1e-10 {
			err = math.Abs(correctorNew-corrector) / math.Abs(correctorNew)
		} else {
			err = math.Abs(correctorNew - corrector)
		}

		if err < tolerance {
			return correctorNew, nil
		}

		corrector = correctorNew
	}
	return xt, fmt.Errorf("implicit step did not converge after %d iterations", maxIter)
}

func (h *HeunsSolver) NextStepImNewton(xt, t float64) (float64, error) {
	halfDt := 0.5 * h.DeltaT
	xPre := xt
	fxt := (h.TdFunc).EvaluateAt(xt, t)
	for i := 0; i < maxIter; i++ {
		fNew := (h.TdFunc).EvaluateAt(xPre, t+h.DeltaT)
		Gx := xPre - xt - halfDt*(fxt+fNew)

		fxPDel := (h.TdFunc).EvaluateAt(xPre+delta, t+h.DeltaT)
		derFxt := (fxPDel - fNew) / delta
		Jacobian := 1 - halfDt*derFxt

		if math.Abs(Jacobian) < 1e-10 {
			return xt, fmt.Errorf("singular Jacobian at iteration %d (Jacobian = %e)", i, Jacobian)
		}
		xNew := xPre - Gx/Jacobian

		if math.IsNaN(xNew) || math.IsInf(xNew, 0) {
			return xt, fmt.Errorf("newton iteration produced invalid value at iteration %d", i)
		}

		if math.Abs(xNew-xPre) < tolerance {
			return xNew, nil
		}

		xPre = xNew
	}
	return xt, fmt.Errorf("implicit step did not converge after %d iterations", maxIter)
}

func (h *HeunsSolver) allocateFixPoint(nPoints int) {

	if nPoints != len(h.corrector) ||
		nPoints != len(h.predictor) ||
		nPoints != len(h.fxt) ||
		nPoints != len(h.diff) {

		h.corrector = make([]float64, nPoints)
		h.predictor = make([]float64, nPoints)
		h.fxt = make([]float64, nPoints)
		h.diff = make([]float64, nPoints)
	}

}

func (h *HeunsSolver) iterate(xt, fxt, predictCorrect []float64, t float64) float64 {
	nPoints := len(xt)

	for i := range predictCorrect {
		h.corrector[i] = fxt[i] + (h.TdFunc).EvaluateAt(predictCorrect[i], t+h.DeltaT)
	}
	trapezoidal := blas64.Vector{N: nPoints, Data: h.corrector, Inc: 1}

	h.diff = slices.Clone(xt)
	xtNew := blas64.Vector{N: nPoints, Data: h.diff, Inc: 1}
	blas64.Axpy(0.5*h.DeltaT, trapezoidal, xtNew)

	h.corrector = slices.Clone(h.diff)

	for i := range h.corrector {
		h.diff[i] -= predictCorrect[i]
	}

	predictCorrect = slices.Clone(h.corrector)

	diffVec := blas64.Vector{N: nPoints, Data: h.diff, Inc: 1}
	correct := blas64.Vector{N: nPoints, Data: predictCorrect, Inc: 1}

	return blas64.Nrm2(diffVec) / blas64.Nrm2(correct)
}

func (h *HeunsSolver) NextStepImOnGrid(xt []float64, t float64) error {
	nPoints := len(xt)
	h.allocateFixPoint(nPoints)
	(h.TdFunc).EvaluateOnRGridInPlace(xt, h.fxt, t)

	h.PredictIni(xt, h.fxt, h.predictor)

	var err float64
	for i := 0; i < maxIter; i++ {
		err = h.iterate(xt, h.fxt, h.predictor, t)
		if err < tolerance {
			xt = slices.Clone(h.predictor)
			return nil
		}
	}
	return fmt.Errorf("implicit step did not converge after %d iterations", maxIter)
}

func (h *HeunsSolver) allocateNewton(nPoints int) {

	if nPoints != len(h.corrector) ||
		nPoints != len(h.predictor) ||
		nPoints != len(h.fxt) ||
		nPoints != len(h.diff) ||
		nPoints != len(h.jacobian) {

		h.corrector = make([]float64, nPoints)
		h.predictor = make([]float64, nPoints)
		h.fxt = make([]float64, nPoints)
		h.jacobian = make([]float64, nPoints)
		h.diff = make([]float64, nPoints)
	}

}

func (h *HeunsSolver) NextStepImOnGridNewton(xt []float64, t float64) error {
	nPoints := len(xt)
	h.allocateNewton(nPoints)

	halfDt := 0.5 * h.DeltaT
	(h.TdFunc).EvaluateOnRGridInPlace(xt, h.fxt, t)
	h.PredictIni(xt, h.fxt, h.predictor)

	for i := 0; i < maxIter; i++ {
		(h.TdFunc).EvaluateOnRGridInPlace(h.predictor, h.diff, t+h.DeltaT)

		for j := range h.predictor {
			h.diff[j] += h.fxt[j]
		}

		h.corrector = slices.Clone(xt)
		trapezoid := blas64.Vector{N: nPoints, Data: h.diff, Inc: 1}
		newXt := blas64.Vector{N: nPoints, Data: h.corrector, Inc: 1}
		blas64.Axpy(halfDt, trapezoid, newXt)

		// Gx := xPre - xt - halfDt*(fxt+fNew)
		for j := range h.predictor {
			h.diff[j] = h.predictor[j] - h.corrector[j]
			fxCenter := (h.TdFunc).EvaluateAt(h.predictor[j], t+h.DeltaT)
			fxPDel := (h.TdFunc).EvaluateAt(h.predictor[j]+delta, t+h.DeltaT)
			derFxt := (fxPDel - fxCenter) / delta
			h.jacobian[j] = 1 - halfDt*derFxt
		}

		singular := false
		for j := range h.jacobian {
			if math.Abs(h.jacobian[j]) < 1e-10 {
				singular = true
				break
			}
		}
		if singular {
			return fmt.Errorf("singular Jacobian at iteration %d", i)
		}

		for j := range h.predictor {
			h.corrector[j] = h.predictor[j] - h.diff[j]/h.jacobian[j]
		}

		xNew := blas64.Nrm2(blas64.Vector{N: nPoints, Data: h.corrector, Inc: 1})
		if math.IsNaN(xNew) || math.IsInf(xNew, 0) {
			return fmt.Errorf("newton iteration produced invalid value at iteration %d", i)
		}

		for k := range h.predictor {
			h.diff[k] = h.corrector[k] - h.predictor[k]
		}

		diffNorm := blas64.Nrm2(blas64.Vector{N: nPoints, Data: h.diff, Inc: 1})

		if diffNorm < tolerance {
			copy(xt, h.corrector)
			return nil
		}
		h.predictor = slices.Clone(h.corrector)
	}

	copy(xt, h.corrector)
	return fmt.Errorf("implicit step did not converge after %d iterations", maxIter)
}

type MidPointSolver struct {
	TdFunc gridData.TDPotentialOp
	DeltaT float64
	xtMid  []float64
	slope  []float64
}

func (m *MidPointSolver) Name() string {
	return "MidPoint Method (Implicit or Explicit)"
}

func (m *MidPointSolver) NextStepEx(xt, t float64) float64 {
	xtMid := xt + (m.TdFunc).EvaluateAt(xt, t)*m.DeltaT/2
	slope := (m.TdFunc).EvaluateAt(xtMid, t+m.DeltaT/2)
	return xt + m.DeltaT*slope
}

func (m *MidPointSolver) NextStepExOnGrid(xt []float64, t float64) {
	nPoints := len(xt)

	(m.TdFunc).EvaluateOnRGridInPlace(xt, m.slope, t)

	m.xtMid = slices.Clone(xt)
	blas64.Axpy(m.DeltaT/2,
		blas64.Vector{N: nPoints, Data: m.slope, Inc: 1},
		blas64.Vector{N: nPoints, Data: m.xtMid, Inc: 1},
	)

	for i := range m.xtMid {
		m.slope[i] = (m.TdFunc).EvaluateOnRGrid(m.xtMid, t+m.DeltaT/2)[i]
	}

	blas64.Axpy(m.DeltaT,
		blas64.Vector{N: nPoints, Data: m.slope, Inc: 1},
		blas64.Vector{N: nPoints, Data: xt, Inc: 1},
	)
}

func (m *MidPointSolver) NextStepIm(xt, t float64) float64 {
	xtMid := xt + (m.TdFunc).EvaluateAt(xt, t)*m.DeltaT/2
	return xt + m.DeltaT*(m.TdFunc).EvaluateAt(xtMid, t+m.DeltaT/2)
}

func (m *MidPointSolver) NextStepImOnGrid(xt []float64, t float64) {}

// RungeKutta4 implements the basic Euler method
type RungeKutta4 struct {
	TdFunc gridData.TDPotentialOp
	DeltaT float64
}

func (rg *RungeKutta4) Name() string {
	return "Runge Kutta order 4"
}

func (rg *RungeKutta4) NextStepEx(xt, t float64) float64 {
	halfDt := 0.5 * rg.DeltaT
	k1 := (rg.TdFunc).EvaluateAt(xt, t)
	k2 := (rg.TdFunc).EvaluateAt(xt+k1*halfDt, t+halfDt)
	k3 := (rg.TdFunc).EvaluateAt(xt+k2*halfDt, t+halfDt)
	k4 := (rg.TdFunc).EvaluateAt(xt+k3*rg.DeltaT, t+rg.DeltaT)
	return xt + rg.DeltaT/6.0*(k1+2*k2+2*k3+k4)
}

func (rg *RungeKutta4) NextStepExOnGrid(xt []float64, t float64) {
	halfDt := 0.5 * rg.DeltaT
	k1 := (rg.TdFunc).EvaluateOnRGrid(xt, t)

	xtPlusDt := slices.Clone(xt)
	blas64.Axpy(0.5*rg.DeltaT,
		blas64.Vector{N: len(k1), Data: k1, Inc: 1},
		blas64.Vector{N: len(xtPlusDt), Data: xtPlusDt, Inc: 1},
	)
	k2 := (rg.TdFunc).EvaluateOnRGrid(xtPlusDt, t+halfDt)
	xtPlusDt = slices.Clone(xt)
	blas64.Axpy(0.5*rg.DeltaT,
		blas64.Vector{N: len(k2), Data: k2, Inc: 1},
		blas64.Vector{N: len(xtPlusDt), Data: xtPlusDt, Inc: 1},
	)

	k3 := (rg.TdFunc).EvaluateOnRGrid(xtPlusDt, t+halfDt)

	xtPlusDt = slices.Clone(xt)
	blas64.Axpy(rg.DeltaT,
		blas64.Vector{N: len(k3), Data: k3, Inc: 1},
		blas64.Vector{N: len(xtPlusDt), Data: xtPlusDt, Inc: 1},
	)

	k4 := (rg.TdFunc).EvaluateOnRGrid(xtPlusDt, t+rg.DeltaT)

	Dtby6 := rg.DeltaT / 6
	for i := range xtPlusDt {
		xt[i] += Dtby6 * (k1[i] + 2*(k2[i]+k3[i]) + k4[i])
	}
}

// PredictorCorrector implements the basic Euler method
type PredictorCorrector struct {
	TdFunc gridData.TDPotentialOp
	DeltaT float64
}
