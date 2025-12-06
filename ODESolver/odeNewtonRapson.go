package EquationSolver

import (
	"GoProject/gridData"
	"fmt"
	"math"
	"slices"

	"gonum.org/v1/gonum/blas/blas64"
)

// EulerImplicitNewton EulerImplicitFixPoint implements the basic Euler method
type EulerImplicitNewton struct {
	TdFunc gridData.TDPotentialOp
	DeltaT float64
	fxt    []float64
	xPre   []float64
	xNew   []float64
	diff   []float64
}

func (eImNw *EulerImplicitNewton) Name() string {
	return "Euler Method"
}

func (eImNw *EulerImplicitNewton) allocate(nPoints int) {
	if nPoints != len(eImNw.fxt) ||
		nPoints != len(eImNw.xNew) ||
		nPoints != len(eImNw.xPre) ||
		nPoints != len(eImNw.diff) {
		eImNw.fxt = make([]float64, nPoints)
		eImNw.xPre = make([]float64, nPoints)
		eImNw.xNew = make([]float64, nPoints)
		eImNw.diff = make([]float64, nPoints)
	}
}

func (eImNw *EulerImplicitNewton) NextStep(xt, t float64) (float64, error) {
	xNext := xt
	for iter := 0; iter < maxIter; iter++ {
		fxt := eImNw.TdFunc.EvaluateAtTime(xNext, t+eImNw.DeltaT)
		Gxt := xNext - (xt + eImNw.DeltaT*fxt)

		if math.Abs(Gxt) < tolerance {
			return xNext, nil
		}

		derGxt := (eImNw.TdFunc.EvaluateAtTime(xNext+delta, t+eImNw.DeltaT) - fxt) / delta
		jacobian := 1 - eImNw.DeltaT*derGxt

		if math.Abs(jacobian) < 1e-14 {
			xNext = xt + eImNw.DeltaT*fxt
			continue
		}

		xNext -= Gxt / jacobian
	}
	return xt, fmt.Errorf("implicit step did not converge after %d iterations", maxIter)
}

func (eImNw *EulerImplicitNewton) NextStepOnGrid(xt []float64, t float64) error {
	nPoints := len(xt)
	eImNw.allocate(nPoints)
	eImNw.xNew = slices.Clone(xt)
	for iter := 0; iter < maxIter; iter++ {
		eImNw.TdFunc.EvaluateOnRGridTimeInPlace(eImNw.xNew, eImNw.fxt, t+eImNw.DeltaT)

		eImNw.xPre = slices.Clone(xt)
		blas64.Axpy(eImNw.DeltaT,
			blas64.Vector{N: nPoints, Data: eImNw.fxt, Inc: 1},
			blas64.Vector{N: nPoints, Data: eImNw.xPre, Inc: 1},
		)

		for i := range eImNw.xNew {
			eImNw.diff[i] = eImNw.xNew[i] - eImNw.xPre[i]
		}

		residualNorm := blas64.Nrm2(blas64.Vector{N: nPoints, Data: eImNw.diff, Inc: 1})
		if residualNorm < tolerance {
			copy(xt, eImNw.xNew)
			return nil
		}

		for i := range eImNw.xNew {
			spaceDerivative := (eImNw.TdFunc.EvaluateAtTime(eImNw.xNew[i]+delta, t+eImNw.DeltaT) - eImNw.fxt[i]) / delta
			eImNw.xPre[i] = 1 - eImNw.DeltaT*spaceDerivative
		}

		JacobianNorm := blas64.Nrm2(blas64.Vector{N: nPoints, Data: eImNw.xPre, Inc: 1})
		if JacobianNorm < 1e-14 {
			eImNw.TdFunc.EvaluateOnRGridTimeInPlace(xt, eImNw.fxt, t+eImNw.DeltaT)
			for i := range eImNw.xNew {
				eImNw.xNew[i] = xt[i] + eImNw.DeltaT*eImNw.fxt[i]
			}
			continue
		}

		for i := range eImNw.xNew {
			eImNw.xNew[i] -= eImNw.diff[i] / eImNw.xPre[i]
		}
	}
	return fmt.Errorf("implicit step did not converge after %d iterations", maxIter)
}

// HeunsImplicitNewton implements the basic Euler method
type HeunsImplicitNewton struct {
	TdFunc    gridData.TDPotentialOp
	DeltaT    float64
	predictor []float64
	corrector []float64
	fxt       []float64
	jacobian  []float64
	diff      []float64
}

func (hIm *HeunsImplicitNewton) Name() string {
	return "Heun's Method (Improved Euler method)"
}

func (hIm *HeunsImplicitNewton) PredictIni(xt, fxt, xP []float64) {
	nPoints := len(xP)
	xP = slices.Clone(xt)
	vecFxt := blas64.Vector{N: nPoints, Data: fxt, Inc: 1}
	xNew := blas64.Vector{N: nPoints, Data: xP, Inc: 1}
	blas64.Axpy(hIm.DeltaT, vecFxt, xNew)
}

func (hIm *HeunsImplicitNewton) NextStep(xt, t float64) (float64, error) {
	halfDt := 0.5 * hIm.DeltaT
	xPre := xt
	fxt := (hIm.TdFunc).EvaluateAtTime(xt, t)
	for i := 0; i < maxIter; i++ {
		fNew := (hIm.TdFunc).EvaluateAtTime(xPre, t+hIm.DeltaT)
		Gx := xPre - xt - halfDt*(fxt+fNew)

		fxPDel := (hIm.TdFunc).EvaluateAtTime(xPre+delta, t+hIm.DeltaT)
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

func (hIm *HeunsImplicitNewton) allocate(nPoints int) {

	if nPoints != len(hIm.corrector) ||
		nPoints != len(hIm.predictor) ||
		nPoints != len(hIm.fxt) ||
		nPoints != len(hIm.diff) ||
		nPoints != len(hIm.jacobian) {

		hIm.corrector = make([]float64, nPoints)
		hIm.predictor = make([]float64, nPoints)
		hIm.fxt = make([]float64, nPoints)
		hIm.jacobian = make([]float64, nPoints)
		hIm.diff = make([]float64, nPoints)
	}

}

func (hIm *HeunsImplicitNewton) NextStepOnGrid(xt []float64, t float64) error {
	nPoints := len(xt)
	hIm.allocate(nPoints)

	halfDt := 0.5 * hIm.DeltaT
	(hIm.TdFunc).EvaluateOnRGridTimeInPlace(xt, hIm.fxt, t)
	hIm.PredictIni(xt, hIm.fxt, hIm.predictor)

	for i := 0; i < maxIter; i++ {
		(hIm.TdFunc).EvaluateOnRGridTimeInPlace(hIm.predictor, hIm.diff, t+hIm.DeltaT)

		for j := range hIm.predictor {
			hIm.diff[j] += hIm.fxt[j]
		}

		hIm.corrector = slices.Clone(xt)
		trapezoid := blas64.Vector{N: nPoints, Data: hIm.diff, Inc: 1}
		newXt := blas64.Vector{N: nPoints, Data: hIm.corrector, Inc: 1}
		blas64.Axpy(halfDt, trapezoid, newXt)

		// Gx := xPre - xt - halfDt*(fxt+fNew)
		for j := range hIm.predictor {
			hIm.diff[j] = hIm.predictor[j] - hIm.corrector[j]
			fxCenter := (hIm.TdFunc).EvaluateAtTime(hIm.predictor[j], t+hIm.DeltaT)
			fxPDel := (hIm.TdFunc).EvaluateAtTime(hIm.predictor[j]+delta, t+hIm.DeltaT)
			derFxt := (fxPDel - fxCenter) / delta
			hIm.jacobian[j] = 1 - halfDt*derFxt
		}

		singular := false
		for j := range hIm.jacobian {
			if math.Abs(hIm.jacobian[j]) < 1e-10 {
				singular = true
				break
			}
		}
		if singular {
			return fmt.Errorf("singular Jacobian at iteration %d", i)
		}

		for j := range hIm.predictor {
			hIm.corrector[j] = hIm.predictor[j] - hIm.diff[j]/hIm.jacobian[j]
		}

		xNew := blas64.Nrm2(blas64.Vector{N: nPoints, Data: hIm.corrector, Inc: 1})
		if math.IsNaN(xNew) || math.IsInf(xNew, 0) {
			return fmt.Errorf("newton iteration produced invalid value at iteration %d", i)
		}

		for k := range hIm.predictor {
			hIm.diff[k] = hIm.corrector[k] - hIm.predictor[k]
		}

		diffNorm := blas64.Nrm2(blas64.Vector{N: nPoints, Data: hIm.diff, Inc: 1})

		if diffNorm < tolerance {
			copy(xt, hIm.corrector)
			return nil
		}
		hIm.predictor = slices.Clone(hIm.corrector)
	}

	copy(xt, hIm.corrector)
	return fmt.Errorf("implicit step did not converge after %d iterations", maxIter)
}

type MidPointNewton struct {
	TdFunc gridData.TDPotentialOp
	DeltaT float64
	xtMid  []float64
	slope  []float64
}

func (m *MidPointNewton) Name() string {
	return "MidPoint Implicit with Newton-Raphson Method"
}

func (m *MidPointNewton) NextStepIm(xt, t float64) (float64, error) {
	halfDt := m.DeltaT / 2
	tMid := t + halfDt
	xNext := xt
	for iter := 0; iter < maxIter; iter++ {
		xMid := (xt + xNext) / 2
		fMid := m.TdFunc.EvaluateAtTime(xMid, tMid)

		Gx := xNext - xt - m.DeltaT*fMid
		if math.Abs(Gx) < tolerance {
			return xNext, nil
		}

		fPlus := m.TdFunc.EvaluateAtTime(xMid+delta, tMid)
		dFbydx := (fPlus - fMid) / delta
		jacobian := 1 - halfDt*dFbydx
		xNext -= Gx / jacobian
	}
	return xt, fmt.Errorf("implicit step did not converge after %d iterations", maxIter)
}
