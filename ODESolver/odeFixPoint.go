package EquationSolver

import (
	"GoProject/gridData"
	"fmt"
	"math"
	"slices"

	"gonum.org/v1/gonum/blas/blas64"
)

// EulerImplicitFixPoint implements the basic Euler method
type EulerImplicitFixPoint struct {
	TdFunc gridData.TDPotentialOp
	DeltaT float64
	fxt    []float64
	xPre   []float64
	xNew   []float64
	diff   []float64
}

func (eImFP *EulerImplicitFixPoint) Name() string {
	return "Euler Method"
}

func (eImFP *EulerImplicitFixPoint) NextStep(xt, t float64) (float64, error) {
	xPre := xt

	for iter := 0; iter < maxIter; iter++ {
		fxt := eImFP.TdFunc.EvaluateAtTime(xPre, t+eImFP.DeltaT)
		xNew := xt + eImFP.DeltaT*fxt

		if math.Abs(xNew-xPre) <= tolerance {
			return xNew, nil
		}
		xPre = xNew
	}

	return xt, fmt.Errorf("implicit step did not converge after %d iterations", maxIter)
}

func (eImFP *EulerImplicitFixPoint) allocate(nPoints int) {
	if nPoints != len(eImFP.fxt) ||
		nPoints != len(eImFP.xNew) ||
		nPoints != len(eImFP.xPre) ||
		nPoints != len(eImFP.diff) {
		eImFP.fxt = make([]float64, nPoints)
		eImFP.xPre = make([]float64, nPoints)
		eImFP.xNew = make([]float64, nPoints)
		eImFP.diff = make([]float64, nPoints)
	}
}

func (eImFP *EulerImplicitFixPoint) NextStepOnGrid(xt []float64, t float64) error {
	nPoints := len(xt)
	eImFP.allocate(nPoints)
	eImFP.xPre = slices.Clone(xt)

	for iter := 0; iter < maxIter; iter++ {
		(eImFP.TdFunc).EvaluateOnRGridTimeInPlace(eImFP.xPre, eImFP.fxt, t+eImFP.DeltaT)

		for i := range eImFP.xNew {
			eImFP.xNew[i] = xt[i] + eImFP.DeltaT*eImFP.fxt[i]
			eImFP.diff[i] = eImFP.xNew[i] - eImFP.xPre[i]
		}

		err := blas64.Nrm2(blas64.Vector{N: nPoints, Data: eImFP.diff, Inc: 1})

		if err <= tolerance {
			copy(xt, eImFP.xNew)
			return nil
		}
		eImFP.xPre = slices.Clone(eImFP.xNew)
	}
	return fmt.Errorf("implicit step did not converge after %d iterations", maxIter)
}

// HeunsImplicitFixPoint implements the basic Euler method
type HeunsImplicitFixPoint struct {
	TdFunc    gridData.TDPotentialOp
	DeltaT    float64
	predictor []float64
	corrector []float64
	fxt       []float64
	jacobian  []float64
	diff      []float64
}

func (hIm *HeunsImplicitFixPoint) Name() string {
	return "Heun's Method (Improved Euler method)"
}

func (hIm *HeunsImplicitFixPoint) PredictIni(xt, fxt, xP []float64) {
	nPoints := len(xP)
	xP = slices.Clone(xt)
	vecFxt := blas64.Vector{N: nPoints, Data: fxt, Inc: 1}
	xNew := blas64.Vector{N: nPoints, Data: xP, Inc: 1}
	blas64.Axpy(hIm.DeltaT, vecFxt, xNew)
}

func (hIm *HeunsImplicitFixPoint) NextStep(xt, t float64) (float64, error) {
	fxt := (hIm.TdFunc).EvaluateAtTime(xt, t)
	predictor := xt + hIm.DeltaT*fxt
	corrector := predictor

	for i := 0; i < maxIter; i++ {
		trapezoid := fxt + (hIm.TdFunc).EvaluateAtTime(corrector, t+hIm.DeltaT)
		correctorNew := xt + 0.5*hIm.DeltaT*trapezoid

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

func (hIm *HeunsImplicitFixPoint) allocate(nPoints int) {

	if nPoints != len(hIm.corrector) ||
		nPoints != len(hIm.predictor) ||
		nPoints != len(hIm.fxt) ||
		nPoints != len(hIm.diff) {

		hIm.corrector = make([]float64, nPoints)
		hIm.predictor = make([]float64, nPoints)
		hIm.fxt = make([]float64, nPoints)
		hIm.diff = make([]float64, nPoints)
	}

}

func (hIm *HeunsImplicitFixPoint) iterate(xt, fxt, predictCorrect []float64, t float64) float64 {
	nPoints := len(xt)

	for i := range predictCorrect {
		hIm.corrector[i] = fxt[i] + (hIm.TdFunc).EvaluateAtTime(predictCorrect[i], t+hIm.DeltaT)
	}
	trapezoidal := blas64.Vector{N: nPoints, Data: hIm.corrector, Inc: 1}

	hIm.diff = slices.Clone(xt)
	xtNew := blas64.Vector{N: nPoints, Data: hIm.diff, Inc: 1}
	blas64.Axpy(0.5*hIm.DeltaT, trapezoidal, xtNew)

	hIm.corrector = slices.Clone(hIm.diff)

	for i := range hIm.corrector {
		hIm.diff[i] -= predictCorrect[i]
	}

	predictCorrect = slices.Clone(hIm.corrector)

	diffVec := blas64.Vector{N: nPoints, Data: hIm.diff, Inc: 1}
	correct := blas64.Vector{N: nPoints, Data: predictCorrect, Inc: 1}

	return blas64.Nrm2(diffVec) / blas64.Nrm2(correct)
}

func (hIm *HeunsImplicitFixPoint) NextStepOnGrid(xt []float64, t float64) error {
	nPoints := len(xt)
	hIm.allocate(nPoints)
	(hIm.TdFunc).EvaluateOnRGridTimeInPlace(xt, hIm.fxt, t)

	hIm.PredictIni(xt, hIm.fxt, hIm.predictor)

	var err float64
	for i := 0; i < maxIter; i++ {
		err = hIm.iterate(xt, hIm.fxt, hIm.predictor, t)
		if err < tolerance {
			xt = slices.Clone(hIm.predictor)
			return nil
		}
	}
	return fmt.Errorf("implicit step did not converge after %d iterations", maxIter)
}

type MidPointFixPoint struct {
	TdFunc gridData.TDPotentialOp
	DeltaT float64
	xMid   []float64
	xPre   []float64
	fxt    []float64
}

func (mIm *MidPointFixPoint) Name() string {
	return "MidPoint Implicit with Fix-Point Method"
}

func (mIm *MidPointFixPoint) NextStep(xt, t float64) (float64, error) {
	tMid := t + mIm.DeltaT/2
	xPre := xt
	var xNext float64
	for iter := 0; iter < maxIter; iter++ {
		xMid := 0.5 * (xt + xPre)
		funcMid := mIm.TdFunc.EvaluateAtTime(xMid, tMid)
		xNext = xt + mIm.DeltaT*funcMid

		if math.Abs(xNext-xPre) < tolerance {
			return xNext, nil
		}
		xPre = xNext
	}
	return xt, fmt.Errorf("implicit fix point step timed out! Maximum Iterations: %d", maxIter)
}

func (mIm *MidPointFixPoint) NextStepOnGrid(xt []float64, t float64) error {
	nPoints := len(xt)

	tMid := t + mIm.DeltaT/2
	copy(mIm.xPre, xt)

	if nPoints != len(mIm.xMid) ||
		nPoints != len(mIm.fxt) ||
		nPoints != len(mIm.xPre) {

		mIm.xPre = make([]float64, nPoints)
		mIm.xMid = make([]float64, nPoints)
		mIm.fxt = make([]float64, nPoints)
	}

	for iter := 0; iter < maxIter; iter++ {

		for i := range xt {
			mIm.xMid[i] = 0.5 * (xt[i] + mIm.xPre[i])
		}

		mIm.TdFunc.EvaluateOnRGridTimeInPlace(mIm.xMid, mIm.fxt, tMid)

		for i := range xt {
			mIm.xMid[i] = xt[i] + mIm.DeltaT*mIm.fxt[i]
			mIm.fxt[i] = mIm.xMid[i] - mIm.xPre[i]
		}

		errDiff := blas64.Vector{
			N:    nPoints,
			Data: mIm.fxt,
			Inc:  1,
		}

		norm := blas64.Nrm2(errDiff)
		if math.Abs(norm) < tolerance {
			copy(xt, mIm.xMid)
			return nil
		}
		copy(mIm.xPre, mIm.xMid)
	}
	return fmt.Errorf("implicit fix point step timed out! Maximum Iterations: %d", maxIter)
}

// PredictorCorrector implements the basic Euler method
type PredictorCorrector struct {
	TdFunc gridData.TDPotentialOp
	DeltaT float64
}
