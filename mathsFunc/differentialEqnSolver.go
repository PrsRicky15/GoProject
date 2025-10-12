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

type EulerExplicit struct {
	TdFunc gridData.TDPotentialOp
	DeltaT float64
	fxt    []float64
}

func (eEx *EulerExplicit) Name() string {
	return "Euler Explicit Method"
}

func (eEx *EulerExplicit) NextStepEx(xt, t float64) float64 {
	fxt := (eEx.TdFunc).EvaluateAt(xt, t)
	return xt + eEx.DeltaT*fxt
}

func (eEx *EulerExplicit) NextStepExOnGrid(xt []float64, t float64) {
	nPoints := len(xt)
	if nPoints != len(eEx.fxt) {
		eEx.fxt = make([]float64, len(xt))
	}

	(eEx.TdFunc).EvaluateOnRGridInPlace(xt, eEx.fxt, t)
	blas64.Axpy(eEx.DeltaT,
		blas64.Vector{N: nPoints, Data: eEx.fxt, Inc: 1},
		blas64.Vector{N: nPoints, Data: xt, Inc: 1},
	)
}

// HeunsExplicit HeunsImplicitFixPoint implements the basic Euler method
type HeunsExplicit struct {
	TdFunc    gridData.TDPotentialOp
	DeltaT    float64
	predictor []float64
	corrector []float64
	fxt       []float64
}

func (hEx *HeunsExplicit) Name() string {
	return "Heun's Explicit Method (Improved Euler method)"
}

func (hEx *HeunsExplicit) NextStepEx(xt, t float64) float64 {
	fxt := (hEx.TdFunc).EvaluateAt(xt, t)
	predictor := xt + fxt*hEx.DeltaT
	slope2 := fxt + (hEx.TdFunc).EvaluateAt(predictor, t+hEx.DeltaT)
	return xt + 0.5*hEx.DeltaT*slope2
}

func (hEx *HeunsExplicit) PredictIni(xt, fxt, xP []float64) {
	nPoints := len(xP)
	xP = slices.Clone(xt)
	vecFxt := blas64.Vector{N: nPoints, Data: fxt, Inc: 1}
	xNew := blas64.Vector{N: nPoints, Data: xP, Inc: 1}
	blas64.Axpy(hEx.DeltaT, vecFxt, xNew)
}

func (hEx *HeunsExplicit) NextStepExOnGrid(xt []float64, t float64) {
	nPoints := len(xt)

	if nPoints != len(hEx.predictor) {
		hEx.predictor = make([]float64, nPoints)
		hEx.corrector = make([]float64, nPoints)
		hEx.fxt = make([]float64, nPoints)
	}

	(hEx.TdFunc).EvaluateOnRGridInPlace(xt, hEx.fxt, t)
	hEx.PredictIni(xt, hEx.fxt, hEx.predictor)
	(hEx.TdFunc).EvaluateOnRGridInPlace(hEx.predictor, hEx.corrector, t)

	for i := range xt {
		hEx.corrector[i] += hEx.fxt[i]
	}

	vecFxt := blas64.Vector{N: nPoints, Data: hEx.corrector, Inc: 1}
	xNew := blas64.Vector{N: nPoints, Data: xt, Inc: 1}
	blas64.Axpy(0.5*hEx.DeltaT, vecFxt, xNew)
}

type MidPointExplicit struct {
	TdFunc gridData.TDPotentialOp
	DeltaT float64
	xtMid  []float64
	slope  []float64
}

func (mEx *MidPointExplicit) Name() string {
	return "MidPoint Explicit Method!"
}

func (mEx *MidPointExplicit) NextStepEx(xt, t float64) float64 {
	xtMid := xt + (mEx.TdFunc).EvaluateAt(xt, t)*mEx.DeltaT/2
	slope := (mEx.TdFunc).EvaluateAt(xtMid, t+mEx.DeltaT/2)
	return xt + mEx.DeltaT*slope
}

func (mEx *MidPointExplicit) NextStepExOnGrid(xt []float64, t float64) {
	nPoints := len(xt)

	(mEx.TdFunc).EvaluateOnRGridInPlace(xt, mEx.slope, t)

	mEx.xtMid = slices.Clone(xt)
	blas64.Axpy(mEx.DeltaT/2,
		blas64.Vector{N: nPoints, Data: mEx.slope, Inc: 1},
		blas64.Vector{N: nPoints, Data: mEx.xtMid, Inc: 1},
	)

	for i := range mEx.xtMid {
		mEx.slope[i] = (mEx.TdFunc).EvaluateOnRGrid(mEx.xtMid, t+mEx.DeltaT/2)[i]
	}

	blas64.Axpy(mEx.DeltaT,
		blas64.Vector{N: nPoints, Data: mEx.slope, Inc: 1},
		blas64.Vector{N: nPoints, Data: xt, Inc: 1},
	)
}

// RungeKutta4Explicit implements the basic Euler method
type RungeKutta4Explicit struct {
	TdFunc gridData.TDPotentialOp
	DeltaT float64
}

func (rgEx *RungeKutta4Explicit) Name() string {
	return "Runge Kutta order 4"
}

func (rgEx *RungeKutta4Explicit) NextStepEx(xt, t float64) float64 {
	halfDt := 0.5 * rgEx.DeltaT
	k1 := (rgEx.TdFunc).EvaluateAt(xt, t)
	k2 := (rgEx.TdFunc).EvaluateAt(xt+k1*halfDt, t+halfDt)
	k3 := (rgEx.TdFunc).EvaluateAt(xt+k2*halfDt, t+halfDt)
	k4 := (rgEx.TdFunc).EvaluateAt(xt+k3*rgEx.DeltaT, t+rgEx.DeltaT)
	return xt + rgEx.DeltaT/6.0*(k1+2*k2+2*k3+k4)
}

func (rgEx *RungeKutta4Explicit) NextStepExOnGrid(xt []float64, t float64) {
	halfDt := 0.5 * rgEx.DeltaT
	k1 := (rgEx.TdFunc).EvaluateOnRGrid(xt, t)

	xtPlusDt := slices.Clone(xt)
	blas64.Axpy(0.5*rgEx.DeltaT,
		blas64.Vector{N: len(k1), Data: k1, Inc: 1},
		blas64.Vector{N: len(xtPlusDt), Data: xtPlusDt, Inc: 1},
	)
	k2 := (rgEx.TdFunc).EvaluateOnRGrid(xtPlusDt, t+halfDt)
	xtPlusDt = slices.Clone(xt)
	blas64.Axpy(0.5*rgEx.DeltaT,
		blas64.Vector{N: len(k2), Data: k2, Inc: 1},
		blas64.Vector{N: len(xtPlusDt), Data: xtPlusDt, Inc: 1},
	)

	k3 := (rgEx.TdFunc).EvaluateOnRGrid(xtPlusDt, t+halfDt)

	xtPlusDt = slices.Clone(xt)
	blas64.Axpy(rgEx.DeltaT,
		blas64.Vector{N: len(k3), Data: k3, Inc: 1},
		blas64.Vector{N: len(xtPlusDt), Data: xtPlusDt, Inc: 1},
	)

	k4 := (rgEx.TdFunc).EvaluateOnRGrid(xtPlusDt, t+rgEx.DeltaT)

	Dtby6 := rgEx.DeltaT / 6
	for i := range xtPlusDt {
		xt[i] += Dtby6 * (k1[i] + 2*(k2[i]+k3[i]) + k4[i])
	}
}

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

func (eImFP *EulerImplicitFixPoint) NextStepIm(xt, t float64) (float64, error) {
	xPre := xt

	for iter := 0; iter < maxIter; iter++ {
		fxt := eImFP.TdFunc.EvaluateAt(xPre, t+eImFP.DeltaT)
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

func (eImFP *EulerImplicitFixPoint) NextStepImOnGrid(xt []float64, t float64) error {
	nPoints := len(xt)
	eImFP.allocate(nPoints)
	eImFP.xPre = slices.Clone(xt)

	for iter := 0; iter < maxIter; iter++ {
		(eImFP.TdFunc).EvaluateOnRGridInPlace(eImFP.xPre, eImFP.fxt, t+eImFP.DeltaT)

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

func (eImNw *EulerImplicitNewton) NextStepIm(xt, t float64) (float64, error) {
	xNext := xt
	for iter := 0; iter < maxIter; iter++ {
		fxt := eImNw.TdFunc.EvaluateAt(xNext, t+eImNw.DeltaT)
		Gxt := xNext - (xt + eImNw.DeltaT*fxt)

		if math.Abs(Gxt) < tolerance {
			return xNext, nil
		}

		derGxt := (eImNw.TdFunc.EvaluateAt(xNext+delta, t+eImNw.DeltaT) - fxt) / delta
		jacobian := 1 - eImNw.DeltaT*derGxt

		if math.Abs(jacobian) < 1e-14 {
			xNext = xt + eImNw.DeltaT*fxt
			continue
		}

		xNext -= Gxt / jacobian
	}
	return xt, fmt.Errorf("implicit step did not converge after %d iterations", maxIter)
}

func (eImNw *EulerImplicitNewton) NextStepImOnGrid(xt []float64, t float64) error {
	nPoints := len(xt)
	eImNw.allocate(nPoints)
	eImNw.xNew = slices.Clone(xt)
	for iter := 0; iter < maxIter; iter++ {
		eImNw.TdFunc.EvaluateOnRGridInPlace(eImNw.xNew, eImNw.fxt, t+eImNw.DeltaT)

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
			spaceDerivative := (eImNw.TdFunc.EvaluateAt(eImNw.xNew[i]+delta, t+eImNw.DeltaT) - eImNw.fxt[i]) / delta
			eImNw.xPre[i] = 1 - eImNw.DeltaT*spaceDerivative
		}

		JacobianNorm := blas64.Nrm2(blas64.Vector{N: nPoints, Data: eImNw.xPre, Inc: 1})
		if JacobianNorm < 1e-14 {
			eImNw.TdFunc.EvaluateOnRGridInPlace(xt, eImNw.fxt, t+eImNw.DeltaT)
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

func (hIm *HeunsImplicitFixPoint) NextStepIm(xt, t float64) (float64, error) {
	fxt := (hIm.TdFunc).EvaluateAt(xt, t)
	predictor := xt + hIm.DeltaT*fxt
	corrector := predictor

	for i := 0; i < maxIter; i++ {
		trapezoid := fxt + (hIm.TdFunc).EvaluateAt(corrector, t+hIm.DeltaT)
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
		hIm.corrector[i] = fxt[i] + (hIm.TdFunc).EvaluateAt(predictCorrect[i], t+hIm.DeltaT)
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

func (hIm *HeunsImplicitFixPoint) NextStepImOnGrid(xt []float64, t float64) error {
	nPoints := len(xt)
	hIm.allocate(nPoints)
	(hIm.TdFunc).EvaluateOnRGridInPlace(xt, hIm.fxt, t)

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

func (hIm *HeunsImplicitNewton) NextStepIm(xt, t float64) (float64, error) {
	halfDt := 0.5 * hIm.DeltaT
	xPre := xt
	fxt := (hIm.TdFunc).EvaluateAt(xt, t)
	for i := 0; i < maxIter; i++ {
		fNew := (hIm.TdFunc).EvaluateAt(xPre, t+hIm.DeltaT)
		Gx := xPre - xt - halfDt*(fxt+fNew)

		fxPDel := (hIm.TdFunc).EvaluateAt(xPre+delta, t+hIm.DeltaT)
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

func (hIm *HeunsImplicitNewton) NextStepImOnGrid(xt []float64, t float64) error {
	nPoints := len(xt)
	hIm.allocate(nPoints)

	halfDt := 0.5 * hIm.DeltaT
	(hIm.TdFunc).EvaluateOnRGridInPlace(xt, hIm.fxt, t)
	hIm.PredictIni(xt, hIm.fxt, hIm.predictor)

	for i := 0; i < maxIter; i++ {
		(hIm.TdFunc).EvaluateOnRGridInPlace(hIm.predictor, hIm.diff, t+hIm.DeltaT)

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
			fxCenter := (hIm.TdFunc).EvaluateAt(hIm.predictor[j], t+hIm.DeltaT)
			fxPDel := (hIm.TdFunc).EvaluateAt(hIm.predictor[j]+delta, t+hIm.DeltaT)
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

type MidPointImplicitFixPoint struct {
	TdFunc gridData.TDPotentialOp
	DeltaT float64
	xtMid  []float64
	slope  []float64
}

func (m *MidPointImplicitFixPoint) Name() string {
	return "MidPoint Implicit with Fix-Point Method"
}

func (m *MidPointImplicitFixPoint) NextStepIm(xt, t float64) (float64, error) {
	tMid := t + m.DeltaT/2
	xPre := xt
	var xNext float64
	for iter := 0; iter < maxIter; iter++ {
		xMid := 0.5 * (xt + xPre)
		funcMid := m.TdFunc.EvaluateAt(xMid, tMid)
		xNext = xt + m.DeltaT*funcMid

		if math.Abs(xNext-xPre) < tolerance {
			return xNext, nil
		}
		xPre = xNext
	}
	return xt, fmt.Errorf("implicit fix point step timed out! Maximum Iterations: %d", maxIter)
}

func (m *MidPointImplicitFixPoint) NextStepImOnGrid(xt []float64, t float64) {}

type MidPointImplicitNewton struct {
	TdFunc gridData.TDPotentialOp
	DeltaT float64
	xtMid  []float64
	slope  []float64
}

func (m *MidPointImplicitNewton) Name() string {
	return "MidPoint Implicit with Newton-Raphson Method"
}

func (m *MidPointImplicitNewton) NextStepIm(xt, t float64) (float64, error) {
	halfDt := m.DeltaT / 2
	tMid := t + halfDt
	xNext := xt
	for iter := 0; iter < maxIter; iter++ {
		xMid := (xt + xNext) / 2
		fMid := m.TdFunc.EvaluateAt(xMid, tMid)

		Gx := xNext - xt - m.DeltaT*fMid
		if math.Abs(Gx) < tolerance {
			return xNext, nil
		}

		fPlus := m.TdFunc.EvaluateAt(xMid+delta, tMid)
		dFbydx := (fPlus - fMid) / delta
		jacobian := 1 - halfDt*dFbydx
		xNext -= Gx / jacobian
	}
	return xt, fmt.Errorf("implicit step did not converge after %d iterations", maxIter)
}

// PredictorCorrector implements the basic Euler method
type PredictorCorrector struct {
	TdFunc gridData.TDPotentialOp
	DeltaT float64
}
