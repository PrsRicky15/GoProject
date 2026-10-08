package OperatorAlgebra

import (
	"GoProject/gridData"
	"fmt"
	"runtime"
	"sync"

	"github.com/jvlmdr/go-fftw/fftw"
)

// FFTW execution is independent per instance, but planning/destruction is global.
var fftwPlannerMu sync.Mutex

// FourierBasis represents the kinetic energy operator in Fourier/DVR basis
type FourierBasis struct {
	grid     *gridData.RadGrid
	mass     float64
	nPoints  int
	fftPlan  fftw.Plan
	ifftPlan fftw.Plan
	Buff     *fftw.Array
	kValues  []float64
	keValues []float64
}

// FFTInit creates a new Fourier struct for fast fourier transform
func FFTInit(grid *gridData.RadGrid, mass float64) *FourierBasis {
	fb := &FourierBasis{}
	fb.initialize(grid, mass)
	runtime.SetFinalizer(fb, (*FourierBasis).destroy)
	return fb
}

func (f *FourierBasis) Redefine(grid *gridData.RadGrid, mass float64) {
	f.destroy()
	f.initialize(grid, mass)
}

func (f *FourierBasis) initialize(grid *gridData.RadGrid, mass float64) {
	fftwPlannerMu.Lock()
	defer fftwPlannerMu.Unlock()
	gridPoints := int(grid.NPoints())
	kVal := grid.KValues()

	kValues := make([]float64, gridPoints)
	keValues := make([]float64, gridPoints)
	invTwoMass := 1.0 / (2.0 * mass)

	for i := 0; i < gridPoints; i++ {
		k := kVal[i]
		kValues[i] = k
		keValues[i] = k * k * invTwoMass
	}

	buff := fftw.NewArray(gridPoints)

	fftPlan := *fftw.NewPlan(buff, buff, fftw.Forward, fftw.Estimate)
	ifftPlan := *fftw.NewPlan(buff, buff, fftw.Backward, fftw.Estimate)

	f.grid = grid
	f.mass = mass
	f.nPoints = gridPoints
	f.fftPlan = fftPlan
	f.ifftPlan = ifftPlan
	f.Buff = buff
	f.kValues = kValues
	f.keValues = keValues
}

func (f *FourierBasis) forwardBuff(in []complex128) {
	copy(f.Buff.Elems, in)
	f.fftPlan.Execute()
}

func (f *FourierBasis) forward(in []complex128, out []complex128) {
	copy(f.Buff.Elems, in)
	f.fftPlan.Execute()
	copy(out, f.Buff.Elems)
}

func (f *FourierBasis) forwardInPlace(InOut []complex128) {
	copy(f.Buff.Elems, InOut)
	f.fftPlan.Execute()
	copy(InOut, f.Buff.Elems)
}

func (f *FourierBasis) backwardBuff(in []complex128) {
	copy(f.Buff.Elems, in)
	f.ifftPlan.Execute()
}

func (f *FourierBasis) backward(in []complex128, out []complex128) {
	copy(f.Buff.Elems, in)
	f.ifftPlan.Execute()
	copy(out, f.Buff.Elems)
}

func (f *FourierBasis) backwardInPlace(InOut []complex128) {
	copy(f.Buff.Elems, InOut)
	f.ifftPlan.Execute()
	copy(InOut, f.Buff.Elems)
}

func (f *FourierBasis) operatorOp(InOut []complex128, Op []float64) {
	if len(InOut) != f.nPoints {
		panic(fmt.Sprintf("FourierBasis: vector length %d, want %d", len(InOut), f.nPoints))
	}
	if len(f.Buff.Elems) != len(Op) {
		panic(fmt.Sprintf("length mismatch: Buff.Elems=%d, kValues=%d",
			len(f.Buff.Elems), len(Op)))
	}

	copy(f.Buff.Elems, InOut)
	f.fftPlan.Execute()

	// FFTW leaves both transforms unnormalized, as in h2p's kinetic operator.
	invN := 1 / float64(f.nPoints)
	for i := range f.Buff.Elems {
		f.Buff.Elems[i] *= complex(Op[i]*invN, 0)
	}

	f.ifftPlan.Execute()
	copy(InOut, f.Buff.Elems)
}

func (f *FourierBasis) MomentumOpInPlace(InOut []complex128) {
	f.operatorOp(InOut, f.kValues)
}

func (f *FourierBasis) MomentumOp(In []complex128, Out []complex128) {
	f.checkVectors(In, Out)
	copy(Out, In)
	f.MomentumOpInPlace(Out)
}

// LaplacianOpInPlace applies kinetic energy -d²/dx²/(2m), retaining its legacy name.
func (f *FourierBasis) LaplacianOpInPlace(InOut []complex128) {
	f.operatorOp(InOut, f.keValues)
}

func (f *FourierBasis) LaplacianOp(In []complex128, Out []complex128) {
	f.checkVectors(In, Out)
	copy(Out, In)
	f.LaplacianOpInPlace(Out)
}

func (f *FourierBasis) checkVectors(in, out []complex128) {
	if len(in) != f.nPoints || len(out) != f.nPoints {
		panic(fmt.Sprintf("FourierBasis: vector lengths %d and %d, want %d", len(in), len(out), f.nPoints))
	}
}

func (f *FourierBasis) destroy() {
	if f != nil {
		fftwPlannerMu.Lock()
		defer fftwPlannerMu.Unlock()
		(*fftw.Plan).Destroy(&f.fftPlan)
		(*fftw.Plan).Destroy(&f.ifftPlan)
	}
}

// Clean cleans up FFTW
func (f *FourierBasis) Clean() {
	runtime.SetFinalizer(f, nil)
	f.destroy()
}
