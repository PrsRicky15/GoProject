package OperatorAlgebra

import "fmt"

// ApplyComplexSpectrum applies an FFT-ordered diagonal multiplier with the
// inverse FFT normalization. It shares the basis's work buffer and lifecycle.
// The caller owns the multiplier; it is neither changed nor retained.
func (f *FourierBasis) ApplyComplexSpectrum(inOut, multiplier []complex128) {
	if len(inOut) != f.nPoints || len(multiplier) != f.nPoints {
		panic(fmt.Sprintf("FourierBasis: spectrum/vector lengths must equal %d", f.nPoints))
	}
	copy(f.Buff.Elems, inOut)
	f.fftPlan.Execute()
	for i := range f.Buff.Elems {
		f.Buff.Elems[i] *= multiplier[i] / complex(float64(f.nPoints), 0)
	}
	f.ifftPlan.Execute()
	copy(inOut, f.Buff.Elems)
}
