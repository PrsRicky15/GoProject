package BasicOneD

import (
	"GoProject/gridData"
	"fmt"
	"math"
	"math/cmplx"
	"os"

	"gonum.org/v1/gonum/mat"
)

type FiniteBarrier struct {
	grid  *gridData.RadGrid
	rFunc gridData.Rfunc
	mass  float64

	nDelta uint
	start  float64
	end    float64
	width  float64
	v0     []float64
}

func NewFiniteBarrier(grid *gridData.RadGrid, rFunc gridData.Rfunc, mass float64) *FiniteBarrier {

	if grid == nil || rFunc == nil {
		panic("grid and rFunc must not be nil")
	}

	return &FiniteBarrier{
		grid:  grid,
		rFunc: rFunc,
		mass:  mass,
	}
}

func (fb *FiniteBarrier) ReDefine(grid *gridData.RadGrid, rFunc gridData.Rfunc) {
	fb.rFunc = rFunc
	fb.grid = grid
}

func (fb *FiniteBarrier) ReDefineFunc(rFunc gridData.Rfunc, nDelta uint) {
	fb.rFunc = rFunc
	fb.nDelta = nDelta
}

func (fb *FiniteBarrier) MaxMinPotent() {
	for i := uint32(0); i < fb.grid.NPoints(); i++ {
		x1 := fb.grid.RMin() + float64(i)*fb.grid.DeltaR()
		x2 := fb.grid.RMax() - float64(i)*fb.grid.DeltaR()
		fx1 := fb.rFunc.EvaluateAt(x1)
		fx2 := fb.rFunc.EvaluateAt(x2)
		if fx1 > 1e-07 && fx2 > 1e-07 {
			fb.start = x2
			fb.end = x1
		}
	}
}

func (fb *FiniteBarrier) DisplayMinMaxBarrier() {
	fmt.Printf("Start: %f\n", fb.start)
	fmt.Printf("End: %f\n", fb.end)
	fmt.Printf("Width: %f\n", fb.width)
	fmt.Printf("N-barrier: %d\n", fb.nDelta)

	fmt.Printf("Barrier v0: %17.7e\n", fb.v0)
}

func (fb *FiniteBarrier) ConvertFuncToBarrier(nDelta uint) {
	fb.nDelta = nDelta
	fb.v0 = make([]float64, nDelta)
	fb.MaxMinPotent()
	fb.width = (fb.end - fb.start) / float64(nDelta)

	for i := uint32(0); i < fb.grid.NPoints(); i++ {
		x := fb.start + float64(i)*fb.width
		for ia := range fb.v0 {
			if math.Abs(x-(fb.start+float64(ia)*fb.width)) < 1e-07 {
				fb.v0[ia] = fb.rFunc.EvaluateAt(x)
			}
		}
	}
}

func (fb *FiniteBarrier) barrierOnGrid(Pot []float64) {
	for i := uint32(0); i < fb.grid.NPoints(); i++ {
		x := fb.grid.RMin() + float64(i)*fb.grid.DeltaR()
		for ia := range fb.v0 {
			barrierStart := fb.start + float64(ia)*fb.width
			barrierEnd := barrierStart + fb.width

			if x >= barrierStart && x < barrierEnd {
				Pot[i] = fb.v0[ia]
				break
			}
		}
	}
}

func (fb *FiniteBarrier) DisplayFuncToDeltaToFile(format string) {
	barPot := make([]float64, fb.grid.NPoints())
	fb.barrierOnGrid(barPot)
	fullFormat := "%14.7e" + "\t" + format + "\t" + format + "\n"
	for i := uint32(0); i < fb.grid.NPoints(); i++ {
		x := fb.grid.RMin() + float64(i)*fb.grid.DeltaR()
		fmt.Printf(fullFormat, x, fb.rFunc.EvaluateAt(x), barPot[i])
	}
}

func (fb *FiniteBarrier) PrintFuncToDeltaToFile(filename string, format string) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer func(file *os.File) {
		err := file.Close()
		if err != nil {

		}
	}(file)

	barPot := make([]float64, fb.grid.NPoints())
	fb.barrierOnGrid(barPot)

	fullFormat := "%14.7e" + "\t" + format + "\t" + format + "\n"
	_, err = fmt.Fprintf(file, "#--------------------------------------------------\n")
	_, err = fmt.Fprintf(file, "#\t\t grid\t\t function value\n")
	_, err = fmt.Fprintf(file, "#--------------------------------------------------\n")
	for i := uint32(0); i < fb.grid.NPoints(); i++ {
		var x = fb.grid.RMin() + float64(i)*fb.grid.DeltaR()
		_, err := fmt.Fprintf(file, fullFormat, x, fb.rFunc.EvaluateAt(x), barPot[i])
		if err != nil {
			return err
		}
	}
	return nil
}

func (fb *FiniteBarrier) kVal(energy, v0 float64) complex128 {
	diff := energy - v0
	if diff >= 0 {
		return complex(math.Sqrt(2*fb.mass*diff), 0)
	}
	return complex(0, math.Sqrt(2*fb.mass*(-diff)))
}

func (fb *FiniteBarrier) mMatrix(kL, kR float64) mat.CMatrix {

	addK := complex((kL+kR)/(2*kR), 0)
	subK := complex((kR-kL)/(2*kR), 0)

	addKLen := complex(0, (kL+kR)*fb.grid.Length())
	subKLen := complex(0, (kL-kR)*fb.grid.Length())

	mMat := mat.NewCDense(2, 2, nil)
	mMat.Set(0, 0, addK*cmplx.Exp(subKLen))
	mMat.Set(0, 1, subK*cmplx.Exp(-addKLen))
	mMat.Set(1, 0, subK*cmplx.Exp(addKLen))
	mMat.Set(1, 1, addK*cmplx.Exp(subKLen))

	return mMat
}

func (fb *FiniteBarrier) scattering(kL, kR float64, left []complex128) []complex128 {
	mMat := fb.mMatrix(kL, kR)

	result := make([]complex128, 2)
	result[0] = mMat.At(0, 0)*left[0] + mMat.At(0, 1)*left[1]
	result[1] = mMat.At(1, 0)*left[0] + mMat.At(1, 1)*left[1]

	return result
}

type TDFiniteBarrier struct {
	barrier *FiniteBarrier
	tdGrid  *gridData.TimeGrid
}
