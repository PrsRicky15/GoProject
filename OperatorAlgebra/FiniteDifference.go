package OperatorAlgebra

import (
	"GoProject/gridData"
	"fmt"

	"gonum.org/v1/gonum/mat"
)

type FiniteDiff struct {
	grid         *gridData.RadGrid
	kineticCoeff float64
	mass         float64
	order        int
	stencil      []float64
}

func (fd *FiniteDiff) RealDiagonalize() (eigenvalues []float64, eigenvectors *mat.Dense, err error) {
	//TODO implement me
	panic("implement me")
}

func (fd *FiniteDiff) ExpDtTo(Dt float64, In []float64, Out []float64) error {
	//TODO implement me
	panic("implement me")
}

func (fd *FiniteDiff) ExpDtInPlace(Dt float64, InOut []float64) error {
	//TODO implement me
	panic("implement me")
}

func (fd *FiniteDiff) ExpIdtTo(Dt float64, In []complex128, Out []complex128) error {
	//TODO implement me
	panic("implement me")
}

func (fd *FiniteDiff) ExpIdtInPlace(Dt float64, InOut []complex128) error {
	//TODO implement me
	panic("implement me")
}

var stencilCoeffs = map[int][]float64{
	3: {-2.0, 1.0},                                                       // 3-point
	5: {-5.0 / 2.0, 4.0 / 3.0, -1.0 / 12.0},                              // 5-point
	7: {-49.0 / 18.0, 3.0 / 2.0, -3.0 / 20.0, 1.0 / 90.0},                // 7-point
	9: {-205.0 / 72.0, 8.0 / 5.0, -1.0 / 5.0, 8.0 / 315.0, -1.0 / 560.0}, // 9-point
}

func NewFiniteDiff(grid *gridData.RadGrid, mass float64, order int) (*FiniteDiff, error) {
	if order != 3 && order != 5 && order != 6 {
		return nil, fmt.Errorf("order must be 3, 5, or 6, got %d", order)
	}

	dxSq := grid.DeltaR() * grid.DeltaR()
	kineticCoeff := -1 / (2 * mass * dxSq)

	fd := &FiniteDiff{
		grid:         grid,
		kineticCoeff: kineticCoeff,
		order:        order,
	}
	return fd, nil
}

func NewFiniteDiff2(grid *gridData.RadGrid, mass float64, order int) (*FiniteDiff, error) {
	stencil, ok := stencilCoeffs[order]
	if !ok {
		return nil, fmt.Errorf("order must be 3, 5, 7, or 9, got %d", order)
	}

	dxSq := grid.DeltaR() * grid.DeltaR()
	kineticCoeff := -1.0 / (2.0 * mass * dxSq)

	fd := &FiniteDiff{
		grid:         grid,
		kineticCoeff: kineticCoeff,
		mass:         mass,
		order:        order,
		stencil:      stencil,
	}
	return fd, nil
}

func (fd *FiniteDiff) ComputeMatrix(keMat mat.Matrix) {
	switch fd.order {
	case 3:
		fd.compute3rdOrderCoefficients(keMat)
	case 5:
		fd.compute5thOrderCoefficients(keMat)
	case 7:
		fd.compute7thOrderCoefficients(keMat)
	case 9:
		fd.compute9thOrderCoefficients(keMat)
	}
}

func (fd *FiniteDiff) compute3rdOrderCoefficients(keMat mat.Matrix) {
	ndims := int(fd.grid.NPoints())
	coeff := 2.0 * fd.kineticCoeff

	for i := 0; i < ndims; i++ {
		keMat.(*mat.Dense).Set(i, i, -coeff)
	}

	for i := 0; i < ndims-1; i++ {
		keMat.(*mat.Dense).Set(i+1, i, fd.kineticCoeff)
		keMat.(*mat.Dense).Set(i, i+1, fd.kineticCoeff)
	}
}

func (fd *FiniteDiff) compute5thOrderCoefficients(keMat mat.Matrix) {
	ndims := int(fd.grid.NPoints())
	k := fd.kineticCoeff

	for i := 0; i < ndims; i++ {
		keMat.(*mat.Dense).Set(i, i, 2.5*k)
	}

	for i := 0; i < ndims-1; i++ {
		keMat.(*mat.Dense).Set(i, i+1, -4.0/3.0*k)
		keMat.(*mat.Dense).Set(i+1, i, -4.0/3.0*k)
	}

	for i := 0; i < ndims-2; i++ {
		keMat.(*mat.Dense).Set(i, i+2, 1.0/12.0*k)
		keMat.(*mat.Dense).Set(i+2, i, 1.0/12.0*k)
	}
}

func (fd *FiniteDiff) compute7thOrderCoefficients(keMat mat.Matrix) {
	ndims := int(fd.grid.NPoints())
	k := fd.kineticCoeff

	for i := 0; i < ndims; i++ {
		keMat.(*mat.Dense).Set(i, i, 49.0/18.0*k)
	}

	for i := 0; i < ndims-1; i++ {
		keMat.(*mat.Dense).Set(i, i+1, -1.5*k)
		keMat.(*mat.Dense).Set(i+1, i, -1.5*k)
	}

	for i := 0; i < ndims-2; i++ {
		keMat.(*mat.Dense).Set(i, i+2, 3.0/20.0*k)
		keMat.(*mat.Dense).Set(i+2, i, 3.0/20.0*k)
	}

	for i := 0; i < ndims-3; i++ {
		keMat.(*mat.Dense).Set(i, i+3, -1.0/90.0*k)
		keMat.(*mat.Dense).Set(i+3, i, -1.0/90.0*k)
	}
}

// 9-point (8th-order)
func (fd *FiniteDiff) compute9thOrderCoefficients(keMat mat.Matrix) {
	ndims := int(fd.grid.NPoints())
	k := fd.kineticCoeff

	for i := 0; i < ndims; i++ {
		keMat.(*mat.Dense).Set(i, i, 205.0/72.0*k)
	}

	for i := 0; i < ndims-1; i++ {
		keMat.(*mat.Dense).Set(i, i+1, -8.0/5.0*k)
		keMat.(*mat.Dense).Set(i+1, i, -8.0/5.0*k)
	}

	for i := 0; i < ndims-2; i++ {
		keMat.(*mat.Dense).Set(i, i+2, 1.0/5.0*k)
		keMat.(*mat.Dense).Set(i+2, i, 1.0/5.0*k)
	}

	for i := 0; i < ndims-3; i++ {
		keMat.(*mat.Dense).Set(i, i+3, -8.0/315.0*k)
		keMat.(*mat.Dense).Set(i+3, i, -8.0/315.0*k)
	}

	for i := 0; i < ndims-4; i++ {
		keMat.(*mat.Dense).Set(i, i+4, 1.0/560.0*k)
		keMat.(*mat.Dense).Set(i+4, i, 1.0/560.0*k)
	}
}

func (fd *FiniteDiff) ComputeMatrixAnotherWay(keMat *mat.Dense) {
	ndims := int(fd.grid.NPoints())
	k := fd.kineticCoeff

	diagCoeff := -fd.stencil[0] * k
	for i := 0; i < ndims; i++ {
		keMat.Set(i, i, diagCoeff)
	}

	for offset := 1; offset < len(fd.stencil); offset++ {
		offDiagCoeff := -fd.stencil[offset] * k
		for i := 0; i < ndims-offset; i++ {
			keMat.Set(i, i+offset, offDiagCoeff)
			keMat.Set(i+offset, i, offDiagCoeff) // Symmetric
		}
	}
}

func (fd *FiniteDiff) ComputeSymMatrix() *mat.SymDense {
	ndims := int(fd.grid.NPoints())
	k := fd.kineticCoeff
	sym := mat.NewSymDense(ndims, nil)

	diagCoeff := -fd.stencil[0] * k
	for i := 0; i < ndims; i++ {
		sym.SetSym(i, i, diagCoeff)
	}

	for offset := 1; offset < len(fd.stencil); offset++ {
		offDiagCoeff := -fd.stencil[offset] * k
		for i := 0; i < ndims-offset; i++ {
			sym.SetSym(i, i+offset, offDiagCoeff)
		}
	}

	return sym
}

func (fd *FiniteDiff) ComputeSymMatrixInPlace(keMat *mat.SymDense) {
	ndims := int(fd.grid.NPoints())
	k := fd.kineticCoeff

	diagCoeff := -fd.stencil[0] * k
	for i := 0; i < ndims; i++ {
		keMat.SetSym(i, i, diagCoeff)
	}

	for offset := 1; offset < len(fd.stencil); offset++ {
		offDiagCoeff := -fd.stencil[offset] * k
		for i := 0; i < ndims-offset; i++ {
			keMat.SetSym(i, i+offset, offDiagCoeff)
		}
	}
}
