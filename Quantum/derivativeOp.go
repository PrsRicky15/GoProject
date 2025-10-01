package Quantum

import (
	"GoProject/gridData"
	"math"

	"golang.org/x/exp/constraints"
	"gonum.org/v1/gonum/blas"
	"gonum.org/v1/gonum/blas/cblas128"

	"gonum.org/v1/gonum/mat"
)

type Number interface {
	constraints.Signed | constraints.Float
}

type MatrixOp interface {
	Mat()
	GetMat() *mat.Dense
	Diagonalize() ([]float64, *mat.Dense, error)
}

type CanTimeSolverOp interface {
	ExpDt(At float64, In []float64) []float64
	ExpDtTo(At float64, In []float64, Out []float64)
	ExpDtInPlace(At float64, InOut []float64)
}

type TimeSolverOp interface {
	ExpDt(Dt float64, In []float64) []float64
	ExpDtTo(Dt float64, In []float64, Out []float64)
	ExpDtInPlace(Dt float64, InOut []float64)

	ExpIdt(Dt float64, In []complex64) []complex64
	ExpIdtTo(Dt float64, In []complex64, Out []complex64)
	ExpIdtInPlace(Dt float64, InOut []complex64)
}

type MomentumOp interface {
	MatrixOp
	TimeSolverOp
}

type KineticOp interface {
	MatrixOp
	TimeSolverOp
}

type CanMomentumOp interface {
	MatrixOp
	CanTimeSolverOp
}

type CanKineticOp interface {
	MatrixOp
	CanTimeSolverOp
}

type CanonicalMomDvrBasis struct {
	grid *gridData.RadGrid
	mass float64
}

type MomDvrBasis struct {
	grid *gridData.RadGrid
	mass float64
}

type CanonicalKeDvrBasis struct {
	grid  *gridData.RadGrid
	mass  float64
	kMat  *mat.Dense
	uMat  *mat.Dense
	keVal []float64
}

type KeDvrBasis struct {
	grid  *gridData.RadGrid
	mass  float64
	kMat  *mat.Dense
	uMat  *mat.Dense
	keVal []float64
}

func NewKeDVR(grid *gridData.RadGrid, mass float64) *KeDvrBasis {
	return &KeDvrBasis{
		grid: grid,
		mass: mass,
	}
}

func (k *KeDvrBasis) mInfinityToInfinity(keMat *mat.Dense) {
	_, dim := keMat.Dims()
	dx2 := k.grid.DeltaR() * k.grid.DeltaR()
	massDx2 := k.mass * dx2
	diagTerm := math.Pi * math.Pi / (6. * massDx2)
	invMassDx2 := 1.0 / massDx2

	keMat = mat.NewDense(dim, dim, nil)

	for i := 0; i < dim; i++ {
		keMat.Set(i, i, diagTerm)
	}

	for i := 1; i < dim; i++ {
		for j := 0; j < i; j++ {
			diff := i - j
			sign := float64(1 - 2*(diff&1))
			diffSq := float64(diff * diff)
			kEval := sign * invMassDx2 / diffSq
			keMat.Set(i, j, kEval)
			keMat.Set(j, i, kEval)
		}
	}
}

func (k *KeDvrBasis) zeroToInfinity(keMat *mat.Dense) {
	_, dim := keMat.Dims()
	dx2 := k.grid.DeltaR() * k.grid.DeltaR()
	massDx2 := k.mass * dx2
	diagTerm := math.Pi * math.Pi / (6. * massDx2)
	invMassDx2 := 1.0 / massDx2

	for i := 0; i < dim; i++ {
		keMat.Set(i, i, diagTerm-0.25/float64(i*i))
	}

	for i := 1; i < dim; i++ {
		for j := 0; j < i; j++ {
			diff := i - j
			add := i + j
			sign := float64(1 - 2*(diff&1))
			diffSq := float64(diff * diff)
			addSq := float64(add * add)
			kEval := sign * invMassDx2 * (1./diffSq - 1./addSq)
			keMat.Set(i, j, kEval)
			keMat.Set(j, i, kEval)
		}
	}
}

func (k *KeDvrBasis) Mat() {
	if math.Abs(k.grid.RMax()-k.grid.RMin()) <= 1 {
		k.zeroToInfinity(k.kMat)
	} else {
		k.mInfinityToInfinity(k.kMat)
	}
}

func (k *KeDvrBasis) GetMat() mat.Matrix { return k.kMat }
func (k *KeDvrBasis) Diagonalize() ([]float64, *mat.Dense, error) {
	return RealDiagonalize(k, int(k.grid.NPoints()))
}

func (k *KeDvrBasis) ExpDt(In []float64) []float64 {
	expMat := mat.NewDense(int(k.grid.NPoints()), int(k.grid.NPoints()), nil)
	expMat.Exp(k.kMat)
	inVec := mat.NewVecDense(len(In), In)
	outVec := mat.NewVecDense(len(In), nil)
	outVec.MulVec(expMat, inVec)

	result := make([]float64, len(In))
	for i := 0; i < len(In); i++ {
		result[i] = outVec.AtVec(i)
	}
	return result
}

func (k *KeDvrBasis) ExpDtTo(In []float64, Out []float64) {
	expMat := mat.NewDense(int(k.grid.NPoints()), int(k.grid.NPoints()), nil)
	expMat.Exp(k.kMat)
	inVec := mat.NewVecDense(len(In), In)
	outVec := mat.NewVecDense(len(In), nil)
	outVec.MulVec(expMat, inVec)

	for i := 0; i < len(In); i++ {
		Out[i] = outVec.AtVec(i)
	}
}

func (k *KeDvrBasis) ExpDtInPlace(InOut []float64) {
	k.ExpDtTo(InOut, InOut)
}

func (k *KeDvrBasis) ExpIdt(In []complex128) []complex128 {
	expMat := mat.NewCDense(int(k.grid.NPoints()), int(k.grid.NPoints()), nil)
	result := make([]complex128, len(In))

	cblas128.Gemv(
		blas.NoTrans,
		1,
		expMat.RawCMatrix(),
		cblas128.Vector{N: len(In), Data: In, Inc: 1},
		0, // beta = 0
		cblas128.Vector{N: len(result), Data: result, Inc: 1},
	)

	return result
}

func (k *KeDvrBasis) ExpIdtTo(In []complex128, Out []complex128) {
	expMat := mat.NewDense(int(k.grid.NPoints()), int(k.grid.NPoints()), nil)

	cblas128.Gemv(
		blas.NoTrans,
		1,
		expMat.RawCMatrix(),
		cblas128.Vector{N: len(In), Data: In, Inc: 1},
		0, // beta = 0
		cblas128.Vector{N: len(Out), Data: Out, Inc: 1},
	)
}

func (k *KeDvrBasis) ExpIdtInPlace(InOut []float64) {
	k.ExpDtTo(InOut, InOut)
}
