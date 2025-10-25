package OperatorAlgebra

import (
	"GoProject/gridData"
	"fmt"
	"math"
	"math/cmplx"

	"gonum.org/v1/gonum/blas"
	"gonum.org/v1/gonum/blas/cblas128"
	"gonum.org/v1/gonum/mat"
)

// KeDvrBasis represents the kinetic energy operator in DVR basis with caching
type KeDvrBasis struct {
	grid        *gridData.RadGrid
	mass        float64
	ndims       int
	dx2         float64
	massDx2     float64
	invMassDx2  float64
	diagTerm    float64
	diagTermZ2I float64
	kMat        *mat.Dense
	kMatCached  bool
}

func (k *KeDvrBasis) ExpDtTo(Dt float64, In []float64, Out []float64) {
	//TODO implement me
	panic("implement me")
}

func (k *KeDvrBasis) ExpIdtTo(Dt float64, In []complex128, Out []complex128) {
	//TODO implement me
	panic("implement me")
}

// NewKeDVR creates and initializes a new kinetic energy DVR basis
func NewKeDVR(grid *gridData.RadGrid, mass float64) *KeDvrBasis {
	ndim := int(grid.NPoints())
	dx2 := grid.DeltaR() * grid.DeltaR()
	massDx2 := mass * dx2
	invMassDx2 := 1.0 / massDx2
	diagTerm := math.Pi * math.Pi / (6.0 * massDx2)

	return &KeDvrBasis{
		grid:        grid,
		mass:        mass,
		ndims:       ndim,
		dx2:         dx2,
		massDx2:     massDx2,
		invMassDx2:  invMassDx2,
		diagTerm:    diagTerm,
		diagTermZ2I: diagTerm - 0.25, // Pre-compute for zero-to-infinity
		kMat:        mat.NewDense(ndim, ndim, nil),
		kMatCached:  false,
	}
}

// isZeroToInfinity determines if the grid domain is [0, infinity)
func (k *KeDvrBasis) isZeroToInfinity() bool {
	return math.Abs(k.grid.RMax()+k.grid.RMin()) >= 1
}

// buildRealMinInftyToInfty constructs the kinetic energy matrix for (-∞, ∞) domain
func (k *KeDvrBasis) buildRealMinInftyToInfty(mat *mat.Dense) {
	for i := 0; i < k.ndims; i++ {
		mat.Set(i, i, k.diagTerm)
	}

	for i := 1; i < k.ndims; i++ {
		for j := 0; j < i; j++ {
			diff := i - j
			sign := float64(1 - 2*(diff&1))
			diffSq := float64(diff * diff)
			val := sign * k.invMassDx2 / diffSq

			mat.Set(i, j, val)
			mat.Set(j, i, val)
		}
	}
}

// buildRealZeroToInfinity constructs the kinetic energy matrix for [0, ∞) domain
func (k *KeDvrBasis) buildRealZeroToInfinity(mat *mat.Dense) {
	for i := 0; i < k.ndims; i++ {
		diagVal := k.diagTerm
		if i > 0 {
			diagVal += k.diagTermZ2I / float64(i*i)
		}
		mat.Set(i, i, diagVal)
	}

	for i := 1; i < k.ndims; i++ {
		for j := 0; j < i; j++ {
			diff := i - j
			add := i + j
			sign := float64(1 - 2*(diff&1))
			diffSq := float64(diff * diff)
			addSq := float64(add * add)

			val := sign * k.invMassDx2 * (1.0/diffSq - 1.0/addSq)

			mat.Set(i, j, val)
			mat.Set(j, i, val)
		}
	}
}

// GetMat returns the kinetic energy matrix, using cache if available
func (k *KeDvrBasis) GetMat() *mat.Dense {
	if !k.kMatCached {
		if k.isZeroToInfinity() {
			k.buildRealZeroToInfinity(k.kMat)
		} else {
			k.buildRealMinInftyToInfty(k.kMat)
		}
		k.kMatCached = true
	}
	return k.kMat
}

// scaleMatrixComplexParam scales a real matrix by exp(-2iθ) and stores in complex matrix
func scaleMatrixComplexParam(realMat *mat.Dense, complexMat *mat.CDense, theta float64) {
	rows, cols := realMat.Dims()
	expFactor := cmplx.Exp(-2i * complex(theta, 0))

	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			val := complex(realMat.At(i, j), 0) * expFactor
			complexMat.Set(i, j, val)
		}
	}
}

// GetComplexMat returns the complex-scaled kinetic energy matrix
func (k *KeDvrBasis) GetComplexMat(theta float64) *mat.CDense {
	realMat := k.GetMat()
	complexMat := mat.NewCDense(k.ndims, k.ndims, nil)
	scaleMatrixComplexParam(realMat, complexMat, theta)
	return complexMat
}

// RealDiagonalize diagonalizes the real kinetic energy matrix
func (k *KeDvrBasis) RealDiagonalize() (eigenvalues []float64, eigenvectors *mat.Dense, err error) {
	eigenvalues = make([]float64, k.ndims)
	eigenvectors = mat.NewDense(k.ndims, k.ndims, nil)
	eigenvectors = k.GetMat()
	err = RealDiagonalizeLapack(eigenvectors, eigenvalues)
	return eigenvalues, eigenvectors, err
}

// ExpDt computes exp(Dt * K) and applies it to input vector
func (k *KeDvrBasis) ExpDt(dt float64, in []float64) ([]float64, error) {
	if len(in) != k.ndims {
		return nil, fmt.Errorf("input vector length %d doesn't match basis dimension %d", len(in), k.ndims)
	}

	kMat := k.GetMat()
	scaledMat := mat.NewDense(k.ndims, k.ndims, nil)
	scaledMat.Scale(dt, kMat)

	expMat := mat.NewDense(k.ndims, k.ndims, nil)
	expMat.Exp(scaledMat)

	result := make([]float64, k.ndims)
	resultVec := mat.NewVecDense(k.ndims, result)
	inVec := mat.NewVecDense(k.ndims, in)
	resultVec.MulVec(expMat, inVec)

	return result, nil
}

// ExpDtInPlace computes exp(Dt * K) and applies it to input vector in-place
func (k *KeDvrBasis) ExpDtInPlace(dt float64, inOut []float64) error {
	if len(inOut) != k.ndims {
		return fmt.Errorf("vector length %d doesn't match basis dimension %d", len(inOut), k.ndims)
	}

	kMat := k.GetMat()
	scaledMat := mat.NewDense(k.ndims, k.ndims, nil)
	scaledMat.Scale(dt, kMat)

	expMat := mat.NewDense(k.ndims, k.ndims, nil)
	expMat.Exp(scaledMat)

	temp := make([]float64, k.ndims)
	tempVec := mat.NewVecDense(k.ndims, temp)
	inOutVec := mat.NewVecDense(k.ndims, inOut)
	tempVec.MulVec(expMat, inOutVec)

	copy(inOut, temp)
	return nil
}

// ExpIdt computes exp(i*Dt*K) and applies it to complex input vector
// This implements the complex-scaled exponential: exp(i*Dt*K)
func (k *KeDvrBasis) ExpIdt(dt float64, in []complex128) ([]complex128, error) {
	if len(in) != k.ndims {
		return nil, fmt.Errorf("input vector length %d doesn't match basis dimension %d", len(in), k.ndims)
	}

	kMat := k.GetMat()
	scaledMat := mat.NewDense(k.ndims, k.ndims, nil)
	scaledMat.Scale(dt, kMat)

	complexScaled := mat.NewCDense(k.ndims, k.ndims, nil)
	for i := 0; i < k.ndims; i++ {
		for j := 0; j < k.ndims; j++ {
			complexScaled.Set(i, j, 1i*complex(scaledMat.At(i, j), 0))
		}
	}

	expMat := mat.NewCDense(k.ndims, k.ndims, nil)
	//	expMat.Exp(complexScaled)

	result := make([]complex128, k.ndims)
	cblas128.Gemv(
		blas.NoTrans,
		complex(1, 0),
		expMat.RawCMatrix(),
		cblas128.Vector{N: len(in), Data: in, Inc: 1},
		complex(0, 0),
		cblas128.Vector{N: len(result), Data: result, Inc: 1},
	)

	return result, nil
}

func (k *KeDvrBasis) ExpIdtInPlace(dt float64, inOut []complex128) error {
	if len(inOut) != k.ndims {
		return fmt.Errorf("vector length %d doesn't match basis dimension %d", len(inOut), k.ndims)
	}

	result, err := k.ExpIdt(dt, inOut)
	if err != nil {
		return err
	}

	copy(inOut, result)
	return nil
}

func (k *KeDvrBasis) Clone() *KeDvrBasis {
	newK := &KeDvrBasis{
		grid:        k.grid,
		mass:        k.mass,
		ndims:       k.ndims,
		dx2:         k.dx2,
		massDx2:     k.massDx2,
		invMassDx2:  k.invMassDx2,
		diagTerm:    k.diagTerm,
		diagTermZ2I: k.diagTermZ2I,
		kMat:        mat.NewDense(k.ndims, k.ndims, nil),
		kMatCached:  false,
	}
	if k.kMatCached {
		newK.kMat.Copy(k.kMat)
		newK.kMatCached = true
	}
	return newK
}

func (k *KeDvrBasis) Clear() {
	k.kMatCached = false
}
