package OperatorAlgebra

import (
	"GoProject/gridData"
	"errors"
	"fmt"
	"math"

	"gonum.org/v1/gonum/blas"
	"gonum.org/v1/gonum/blas/blas64"
	"gonum.org/v1/gonum/mat"
)

// KeDvrBasis represents the kinetic energy operator
type KeDvrBasis struct {
	grid       *gridData.RadGrid
	mass       float64
	ndims      int
	dx2        float64
	invMassDx2 float64
	diagTerm   float64
	kMat       *mat.Dense

	kMatCached bool
	eigVals    []float64
	eigVecs    *mat.Dense
}

// NewKeDVR creates and initializes DVR basis
func NewKeDVR(grid *gridData.RadGrid, mass float64) *KeDvrBasis {
	ndim := int(grid.NPoints())
	dx2 := grid.DeltaR() * grid.DeltaR()
	invMassDx2 := 1.0 / (mass * dx2)
	diagTerm := math.Pi * math.Pi / (6.0 * mass * dx2)

	return &KeDvrBasis{
		grid:       grid,
		mass:       mass,
		ndims:      ndim,
		dx2:        dx2,
		invMassDx2: invMassDx2,
		diagTerm:   diagTerm,
		kMat:       mat.NewDense(ndim, ndim, nil),
		kMatCached: false,
	}
}

func (k *KeDvrBasis) String() string {
	return fmt.Sprintf("KeDVR(xmin:\t%14.7f\n, "+
		"xmax:\t %14.7f\n, "+
		"mass:\t %14.7f\n,"+
		" Dim:%d\n)", k.grid.RMin(), k.grid.RMax(), k.mass, k.ndims)
}

func (k *KeDvrBasis) Ndims() int {
	return k.ndims
}

func (k *KeDvrBasis) Mass() float64 {
	return k.mass
}

func (k *KeDvrBasis) Clear() {
	k.kMatCached = false
	k.eigVals = nil
	k.eigVecs = nil
}

func (k *KeDvrBasis) isZeroToInfinity() bool {
	return math.Abs(k.grid.RMax()+k.grid.RMin()) >= 1
}

func (k *KeDvrBasis) MinInftyToInfty(m *mat.Dense) {
	for i := 0; i < k.ndims; i++ {
		m.Set(i, i, k.diagTerm)
	}
	for i := 1; i < k.ndims; i++ {
		for j := 0; j < i; j++ {
			diff := i - j
			sign := float64(1 - 2*(diff&1))
			val := sign * k.invMassDx2 / float64(diff*diff)
			m.Set(i, j, val)
			m.Set(j, i, val)
		}
	}
}

func (k *KeDvrBasis) ZeroToInfinity(m *mat.Dense) {
	constVal := 0.25 * k.invMassDx2
	for i := 0; i < k.ndims; i++ {
		m.Set(i, i, k.diagTerm-constVal/float64((i+1)*(i+1)))
	}

	for i := 1; i < k.ndims; i++ {
		for j := 0; j < i; j++ {
			diff := i - j
			add := i + j
			sign := float64(1 - 2*(diff&1))
			diffSq := float64(diff * diff)
			addSq := float64(add * add)
			val := sign * k.invMassDx2 * (1.0/diffSq - 1.0/addSq)
			m.Set(i, j, val)
			m.Set(j, i, val)
		}
	}
}

func (k *KeDvrBasis) SetMat() {
	if !k.kMatCached {
		if k.isZeroToInfinity() {
			k.ZeroToInfinity(k.kMat)
		} else {
			k.MinInftyToInfty(k.kMat)
		}
		k.kMatCached = true
	}
}

func (k *KeDvrBasis) GetMat() *mat.Dense {
	k.SetMat()
	return k.kMat
}

func (k *KeDvrBasis) EvalEigenVectors() error {
	if k.eigVals != nil && k.eigVecs != nil {
		return nil
	}

	k.eigVals = make([]float64, k.ndims)
	k.eigVecs = mat.NewDense(k.ndims, k.ndims, nil)
	KM := k.GetMat()

	sym := mat.NewSymDense(k.ndims, nil)
	for i := 0; i < k.ndims; i++ {
		for j := i; j < k.ndims; j++ {
			sym.SetSym(i, j, KM.At(i, j))
		}
	}

	var es mat.EigenSym
	ok := es.Factorize(sym, true)
	if !ok {
		return errors.New("eigendecomposition failed (EigenSym.Factorize returned false)")
	}

	k.eigVals = es.Values(k.eigVals)
	es.VectorsTo(k.eigVecs)

	return nil
}

func (k *KeDvrBasis) GetEigenValuesVectors(eigenvalues []float64, eigenvectors *mat.Dense) error {
	if err := k.EvalEigenVectors(); err != nil {
		return err
	}

	if eigenvalues != nil && len(eigenvalues) != k.ndims {
		return errors.New("eigenvalues array length mismatch")
	}
	if eigenvectors != nil {
		r, c := eigenvectors.Dims()
		if r != k.ndims || c != k.ndims {
			return errors.New("eigenvectors matrix dimensions mismatch")
		}
	}

	if eigenvalues != nil {
		copy(eigenvalues, k.eigVals)
	}
	if eigenvectors != nil {
		eigenvectors.Copy(k.eigVecs)
	}

	return nil
}

func (k *KeDvrBasis) ExpDtTo(Dt float64, In []float64, Out []float64) error {
	if len(In) != k.ndims || len(Out) != k.ndims {
		return fmt.Errorf("vector length mismatch: got In=%d Out=%d, expected %d", len(In), len(Out), k.ndims)
	}

	if err := k.EvalEigenVectors(); err != nil {
		return err
	}

	V := k.eigVecs
	lams := k.eigVals

	tmp := make([]float64, k.ndims)

	vData := V.RawMatrix()
	blas64.Gemv(blas.Trans, 1.0, vData, blas64.Vector{N: k.ndims, Data: In, Inc: 1},
		0.0, blas64.Vector{N: k.ndims, Data: tmp, Inc: 1})

	for i := 0; i < k.ndims; i++ {
		tmp[i] *= math.Exp(Dt * lams[i])
	}

	blas64.Gemv(blas.NoTrans, 1.0, vData, blas64.Vector{N: k.ndims, Data: tmp, Inc: 1},
		0.0, blas64.Vector{N: k.ndims, Data: Out, Inc: 1})

	return nil
}

func (k *KeDvrBasis) ExpDtToMatWrapper(Dt float64, In []float64, Out []float64) error {
	if len(In) != k.ndims || len(Out) != k.ndims {
		return fmt.Errorf("vector length mismatch: got In=%d Out=%d, expected %d", len(In), len(Out), k.ndims)
	}

	if err := k.EvalEigenVectors(); err != nil {
		return err
	}

	V := k.eigVecs
	lams := k.eigVals

	inVec := mat.NewVecDense(k.ndims, In)
	tmpVec := mat.NewVecDense(k.ndims, nil)
	outVec := mat.NewVecDense(k.ndims, Out)

	tmpVec.MulVec(V.T(), inVec)

	tmpData := tmpVec.RawVector().Data
	for i := 0; i < k.ndims; i++ {
		tmpData[i] *= math.Exp(Dt * lams[i])
	}

	outVec.MulVec(V, tmpVec)

	return nil
}

func (k *KeDvrBasis) ExpDt(dt float64, in []float64) ([]float64, error) {
	if len(in) != k.ndims {
		return nil, fmt.Errorf("input vector length %d doesn't match basis dimension %d", len(in), k.ndims)
	}
	out := make([]float64, k.ndims)
	if err := k.ExpDtTo(dt, in, out); err != nil {
		return nil, err
	}
	return out, nil
}

func (k *KeDvrBasis) ExpDtInPlace(dt float64, inOut []float64) error {
	if len(inOut) != k.ndims {
		return fmt.Errorf("vector length %d doesn't match basis dimension %d", len(inOut), k.ndims)
	}
	tmp := make([]float64, k.ndims)
	if err := k.ExpDtTo(dt, inOut, tmp); err != nil {
		return err
	}
	copy(inOut, tmp)
	return nil
}
