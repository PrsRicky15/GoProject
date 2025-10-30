package NumericalRecipies

import "gonum.org/v1/gonum/mat"

type dataType interface {
	float64 | complex128 | []float64 | []complex128 | *mat.Dense | *mat.CDense
}

type ChebyshevFirstKind[T dataType] struct{}

func (ch ChebyshevFirstKind[T]) Compute(x T, order uint32) T {
	if order == 0 {
		return ch.initPhi0(x)
	}
	if order == 1 {
		return x
	}

	switch v := any(x).(type) {
	case float64:
		return ch.computeFloat(v, order)
	case complex128:
		return ch.computeComplex(v, order)
	}

	return x
}

func (ch ChebyshevFirstKind[T]) initPhi0(x T) T {
	switch v := any(x).(type) {
	case float64:
		return any(1.0).(T)
	case complex128:
		return any(complex(1.0, 0.0)).(T)
	case []float64:
		result := make([]float64, len(v))
		for i := range result {
			result[i] = 1.0
		}
		return any(result).(T)
	case []complex128:
		result := make([]complex128, len(v))
		for i := range result {
			result[i] = complex(1.0, 0.0)
		}
		return any(result).(T)
	case *mat.Dense:
		rows, cols := v.Dims()
		result := mat.NewDense(rows, cols, nil)
		for i := 0; i < rows && i < cols; i++ {
			result.Set(i, i, 1.0)
		}
		return any(result).(T)
	case *mat.CDense:
		rows, cols := v.Dims()
		result := mat.NewCDense(rows, cols, nil)
		for i := 0; i < rows && i < cols; i++ {
			result.Set(i, i, complex(1.0, 0.0))
		}
		return any(result).(T)
	}
	return x
}

func (ch ChebyshevFirstKind[T]) computeFloat(x float64, order uint32) T {
	if order == 0 {
		return any(1.0).(T)
	}
	if order == 1 {
		return any(x).(T)
	}

	phi0 := 1.0
	phi1 := x
	var phi2 float64

	for i := uint32(2); i <= order; i++ {
		phi2 = 2.0*phi1*x - phi0
		phi0 = phi1
		phi1 = phi2
	}

	return any(phi2).(T)
}

func (ch ChebyshevFirstKind[T]) computeComplex(x complex128, order uint32) T {
	if order == 0 {
		return any(complex(1.0, 0.0)).(T)
	}
	if order == 1 {
		return any(x).(T)
	}

	phi0 := complex(1.0, 0.0)
	phi1 := x
	var phi2 complex128

	for i := uint32(2); i <= order; i++ {
		phi2 = 2*x*phi1 - phi0
		phi0 = phi1
		phi1 = phi2
	}

	return any(phi2).(T)
}

type ChebyshevSecondKind[T dataType] struct{}

type ChebyshevComplex[T dataType] struct{}
