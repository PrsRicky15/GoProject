package gridData

import "math/cmplx"

// NDimEvaluateAt defines the Cartesian Morse extension as a sum of independent
// one-coordinate Morse potentials. Cross-coordinate couplings can be added with
// NDimPotentialFunc. Existing scalar evaluation and force APIs are unchanged.
func (m Morse[T]) NDimEvaluateAt(x []T) T {
	var v T
	for _, xi := range x {
		v += m.EvaluateAt(xi)
	}
	return v
}

// NDimEvaluateAt sums the polynomial along each coordinate, including its constant.
func (p Polynomial[T]) NDimEvaluateAt(x []T) T {
	var v T
	for _, xi := range x {
		v += p.EvaluateAt(xi)
	}
	return v
}

// NDimEvaluateAt is an isotropic Gaussian. Cen is the center of EVERY coordinate;
// use a translated/scaled coordinate callback for an arbitrary center/ellipsoid.
func (g Gaussian[T]) NDimEvaluateAt(x []T) T {
	r2 := ndRadiusSquared(x, g.Cen) / complex(g.Sigma*g.Sigma, 0)
	return ndScalar[T](complex(g.Strength, 0) * cmplx.Exp(-r2/2))
}

// NDimEvaluateAt is Charge/sqrt(sum_d(x[d]-Centre)^2 + SoftParam^2).
// Complex coordinates use the analytic squared radius, not absolute squares.
func (s SoftCore[T]) NDimEvaluateAt(x []T) T {
	r2 := ndRadiusSquared(x, s.Centre) + complex(s.SoftParam*s.SoftParam, 0)
	return ndScalar[T](complex(s.Charge, 0) / cmplx.Sqrt(r2))
}

// NDimEvaluateAt is the radial super-Gaussian Strength*exp(-(r/Sigma)^Order).
// For even Order this agrees with the scalar definition in one dimension.
func (s SuperGaussian[T]) NDimEvaluateAt(x []T) T {
	r2 := ndRadiusSquared(x, s.Cen) / complex(s.Sigma*s.Sigma, 0)
	return ndScalar[T](complex(s.Strength, 0) * cmplx.Exp(-cmplx.Pow(r2, complex(float64(s.Order)/2, 0))))
}

// NDimEvaluateAt extends the original chain of Gaussians additively along axes.
// It is distinct from a pair of localized n-dimensional Gaussian wells.
func (m MultiGaussian[T]) NDimEvaluateAt(x []T) T {
	var v T
	for _, xi := range x {
		v += m.EvaluateAt(xi)
	}
	return v
}

func ndRadiusSquared[T VarType](x []T, center float64) complex128 {
	var r2 complex128
	for _, xi := range x {
		var z complex128
		switch v := any(xi).(type) {
		case float64:
			z = complex(v, 0)
		case complex128:
			z = v
		}
		z -= complex(center, 0)
		r2 += z * z
	}
	return r2
}
func ndScalar[T VarType](z complex128) T {
	var zero T
	if _, ok := any(zero).(float64); ok {
		return any(real(z)).(T)
	}
	return any(z).(T)
}
