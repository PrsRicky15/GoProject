package BasicNumerics

type Integrator1D interface {
	Integrate(f Func1D, bounds Bounds1D, nPoints int) float64
	Name() string
}

type Integrator2D interface {
	Integrate(f Func2D, bounds Bounds2D, nPoints int) float64
	Name() string
}

type Integrator3D interface {
	Integrate(f Func3D, bounds Bounds3D, nPoints int) float64
	Name() string
}

type Func1D func(x float64) float64
type Func2D func(x, y float64) float64
type Func3D func(x, y, z float64) float64

type Bounds1D struct {
	A, B float64 // [a, b]
}

type Bounds2D struct {
	X Bounds1D // [x_min, x_max]
	Y Bounds1D // [y_min, y_max]
}

type Bounds3D struct {
	X Bounds1D // [x_min, x_max]
	Y Bounds1D // [y_min, y_max]
	Z Bounds1D // [z_min, z_max]
}

type RectangleIntegrator1D struct{}

func (r *RectangleIntegrator1D) Name() string {
	return "Rectangle (1D)"
}

func (r *RectangleIntegrator1D) Integrate(f Func1D, bounds Bounds1D, nPoints int) float64 {
	dx := (bounds.B - bounds.A) / float64(nPoints)
	sum := 0.0
	x := bounds.A

	for i := 0; i < nPoints; i++ {
		sum += f(x)
		x += dx
	}

	return sum * dx
}

type TrapezoidalIntegrator1D struct {
	Bounds1D
}

func (t *TrapezoidalIntegrator1D) Name() string {
	return "Trapezoidal (1D)"
}

func (t *TrapezoidalIntegrator1D) SetBounds(bounds Bounds1D) {
	t.Bounds1D = bounds
}

func (t *TrapezoidalIntegrator1D) Integrate(f Func1D, nPoints int) float64 {
	dx := (t.Bounds1D.B - t.Bounds1D.A) / float64(nPoints)
	sum := 0.0

	for i := 1; i < nPoints; i++ {
		x := t.Bounds1D.A + float64(i)*dx
		sum += f(x)
	}

	return (sum + 0.5*(f(t.Bounds1D.A)+f(t.Bounds1D.B))) * dx
}

// SimpsonIntegrator1D implements Simpson's rule
type SimpsonIntegrator1D struct{}

func (s *SimpsonIntegrator1D) Name() string {
	return "Simpson (1D)"
}

func (s *SimpsonIntegrator1D) Integrate(f Func1D, bounds Bounds1D, nPoints int) float64 {
	dx := (bounds.B - bounds.A) / float64(nPoints)
	sumOdd := 0.0
	sumEven := 0.0

	for i := 1; i < nPoints; i += 2 {
		x := bounds.A + float64(i)*dx
		sumOdd += f(x)
	}

	for i := 2; i < nPoints; i += 2 {
		x := bounds.A + float64(i)*dx
		sumEven += f(x)
	}

	return dx / 3.0 * (f(bounds.A) + 4*sumOdd + 2*sumEven + f(bounds.B))
}

// ============================================================================
// 2D Integration Methods
// ============================================================================

// RectangleIntegrator2D implements the rectangle rule for 2D
type RectangleIntegrator2D struct{}

func (r *RectangleIntegrator2D) Name() string {
	return "Rectangle (2D)"
}

func (r *RectangleIntegrator2D) Integrate(f Func2D, bounds Bounds2D, nPoints int) float64 {
	dx := (bounds.X.B - bounds.X.A) / float64(nPoints)
	dy := (bounds.Y.B - bounds.Y.A) / float64(nPoints)
	sum := 0.0

	x := bounds.X.A
	for i := 0; i < nPoints; i++ {
		y := bounds.Y.A
		for j := 0; j < nPoints; j++ {
			sum += f(x, y)
			y += dy
		}
		x += dx
	}

	return sum * dx * dy
}

// TrapezoidalIntegrator2D implements the trapezoidal rule for 2D
type TrapezoidalIntegrator2D struct{}

func (t *TrapezoidalIntegrator2D) Name() string {
	return "Trapezoidal (2D)"
}

func (t *TrapezoidalIntegrator2D) Integrate(f Func2D, bounds Bounds2D, nPoints int) float64 {
	dx := (bounds.X.B - bounds.X.A) / float64(nPoints)
	dy := (bounds.Y.B - bounds.Y.A) / float64(nPoints)
	sum := 0.0

	// Interior points
	for i := 1; i < nPoints; i++ {
		x := bounds.X.A + float64(i)*dx
		for j := 1; j < nPoints; j++ {
			y := bounds.Y.A + float64(j)*dy
			sum += f(x, y)
		}
	}

	// Edge points (weight 0.5)
	edgeSum := 0.0
	for i := 1; i < nPoints; i++ {
		x := bounds.X.A + float64(i)*dx
		edgeSum += f(x, bounds.Y.A) + f(x, bounds.Y.B)
	}
	for j := 1; j < nPoints; j++ {
		y := bounds.Y.A + float64(j)*dy
		edgeSum += f(bounds.X.A, y) + f(bounds.X.B, y)
	}

	// Corner points (weight 0.25)
	cornerSum := f(bounds.X.A, bounds.Y.A) + f(bounds.X.A, bounds.Y.B) +
		f(bounds.X.B, bounds.Y.A) + f(bounds.X.B, bounds.Y.B)

	return (sum + 0.5*edgeSum + 0.25*cornerSum) * dx * dy
}

// SimpsonIntegrator2D implements Simpson's rule for 2D
type SimpsonIntegrator2D struct{}

func (s *SimpsonIntegrator2D) Name() string {
	return "Simpson (2D)"
}

func (s *SimpsonIntegrator2D) Integrate(f Func2D, bounds Bounds2D, nPoints int) float64 {
	dx := (bounds.X.B - bounds.X.A) / float64(nPoints)
	dy := (bounds.Y.B - bounds.Y.A) / float64(nPoints)

	sum := 0.0

	for i := 0; i <= nPoints; i++ {
		x := bounds.X.A + float64(i)*dx
		weightX := s.getWeight(i, nPoints)

		for j := 0; j <= nPoints; j++ {
			y := bounds.Y.A + float64(j)*dy
			weightY := s.getWeight(j, nPoints)

			sum += weightX * weightY * f(x, y)
		}
	}

	return sum * dx * dy / 9.0
}

func (s *SimpsonIntegrator2D) getWeight(index, nPoints int) float64 {
	if index == 0 || index == nPoints {
		return 1.0
	} else if index%2 == 1 {
		return 4.0
	}
	return 2.0
}

// ============================================================================
// 3D Integration Methods
// ============================================================================

// RectangleIntegrator3D implements the rectangle rule for 3D
type RectangleIntegrator3D struct{}

func (r *RectangleIntegrator3D) Name() string {
	return "Rectangle (3D)"
}

func (r *RectangleIntegrator3D) Integrate(f Func3D, bounds Bounds3D, nPoints int) float64 {
	dx := (bounds.X.B - bounds.X.A) / float64(nPoints)
	dy := (bounds.Y.B - bounds.Y.A) / float64(nPoints)
	dz := (bounds.Z.B - bounds.Z.A) / float64(nPoints)
	sum := 0.0

	x := bounds.X.A
	for i := 0; i < nPoints; i++ {
		y := bounds.Y.A
		for j := 0; j < nPoints; j++ {
			z := bounds.Z.A
			for k := 0; k < nPoints; k++ {
				sum += f(x, y, z)
				z += dz
			}
			y += dy
		}
		x += dx
	}

	return sum * dx * dy * dz
}

// TrapezoidalIntegrator3D implements the trapezoidal rule for 3D
type TrapezoidalIntegrator3D struct{}

func (t *TrapezoidalIntegrator3D) Name() string {
	return "Trapezoidal (3D)"
}

func (t *TrapezoidalIntegrator3D) Integrate(f Func3D, bounds Bounds3D, nPoints int) float64 {
	dx := (bounds.X.B - bounds.X.A) / float64(nPoints)
	dy := (bounds.Y.B - bounds.Y.A) / float64(nPoints)
	dz := (bounds.Z.B - bounds.Z.A) / float64(nPoints)

	sum := 0.0

	for i := 0; i <= nPoints; i++ {
		x := bounds.X.A + float64(i)*dx
		weightX := t.getWeight(i, nPoints)

		for j := 0; j <= nPoints; j++ {
			y := bounds.Y.A + float64(j)*dy
			weightY := t.getWeight(j, nPoints)

			for k := 0; k <= nPoints; k++ {
				z := bounds.Z.A + float64(k)*dz
				weightZ := t.getWeight(k, nPoints)

				sum += weightX * weightY * weightZ * f(x, y, z)
			}
		}
	}

	return sum * dx * dy * dz / 8.0
}

func (t *TrapezoidalIntegrator3D) getWeight(index, nPoints int) float64 {
	if index == 0 || index == nPoints {
		return 1.0
	}
	return 2.0
}

// SimpsonIntegrator3D implements Simpson's rule for 3D
type SimpsonIntegrator3D struct{}

func (s *SimpsonIntegrator3D) Name() string {
	return "Simpson (3D)"
}

func (s *SimpsonIntegrator3D) Integrate(f Func3D, bounds Bounds3D, nPoints int) float64 {
	dx := (bounds.X.B - bounds.X.A) / float64(nPoints)
	dy := (bounds.Y.B - bounds.Y.A) / float64(nPoints)
	dz := (bounds.Z.B - bounds.Z.A) / float64(nPoints)

	sum := 0.0

	for i := 0; i <= nPoints; i++ {
		x := bounds.X.A + float64(i)*dx
		weightX := s.getWeight(i, nPoints)

		for j := 0; j <= nPoints; j++ {
			y := bounds.Y.A + float64(j)*dy
			weightY := s.getWeight(j, nPoints)

			for k := 0; k <= nPoints; k++ {
				z := bounds.Z.A + float64(k)*dz
				weightZ := s.getWeight(k, nPoints)

				sum += weightX * weightY * weightZ * f(x, y, z)
			}
		}
	}

	return sum * dx * dy * dz / 27.0
}

func (s *SimpsonIntegrator3D) getWeight(index, nPoints int) float64 {
	if index == 0 || index == nPoints {
		return 1.0
	} else if index%2 == 1 {
		return 4.0
	}
	return 2.0
}
