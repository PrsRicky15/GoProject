package OperatorAlgebra

import (
	"GoProject/gridData"
	"fmt"
	"math"
	"math/cmplx"
	"strings"
	"testing"

	"gonum.org/v1/gonum/mat"
)

const tolerance = 1e-10

func TestKeDVR_Evaluate(t *testing.T) {
	rgrid, err := gridData.NewFromLength(10., 30)
	if err != nil {
		panic(err)
	}
	kinE := NewKeDVR(rgrid, 1.)
	kinE.GetMat()
}

// Helper function to check if two floats are approximately equal
func floatEqual(a, b, tol float64) bool {
	return math.Abs(a-b) < tol
}

// Helper function to check if two complex numbers are approximately equal
func complexEqual(a, b complex128, tol float64) bool {
	return cmplx.Abs(a-b) < tol
}

func TestNewKeDVR(t *testing.T) {
	tests := []struct {
		name    string
		rMin    float64
		rMax    float64
		nPoints int
		mass    float64
	}{
		{"small_symmetric_grid", -5.0, 5.0, 10, 1.0},
		{"zero_to_infinity_grid", 0.0, 10.0, 20, 1.0},
		{"heavy_mass", -10.0, 10.0, 15, 10.0},
		{"light_mass", 0.0, 5.0, 8, 0.5},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			grid, _ := gridData.NewRGrid(tt.rMin, tt.rMax, uint32(tt.nPoints))
			ke := NewKeDVR(grid, tt.mass)

			if ke == nil {
				t.Fatal("NewKeDVR returned nil")
			}

			if ke.Ndims() != tt.nPoints {
				t.Errorf("Ndims() = %d, want %d", ke.Ndims(), tt.nPoints)
			}

			if ke.Mass() != tt.mass {
				t.Errorf("Mass() = %f, want %f", ke.Mass(), tt.mass)
			}
		})
	}
}

func TestIsZeroToInfinity(t *testing.T) {
	tests := []struct {
		name     string
		rMin     float64
		rMax     float64
		expected bool
	}{
		{"symmetric_grid", -5.0, 5.0, false},
		{"zero_to_positive", 0.0, 10.0, true},
		{"small_asymmetry", -0.4, 0.5, false},
		{"large_positive", 1.0, 10.0, true},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			grid, _ := gridData.NewRGrid(tt.rMin, tt.rMax, 10)
			ke := NewKeDVR(grid, 1.0)

			result := ke.isZeroToInfinity()
			if result != tt.expected {
				t.Errorf("isZeroToInfinity() = %v, want %v (rMin=%f, rMax=%f, sum=%f)",
					result, tt.expected, tt.rMin, tt.rMax, tt.rMin+tt.rMax)
			}
		})
	}
}

func TestGetMatSymmetry(t *testing.T) {
	tests := []struct {
		name    string
		rMin    float64
		rMax    float64
		nPoints int
	}{
		{"symmetric_domain", -5.0, 5.0, 10},
		{"zero_to_infinity", 0.0, 10.0, 15},
		{"large_grid", -20.0, 20.0, 50},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			grid, _ := gridData.NewRGrid(tt.rMin, tt.rMax, uint32(tt.nPoints))
			ke := NewKeDVR(grid, 1.0)
			kMat := ke.GetMat()

			if !isSymmetric(kMat) {
				t.Error("Kinetic energy matrix is not symmetric")
			}
		})
	}
}

func TestGetMatCaching(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 10)
	ke := NewKeDVR(grid, 1.0)
	mat1 := ke.GetMat()
	mat2 := ke.GetMat()

	if mat1 != mat2 {
		t.Error("GetMat() should return cached matrix on second call")
	}
}

func TestClear(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 10)
	ke := NewKeDVR(grid, 1.0)
	ke.GetMat()
	err := ke.EvalEigenVectors()
	if err != nil {
		return
	}
	ke.Clear()

	mat1 := ke.GetMat()
	if mat1 == nil {
		t.Error("GetMat() should work after Clear()")
	}
	return
}

func TestMinInftyToInfty(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 5)
	ke := NewKeDVR(grid, 1.0)

	m := mat.NewDense(5, 5, nil)
	ke.MinInftyToInfty(m)

	for i := 0; i < 5; i++ {
		if m.At(i, i) <= 0 {
			t.Errorf("Diagonal element [%d,%d] = %f should be positive", i, i, m.At(i, i))
		}
	}

	if !isSymmetric(m) {
		t.Error("MinInftyToInfty matrix is not symmetric")
	}

	for i := 1; i < 5; i++ {
		diff := i - 0
		expectedSign := float64(1 - 2*(diff&1))
		actualSign := math.Copysign(1.0, m.At(i, 0))
		if expectedSign != actualSign {
			t.Errorf("Off-diagonal sign at [%d,0] incorrect: expected %f, got %f",
				i, expectedSign, actualSign)
		}
	}
}

func TestZeroToInfinity(t *testing.T) {
	grid, _ := gridData.NewRGrid(0.0, 10.0, 5)
	ke := NewKeDVR(grid, 1.0)

	m := mat.NewDense(5, 5, nil)
	ke.ZeroToInfinity(m)

	if !isSymmetric(m) {
		t.Error("ZeroToInfinity matrix is not symmetric")
	}

	for i := 0; i < 5; i++ {
		if m.At(i, i) <= 0 {
			t.Errorf("Diagonal element [%d,%d] = %f should be positive", i, i, m.At(i, i))
		}
	}
}

func TestRealDiagonalize(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 10)
	ke := NewKeDVR(grid, 1.0)

	eigenvalues, eigenvectors, err := ke.GetEigenValuesVectors()
	if err != nil {
		t.Fatalf("Get_Eigen_Values_Vectors failed: %v", err)
	}

	if len(eigenvalues) != 10 {
		t.Errorf("Expected 10 eigenvalues, got %d", len(eigenvalues))
	}

	r, c := eigenvectors.Dims()
	if r != 10 || c != 10 {
		t.Errorf("Eigenvector matrix dimensions = (%d,%d), want (10,10)", r, c)
	}

	// Eigenvalues should be non-negative for kinetic energy
	for i, ev := range eigenvalues {
		if ev < -tolerance {
			t.Errorf("Eigenvalue[%d] = %f is negative", i, ev)
		}
	}

	// Verify V^T * K * V ≈ diagonal(eigenvalues)
	K := ke.GetMat()
	V := eigenvectors

	// Compute K * V
	KV := mat.NewDense(10, 10, nil)
	KV.Mul(K, V)

	// Compute V^T * K * V
	VtKV := mat.NewDense(10, 10, nil)
	VtKV.Mul(V.T(), KV)

	// Check that it's approximately diagonal with eigenvalues
	for i := 0; i < 10; i++ {
		for j := 0; j < 10; j++ {
			expected := 0.0
			if i == j {
				expected = eigenvalues[i]
			}
			if !floatEqual(VtKV.At(i, j), expected, 1e-8) {
				t.Errorf("V^T*K*V[%d,%d] = %f, expected %f", i, j, VtKV.At(i, j), expected)
			}
		}
	}
}

func TestExpDtTo(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 5)
	ke := NewKeDVR(grid, 1.0)

	in := []float64{1.0, 0.0, 0.0, 0.0, 0.0}
	out := make([]float64, 5)

	err := ke.ExpDtTo(0.0, in, out)
	if err != nil {
		t.Fatalf("ExpDtTo failed: %v", err)
	}

	// With dt=0, exp(0*K) = I, so out should equal in
	for i := range in {
		if !floatEqual(out[i], in[i], tolerance) {
			t.Errorf("ExpDtTo with dt=0: out[%d] = %f, want %f", i, out[i], in[i])
		}
	}
}

func TestExpDtToVectorLengthMismatch(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 5)
	ke := NewKeDVR(grid, 1.0)

	in := []float64{1.0, 0.0, 0.0} // Wrong size
	out := make([]float64, 5)

	err := ke.ExpDtTo(0.1, in, out)
	if err == nil {
		t.Error("ExpDtTo should return error for mismatched vector length")
	}
}

func TestExpIdtTo(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 5)
	ke := NewKeDVR(grid, 1.0)

	in := []complex128{1.0, 0.0, 0.0, 0.0, 0.0}
	out := make([]complex128, 5)

	err := ke.ExpIdtTo(0.0, in, out)
	if err != nil {
		t.Fatalf("ExpIdtTo failed: %v", err)
	}

	// With dt=0, exp(i*0*K) = I, so out should equal in
	for i := range in {
		if !complexEqual(out[i], in[i], tolerance) {
			t.Errorf("ExpIdtTo with dt=0: out[%d] = %v, want %v", i, out[i], in[i])
		}
	}
}

func TestExpIdtUnitarity(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 8)
	ke := NewKeDVR(grid, 1.0)

	// Create a normalized input vector
	in := make([]complex128, 8)
	norm := 0.0
	for i := range in {
		in[i] = complex(float64(i+1), float64(i))
		norm += real(in[i])*real(in[i]) + imag(in[i])*imag(in[i])
	}
	norm = math.Sqrt(norm)
	for i := range in {
		in[i] /= complex(norm, 0)
	}

	out := make([]complex128, 8)
	err := ke.ExpIdtTo(0.1, in, out)
	if err != nil {
		t.Fatalf("ExpIdtTo failed: %v", err)
	}

	// Check that output is also normalized (unitary evolution preserves norm)
	outNorm := 0.0
	for _, v := range out {
		outNorm += real(v)*real(v) + imag(v)*imag(v)
	}
	outNorm = math.Sqrt(outNorm)

	if !floatEqual(outNorm, 1.0, 1e-8) {
		t.Errorf("Unitary evolution did not preserve norm: got %f, want 1.0", outNorm)
	}
}

func TestExpDt(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 5)
	ke := NewKeDVR(grid, 1.0)

	in := []float64{1.0, 0.5, 0.0, -0.5, -1.0}
	out, err := ke.ExpDt(0.01, in)
	if err != nil {
		t.Fatalf("ExpDt failed: %v", err)
	}

	if len(out) != 5 {
		t.Errorf("ExpDt output length = %d, want 5", len(out))
	}
}

func TestExpDtInPlace(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 5)
	ke := NewKeDVR(grid, 1.0)

	in := []float64{1.0, 0.5, 0.0, -0.5, -1.0}
	expected, _ := ke.ExpDt(0.01, in)

	inOut := []float64{1.0, 0.5, 0.0, -0.5, -1.0}
	err := ke.ExpDtInPlace(0.01, inOut)
	if err != nil {
		t.Fatalf("ExpDtInPlace failed: %v", err)
	}

	for i := range expected {
		if !floatEqual(inOut[i], expected[i], tolerance) {
			t.Errorf("ExpDtInPlace[%d] = %f, want %f", i, inOut[i], expected[i])
		}
	}
}

func TestExpIdt(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 5)
	ke := NewKeDVR(grid, 1.0)

	in := []complex128{1.0, 0.5i, 0.0, -0.5i, -1.0}
	out, err := ke.ExpIdt(0.01, in)
	if err != nil {
		t.Fatalf("ExpIdt failed: %v", err)
	}

	if len(out) != 5 {
		t.Errorf("ExpIdt output length = %d, want 5", len(out))
	}
}

func TestExpIdtInPlace(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 5)
	ke := NewKeDVR(grid, 1.0)

	in := []complex128{1.0, 0.5i, 0.0, -0.5i, -1.0}
	expected, _ := ke.ExpIdt(0.01, in)

	inOut := []complex128{1.0, 0.5i, 0.0, -0.5i, -1.0}
	err := ke.ExpIdtInPlace(0.01, inOut)
	if err != nil {
		t.Fatalf("ExpIdtInPlace failed: %v", err)
	}

	// Results should match
	for i := range expected {
		if !complexEqual(inOut[i], expected[i], tolerance) {
			t.Errorf("ExpIdtInPlace[%d] = %v, want %v", i, inOut[i], expected[i])
		}
	}
}

func TestClone(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 10)
	ke := NewKeDVR(grid, 2.0)

	ke.GetMat()
	err := ke.EvalEigenVectors()
	if err != nil {
		return
	}

	keClone := ke.Clone()

	if keClone.Ndims() != ke.Ndims() {
		t.Errorf("Clone Ndims = %d, want %d", keClone.Ndims(), ke.Ndims())
	}

	if keClone.Mass() != ke.Mass() {
		t.Errorf("Clone Mass = %f, want %f", keClone.Mass(), ke.Mass())
	}

	ke.Clear()
	if keClone.eigVals == nil {
		t.Error("Clone should have independent eigen cache")
	}
}

func TestCloneWithoutCache(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 10)
	ke := NewKeDVR(grid, 2.0)
	keClone := ke.Clone()

	_, err, _ := keClone.GetEigenValuesVectors()
	if err != nil {
		t.Fatalf("Clone Get_Eigen_Values_Vectors failed: %v", err)
	}
}

func TestEigenvaluesOrdering(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 15)
	ke := NewKeDVR(grid, 1.0)

	eigenvalues, _, err := ke.GetEigenValuesVectors()
	if err != nil {
		t.Fatalf("Get_Eigen_Values_Vectors failed: %v", err)
	}

	// Check eigenvalues are in ascending order (gonum returns them sorted)
	for i := 1; i < len(eigenvalues); i++ {
		if eigenvalues[i] < eigenvalues[i-1]-tolerance {
			t.Errorf("Eigenvalues not sorted: eigenvalues[%d]=%f < eigenvalues[%d]=%f",
				i, eigenvalues[i], i-1, eigenvalues[i-1])
		}
	}
}

func TestExpDtConsistency(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 8)
	ke := NewKeDVR(grid, 1.0)

	in := make([]float64, 8)
	for i := range in {
		in[i] = math.Sin(float64(i) * 0.5)
	}

	// Apply exp(dt*K) twice with dt=0.1
	out1, _ := ke.ExpDt(0.1, in)
	out2, _ := ke.ExpDt(0.1, out1)

	// Apply exp(dt*K) once with dt=0.2
	out3, _ := ke.ExpDt(0.2, in)

	// Results should match: exp(0.2*K) = exp(0.1*K) * exp(0.1*K)
	for i := range out2 {
		if !floatEqual(out2[i], out3[i], 1e-8) {
			t.Errorf("Consistency check failed at [%d]: exp(0.1)^2 = %f, exp(0.2) = %f",
				i, out2[i], out3[i])
		}
	}
}

func TestLargeDimension(t *testing.T) {
	// Test with a larger grid to ensure no issues with larger matrices
	grid, _ := gridData.NewRGrid(-10.0, 10.0, 100)
	ke := NewKeDVR(grid, 1.0)

	if ke.Ndims() != 100 {
		t.Errorf("Ndims = %d, want 100", ke.Ndims())
	}

	// Should be able to compute matrix and eigenvalues
	eigenvalues, _, err := ke.GetEigenValuesVectors()
	if err != nil {
		t.Fatalf("Get_Eigen_Values_Vectors failed for large matrix: %v", err)
	}

	if len(eigenvalues) != 100 {
		t.Errorf("Expected 100 eigenvalues, got %d", len(eigenvalues))
	}
}

func TestDifferentMasses(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 10)

	// Heavier mass should give smaller kinetic energy eigenvalues
	keLightMass := NewKeDVR(grid, 1.0)
	keHeavyMass := NewKeDVR(grid, 10.0)

	evLight, _, _ := keLightMass.GetEigenValuesVectors()
	evHeavy, _, _ := keHeavyMass.GetEigenValuesVectors()

	// Each eigenvalue for heavy mass should be smaller (by factor of a mass ratio)
	for i := range evLight {
		ratio := evLight[i] / evHeavy[i]
		expectedRatio := 10.0
		if !floatEqual(ratio, expectedRatio, 0.01) {
			t.Errorf("Mass scaling incorrect at eigenvalue %d: ratio=%f, expected=%f",
				i, ratio, expectedRatio)
		}
	}
}

func BenchmarkGetMat(b *testing.B) {
	grid, _ := gridData.NewRGrid(-10.0, 10.0, 50)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ke := NewKeDVR(grid, 1.0)
		ke.GetMat()
	}
}

func BenchmarkRealDiagonalize(b *testing.B) {
	grid, _ := gridData.NewRGrid(-10.0, 10.0, 50)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ke := NewKeDVR(grid, 1.0)
		_, _, err := ke.GetEigenValuesVectors()
		if err != nil {
			return
		}
	}
}

func BenchmarkExpDt(b *testing.B) {
	grid, _ := gridData.NewRGrid(-10.0, 10.0, 50)
	ke := NewKeDVR(grid, 1.0)

	in := make([]float64, 50)
	for i := range in {
		in[i] = math.Sin(float64(i) * 0.1)
	}

	err := ke.EvalEigenVectors()
	if err != nil {
		return
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, err := ke.ExpDt(0.01, in)
		if err != nil {
			return
		}
	}
}

func BenchmarkExpIdt(b *testing.B) {
	grid, _ := gridData.NewRGrid(-10.0, 10.0, 50)
	ke := NewKeDVR(grid, 1.0)

	in := make([]complex128, 50)
	for i := range in {
		in[i] = complex(math.Sin(float64(i)*0.1), math.Cos(float64(i)*0.1))
	}

	err := ke.EvalEigenVectors()
	if err != nil {
		return
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, err2 := ke.ExpIdt(0.01, in)
		if err2 != nil {
			return
		}
	}
}

// TestNewFiniteDiff tests the constructor
func TestNewFiniteDiff(t *testing.T) {
	tests := []struct {
		name      string
		mass      float64
		order     int
		wantError bool
	}{
		{"Valid order 3", 1.0, 3, false},
		{"Valid order 5", 1.0, 5, false},
		{"Valid order 6", 1.0, 6, false},
		{"Invalid order 2", 1.0, 2, true},
		{"Invalid order 4", 1.0, 4, true},
		{"Invalid order 7", 1.0, 7, true},
	}

	grid, err := gridData.NewFromLength(20, 100)
	if err != nil {
		t.Errorf("NewFromLength failed with error: %v", err)
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			fd, err := NewFiniteDiff(grid, tt.mass, tt.order)

			if tt.wantError {
				if err == nil {
					t.Errorf("NewFiniteDiff() expected error but got none")
				}
				if fd != nil {
					t.Errorf("NewFiniteDiff() expected nil FiniteDiff on error")
				}
			} else {
				if err != nil {
					t.Errorf("NewFiniteDiff() unexpected error: %v", err)
				}
				if fd == nil {
					t.Errorf("NewFiniteDiff() expected non-nil FiniteDiff")
				}
			}
		})
	}
}

// TestKineticCoefficient verifies the kinetic coefficient calculation
func TestKineticCoefficient(t *testing.T) {
	grid, err := gridData.NewFromLength(20, 100)
	if err != nil {
		t.Errorf("NewFromLength failed with error: %v", err)
	}
	mass := 2.0

	fd, err := NewFiniteDiff(grid, mass, 3)
	if err != nil {
		t.Fatalf("NewFiniteDiff() error: %v", err)
	}

	dxSq := grid.DeltaR() * grid.DeltaR()
	expected := -1.0 / (2.0 * mass * dxSq)

	if math.Abs(fd.kineticCoeff-expected) > 1e-10 {
		t.Errorf("kineticCoeff = %v, want %v", fd.kineticCoeff, expected)
	}
}

// TestCompute3rdOrderCoefficients tests 3-point stencil
func TestCompute3rdOrderCoefficients(t *testing.T) {
	grid, err := gridData.NewFromLength(20, 100)
	if err != nil {
		t.Errorf("NewFromLength failed with error: %v", err)
	}
	fd, _ := NewFiniteDiff(grid, 1.0, 3)

	keMat := mat.NewDense(5, 5, nil)
	fd.compute3rdOrderCoefficients(keMat)

	expectedDiag := -2.0 * fd.kineticCoeff
	for i := 0; i < 5; i++ {
		if math.Abs(keMat.At(i, i)-expectedDiag) > 1e-10 {
			t.Errorf("Diagonal[%d] = %v, want %v", i, keMat.At(i, i), expectedDiag)
		}
	}

	for i := 0; i < 4; i++ {
		if math.Abs(keMat.At(i, i+1)-fd.kineticCoeff) > 1e-10 {
			t.Errorf("Upper diagonal[%d] = %v, want %v", i, keMat.At(i, i+1), fd.kineticCoeff)
		}
		if math.Abs(keMat.At(i+1, i)-fd.kineticCoeff) > 1e-10 {
			t.Errorf("Lower diagonal[%d] = %v, want %v", i, keMat.At(i+1, i), fd.kineticCoeff)
		}
	}

	if !isSymmetric(keMat) {
		t.Error("Matrix is not symmetric")
	}
}

// TestCompute5thOrderCoefficients tests 5-point stencil
func TestCompute5thOrderCoefficients(t *testing.T) {
	grid, err := gridData.NewFromLength(20, 100)
	if err != nil {
		t.Errorf("NewFromLength failed with error: %v", err)
	}
	fd, _ := NewFiniteDiff(grid, 1.0, 5)

	keMat := mat.NewDense(6, 6, nil)
	fd.compute5thOrderCoefficients(keMat)

	k := fd.kineticCoeff

	for i := 0; i < 6; i++ {
		expected := 2.5 * k
		if math.Abs(keMat.At(i, i)-expected) > 1e-10 {
			t.Errorf("Diagonal[%d] = %v, want %v", i, keMat.At(i, i), expected)
		}
	}

	for i := 0; i < 5; i++ {
		expected := -4.0 / 3.0 * k
		if math.Abs(keMat.At(i, i+1)-expected) > 1e-10 {
			t.Errorf("First off-diag[%d] = %v, want %v", i, keMat.At(i, i+1), expected)
		}
	}

	for i := 0; i < 4; i++ {
		expected := 1.0 / 12.0 * k
		if math.Abs(keMat.At(i, i+2)-expected) > 1e-10 {
			t.Errorf("Second off-diag[%d] = %v, want %v", i, keMat.At(i, i+2), expected)
		}
	}

	if !isSymmetric(keMat) {
		t.Error("Matrix is not symmetric")
	}
}

// TestCompute7thOrderCoefficients tests 7-point stencil
func TestCompute7thOrderCoefficients(t *testing.T) {
	grid, err := gridData.NewFromLength(20, 100)
	if err != nil {
		t.Errorf("NewFromLength failed with error: %v", err)
	}
	fd, _ := NewFiniteDiff(grid, 1.0, 6) // Note: using order 6 to call 7th order

	keMat := mat.NewDense(8, 8, nil)
	fd.compute7thOrderCoefficients(keMat)

	k := fd.kineticCoeff
	expectedDiag := 49.0 / 18.0 * k
	if math.Abs(keMat.At(0, 0)-expectedDiag) > 1e-10 {
		t.Errorf("Diagonal = %v, want %v", keMat.At(0, 0), expectedDiag)
	}

	expectedOff1 := -1.5 * k
	if math.Abs(keMat.At(0, 1)-expectedOff1) > 1e-10 {
		t.Errorf("First off-diag = %v, want %v", keMat.At(0, 1), expectedOff1)
	}

	expectedOff2 := 3.0 / 20.0 * k
	if math.Abs(keMat.At(0, 2)-expectedOff2) > 1e-10 {
		t.Errorf("Second off-diag = %v, want %v", keMat.At(0, 2), expectedOff2)
	}

	expectedOff3 := -1.0 / 90.0 * k
	if math.Abs(keMat.At(0, 3)-expectedOff3) > 1e-10 {
		t.Errorf("Third off-diag = %v, want %v", keMat.At(0, 3), expectedOff3)
	}

	if !isSymmetric(keMat) {
		t.Error("Matrix is not symmetric")
	}
}

// TestCompute9thOrderCoefficients tests 9-point stencil
func TestCompute9thOrderCoefficients(t *testing.T) {
	grid, err := gridData.NewFromLength(20, 100)
	if err != nil {
		t.Errorf("NewFromLength failed with error: %v", err)
	}
	fd, _ := NewFiniteDiff(grid, 1.0, 5)

	keMat := mat.NewDense(10, 10, nil)
	fd.compute9thOrderCoefficients(keMat)

	k := fd.kineticCoeff
	coeffs := []struct {
		offset int
		value  float64
	}{
		{0, 205.0 / 72.0},
		{1, -8.0 / 5.0},
		{2, 1.0 / 5.0},
		{3, -8.0 / 315.0},
		{4, 1.0 / 560.0},
	}

	for _, c := range coeffs {
		expected := c.value * k
		for i := 0; i < 10-c.offset; i++ {
			if math.Abs(keMat.At(i, i+c.offset)-expected) > 1e-10 {
				t.Errorf("Offset %d at [%d,%d] = %v, want %v",
					c.offset, i, i+c.offset, keMat.At(i, i+c.offset), expected)
			}
		}
	}

	if !isSymmetric(keMat) {
		t.Error("Matrix is not symmetric")
	}
}

// TestComputeMatrix tests the dispatcher method
func TestComputeMatrix(t *testing.T) {
	grid, err := gridData.NewFromLength(20, 100)
	if err != nil {
		t.Errorf("NewFromLength failed with error: %v", err)
	}

	orders := []int{3, 5}
	for _, order := range orders {
		t.Run("Order"+string(rune(order+'0')), func(t *testing.T) {
			fd, err := NewFiniteDiff(grid, 1.0, order)
			if err != nil {
				t.Fatalf("NewFiniteDiff() error: %v", err)
			}

			keMat := mat.NewDense(5, 5, nil)
			fd.ComputeMatrix(keMat)

			hasNonZero := false
			for i := 0; i < 5; i++ {
				for j := 0; j < 5; j++ {
					if math.Abs(keMat.At(i, j)) > 1e-10 {
						hasNonZero = true
						break
					}
				}
			}

			if !hasNonZero {
				t.Error("Matrix is all zeros")
			}
		})
	}
}

// Helper function to check matrix symmetry
func isSymmetric(m *mat.Dense) bool {
	r, c := m.Dims()
	if r != c {
		return false
	}

	for i := 0; i < r; i++ {
		for j := i + 1; j < c; j++ {
			if math.Abs(m.At(i, j)-m.At(j, i)) > 1e-10 {
				return false
			}
		}
	}
	return true
}

// Benchmark tests
func BenchmarkCompute3rdOrder(b *testing.B) {
	grid, err := gridData.NewFromLength(20, 100)
	if err != nil {
		b.Errorf("NewFromLength failed with error: %v", err)
	}

	fd, _ := NewFiniteDiff(grid, 1.0, 3)
	keMat := mat.NewDense(100, 100, nil)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		fd.compute3rdOrderCoefficients(keMat)
	}
}

func BenchmarkCompute9thOrder(b *testing.B) {
	grid, err := gridData.NewFromLength(20, 100)
	if err != nil {
		b.Errorf("NewFromLength failed with error: %v", err)
	}

	fd, _ := NewFiniteDiff(grid, 1.0, 5)
	keMat := mat.NewDense(100, 100, nil)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		fd.compute9thOrderCoefficients(keMat)
	}
}

// Fourier Basis test
// Helper function to check if two complex slices are approximately equal
func complexSlicesEqual(a, b []complex128, tolerance float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if cmplx.Abs(a[i]-b[i]) > tolerance {
			return false
		}
	}
	return true
}

// Helper function to compute L2 norm
func l2Norm(a []complex128) float64 {
	sum := 0.0
	for _, v := range a {
		sum += real(v)*real(v) + imag(v)*imag(v)
	}
	return math.Sqrt(sum)
}

// Helper function to compute relative error
func relativeError(expected, actual []complex128) float64 {
	diff := make([]complex128, len(expected))
	for i := range expected {
		diff[i] = expected[i] - actual[i]
	}
	normDiff := l2Norm(diff)
	normExpected := l2Norm(expected)
	if normExpected == 0 {
		return normDiff
	}
	return normDiff / normExpected
}

// ============================================================================
// Test Initialization
// ============================================================================

func TestFFTInit(t *testing.T) {
	grid, _ := gridData.NewFromLength(10, 100)
	mass := 1.0

	fb := FFTInit(grid, mass)
	if fb == nil {
		t.Fatal("FFTInit returned nil")
	}
	defer fb.Clean()

	if fb.nPoints != 64 {
		t.Errorf("Expected nPoints=64, got %d", fb.nPoints)
	}

	if fb.mass != mass {
		t.Errorf("Expected mass=%f, got %f", mass, fb.mass)
	}

	if len(fb.kValues) != 64 {
		t.Errorf("Expected kValues length=64, got %d", len(fb.kValues))
	}

	if len(fb.keValues) != 64 {
		t.Errorf("Expected keValues length=64, got %d", len(fb.keValues))
	}

	if fb.Buff == nil {
		t.Error("Buff should not be nil")
	}
}

func TestFFTInitWithDifferentSizes(t *testing.T) {
	testSizes := []int{8, 16, 32, 64, 128, 256}

	for _, size := range testSizes {
		t.Run(fmt.Sprintf("size_%d", size), func(t *testing.T) {
			grid, _ := gridData.NewFromLength(6.4, uint32(size))
			fb := FFTInit(grid, 1.0)
			if fb == nil {
				t.Fatal("FFTInit returned nil")
			}
			defer fb.Clean()

			if fb.nPoints != size {
				t.Errorf("Expected nPoints=%d, got %d", size, fb.nPoints)
			}
		})
	}
}

func TestRedefine(t *testing.T) {
	grid1, _ := gridData.NewFromLength(6.4, 64)
	grid2, _ := gridData.NewFromLength(6.4, 128)

	fb := FFTInit(grid1, 1.0)
	defer fb.Clean()

	if fb.nPoints != 64 {
		t.Errorf("Initial nPoints should be 64, got %d", fb.nPoints)
	}

	fb.Redefine(grid2, 2.0)

	if fb.nPoints != 128 {
		t.Errorf("After Redefine nPoints should be 128, got %d", fb.nPoints)
	}

	if fb.mass != 2.0 {
		t.Errorf("After Redefine mass should be 2.0, got %f", fb.mass)
	}
}

// ============================================================================
// Test Forward and Backward Transforms
// ============================================================================

func TestForwardBackwardRoundTrip(t *testing.T) {
	grid, _ := gridData.NewFromLength(6.4, 64)
	fb := FFTInit(grid, 1.0)
	defer fb.Clean()

	input := make([]complex128, 64)
	for i := range input {
		input[i] = complex(float64(i), float64(-i))
	}

	original := make([]complex128, 64)
	copy(original, input)

	output := make([]complex128, 64)
	temp := make([]complex128, 64)

	fb.forward(input, temp)
	fb.backward(temp, output)

	nPoints := float64(fb.nPoints)
	for i := range output {
		output[i] /= complex(nPoints, 0)
	}

	if !complexSlicesEqual(original, output, 1e-10) {
		t.Error("Forward-Backward round trip failed")
		t.Logf("Expected first element: %v", original[0])
		t.Logf("Got first element: %v", output[0])
	}
}

func TestForwardInPlace(t *testing.T) {
	grid, _ := gridData.NewFromLength(6.4, 64)
	fb := FFTInit(grid, 1.0)
	defer fb.Clean()

	input1 := make([]complex128, 64)
	input2 := make([]complex128, 64)
	for i := range input1 {
		input1[i] = complex(float64(i), 0)
		input2[i] = complex(float64(i), 0)
	}

	output := make([]complex128, 64)
	fb.forward(input1, output)
	fb.forwardInPlace(input2)

	if !complexSlicesEqual(output, input2, 1e-10) {
		t.Error("forwardInPlace doesn't match forward")
	}
}

func TestBackwardInPlace(t *testing.T) {
	grid, _ := gridData.NewFromLength(6.4, 64)
	fb := FFTInit(grid, 1.0)
	defer fb.Clean()

	input1 := make([]complex128, 64)
	input2 := make([]complex128, 64)
	for i := range input1 {
		input1[i] = complex(float64(i), float64(i*2))
		input2[i] = complex(float64(i), float64(i*2))
	}

	output := make([]complex128, 64)
	fb.backward(input1, output)
	fb.backwardInPlace(input2)

	if !complexSlicesEqual(output, input2, 1e-10) {
		t.Error("backwardInPlace doesn't match backward")
	}
}

// ============================================================================
// Test Parseval's Theorem (energy conservation)
// ============================================================================

func TestParsevalTheorem(t *testing.T) {
	grid, _ := gridData.NewFromLength(6.4, 64)
	fb := FFTInit(grid, 1.0)
	defer fb.Clean()

	input := make([]complex128, 64)
	for i := range input {
		input[i] = complex(math.Sin(float64(i)*0.1), math.Cos(float64(i)*0.2))
	}

	output := make([]complex128, 64)
	fb.forward(input, output)

	inputEnergy := 0.0
	outputEnergy := 0.0
	for i := range input {
		inputEnergy += real(input[i])*real(input[i]) + imag(input[i])*imag(input[i])
		outputEnergy += real(output[i])*real(output[i]) + imag(output[i])*imag(output[i])
	}

	expectedRatio := float64(fb.nPoints)
	actualRatio := outputEnergy / inputEnergy

	if math.Abs(actualRatio-expectedRatio) > 1e-10 {
		t.Errorf("Parseval's theorem violated: expected ratio %f, got %f", expectedRatio, actualRatio)
	}
}

// ============================================================================
// Test Operator Operations
// ============================================================================

func TestMomentumOpInPlace(t *testing.T) {
	grid, _ := gridData.NewFromLength(6.4, 64)
	fb := FFTInit(grid, 1.0)
	defer fb.Clean()

	input := make([]complex128, 64)
	for i := range input {
		input[i] = complex(1.0, 0)
	}

	fb.MomentumOpInPlace(input)

	// For a constant function, momentum operator should give zero (except at k=0)
	// This is because the FT of a constant is a delta at k=0
	// Verify the operation completed without panic
	if len(input) != 64 {
		t.Error("MomentumOpInPlace changed array length")
	}
}

func TestLaplacianOpInPlace(t *testing.T) {
	grid, _ := gridData.NewFromLength(6.4, 64)
	fb := FFTInit(grid, 1.0)
	defer fb.Clean()

	input := make([]complex128, 64)
	for i := range input {
		input[i] = complex(1.0, 0)
	}

	fb.LaplacianOpInPlace(input)

	// Verify the operation completed without panic
	if len(input) != 64 {
		t.Error("LaplacianOpInPlace changed array length")
	}
}

// BUG TEST: MomentumOp has redundant copy
func TestMomentumOpRedundantCopy(t *testing.T) {
	grid, _ := gridData.NewFromLength(6.4, 64)
	fb := FFTInit(grid, 1.0)
	defer fb.Clean()

	input := make([]complex128, 64)
	output1 := make([]complex128, 64)
	output2 := make([]complex128, 64)

	for i := range input {
		input[i] = complex(math.Sin(float64(i)*0.1), 0)
	}

	fb.MomentumOp(input, output1)
	copy(output2, input)
	fb.MomentumOpInPlace(output2)

	if !complexSlicesEqual(output1, output2, 1e-10) {
		t.Error("MomentumOp and MomentumOpInPlace give different results")
	}
}

// BUG TEST: LaplacianOp has redundant copy
func TestLaplacianOpRedundantCopy(t *testing.T) {
	grid, _ := gridData.NewFromLength(6.4, 64)
	fb := FFTInit(grid, 1.0)
	defer fb.Clean()

	input := make([]complex128, 64)
	output1 := make([]complex128, 64)
	output2 := make([]complex128, 64)

	for i := range input {
		input[i] = complex(math.Sin(float64(i)*0.1), 0)
	}

	fb.LaplacianOp(input, output1)
	copy(output2, input)
	fb.LaplacianOpInPlace(output2)

	if !complexSlicesEqual(output1, output2, 1e-10) {
		t.Error("LaplacianOp and LaplacianOpInPlace give different results")
	}
}

// ============================================================================
// Test Edge Cases and Error Handling
// ============================================================================

func TestOperatorOpLengthMismatch(t *testing.T) {
	grid, _ := gridData.NewFromLength(6.4, 64)
	fb := FFTInit(grid, 1.0)
	defer fb.Clean()

	defer func() {
		if r := recover(); r != nil {
			panicMsg := fmt.Sprintf("%v", r)
			if !strings.Contains(panicMsg, "length mismatch") {
				t.Errorf("Expected 'length mismatch' in panic, got: %s", panicMsg)
			}
		}
	}()

	input := make([]complex128, 64)
	wrongSizeOp := make([]float64, 32) // Wrong size

	fb.operatorOp(input, wrongSizeOp)
	t.Error("Expected panic for length mismatch")
}

func TestEmptyInput(t *testing.T) {
	grid, _ := gridData.NewFromLength(6.4, 64)
	fb := FFTInit(grid, 1.0)
	defer fb.Clean()

	input := make([]complex128, 1)
	input[0] = complex(1.0, 0)

	output := make([]complex128, 1)
	fb.forward(input, output)
	if cmplx.Abs(input[0]-output[0]) > 1e-10 {
		t.Errorf("FFT of single element failed: expected %v, got %v", input[0], output[0])
	}
}

// ============================================================================
// Test keValues Calculation
// ============================================================================

func TestKEValuesCalculation(t *testing.T) {
	grid, _ := gridData.NewFromLength(6.4, 64)
	mass := 2.0
	fb := FFTInit(grid, mass)
	defer fb.Clean()
	for i := 0; i < fb.nPoints; i++ {
		expected := fb.kValues[i] * fb.kValues[i] / (2.0 * mass)
		if math.Abs(fb.keValues[i]-expected) > 1e-14 {
			t.Errorf("keValues[%d] incorrect: expected %f, got %f", i, expected, fb.keValues[i])
		}
	}
}

func TestPlaneWaveMomentum(t *testing.T) {
	n := 128
	L := 10.0
	grid, _ := gridData.NewFromLength(L, uint32(n))
	fb := FFTInit(grid, 1.0)
	defer fb.Clean()
	k0 := 2.0 * math.Pi / L * 3.0 // k = 3 * 2π/L

	psi := make([]complex128, n)
	for i := 0; i < n; i++ {
		x := float64(i) * grid.DeltaR()
		psi[i] = cmplx.Exp(complex(0, k0*x))
	}

	result := make([]complex128, n)
	copy(result, psi)
	fb.MomentumOpInPlace(result)

	norm := complex(float64(n), 0)
	for i := range result {
		result[i] /= norm
	}

	avgMomentum := complex(0, 0)
	for i := range result {
		if cmplx.Abs(psi[i]) > 1e-10 {
			avgMomentum += result[i] / psi[i]
		}
	}
	avgMomentum /= complex(float64(n), 0)
	t.Logf("Expected momentum: %f, Computed average: %v", k0, avgMomentum)
}

// ============================================================================
// Benchmarks
// ============================================================================

func BenchmarkForward(b *testing.B) {
	grid, _ := gridData.NewFromLength(102.4, 1024)
	fb := FFTInit(grid, 1.0)
	defer fb.Clean()

	input := make([]complex128, 1024)
	output := make([]complex128, 1024)
	for i := range input {
		input[i] = complex(float64(i), 0)
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		fb.forward(input, output)
	}
}

func BenchmarkForwardInPlace(b *testing.B) {
	grid, _ := gridData.NewFromLength(102.4, 1024)
	fb := FFTInit(grid, 1.0)
	defer fb.Clean()

	input := make([]complex128, 1024)
	for i := range input {
		input[i] = complex(float64(i), 0)
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		fb.forwardInPlace(input)
	}
}

func BenchmarkMomentumOp(b *testing.B) {
	grid, _ := gridData.NewFromLength(102.4, 1024)
	fb := FFTInit(grid, 1.0)
	defer fb.Clean()

	input := make([]complex128, 1024)
	output := make([]complex128, 1024)
	for i := range input {
		input[i] = complex(math.Sin(float64(i)*0.01), 0)
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		fb.MomentumOp(input, output)
	}
}

// This test ensures fmt is imported (needed for the panic message fix)
func TestFmtImported(t *testing.T) {
	_ = fmt.Sprintf("test")
}
