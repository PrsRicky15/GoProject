package model

import (
	"math"
	"testing"
)

func testFunc2(x float64) float64 {
	return math.Sin(x) // integral from 0 to π = 2
}

func testFunc3(x float64) float64 {
	return math.Exp(x) // integral from 0 to 1 = e - 1
}

func testFunc4(x float64) float64 {
	return 1.0 / (1.0 + x*x) // integral from 0 to 1 = π/4
}

// TestRectangleMethod tests the Rectangle method
func TestRectangleMethod(t *testing.T) {
	tests := []struct {
		name     string
		f        func(float64) float64
		a, b     float64
		expected float64
		nPoints  int
		tol      float64
	}{
		{"x^2 from 0 to 1", testFunc1, 0, 1, 1.0 / 3.0, 1000, 0.001},
		{"sin(x) from 0 to π", testFunc2, 0, math.Pi, 2.0, 1000, 0.01},
		{"e^x from 0 to 1", testFunc3, 0, 1, math.E - 1, 1000, 0.002},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := Rectangle(tt.f, tt.a, tt.b, tt.nPoints)
			error := math.Abs(result - tt.expected)
			if error > tt.tol {
				t.Errorf("Rectangle(%s) = %v, want %v (error: %v > tolerance: %v)",
					tt.name, result, tt.expected, error, tt.tol)
			}
		})
	}
}

// TestTrapezoidalMethod tests the Trapezoidal method
func TestTrapezoidalMethod(t *testing.T) {
	tests := []struct {
		name     string
		f        func(float64) float64
		a, b     float64
		expected float64
		nPoints  int
		tol      float64
	}{
		{"x^2 from 0 to 1", testFunc1, 0, 1, 1.0 / 3.0, 1000, 0.0001},
		{"sin(x) from 0 to π", testFunc2, 0, math.Pi, 2.0, 1000, 0.0001},
		{"e^x from 0 to 1", testFunc3, 0, 1, math.E - 1, 1000, 0.0001},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := Trapezoidal(tt.f, tt.a, tt.b, tt.nPoints)
			error := math.Abs(result - tt.expected)
			if error > tt.tol {
				t.Errorf("Trapezoidal(%s) = %v, want %v (error: %v > tolerance: %v)",
					tt.name, result, tt.expected, error, tt.tol)
			}
		})
	}
}

// TestSimpsonMethod tests the Simpson method
func TestSimpsonMethod(t *testing.T) {
	tests := []struct {
		name     string
		f        func(float64) float64
		a, b     float64
		expected float64
		nPoints  int
		tol      float64
	}{
		{"x^2 from 0 to 1", testFunc1, 0, 1, 1.0 / 3.0, 100, 0.00001},
		{"sin(x) from 0 to π", testFunc2, 0, math.Pi, 2.0, 100, 0.00001},
		{"e^x from 0 to 1", testFunc3, 0, 1, math.E - 1, 100, 0.00001},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := Simpson(tt.f, tt.a, tt.b, tt.nPoints)
			error := math.Abs(result - tt.expected)
			if error > tt.tol {
				t.Errorf("Simpson(%s) = %v, want %v (error: %v > tolerance: %v)",
					tt.name, result, tt.expected, error, tt.tol)
			}
		})
	}
}
