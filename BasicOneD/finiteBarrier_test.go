package BasicOneD

import (
	"math"
	"testing"
)

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
		{"1/(1+x^2) from 0 to 1", testFunc4, 0, 1, math.Pi / 4., 1000, 0.002},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := Rectangle(tt.f, tt.a, tt.b, tt.nPoints)
			err := math.Abs(result - tt.expected)
			if err > tt.tol {
				t.Errorf("Rectangle(%s) = %v, want %v (err: %v > tolerance: %v)",
					tt.name, result, tt.expected, err, tt.tol)
			}
		})
	}
}

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
		{"1/(1+x^2) from 0 to 1", testFunc4, 0, 1, math.Pi / 4., 100, 0.0001},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := Trapezoidal(tt.f, tt.a, tt.b, tt.nPoints)
			err := math.Abs(result - tt.expected)
			if err > tt.tol {
				t.Errorf("Trapezoidal(%s) = %v, want %v (err: %v > tolerance: %v)",
					tt.name, result, tt.expected, err, tt.tol)
			}
		})
	}
}

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
		{"1/(1+x^2) from 0 to 1", testFunc4, 0, 1, math.Pi / 4., 100, 0.00001},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := Simpson(tt.f, tt.a, tt.b, tt.nPoints)
			err := math.Abs(result - tt.expected)
			if err > tt.tol {
				t.Errorf("Simpson(%s) = %v, want %v (err: %v > tolerance: %v)",
					tt.name, result, tt.expected, err, tt.tol)
			}
		})
	}
}
