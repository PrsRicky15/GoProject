package EquationSolver

// a simple decay equation: dx/dt = -x
type MockPotential struct{}

// EvaluateOnRGridTime implements [gridData.TDPotentialOp].
func (m *MockPotential) EvaluateOnRGridTime(x []float64, t float64) []float64 {
	panic("unimplemented")
}

func (m *MockPotential) EvaluateAtTime(x, t float64) float64 {
	// dx/dt = -x
	return -x
}

func (m *MockPotential) EvaluateOnRGridTimeInPlace(x []float64, buff []float64, t float64) {
	for i, v := range x {
		buff[i] = -v
	}
}

type MockGaussianPotential struct{}

func (m *MockGaussianPotential) EvaluateOnRGridTime(x []float64, t float64) []float64 {
	panic("unimplemented")
}

func (m *MockGaussianPotential) EvaluateAtTime(x, t float64) float64 {
	// dx/dt = -2tx => x(t) = x0 * exp(-t^2)
	return -2 * t * x
}

func (m *MockGaussianPotential) EvaluateOnRGridTimeInPlace(x []float64, buff []float64, t float64) {
	for i, v := range x {
		buff[i] = -2 * t * v
	}
}