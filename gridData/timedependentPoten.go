package gridData

// TDPotentialOp General interface for the evaluating the potential on a grid
type TDPotentialOp interface {
	EvaluateAt(x float64, t float64) float64
	EvaluateOnRGrid(x []float64, t float64) []float64
	EvaluateOnRGridInPlace(x, res []float64, t float64)
}
