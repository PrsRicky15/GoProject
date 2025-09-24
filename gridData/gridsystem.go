package gridData

type SystemState struct {
	Value []float64
}

type ModelSystem struct {
	state SystemState
	mass  float64
	rgrid RadGrid
	tgrid TimeGrid
}
