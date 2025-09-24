package gridData

type SystemInfo struct {
	position float64
	momentum float64
	kinE     float64
	potenE   float64
}

type ModelSystem struct {
	state SystemInfo
	mass  float64
	rgrid RadGrid
	tgrid TimeGrid
}

func NewModelSystem(state SystemInfo, mass float64, rgrid RadGrid, tgrid TimeGrid) *ModelSystem {
	return &ModelSystem{
		state: state,
		mass:  mass,
		rgrid: rgrid,
		tgrid: tgrid,
	}
}
