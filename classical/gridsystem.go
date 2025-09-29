package classical

import "GoProject/gridData"

type SystemInfo struct {
	position float64
	momentum float64
	kinE     float64
	potenE   float64
}

type ModelSystem struct {
	state SystemInfo
	mass  float64
	rgrid gridData.RadGrid
	tgrid gridData.TimeGrid
}

func NewModelSystem(state SystemInfo, mass float64, rgrid gridData.RadGrid, tgrid gridData.TimeGrid) *ModelSystem {
	return &ModelSystem{
		state: state,
		mass:  mass,
		rgrid: rgrid,
		tgrid: tgrid,
	}
}
