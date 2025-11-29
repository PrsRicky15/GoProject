package BasicOneD

import (
	"GoProject/gridData"
)

type FiniteBarrier struct {
	grid  *gridData.RadGrid
	rFunc *gridData.Rfunc
}

func NewFiniteBarrier(grid *gridData.RadGrid, rFunc *gridData.Rfunc) *FiniteBarrier {

	if grid == nil || rFunc == nil {
		panic("grid and rFunc must not be nil")
	}

	return &FiniteBarrier{
		grid:  grid,
		rFunc: rFunc,
	}
}

func (fb *FiniteBarrier) ReDefine(grid *gridData.RadGrid, rFunc *gridData.Rfunc) {
	fb.rFunc = rFunc
	fb.grid = grid
}

func (fb *FiniteBarrier) ReDefineFunc(rFunc *gridData.Rfunc) {
	fb.rFunc = rFunc
}

type TDFiniteBarrier struct {
	barrier *FiniteBarrier
	tdGrid  *gridData.TimeGrid
}
