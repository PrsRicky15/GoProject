package input

import (
	"fmt"
	"math"
)

type RGrid struct {
	rmin    float64
	rmax    float64
	nr      uint32
	dR      float64
	len     float64
	dK      float64
	kmin    float64
	kmax    float64
	cutOffe float64
}

func NewRGrid(rmin float64, rmax float64, ngrid uint32) *RGrid {
	dR := (rmax - rmin) / float64(ngrid)
	length := rmax - rmin
	dK := 2 * math.Pi / length
	kmin := -math.Pi / dR
	kmax := math.Pi / dR
	cutOffe := math.Pow(kmax, 2) / 2

	return &RGrid{
		rmin:    rmin,
		rmax:    rmax,
		nr:      ngrid,
		dR:      dR,
		len:     length,
		dK:      dK,
		kmin:    kmin,
		kmax:    kmax,
		cutOffe: cutOffe,
	}
}

func FromLength(length float64, ngrid uint32) *RGrid {
	grdMin := -length / 2
	grdMax := length / 2
	return NewRGrid(grdMin, grdMax, ngrid)
}

func (grid RGrid) DisplayInfo() {
	fmt.Println("Rmin: {} Rmax: {} Nr: {}", grid.rmin, grid.rmax, grid.dR)
	fmt.Println(grid.kmin, grid.kmax, grid.dK)
	fmt.Println(grid.cutOffe)
	fmt.Println(grid.len)
	fmt.Println(grid.nr)
}
