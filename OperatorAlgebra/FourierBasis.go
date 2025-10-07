package OperatorAlgebra

import (
	"GoProject/gridData"
)

// FourierBasis represents the kinetic energy operator in DVR basis with caching
type FourierBasis struct {
	grid  *gridData.RadGrid
	mass  float64
	ndims int
}
