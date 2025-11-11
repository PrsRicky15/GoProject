package OperatorAlgebra

import (
	"GoProject/gridData"
	"testing"
)

func TestKeDVR_Evaluate(t *testing.T) {
	rgrid, err := gridData.NewFromLength(10., 30)
	if err != nil {
		panic(err)
	}
	kinE := NewKeDVR(rgrid, 1.)
	kinE.GetMat()
}
