package main

import (
	"GoProject/gridData"
	"GoProject/model"
)

func main() {
	grid, err := gridData.NewRGrid(-3, 8, 100)
	if err != nil {
		panic(err)
	}
	v0x := make([]float64, 3)
	v0x[0] = 0
	v0x[1] = 1
	v0x[2] = -0.5

	aPoint := make([]float64, 3)
	aPoint[0] = 0
	aPoint[1] = 5
	aPoint[2] = 8

	myfunc := gridData.NewRectBarrier(v0x, aPoint)
	model.Barrier(*grid, myfunc.EvaluateAt)

}
