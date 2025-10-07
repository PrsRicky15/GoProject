package main

import (
	"GoProject/gridData"
	"fmt"
	"os"
)

func main() {

	if len(os.Args) < 2 {
		fmt.Println("Usage: go run main.go <file-path>")
		os.Exit(1)
	}

	filePath := os.Args[1]

	s := "Ricky"
	fmt.Printf("Hello and welcome, %s!\n", s)

	potent := gridData.SoftCore[float64]{Charge: 1., SoftParam: 0.5}
	grid, _ := gridData.NewRGridFromFile(filePath)
	fmt.Println(grid)
	grid.DisplayInfo()
	err := grid.PrintPotentToFileRe(potent, "softcore.dat", "%21.14e")
	err = grid.PrintForceToFileRe(potent, "forceSoftcore.dat", "%21.14e")
	if err != nil {
		panic(err)
	}
}
