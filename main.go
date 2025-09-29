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

	potent := gridData.SoftCore{Charge: 1., SoftParam: 0.5}
	grid, _ := gridData.NewRGridFromFile(filePath)
	err := grid.PrintToFile(potent, "softcore.dat", "%21.14e")
	if err != nil {
		return
	}
	grid.DisplayInfo()
	fmt.Println("kGrid points:")
	grid.DisplayRgrid()
	cards := newDeck()

	fmt.Println("")
	fmt.Println("Test card deck:")
	hand, _ := deal(cards, 5)
	hand.print()
}
