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

	grid, _ := gridData.NewRGridFromFile(filePath)
	grid.DisplayInfo()
	fmt.Println("kGrid points:")
	grid.DisplayRgrid()
	cards := newDeck()

	fmt.Println("")
	fmt.Println("Test card deck:")
	hand, _ := deal(cards, 5)
	hand.print()
}
