package main

import (
	"GoProject/gridData"
	"fmt"
)

func main() {
	s := "Ricky"
	fmt.Printf("Hello and welcome, %s!\n", s)

	grid, _ := gridData.NewFromLength(20., 30)
	grid.DisplayInfo()
	fmt.Println("kGrid points:")
	grid.DisplayRgrid()
	cards := newDeck()

	fmt.Println("")
	fmt.Println("Test card deck:")
	hand, _ := deal(cards, 5)
	hand.print()
}
