package main

import (
	"GoProject/input"
	"fmt"
)

func main() {
	s := "Ricky"
	fmt.Printf("Hello and welcome, %s!\n", s)

	grid, _ := input.NewFromLength(20., 30)
	fmt.Println("kGrid points:")
	grid.DisplayKgrid()
	cards := newDeck()

	fmt.Println("")
	fmt.Println("Test card deck:")
	hand, _ := deal(cards, 5)
	hand.print()
}
