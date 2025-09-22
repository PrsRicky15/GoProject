package main

import (
	"GoProject/input"
	"fmt"
)

func main() {
	s := "Ricky"
	fmt.Printf("Hello and welcome, %s!\n", s)

	grid := input.FromLength(20., 60)
	grid.Display_info()
	cards := newDeck()
	hand, _ := deal(cards, 5)

	hand.print()

	for i := 0; i <= 5; i++ {
		fmt.Println("i =", i)
	}
}
