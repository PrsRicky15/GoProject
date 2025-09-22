package main

import (
	"fmt"
)

func main() {
	s := "Ricky"
	fmt.Printf("Hello and welcome, %s!\n", s)

	cards := newDeck()
	cards.printWithoutIndex()

	for i := 0; i <= 5; i++ {
		fmt.Println("i =", i)
	}
}
