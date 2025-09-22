package main

import "fmt"

type cardDeck []string

func newDeck() cardDeck {
	nDeck := cardDeck{}
	cardSuit := []string{"Spades", "Hearts", "Clubs", "Diamonds"}
	cardValues := []string{"Ace", "Two", "Three", "Four", "Five", "Six",
		"Seven", "Eight", "Nine", "Ten", "Jack", "Queen", "King"}

	for _, suit := range cardSuit {
		for _, value := range cardValues {
			nDeck = append(nDeck, suit+"_of_"+value)
		}
	}
	return nDeck
}

func (d cardDeck) print() {
	for i, card := range d {
		fmt.Println(i, card)
	}
}

func (d cardDeck) printWithoutIndex() {
	for _, card := range d {
		fmt.Println(card)
	}
}
