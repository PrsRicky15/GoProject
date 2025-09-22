package main

import "fmt"

func basic() {
	var MyName string
	fmt.Println("Enter your name\n")
	_, _ = fmt.Scan(&MyName)
	fmt.Println("My name is ", MyName, " and I am an Engineer")

	var age int
	var months int
	age = 20
	months = 4
	fmt.Printf("My age is %v years and %v months\n", age, months)
}
