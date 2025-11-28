package NumericalRecipies

import "math/rand"

func RandomNumberInt(n int) int {
	value := rand.Intn(n)
	return value
}
