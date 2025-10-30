package NumericalRecipies

type Integration interface {
	Trapazoidal(dx float64) float64
	LogTrapazoidal(dx float64) float64
	Simpsons(dx float64) float64
	LogSimpsons(dx float64) float64
}
