package NumericalRecipies

type Integration interface {
	Integrate(min, max float64, f func(x float64) float64) float64
	LogIntegrate(min, max float64, f func(x float64) float64) float64
}

type SimpsonIntegrator struct{}
type TrapezoidIntegrator struct{}
type GaussIntegrator struct{}
