package classical

import (
	"fmt"
	"math"
	"math/rand/v2"
)

// LJSystem is a 2D periodic, monatomic model in reduced units: m=sigma=epsilon=1.
// It uses a force-shifted cutoff and velocity-Verlet. It is not a chemical force field.
type LJSystem struct {
	Box, Cutoff  float64
	X, Y, VX, VY []float64
	fx, fy       []float64
	Potential    float64
}

func NewLJSystem(n int, box, temperature float64, seed uint64) (*LJSystem, error) {
	if n < 4 || n > 64 || !finiteMD(box) || box < 5 || box > 30 || !finiteMD(temperature) || temperature < 0 || temperature > 2 {
		return nil, fmt.Errorf("use 4–64 particles, box length 5–30, and temperature 0–2")
	}
	side := int(math.Ceil(math.Sqrt(float64(n))))
	spacing := box / float64(side)
	if spacing < 1.1 {
		return nil, fmt.Errorf("box is too small for the initial lattice; increase box size")
	}
	s := &LJSystem{Box: box, Cutoff: math.Min(2.5, box/2), X: make([]float64, n), Y: make([]float64, n), VX: make([]float64, n), VY: make([]float64, n), fx: make([]float64, n), fy: make([]float64, n)}
	rng := rand.New(rand.NewPCG(seed, seed+1))
	mx, my := 0., 0.
	for i := range n {
		s.X[i] = (float64(i%side) + .5) * spacing
		s.Y[i] = (float64(i/side) + .5) * spacing
		s.VX[i] = rng.NormFloat64()
		s.VY[i] = rng.NormFloat64()
		mx += s.VX[i]
		my += s.VY[i]
	}
	for i := range n {
		s.VX[i] -= mx / float64(n)
		s.VY[i] -= my / float64(n)
	}
	scale := math.Sqrt(float64(n-1) * temperature / s.Kinetic())
	for i := range n {
		s.VX[i] *= scale
		s.VY[i] *= scale
	}
	if err := s.Forces(); err != nil {
		return nil, err
	}
	return s, nil
}

func finiteMD(x float64) bool { return !math.IsNaN(x) && !math.IsInf(x, 0) }
func rawLJ(r float64) (energy, force float64) {
	q := math.Pow(1/r, 6)
	return 4 * q * (q - 1), 24 * q * (2*q - 1) / r
}

func (s *LJSystem) Forces() error {
	clear(s.fx)
	clear(s.fy)
	s.Potential = 0
	uc, fc := rawLJ(s.Cutoff)
	for i := range s.X {
		for j := 0; j < i; j++ {
			dx, dy := s.X[i]-s.X[j], s.Y[i]-s.Y[j]
			dx -= s.Box * math.Round(dx/s.Box)
			dy -= s.Box * math.Round(dy/s.Box)
			r2 := dx*dx + dy*dy
			if !finiteMD(r2) || r2 < .36 {
				return fmt.Errorf("particles became too close; reduce dt or temperature")
			}
			if r2 >= s.Cutoff*s.Cutoff {
				continue
			}
			r := math.Sqrt(r2)
			u, f := rawLJ(r)
			f = (f - fc) / r
			s.Potential += u - uc + (r-s.Cutoff)*fc
			s.fx[i] += f * dx
			s.fy[i] += f * dy
			s.fx[j] -= f * dx
			s.fy[j] -= f * dy
		}
	}
	return nil
}
func (s *LJSystem) Kinetic() float64 {
	k := 0.
	for i := range s.VX {
		k += .5 * (s.VX[i]*s.VX[i] + s.VY[i]*s.VY[i])
	}
	return k
}
func (s *LJSystem) Temperature() float64 { return s.Kinetic() / float64(len(s.X)-1) }
func (s *LJSystem) Step(dt float64) error {
	if !finiteMD(dt) || dt <= 0 || dt > .01 {
		return fmt.Errorf("Lennard-Jones dt must be positive and <= 0.01")
	}
	for i := range s.X {
		s.VX[i] += .5 * dt * s.fx[i]
		s.VY[i] += .5 * dt * s.fy[i]
		s.X[i] += dt * s.VX[i]
		s.Y[i] += dt * s.VY[i]
		s.X[i] -= s.Box * math.Floor(s.X[i]/s.Box)
		s.Y[i] -= s.Box * math.Floor(s.Y[i]/s.Box)
	}
	if err := s.Forces(); err != nil {
		return err
	}
	for i := range s.X {
		s.VX[i] += .5 * dt * s.fx[i]
		s.VY[i] += .5 * dt * s.fy[i]
	}
	return nil
}
