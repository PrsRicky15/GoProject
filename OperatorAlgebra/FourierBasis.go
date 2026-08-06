package OperatorAlgebra

import (
	"GoProject/gridData"
	"fmt"
	"math"
	"math/rand/v2"

	"github.com/jvlmdr/go-fftw/fftw"
	"gonum.org/v1/gonum/mat"
)

type FFTCalculator struct {
	grid     *gridData.RadGrid
	mass     float64
	nPoints  int
	fftPlan  fftw.Plan
	ifftPlan fftw.Plan
	buffer   *fftw.Array
	kValues  []float64
	keValues []float64
}

func NewFFTCalculator(grid *gridData.RadGrid, mass float64) *FFTCalculator {
	calc := &FFTCalculator{}
	calc.initialize(grid, mass)
	return calc
}

func (c *FFTCalculator) Reconfigure(grid *gridData.RadGrid, mass float64) {
	c.destroyPlans()
	c.initialize(grid, mass)
}

// initialize sets up the FFT calculator with the given grid and mass.
func (c *FFTCalculator) initialize(grid *gridData.RadGrid, mass float64) {
	gridPoints := int(grid.NPoints())
	kVals := grid.KValues()

	kValues := make([]float64, gridPoints)
	keValues := make([]float64, gridPoints)
	invTwoMass := 1.0 / (2.0 * mass)

	for i := 0; i < gridPoints; i++ {
		k := kVals[i]
		kValues[i] = k
		keValues[i] = k * k * invTwoMass
	}

	buffer := fftw.NewArray(gridPoints)
	for i := 0; i < gridPoints; i++ {
		re := rand.NormFloat64()
		im := rand.NormFloat64()
		mag := math.Sqrt(re*re + im*im)
		buffer.Set(i, complex(re/mag, im/mag))
	}

	// Create FFTW plans
	fftPlan := *fftw.NewPlan(buffer, buffer, fftw.Forward, fftw.Estimate)
	ifftPlan := *fftw.NewPlan(buffer, buffer, fftw.Backward, fftw.Estimate)

	c.grid = grid
	c.mass = mass
	c.nPoints = gridPoints
	c.fftPlan = fftPlan
	c.ifftPlan = ifftPlan
	c.buffer = buffer
	c.kValues = kValues
	c.keValues = keValues
}

// Forward performs a forward FFT on the input and stores the result in output.
func (c *FFTCalculator) Forward(input, output []complex128) {
	if len(input) != c.nPoints || len(output) != c.nPoints {
		panic(fmt.Sprintf("invalid slice lengths: input=%d, output=%d, expected=%d",
			len(input), len(output), c.nPoints))
	}
	copy(c.buffer.Elems, input)
	c.fftPlan.Execute()
	copy(output, c.buffer.Elems)
}

// ForwardInPlace performs a forward FFT in place on the input slice.
func (c *FFTCalculator) ForwardInPlace(data []complex128) {
	if len(data) != c.nPoints {
		panic(fmt.Sprintf("invalid slice length: got=%d, expected=%d",
			len(data), c.nPoints))
	}
	copy(c.buffer.Elems, data)
	c.fftPlan.Execute()
	copy(data, c.buffer.Elems)
}

// Backward performs a backward FFT on the input and stores the result in output.
func (c *FFTCalculator) Backward(input, output []complex128) {
	if len(input) != c.nPoints || len(output) != c.nPoints {
		panic(fmt.Sprintf("invalid slice lengths: input=%d, output=%d, expected=%d",
			len(input), len(output), c.nPoints))
	}
	copy(c.buffer.Elems, input)
	c.ifftPlan.Execute()
	copy(output, c.buffer.Elems)
}

// BackwardInPlace performs a backward FFT in place on the input slice.
func (c *FFTCalculator) BackwardInPlace(data []complex128) {
	if len(data) != c.nPoints {
		panic(fmt.Sprintf("invalid slice length: got=%d, expected=%d",
			len(data), c.nPoints))
	}
	copy(c.buffer.Elems, data)
	c.ifftPlan.Execute()
	copy(data, c.buffer.Elems)
}

// ApplyOperator applies a diagonal operator in Fourier space to the input.
// It performs: output = IFFT(operator * FFT(input))
func (c *FFTCalculator) ApplyOperator(data []complex128, operator []float64) {
	if len(data) != c.nPoints || len(operator) != c.nPoints {
		panic(fmt.Sprintf("length mismatch: data=%d, operator=%d, expected=%d",
			len(data), len(operator), c.nPoints))
	}

	// Forward transform
	copy(c.buffer.Elems, data)
	c.fftPlan.Execute()

	// Apply operator in Fourier space
	for i := range c.buffer.Elems {
		c.buffer.Elems[i] *= complex(operator[i], 0)
	}

	// Backward transform
	c.ifftPlan.Execute()
	copy(data, c.buffer.Elems)
}

// ApplyMomentum applies the momentum operator to the input in place.
func (c *FFTCalculator) ApplyMomentum(data []complex128) {
	c.ApplyOperator(data, c.kValues)
}

// ApplyLaplacian applies the Laplacian operator (kinetic energy) to the input in place.
func (c *FFTCalculator) ApplyLaplacian(data []complex128) {
	c.ApplyOperator(data, c.keValues)
}

// BuildHamiltonian constructs the Hamiltonian matrix using the Fourier Grid Hamiltonian method.
// Returns a symmetric dense matrix representing the kinetic energy operator.
func (c *FFTCalculator) BuildHamiltonian() *mat.SymDense {
	// Compute kinetic energy values for each k mode
	kineticEnergy := make([]float64, c.nPoints)
	prefactor := 1.0 / (2.0 * c.mass)

	for r := 1; r <= c.nPoints; r++ {
		k := 2.0 * math.Pi * float64(r) / c.grid.Length()
		kineticEnergy[r-1] = prefactor * k * k
	}

	// Build Hamiltonian matrix elements
	hamiltonian := mat.NewSymDense(c.nPoints, nil)
	for i := 0; i < c.nPoints; i++ {
		for j := 0; j <= i; j++ {
			delta := i - j
			sum := 0.0

			// Sum over all k modes
			for r := 1; r <= c.nPoints; r++ {
				angle := 2.0 * math.Pi * float64(r*delta) / float64(c.nPoints)
				sum += math.Cos(angle) * kineticEnergy[r-1]
			}

			hij := (2.0 / float64(c.nPoints)) * sum
			hamiltonian.SetSym(i, j, hij)
		}
	}

	return hamiltonian
}

// Diagonalize computes the eigenvalues and eigenvectors of the Hamiltonian.
// Returns eigenvalues as a slice and eigenvectors as columns of a dense matrix.
func (c *FFTCalculator) Diagonalize() (eigenvalues []float64, eigenvectors *mat.Dense, err error) {
	hamiltonian := c.BuildHamiltonian()

	var eigenSolver mat.EigenSym
	ok := eigenSolver.Factorize(hamiltonian, true)
	if !ok {
		return nil, nil, fmt.Errorf("eigenvalue decomposition failed")
	}

	eigenvalues = eigenSolver.Values(nil)
	eigenvectors = mat.NewDense(len(eigenvalues), len(eigenvalues), nil)
	eigenSolver.VectorsTo(eigenvectors)

	return eigenvalues, eigenvectors, nil
}

// destroyPlans releases the FFTW plan resources.
func (c *FFTCalculator) destroyPlans() {
	c.fftPlan.Destroy()
	c.ifftPlan.Destroy()
}

// cleanup is called by the finalizer to release resources.
func (c *FFTCalculator) cleanup() {
	c.destroyPlans()
}

func (c *FFTCalculator) Close() error {
	c.destroyPlans()
	return nil
}

// FourierGridHamiltonian represents a Hamiltonian operator in the Fourier grid basis.
type FourierGridHamiltonian struct {
	calculator *FFTCalculator
}

// NewFourierGridHamiltonian creates a new Fourier grid Hamiltonian.
func NewFourierGridHamiltonian(grid *gridData.RadGrid, mass float64) *FourierGridHamiltonian {
	return &FourierGridHamiltonian{
		calculator: NewFFTCalculator(grid, mass),
	}
}

// ApplyMomentum applies the momentum operator to the wavefunction in place.
func (h *FourierGridHamiltonian) ApplyMomentum(wavefunction []complex128) {
	h.calculator.ApplyMomentum(wavefunction)
}

// ApplyKineticEnergy applies the kinetic energy operator to the wavefunction in place.
func (h *FourierGridHamiltonian) ApplyKineticEnergy(wavefunction []complex128) {
	h.calculator.ApplyLaplacian(wavefunction)
}

// Close releases all resources associated with the Hamiltonian.
func (h *FourierGridHamiltonian) Close() error {
	return h.calculator.Close()
}
