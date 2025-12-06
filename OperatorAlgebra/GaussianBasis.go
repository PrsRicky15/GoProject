package OperatorAlgebra

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"os"
	"strconv"
	"strings"
	"time"

	"gonum.org/v1/gonum/mat"
)

type Config struct {
	Norb      int
	Nelec     int
	Threshold float64
}

type IntegralData struct {
	Hcore []float64 // Core Hamiltonian (lower triangular)
	S     []float64 // Overlap matrix (lower triangular)
	TwoEl []float64 // Two-electron integrals
}

// HartreeFock represents the main HF calculation engine
type HartreeFock struct {
	config    Config
	integrals *IntegralData
	xMatrix   *mat.Dense
	coef      *mat.Dense
	density   *mat.Dense
	fock      *mat.Dense
	energy    float64
	iteration int
}

// NewHartreeFock creates a new HF calculator
func NewHartreeFock(config Config, integrals *IntegralData) *HartreeFock {
	norb := config.Norb
	return &HartreeFock{
		config:    config,
		integrals: integrals,
		coef:      mat.NewDense(norb, norb, nil),
		density:   mat.NewDense(norb, norb, nil),
		fock:      mat.NewDense(norb, norb, nil),
	}
}

// Run performs the SCF calculation
func (hf *HartreeFock) Run() error {
	start := time.Now()
	hf.xMatrix = hf.calculateXMatrix()
	hf.initialGuess()
	if err := hf.scfLoop(); err != nil {
		return err
	}

	elapsed := time.Since(start)
	fmt.Printf("SCF converged in %d iterations\n", hf.iteration)
	fmt.Printf("Final energy: %.8f\n", hf.energy)
	fmt.Printf("Time taken: %.5f seconds\n", elapsed.Seconds())

	return nil
}

// initialGuess generates the initial coefficient matrix
func (hf *HartreeFock) initialGuess() {
	norb := hf.config.Norb
	hmat := mat.NewDense(norb, norb, nil)
	idx := 0
	for i := 0; i < norb; i++ {
		for j := 0; j <= i; j++ {
			hmat.Set(i, j, hf.integrals.Hcore[idx])
			hmat.Set(j, i, hf.integrals.Hcore[idx])
			idx++
		}
	}

	temp := mat.NewDense(norb, norb, nil)
	temp.Mul(hmat, hf.xMatrix)

	fPrime := mat.NewDense(norb, norb, nil)
	fPrime.Mul(hf.xMatrix.T(), temp)
	hf.diagonalize(fPrime)
	hf.coef.Mul(hf.xMatrix, fPrime)
}

// scfLoop performs the self-consistent field iterations
func (hf *HartreeFock) scfLoop() error {
	norb := hf.config.Norb
	maxIter := 200

	errFile, err := os.Create("hferr.txt")
	if err != nil {
		return err
	}
	defer func(errFile *os.File) {
		err := errFile.Close()
		if err != nil {
			log.Fatal(err)
		}
	}(errFile)

	err = hf.writeErrorHeader(errFile)
	if err != nil {
		return err
	}

	prevCoef := mat.NewDense(norb, norb, nil)
	prevDensity := mat.NewDense(norb, norb, nil)
	prevEnergy := 0.

	for hf.iteration = 1; hf.iteration <= maxIter; hf.iteration++ {
		prevCoef.Copy(hf.coef)
		prevDensity.Copy(hf.density)
		prevEnergy = hf.energy

		hf.buildDensityMatrix()
		hf.buildFockMatrix()
		hf.solveFockEquation()
		hf.calculateEnergy()
		energyErr := math.Abs(prevEnergy - hf.energy)
		densityErr := hf.matrixDifference(prevDensity, hf.density)
		coefErr := hf.matrixDifference(prevCoef, hf.coef)

		_, err2 := fmt.Fprintf(errFile, "%-5d%-17.8f%-17.8f%-17.8f%-15.5f\n",
			hf.iteration, energyErr, densityErr, coefErr, hf.energy)
		if err2 != nil {
			return err2
		}

		if densityErr < hf.config.Threshold {
			return nil
		}
	}

	return fmt.Errorf("SCF did not converge in %d iterations", maxIter)
}

// buildDensityMatrix constructs the density matrix from coefficients
func (hf *HartreeFock) buildDensityMatrix() {
	norb := hf.config.Norb
	nocc := hf.config.Nelec / 2

	hf.density.Zero()
	for i := 0; i < norb; i++ {
		for j := 0; j < norb; j++ {
			sum := 0.
			for k := 0; k < nocc; k++ {
				sum += hf.coef.At(i, k) * hf.coef.At(j, k)
			}
			hf.density.Set(i, j, sum)
		}
	}
}

// buildFockMatrix constructs the Fock matrix
func (hf *HartreeFock) buildFockMatrix() {
	norb := hf.config.Norb
	gmat := hf.buildGMatrix()

	idx := 0
	for i := 0; i < norb; i++ {
		for j := 0; j <= i; j++ {
			val := hf.integrals.Hcore[idx] + gmat.At(i, j)
			hf.fock.Set(i, j, val)
			hf.fock.Set(j, i, val)
			idx++
		}
	}
}

// buildGMatrix constructs the two-electron contribution matrix
func (hf *HartreeFock) buildGMatrix() *mat.Dense {
	norb := hf.config.Norb
	gmat := mat.NewDense(norb, norb, nil)
	idx := 0

	for i := 0; i < norb; i++ {
		for j := 0; j < norb; j++ {
			for k := 0; k < norb; k++ {
				for l := 0; l < norb; l++ {
					val := hf.integrals.TwoEl[idx]
					idx++

					if math.Abs(val) > threshold {
						gmat.Set(i, j, gmat.At(i, j)+val*hf.density.At(k, l))
						gmat.Set(k, l, gmat.At(k, l)+val*hf.density.At(i, j))

						qval := 0.25 * val
						gmat.Set(i, k, gmat.At(i, k)-qval*hf.density.At(j, l))
						gmat.Set(j, l, gmat.At(j, l)-qval*hf.density.At(i, k))
						gmat.Set(i, l, gmat.At(i, l)-qval*hf.density.At(j, k))
						gmat.Set(j, k, gmat.At(j, k)-qval*hf.density.At(i, l))
					}
				}
			}
		}
	}

	return gmat
}

// solveFockEquation diagonalizes the Fock matrix in orthogonal basis
func (hf *HartreeFock) solveFockEquation() {
	norb := hf.config.Norb
	temp := mat.NewDense(norb, norb, nil)
	temp.Mul(hf.fock, hf.xMatrix)
	fPrime := mat.NewDense(norb, norb, nil)
	fPrime.Mul(hf.xMatrix.T(), temp)
	hf.diagonalize(fPrime)
	hf.coef.Mul(hf.xMatrix, fPrime)
}

// calculateEnergy computes the total electronic energy
func (hf *HartreeFock) calculateEnergy() {
	norb := hf.config.Norb
	hcore := mat.NewDense(norb, norb, nil)
	idx := 0
	for i := 0; i < norb; i++ {
		for j := 0; j <= i; j++ {
			hcore.Set(i, j, hf.integrals.Hcore[idx])
			hcore.Set(j, i, hf.integrals.Hcore[idx])
			idx++
		}
	}

	hf.energy = 0.
	for i := 0; i < norb; i++ {
		for j := 0; j < norb; j++ {
			hf.energy += hf.density.At(i, j) * (hcore.At(i, j) + hf.fock.At(i, j))
		}
	}
}

// calculateXMatrix computes the transformation matrix X = U * s^(-1/2)
func (hf *HartreeFock) calculateXMatrix() *mat.Dense {
	norb := hf.config.Norb

	// Build overlap matrix
	smat := mat.NewDense(norb, norb, nil)
	idx := 0
	for i := 0; i < norb; i++ {
		for j := 0; j <= i; j++ {
			smat.Set(i, j, hf.integrals.S[idx])
			smat.Set(j, i, hf.integrals.S[idx])
			idx++
		}
	}

	eigvals := hf.diagonalize(smat)

	sPowNegHalf := mat.NewDense(norb, norb, nil)
	for i := 0; i < norb; i++ {
		sPowNegHalf.Set(i, i, 1./math.Sqrt(eigvals[i]))
	}

	xTemp := mat.NewDense(norb, norb, nil)
	for i := 0; i < norb; i++ {
		for j := 0; j < norb; j++ {
			xTemp.Set(i, j, sPowNegHalf.At(i, i)*smat.At(j, i))
		}
	}

	xmat := mat.NewDense(norb, norb, nil)
	xmat.Mul(smat, xTemp)

	return xmat
}

// diagonalize performs eigenvalue decomposition
func (hf *HartreeFock) diagonalize(m *mat.Dense) []float64 {
	var eig mat.EigenSym
	ok := eig.Factorize(m, true)
	if !ok {
		panic("eigenvalue decomposition failed")
	}

	norb := hf.config.Norb
	eigvals := make([]float64, norb)
	eig.Values(eigvals)

	eigvecs := mat.NewDense(norb, norb, nil)
	eig.VectorsTo(eigvecs)

	// Copy eigenvectors back to m
	for i := 0; i < norb; i++ {
		for j := 0; j < norb; j++ {
			m.Set(i, j, eigvecs.At(i, j))
		}
	}

	return eigvals
}

// matrixDifference calculates the Frobenius norm of the difference
func (hf *HartreeFock) matrixDifference(a, b *mat.Dense) float64 {
	r, c := a.Dims()
	sum := 0.
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			diff := a.At(i, j) - b.At(i, j)
			sum += diff * diff
		}
	}
	return math.Sqrt(sum)
}

// writeErrorHeader writes the convergence tracking header
func (hf *HartreeFock) writeErrorHeader(f *os.File) error {
	dash := "---------------------"
	_, err := fmt.Fprintf(f, "%s%s%s%s%s\n", dash, dash, dash, dash, dash)
	if err != nil {
		return err
	}
	_, err = fmt.Fprintf(f, "%-5s%-17s%-17s%-17s%-15s\n",
		"S.NO.", "Energy Err.", "Density Err.", "Coef. Err.", "Energy")
	if err != nil {
		return err
	}
	_, err = fmt.Fprintf(f, "%s%s%s%s%s\n", dash, dash, dash, dash, dash)
	if err != nil {
		return err
	}
	return nil
}

// WriteCoefficients exports the coefficient matrix to a file
func (hf *HartreeFock) WriteCoefficients(filename string) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer func(file *os.File) {
		err := file.Close()
		if err != nil {
			log.Fatal(err)
		}
	}(file)

	norb := hf.config.Norb
	for i := 0; i < norb; i++ {
		for j := 0; j < norb; j++ {
			_, err := fmt.Fprintf(file, "%.15e ", hf.coef.At(i, j))
			if err != nil {
				return err
			}
		}
		_, err := fmt.Fprintln(file)
		if err != nil {
			return err
		}
	}

	return nil
}

// ConfigReader handles reading configuration files
type ConfigReader struct{}

func (cr *ConfigReader) ReadConfig(filename string) (Config, error) {
	file, err := os.Open(filename)
	if err != nil {
		return Config{}, err
	}
	defer func(file *os.File) {
		err := file.Close()
		if err != nil {
			log.Fatal(err)
		}
	}(file)

	scanner := bufio.NewScanner(file)
	var config Config

	scanner.Scan()
	config.Norb, _ = strconv.Atoi(strings.TrimSpace(scanner.Text()))
	scanner.Scan()
	config.Nelec, _ = strconv.Atoi(strings.TrimSpace(scanner.Text()))
	scanner.Scan()
	config.Threshold, _ = strconv.ParseFloat(strings.TrimSpace(scanner.Text()), 64)

	return config, nil
}

// IntegralReader handles reading integral files
type IntegralReader struct{}

func (ir *IntegralReader) ReadIntegrals(norb int) (*IntegralData, error) {
	nlowr := (norb * (norb + 1)) / 2
	norb4 := norb * norb * norb * norb

	integrals := &IntegralData{
		Hcore: make([]float64, nlowr),
		S:     make([]float64, nlowr),
		TwoEl: make([]float64, norb4),
	}

	// Read one-electron integrals
	if err := ir.readOneElectron("../matrx/onelec.txt", norb, integrals); err != nil {
		return nil, err
	}

	// Read two-electron integrals
	if err := ir.readTwoElectron("../matrx/twelc.txt", norb, integrals); err != nil {
		return nil, err
	}

	return integrals, nil
}

func (ir *IntegralReader) readOneElectron(filename string, norb int,
	integrals *IntegralData) error {
	file, err := os.Open(filename)
	if err != nil {
		return err
	}
	defer func(file *os.File) {
		err := file.Close()
		if err != nil {
			log.Fatal(err)
		}
	}(file)

	scanner := bufio.NewScanner(file)
	idx := 0

	for i := 1; i <= norb; i++ {
		for j := 1; j <= i; j++ {
			scanner.Scan()
			fields := strings.Fields(scanner.Text())

			tval, _ := strconv.ParseFloat(fields[0], 64)
			vval, _ := strconv.ParseFloat(fields[1], 64)
			sval, _ := strconv.ParseFloat(fields[3], 64)

			integrals.S[idx] = sval
			integrals.Hcore[idx] = tval + vval
			idx++
		}
	}

	return nil
}

func (ir *IntegralReader) readTwoElectron(filename string, norb int,
	integrals *IntegralData) error {
	file, err := os.Open(filename)
	if err != nil {
		return err
	}
	defer func(file *os.File) {
		err := file.Close()
		if err != nil {
			log.Fatal(err)
		}
	}(file)

	scanner := bufio.NewScanner(file)
	scanner.Scan()
	num, _ := strconv.Atoi(strings.TrimSpace(scanner.Text()))

	norb2 := norb * norb
	norb3 := norb * norb * norb

	for i := 0; i < num; i++ {
		scanner.Scan()
		fields := strings.Fields(scanner.Text())

		iorb, _ := strconv.Atoi(fields[0])
		jorb, _ := strconv.Atoi(fields[1])
		korb, _ := strconv.Atoi(fields[2])
		lorb, _ := strconv.Atoi(fields[3])
		vval, _ := strconv.ParseFloat(fields[5], 64)

		// Store all 8 permutations
		indices := []int{
			norb3*(iorb-1) + norb2*(jorb-1) + norb*(korb-1) + lorb - 1,
			norb3*(iorb-1) + norb2*(jorb-1) + norb*(lorb-1) + korb - 1,
			norb3*(jorb-1) + norb2*(iorb-1) + norb*(korb-1) + lorb - 1,
			norb3*(jorb-1) + norb2*(iorb-1) + norb*(lorb-1) + korb - 1,
			norb3*(korb-1) + norb2*(lorb-1) + norb*(iorb-1) + jorb - 1,
			norb3*(korb-1) + norb2*(lorb-1) + norb*(jorb-1) + iorb - 1,
			norb3*(lorb-1) + norb2*(korb-1) + norb*(iorb-1) + jorb - 1,
			norb3*(lorb-1) + norb2*(korb-1) + norb*(jorb-1) + iorb - 1,
		}

		for _, idx := range indices {
			integrals.TwoEl[idx] = vval
		}
	}

	return nil
}

func CheckCode() {
	configReader := &ConfigReader{}
	config, err := configReader.ReadConfig("inpar.inp")
	if err != nil {
		panic(err)
	}

	integralReader := &IntegralReader{}
	integrals, err := integralReader.ReadIntegrals(config.Norb)
	if err != nil {
		panic(err)
	}

	hf := NewHartreeFock(config, integrals)
	if err := hf.Run(); err != nil {
		panic(err)
	}

	if err := hf.WriteCoefficients("coef.txt"); err != nil {
		panic(err)
	}
}
