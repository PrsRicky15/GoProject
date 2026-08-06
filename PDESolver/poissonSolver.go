package PDESolver

import (
	"GoProject/gridData"
	"errors"
	"fmt"
	"math"
)

type PoissonSolver1D struct {
	grid *gridData.RadGrid
	fx   gridData.Rfunc

	maxIter int
	tol     float64

	nMin1  int
	dx2    float64
	rho    []float64
	phiPre []float64
	phiNew []float64
}

func NewOneDimPoissonSolver(grid *gridData.RadGrid, fx gridData.Rfunc) (*PoissonSolver1D, error) {
	if grid == nil || fx == nil {
		return nil, errors.New("grid or function is not defined")
	}

	return &PoissonSolver1D{
		grid: grid,
		fx:   fx,
	}, nil
}

func (PS1D *PoissonSolver1D) Redefine(grid *gridData.RadGrid, fx gridData.Rfunc) error {
	if grid == nil {
		return errors.New("grid cannot be nil")
	}
	if fx == nil {
		return errors.New("function cannot be nil")
	}

	PS1D.grid = grid
	PS1D.fx = fx
	return nil
}

func (PS1D *PoissonSolver1D) RedefineFunc(fx gridData.Rfunc) error {
	if fx == nil {
		return errors.New("function cannot be nil")
	}

	PS1D.fx = fx
	return nil
}

func (PS1D *PoissonSolver1D) Print() {
	fmt.Println(PS1D.fx)
	fmt.Println(PS1D.grid)
}

func (PS1D *PoissonSolver1D) Initialize(maxIter int, tolerance, phiA, phiB float64) error {
	if maxIter < 1 {
		return errors.New("maxIter must be greater than zero")
	}
	if tolerance <= 0 || tolerance > 1 {
		return errors.New("tolerance must be between 0 and 1")
	}
	if math.IsNaN(phiA) || math.IsInf(phiA, 0) {
		return errors.New("phiA must be a finite number")
	}
	if math.IsNaN(phiB) || math.IsInf(phiB, 0) {
		return errors.New("phiB must be a finite number")
	}

	PS1D.maxIter = maxIter
	PS1D.tol = tolerance
	PS1D.dx2 = PS1D.grid.DeltaR() * PS1D.grid.DeltaR()
	PS1D.nMin1 = int(PS1D.grid.NPoints() - 1)

	// Boundary condition
	PS1D.phiPre = make([]float64, PS1D.grid.NPoints())
	PS1D.phiNew = make([]float64, PS1D.grid.NPoints())
	PS1D.phiPre[0] = phiA
	PS1D.phiPre[PS1D.nMin1] = phiB // Use nMin1 for clarity

	copy(PS1D.phiNew, PS1D.phiPre)

	PS1D.rho = make([]float64, PS1D.grid.NPoints())
	PS1D.grid.FunctionOnGridInPlace(PS1D.fx, PS1D.rho)

	return nil
}

// IteratingStepMethod Jacobi iteration method
func (PS1D *PoissonSolver1D) IteratingStepMethod() ([]float64, int, error) {
	for iter := 0; iter < PS1D.maxIter; iter++ {
		errL2Norm := 0.0

		// Update all interior points
		for i := 1; i < PS1D.nMin1; i++ {
			PS1D.phiNew[i] = 0.5 * (PS1D.phiPre[i+1] + PS1D.phiPre[i-1] - PS1D.rho[i]*PS1D.dx2)
			diff := PS1D.phiNew[i] - PS1D.phiPre[i]
			errL2Norm += diff * diff
		}

		// Copy new solution to previous
		copy(PS1D.phiPre, PS1D.phiNew)

		// Check convergence
		if math.Sqrt(errL2Norm) < PS1D.tol {
			return PS1D.phiNew, iter + 1, nil
		}
	}

	return nil, PS1D.maxIter, fmt.Errorf("did not converge in %d iterations", PS1D.maxIter)
}

// IteratingStepWithOmega
//
//	omega = 1.0 → Gauss-Seidel
//	omega > 1.0 → over-relaxation
//	omega < 1.0 → under-relaxation
func (PS1D *PoissonSolver1D) IteratingStepWithOmega(omega float64) ([]float64, int, error) {
	if omega <= 0 || omega >= 2 {
		return nil, 0, errors.New("omega must be in range (0, 2) for stability")
	}

	for iter := 0; iter < PS1D.maxIter; iter++ {
		errL2Norm := 0.0

		// SOR update: use in-place updates (Gauss-Seidel style)
		for i := 1; i < PS1D.nMin1; i++ {
			// Standard Jacobi update
			phiJacobi := 0.5 * (PS1D.phiNew[i+1] + PS1D.phiNew[i-1] - PS1D.rho[i]*PS1D.dx2)

			// SOR: weighted combination of old and new values
			phiOld := PS1D.phiNew[i]
			PS1D.phiNew[i] = (1.0-omega)*phiOld + omega*phiJacobi

			diff := PS1D.phiNew[i] - phiOld
			errL2Norm += diff * diff
		}

		// Check convergence
		if math.Sqrt(errL2Norm) < PS1D.tol {
			// Copy the final result to phiPre for consistency
			copy(PS1D.phiPre, PS1D.phiNew)
			return PS1D.phiNew, iter + 1, nil
		}
	}

	return nil, PS1D.maxIter, fmt.Errorf("did not converge in %d iterations", PS1D.maxIter)
}

func (PS1D *PoissonSolver1D) CheckError() error {
	// Test without an omega (Jacobi method)
	err := PS1D.Initialize(5000, 5.0e-3, 0., 0.)
	if err != nil {
		return fmt.Errorf("initialization failed: %w", err)
	}

	phi, iter, errSolve := PS1D.IteratingStepMethod()
	fmt.Printf("Jacobi Method:\t iterations=%d, converged=%v\n", iter, errSolve == nil)
	if errSolve != nil {
		fmt.Printf("  Error: %v\n", errSolve)
	}

	// Test with omega (SOR method) - use optimal omega around 1.5-1.8
	err = PS1D.Initialize(5000, 5.0e-3, 0., 0.)
	if err != nil {
		return fmt.Errorf("re-initialization failed: %w", err)
	}

	phiOmg, iterOmg, errOmg := PS1D.IteratingStepWithOmega(1.5) // Changed from 2 to 1.5
	fmt.Printf("SOR Method (ω=1.5):\t iterations=%d, converged=%v\n", iterOmg, errOmg == nil)
	if errOmg != nil {
		fmt.Printf("  Error: %v\n", errOmg)
	}

	// Save results
	if phi != nil {
		err = PS1D.grid.PrintVectorToFileRe(phi, "phiWO_omega.dat", "%21.14e")
		if err != nil {
			return fmt.Errorf("failed to write Jacobi results: %w", err)
		}
	}

	if phiOmg != nil {
		err = PS1D.grid.PrintVectorToFileRe(phiOmg, "phi_with_omega.dat", "%21.14e")
		if err != nil {
			return fmt.Errorf("failed to write SOR results: %w", err)
		}
	}

	return nil
}

type PoissonSolver2D struct {
	gridX *gridData.RadGrid
	gridY *gridData.RadGrid
	fxy   func(x, y float64) float64

	maxIter int
	tol     float64

	nx, ny   int
	dx2, dy2 float64
	rho      [][]float64
	phiPre   [][]float64
	phiNew   [][]float64
}

func NewTwoDimPoissonSolver(gridX, gridY *gridData.RadGrid, fxy func(x, y float64) float64) (*PoissonSolver2D, error) {
	if gridX == nil || gridY == nil || fxy == nil {
		return nil, errors.New("grids or function is not defined")
	}

	return &PoissonSolver2D{
		gridX: gridX,
		gridY: gridY,
		fxy:   fxy,
	}, nil
}

func (PS2D *PoissonSolver2D) Redefine(gridX, gridY *gridData.RadGrid, fxy func(x, y float64) float64) error {
	if gridX == nil || gridY == nil {
		return errors.New("grids cannot be nil")
	}
	if fxy == nil {
		return errors.New("function cannot be nil")
	}

	PS2D.gridX = gridX
	PS2D.gridY = gridY
	PS2D.fxy = fxy
	return nil
}

func (PS2D *PoissonSolver2D) RedefineFunc(fxy func(x, y float64) float64) error {
	if fxy == nil {
		return errors.New("function cannot be nil")
	}
	PS2D.fxy = fxy
	return nil
}

func (PS2D *PoissonSolver2D) Initialize(maxIter int, tolerance float64, boundaryCondition func(i, j, nx, ny int) float64) error {
	if maxIter < 1 {
		return errors.New("maxIter must be greater than zero")
	}
	if tolerance <= 0 || tolerance > 1 {
		return errors.New("tolerance must be between 0 and 1")
	}

	PS2D.maxIter = maxIter
	PS2D.tol = tolerance
	PS2D.nx = int(PS2D.gridX.NPoints())
	PS2D.ny = int(PS2D.gridY.NPoints())
	PS2D.dx2 = PS2D.gridX.DeltaR() * PS2D.gridX.DeltaR()
	PS2D.dy2 = PS2D.gridY.DeltaR() * PS2D.gridY.DeltaR()

	// Allocate arrays
	PS2D.phiPre = make([][]float64, PS2D.nx)
	PS2D.phiNew = make([][]float64, PS2D.nx)
	PS2D.rho = make([][]float64, PS2D.nx)

	for i := 0; i < PS2D.nx; i++ {
		PS2D.phiPre[i] = make([]float64, PS2D.ny)
		PS2D.phiNew[i] = make([]float64, PS2D.ny)
		PS2D.rho[i] = make([]float64, PS2D.ny)
	}

	// Set boundary conditions
	for i := 0; i < PS2D.nx; i++ {
		for j := 0; j < PS2D.ny; j++ {
			if i == 0 || i == PS2D.nx-1 || j == 0 || j == PS2D.ny-1 {
				PS2D.phiPre[i][j] = boundaryCondition(i, j, PS2D.nx, PS2D.ny)
				PS2D.phiNew[i][j] = PS2D.phiPre[i][j]
			}
		}
	}

	// Compute source term on grid
	for i := 0; i < PS2D.nx; i++ {
		for j := 0; j < PS2D.ny; j++ {
			x := PS2D.gridX.RMin() + float64(i)*PS2D.gridX.DeltaR()
			y := PS2D.gridY.RMin() + float64(j)*PS2D.gridY.DeltaR()
			PS2D.rho[i][j] = PS2D.fxy(x, y)
		}
	}

	return nil
}

// IteratingStepMethod Jacobi iteration method for 2D
func (PS2D *PoissonSolver2D) IteratingStepMethod() ([][]float64, int, error) {
	invDenom := 1.0 / (2.0/PS2D.dx2 + 2.0/PS2D.dy2)

	for iter := 0; iter < PS2D.maxIter; iter++ {
		errL2Norm := 0.0

		// Update interior points
		for i := 1; i < PS2D.nx-1; i++ {
			for j := 1; j < PS2D.ny-1; j++ {
				PS2D.phiNew[i][j] = invDenom * ((PS2D.phiPre[i+1][j]+PS2D.phiPre[i-1][j])/PS2D.dx2 +
					(PS2D.phiPre[i][j+1]+PS2D.phiPre[i][j-1])/PS2D.dy2 -
					PS2D.rho[i][j])

				diff := PS2D.phiNew[i][j] - PS2D.phiPre[i][j]
				errL2Norm += diff * diff
			}
		}

		// Copy new to previous
		for i := 1; i < PS2D.nx-1; i++ {
			copy(PS2D.phiPre[i], PS2D.phiNew[i])
		}

		// Check convergence
		if math.Sqrt(errL2Norm) < PS2D.tol {
			return PS2D.phiNew, iter + 1, nil
		}
	}

	return nil, PS2D.maxIter, fmt.Errorf("did not converge in %d iterations", PS2D.maxIter)
}

// IteratingStepWithOmega SOR method for 2D
func (PS2D *PoissonSolver2D) IteratingStepWithOmega(omega float64) ([][]float64, int, error) {
	if omega <= 0 || omega >= 2 {
		return nil, 0, errors.New("omega must be in range (0, 2) for stability")
	}

	invDenom := 1.0 / (2.0/PS2D.dx2 + 2.0/PS2D.dy2)

	for iter := 0; iter < PS2D.maxIter; iter++ {
		errL2Norm := 0.0

		// SOR update with Red-Black ordering for better convergence
		for i := 1; i < PS2D.nx-1; i++ {
			for j := 1; j < PS2D.ny-1; j++ {
				phiJacobi := invDenom * ((PS2D.phiNew[i+1][j]+PS2D.phiNew[i-1][j])/PS2D.dx2 +
					(PS2D.phiNew[i][j+1]+PS2D.phiNew[i][j-1])/PS2D.dy2 -
					PS2D.rho[i][j])

				phiOld := PS2D.phiNew[i][j]
				PS2D.phiNew[i][j] = (1.0-omega)*phiOld + omega*phiJacobi

				diff := PS2D.phiNew[i][j] - phiOld
				errL2Norm += diff * diff
			}
		}

		// Check convergence
		if math.Sqrt(errL2Norm) < PS2D.tol {
			return PS2D.phiNew, iter + 1, nil
		}
	}

	return nil, PS2D.maxIter, fmt.Errorf("did not converge in %d iterations", PS2D.maxIter)
}

// GaussSeidelRedBlack Gauss-Seidel method (Red-Black ordering)
func (PS2D *PoissonSolver2D) GaussSeidelRedBlack() ([][]float64, int, error) {
	invDenom := 1.0 / (2.0/PS2D.dx2 + 2.0/PS2D.dy2)

	for iter := 0; iter < PS2D.maxIter; iter++ {
		errL2Norm := 0.0

		// Red points (i+j even)
		for i := 1; i < PS2D.nx-1; i++ {
			for j := 1; j < PS2D.ny-1; j++ {
				if (i+j)%2 == 0 {
					phiOld := PS2D.phiNew[i][j]
					PS2D.phiNew[i][j] = invDenom * ((PS2D.phiNew[i+1][j]+PS2D.phiNew[i-1][j])/PS2D.dx2 +
						(PS2D.phiNew[i][j+1]+PS2D.phiNew[i][j-1])/PS2D.dy2 -
						PS2D.rho[i][j])

					diff := PS2D.phiNew[i][j] - phiOld
					errL2Norm += diff * diff
				}
			}
		}

		// Black points (i+j odd)
		for i := 1; i < PS2D.nx-1; i++ {
			for j := 1; j < PS2D.ny-1; j++ {
				if (i+j)%2 == 1 {
					phiOld := PS2D.phiNew[i][j]
					PS2D.phiNew[i][j] = invDenom * ((PS2D.phiNew[i+1][j]+PS2D.phiNew[i-1][j])/PS2D.dx2 +
						(PS2D.phiNew[i][j+1]+PS2D.phiNew[i][j-1])/PS2D.dy2 -
						PS2D.rho[i][j])

					diff := PS2D.phiNew[i][j] - phiOld
					errL2Norm += diff * diff
				}
			}
		}

		// Check convergence
		if math.Sqrt(errL2Norm) < PS2D.tol {
			return PS2D.phiNew, iter + 1, nil
		}
	}

	return nil, PS2D.maxIter, fmt.Errorf("did not converge in %d iterations", PS2D.maxIter)
}

func (PS2D *PoissonSolver2D) SaveToFile(phi [][]float64, filename string) error {
	// Implementation depends on your file writing utilities
	// This is a placeholder
	return nil
}

func (PS2D *PoissonSolver2D) CheckError() error {
	// Zero boundary conditions
	zeroBoundary := func(i, j, nx, ny int) float64 {
		return 0.0
	}

	// Test Jacobi
	err := PS2D.Initialize(10000, 1.0e-4, zeroBoundary)
	if err != nil {
		return fmt.Errorf("initialization failed: %w", err)
	}

	phi, iter, errSolve := PS2D.IteratingStepMethod()
	fmt.Printf("2D Jacobi Method:\t iterations=%d, converged=%v\n", iter, errSolve == nil)

	// Test SOR
	err = PS2D.Initialize(10000, 1.0e-4, zeroBoundary)
	if err != nil {
		return fmt.Errorf("re-initialization failed: %w", err)
	}

	phiOmg, iterOmg, errOmg := PS2D.IteratingStepWithOmega(1.7)
	fmt.Printf("2D SOR Method (ω=1.7):\t iterations=%d, converged=%v\n", iterOmg, errOmg == nil)

	// Test Gauss-Seidel
	err = PS2D.Initialize(10000, 1.0e-4, zeroBoundary)
	if err != nil {
		return fmt.Errorf("re-initialization failed: %w", err)
	}

	phiGS, iterGS, errGS := PS2D.GaussSeidelRedBlack()
	fmt.Printf("2D Gauss-Seidel:\t iterations=%d, converged=%v\n", iterGS, errGS == nil)

	// Save results if needed
	if phi != nil {
		err := PS2D.SaveToFile(phi, "phi2D_jacobi.dat")
		if err != nil {
			return err
		}
	}
	if phiOmg != nil {
		err := PS2D.SaveToFile(phiOmg, "phi2D_sor.dat")
		if err != nil {
			return err
		}
	}
	if phiGS != nil {
		err := PS2D.SaveToFile(phiGS, "phi2D_gs.dat")
		if err != nil {
			return err
		}
	}

	return nil
}

type PoissonSolver3D struct {
	gridX *gridData.RadGrid
	gridY *gridData.RadGrid
	gridZ *gridData.RadGrid
	fxyz  func(x, y, z float64) float64

	maxIter int
	tol     float64

	nx, ny, nz    int
	dx2, dy2, dz2 float64
	rho           [][][]float64
	phiPre        [][][]float64
	phiNew        [][][]float64
}

func NewThreeDimPoissonSolver(gridX, gridY, gridZ *gridData.RadGrid, fxyz func(x, y, z float64) float64) (*PoissonSolver3D, error) {
	if gridX == nil || gridY == nil || gridZ == nil || fxyz == nil {
		return nil, errors.New("grids or function is not defined")
	}

	return &PoissonSolver3D{
		gridX: gridX,
		gridY: gridY,
		gridZ: gridZ,
		fxyz:  fxyz,
	}, nil
}

func (PS3D *PoissonSolver3D) Redefine(gridX, gridY, gridZ *gridData.RadGrid, fxyz func(x, y, z float64) float64) error {
	if gridX == nil || gridY == nil || gridZ == nil {
		return errors.New("grids cannot be nil")
	}
	if fxyz == nil {
		return errors.New("function cannot be nil")
	}

	PS3D.gridX = gridX
	PS3D.gridY = gridY
	PS3D.gridZ = gridZ
	PS3D.fxyz = fxyz
	return nil
}

func (PS3D *PoissonSolver3D) RedefineFunc(fxyz func(x, y, z float64) float64) error {
	if fxyz == nil {
		return errors.New("function cannot be nil")
	}
	PS3D.fxyz = fxyz
	return nil
}

func (PS3D *PoissonSolver3D) Initialize(maxIter int, tolerance float64, boundaryCondition func(i, j, k, nx, ny, nz int) float64) error {
	if maxIter < 1 {
		return errors.New("maxIter must be greater than zero")
	}
	if tolerance <= 0 || tolerance > 1 {
		return errors.New("tolerance must be between 0 and 1")
	}

	PS3D.maxIter = maxIter
	PS3D.tol = tolerance
	PS3D.nx = int(PS3D.gridX.NPoints())
	PS3D.ny = int(PS3D.gridY.NPoints())
	PS3D.nz = int(PS3D.gridZ.NPoints())
	PS3D.dx2 = PS3D.gridX.DeltaR() * PS3D.gridX.DeltaR()
	PS3D.dy2 = PS3D.gridY.DeltaR() * PS3D.gridY.DeltaR()
	PS3D.dz2 = PS3D.gridZ.DeltaR() * PS3D.gridZ.DeltaR()

	// Allocate 3D arrays
	PS3D.phiPre = make([][][]float64, PS3D.nx)
	PS3D.phiNew = make([][][]float64, PS3D.nx)
	PS3D.rho = make([][][]float64, PS3D.nx)

	for i := 0; i < PS3D.nx; i++ {
		PS3D.phiPre[i] = make([][]float64, PS3D.ny)
		PS3D.phiNew[i] = make([][]float64, PS3D.ny)
		PS3D.rho[i] = make([][]float64, PS3D.ny)

		for j := 0; j < PS3D.ny; j++ {
			PS3D.phiPre[i][j] = make([]float64, PS3D.nz)
			PS3D.phiNew[i][j] = make([]float64, PS3D.nz)
			PS3D.rho[i][j] = make([]float64, PS3D.nz)
		}
	}

	// Set boundary conditions
	for i := 0; i < PS3D.nx; i++ {
		for j := 0; j < PS3D.ny; j++ {
			for k := 0; k < PS3D.nz; k++ {
				if i == 0 || i == PS3D.nx-1 ||
					j == 0 || j == PS3D.ny-1 ||
					k == 0 || k == PS3D.nz-1 {
					PS3D.phiPre[i][j][k] = boundaryCondition(i, j, k, PS3D.nx, PS3D.ny, PS3D.nz)
					PS3D.phiNew[i][j][k] = PS3D.phiPre[i][j][k]
				}
			}
		}
	}

	// Compute source term on grid
	for i := 0; i < PS3D.nx; i++ {
		for j := 0; j < PS3D.ny; j++ {
			for k := 0; k < PS3D.nz; k++ {
				x := PS3D.gridX.RMin() + float64(i)*PS3D.gridX.DeltaR()
				y := PS3D.gridY.RMin() + float64(j)*PS3D.gridY.DeltaR()
				z := PS3D.gridZ.RMin() + float64(k)*PS3D.gridZ.DeltaR()
				PS3D.rho[i][j][k] = PS3D.fxyz(x, y, z)
			}
		}
	}

	return nil
}

// IteratingStepMethod Jacobi iteration method for 3D
func (PS3D *PoissonSolver3D) IteratingStepMethod() ([][][]float64, int, error) {
	invDenom := 1.0 / (2.0/PS3D.dx2 + 2.0/PS3D.dy2 + 2.0/PS3D.dz2)

	for iter := 0; iter < PS3D.maxIter; iter++ {
		errL2Norm := 0.0

		// Update interior points
		for i := 1; i < PS3D.nx-1; i++ {
			for j := 1; j < PS3D.ny-1; j++ {
				for k := 1; k < PS3D.nz-1; k++ {
					PS3D.phiNew[i][j][k] = invDenom * ((PS3D.phiPre[i+1][j][k]+PS3D.phiPre[i-1][j][k])/PS3D.dx2 +
						(PS3D.phiPre[i][j+1][k]+PS3D.phiPre[i][j-1][k])/PS3D.dy2 +
						(PS3D.phiPre[i][j][k+1]+PS3D.phiPre[i][j][k-1])/PS3D.dz2 -
						PS3D.rho[i][j][k])

					diff := PS3D.phiNew[i][j][k] - PS3D.phiPre[i][j][k]
					errL2Norm += diff * diff
				}
			}
		}

		// Copy new to previous
		for i := 1; i < PS3D.nx-1; i++ {
			for j := 1; j < PS3D.ny-1; j++ {
				copy(PS3D.phiPre[i][j], PS3D.phiNew[i][j])
			}
		}

		// Check convergence
		if math.Sqrt(errL2Norm) < PS3D.tol {
			return PS3D.phiNew, iter + 1, nil
		}
	}

	return nil, PS3D.maxIter, fmt.Errorf("did not converge in %d iterations", PS3D.maxIter)
}

// IteratingStepWithOmega SOR method for 3D
func (PS3D *PoissonSolver3D) IteratingStepWithOmega(omega float64) ([][][]float64, int, error) {
	if omega <= 0 || omega >= 2 {
		return nil, 0, errors.New("omega must be in range (0, 2) for stability")
	}

	invDenom := 1.0 / (2.0/PS3D.dx2 + 2.0/PS3D.dy2 + 2.0/PS3D.dz2)

	for iter := 0; iter < PS3D.maxIter; iter++ {
		errL2Norm := 0.0

		// SOR update
		for i := 1; i < PS3D.nx-1; i++ {
			for j := 1; j < PS3D.ny-1; j++ {
				for k := 1; k < PS3D.nz-1; k++ {
					phiJacobi := invDenom * ((PS3D.phiNew[i+1][j][k]+PS3D.phiNew[i-1][j][k])/PS3D.dx2 +
						(PS3D.phiNew[i][j+1][k]+PS3D.phiNew[i][j-1][k])/PS3D.dy2 +
						(PS3D.phiNew[i][j][k+1]+PS3D.phiNew[i][j][k-1])/PS3D.dz2 -
						PS3D.rho[i][j][k])

					phiOld := PS3D.phiNew[i][j][k]
					PS3D.phiNew[i][j][k] = (1.0-omega)*phiOld + omega*phiJacobi

					diff := PS3D.phiNew[i][j][k] - phiOld
					errL2Norm += diff * diff
				}
			}
		}

		// Check convergence
		if math.Sqrt(errL2Norm) < PS3D.tol {
			return PS3D.phiNew, iter + 1, nil
		}
	}

	return nil, PS3D.maxIter, fmt.Errorf("did not converge in %d iterations", PS3D.maxIter)
}

// GaussSeidelRedBlack Gauss-Seidel with Red-Black ordering for 3D
func (PS3D *PoissonSolver3D) GaussSeidelRedBlack() ([][][]float64, int, error) {
	invDenom := 1.0 / (2.0/PS3D.dx2 + 2.0/PS3D.dy2 + 2.0/PS3D.dz2)

	for iter := 0; iter < PS3D.maxIter; iter++ {
		errL2Norm := 0.0

		// Red-Black ordering: red points have (i+j+k) even, black have (i+j+k) odd
		for color := 0; color < 2; color++ {
			for i := 1; i < PS3D.nx-1; i++ {
				for j := 1; j < PS3D.ny-1; j++ {
					for k := 1; k < PS3D.nz-1; k++ {
						if (i+j+k)%2 == color {
							phiOld := PS3D.phiNew[i][j][k]
							PS3D.phiNew[i][j][k] = invDenom * ((PS3D.phiNew[i+1][j][k]+PS3D.phiNew[i-1][j][k])/PS3D.dx2 +
								(PS3D.phiNew[i][j+1][k]+PS3D.phiNew[i][j-1][k])/PS3D.dy2 +
								(PS3D.phiNew[i][j][k+1]+PS3D.phiNew[i][j][k-1])/PS3D.dz2 -
								PS3D.rho[i][j][k])

							diff := PS3D.phiNew[i][j][k] - phiOld
							errL2Norm += diff * diff
						}
					}
				}
			}
		}

		// Check convergence
		if math.Sqrt(errL2Norm) < PS3D.tol {
			return PS3D.phiNew, iter + 1, nil
		}
	}

	return nil, PS3D.maxIter, fmt.Errorf("did not converge in %d iterations", PS3D.maxIter)
}

// MultigridVCycle Multigrid V-cycle (basic implementation)
func (PS3D *PoissonSolver3D) MultigridVCycle(levels int) ([][][]float64, int, error) {
	// This is a placeholder for a more advanced multigrid implementation
	// Multigrid is the most efficient method for 3D Poisson problems
	// but requires careful implementation of restriction and prolongation operators
	return nil, 0, errors.New("multigrid not yet implemented")
}

func (PS3D *PoissonSolver3D) SaveToFile(phi [][][]float64, filename string) error {
	// Implementation depends on your file writing utilities
	// This is a placeholder
	return nil
}

func (PS3D *PoissonSolver3D) CheckError() error {
	// Zero boundary conditions
	zeroBoundary := func(i, j, k, nx, ny, nz int) float64 {
		return 0.0
	}

	// Test Jacobi
	err := PS3D.Initialize(5000, 1.0e-3, zeroBoundary)
	if err != nil {
		return fmt.Errorf("initialization failed: %w", err)
	}

	phi, iter, errSolve := PS3D.IteratingStepMethod()
	fmt.Printf("3D Jacobi Method:\t iterations=%d, converged=%v\n", iter, errSolve == nil)

	// Test SOR
	err = PS3D.Initialize(5000, 1.0e-3, zeroBoundary)
	if err != nil {
		return fmt.Errorf("re-initialization failed: %w", err)
	}

	phiOmg, iterOmg, errOmg := PS3D.IteratingStepWithOmega(1.6)
	fmt.Printf("3D SOR Method (ω=1.6):\t iterations=%d, converged=%v\n", iterOmg, errOmg == nil)

	// Test Gauss-Seidel
	err = PS3D.Initialize(5000, 1.0e-3, zeroBoundary)
	if err != nil {
		return fmt.Errorf("re-initialization failed: %w", err)
	}

	phiGS, iterGS, errGS := PS3D.GaussSeidelRedBlack()
	fmt.Printf("3D Gauss-Seidel:\t iterations=%d, converged=%v\n", iterGS, errGS == nil)

	// Save results if needed
	if phi != nil {
		err := PS3D.SaveToFile(phi, "phi3D_jacobi.dat")
		if err != nil {
			return err
		}
	}
	if phiOmg != nil {
		err := PS3D.SaveToFile(phiOmg, "phi3D_sor.dat")
		if err != nil {
			return err
		}
	}
	if phiGS != nil {
		err := PS3D.SaveToFile(phiGS, "phi3D_gs.dat")
		if err != nil {
			return err
		}
	}

	return nil
}

// ComputeResidual Helper function to compute residual (useful for debugging)
func (PS3D *PoissonSolver3D) ComputeResidual(phi [][][]float64) float64 {
	residual := 0.0
	invDenom := 1.0 / (2.0/PS3D.dx2 + 2.0/PS3D.dy2 + 2.0/PS3D.dz2)

	for i := 1; i < PS3D.nx-1; i++ {
		for j := 1; j < PS3D.ny-1; j++ {
			for k := 1; k < PS3D.nz-1; k++ {
				expected := invDenom * ((phi[i+1][j][k]+phi[i-1][j][k])/PS3D.dx2 +
					(phi[i][j+1][k]+phi[i][j-1][k])/PS3D.dy2 +
					(phi[i][j][k+1]+phi[i][j][k-1])/PS3D.dz2 -
					PS3D.rho[i][j][k])

				diff := phi[i][j][k] - expected
				residual += diff * diff
			}
		}
	}

	return math.Sqrt(residual)
}
