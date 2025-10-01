package gridData

import (
	"fmt"
	"math"
	"os"
)

// setGrid Setup some parameters for time- and space- grid
type setGrid struct {
	length  float64
	deltaS  float64
	deltaCS float64
	cMin    float64
	cMax    float64
}

// getGridData get the basic input and setup time- and space- grid
type getGridData interface {
	getMin() float64
	getdS() float64
	getNgrid() uint32
	getdCS() float64
}

// GridFunctionality Functionality used in both R-Grid and T-Grid
type GridFunctionality interface {
	PotentialAt(pot PotentialOp, x float64) float64
	DisplayInfo()
	PrintPotentToFile(Pot PotentialOp, filename string, format string) error
}

/*
 Some Important function which defines R-Grid and T-grid properly without duplication
*/

// displayFunc Display function on a grid
func displayFunc(g getGridData, Pot PotentialOp, format string,
	f func(evaluate PotentialOp, x float64) float64) {
	fullFormat := "%14.7e" + "\t" + format + "\n"
	fmt.Printf("#--------------------------------------------------\n")
	fmt.Printf("#\t\t grid\t\t function value\n")
	fmt.Printf("#--------------------------------------------------\n")
	for i := uint32(0); i < g.getNgrid(); i++ {
		var x = g.getMin() + float64(i)*g.getdS()
		fmt.Printf(fullFormat, x, f(Pot, x))
	}
}

// functionToFile print function on a grid to a File
func functionToFile(g getGridData, Pot PotentialOp, filename string, format string,
	f func(evaluate PotentialOp, x float64) float64) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer func(file *os.File) {
		err := file.Close()
		if err != nil {

		}
	}(file)

	fullFormat := "%14.7e" + "\t" + format + "\n"
	_, err = fmt.Fprintf(file, "#--------------------------------------------------\n")
	_, err = fmt.Fprintf(file, "#\t\t grid\t\t function value\n")
	_, err = fmt.Fprintf(file, "#--------------------------------------------------\n")
	for i := uint32(0); i < g.getNgrid(); i++ {
		var x = g.getMin() + float64(i)*g.getdS()
		_, err := fmt.Fprintf(file, fullFormat, x, f(Pot, x))
		if err != nil {
			return err
		}
	}
	return nil
}

func vectorToFile(g getGridData, vec []float64, filename string, format string) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer func(file *os.File) {
		err := file.Close()
		if err != nil {
			panic(err)
		}
	}(file)

	fullFormat := "%14.7e" + "\t" + format + "\n"
	_, err = fmt.Fprintf(file, "#--------------------------------------------------\n")
	_, err = fmt.Fprintf(file, "#\t\t grid\t\t vectorPoints\n")
	_, err = fmt.Fprintf(file, "#--------------------------------------------------\n")
	for i := uint32(0); i < g.getNgrid(); i++ {
		x := g.getMin() + float64(i)*g.getdS()
		if _, err := fmt.Fprintf(file, fullFormat, x, vec[i]); err != nil {
			return err
		}
	}
	return nil
}

func displayGrid(space func() []float64) {
	points := space()
	for i, val := range points {
		fmt.Printf("r-%d %14.7e\n", i, val)
	}
}

func generatePoints(g getGridData) []float64 {
	values := make([]float64, g.getNgrid())
	for i := uint32(0); i < g.getNgrid(); i++ {
		values[i] = g.getMin() + float64(i)*g.getdS()
	}
	return values
}

func generateConjugatePoints(g getGridData) []float64 {
	nby2 := g.getNgrid() / 2
	values := make([]float64, g.getNgrid())
	values[0] = 0.
	values[g.getNgrid()/2] = -float64(nby2) * g.getdCS()
	for i := uint32(1); i < nby2; i++ {
		values[i] = -float64(i) * g.getdCS()
		values[i+nby2] = float64(nby2-i) * g.getdCS()
	}
	return values
}

func createGrid(min, max float64, nPoints uint32, paramName string) (setGrid, error) {
	if nPoints == 0 {
		return setGrid{0., 0., 0., 0., 0.},
			fmt.Errorf("number of grid points must be positive")
	}
	if max <= min {
		return setGrid{0., 0., 0., 0., 0.},
			fmt.Errorf("%sMax (%g) must be greater than %sMin (%g)", paramName, max, paramName, min)
	}
	dSpace := (max - min) / float64(nPoints)
	length := max - min
	dConjugate := 2 * math.Pi / length
	cMin := -math.Pi / dSpace
	cMax := math.Pi / dSpace
	return setGrid{length, dSpace, dConjugate, cMin, cMax}, nil
}
