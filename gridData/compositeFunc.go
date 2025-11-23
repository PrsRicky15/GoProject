package gridData

type ProductFunc struct {
	func1 Rfunc
	func2 Rfunc
}

func NewProductFunc(func1 Rfunc, func2 Rfunc) *ProductFunc {
	return &ProductFunc{func1, func2}
}

func (PF *ProductFunc) Redefine(func1 Rfunc, func2 Rfunc) {
	PF.func1 = func1
	PF.func2 = func2
}

func (PF *ProductFunc) EvaluateAt(x float64) float64 {
	return PF.func1.EvaluateAt(x) * PF.func2.EvaluateAt(x)
}

type SumOfFunc struct {
	func1 Rfunc
	func2 Rfunc
}

func (SF SumOfFunc) EvaluateAt(x float64) float64 {
	return SF.func1.EvaluateAt(x) * SF.func2.EvaluateAt(x)
}

// CompositePotential multiplication of two function
type CompositePotential[T VarType] struct {
	func1 PotentialOp[T]
	func2 PotentialOp[T]
}

// compositeOnGrid Made generic to work with VarType
func compositeOnGrid[T VarType](f func(T) T, f2 func(T) T, x []T) []T {
	results := make([]T, len(x))
	for i, val := range x {
		results[i] = f(val) * f2(val)
	}
	return results
}

// compositeOnGrid Made generic to work with VarType
func compositeOnGridInPlace[T VarType](f func(T) T, f2 func(T) T, fn, x []T) {
	for i, val := range x {
		fn[i] = f(val) * f2(val)
	}
}

func (cF CompositePotential[T]) EvaluateAt(x T) T {
	return cF.func1.EvaluateAt(x) * cF.func2.EvaluateAt(x)
}

func (cF CompositePotential[T]) ForceAt(x T) T {
	return cF.func1.ForceAt(x) * cF.func2.ForceAt(x)
}

func (cF CompositePotential[T]) EvaluateOnGrid(x []T) []T {
	return compositeOnGrid(cF.func1.EvaluateAt, cF.func2.EvaluateAt, x)
}

func (cF CompositePotential[T]) ForceOnGrid(x []T) []T {
	return compositeOnGrid(cF.func1.ForceAt, cF.func2.ForceAt, x)
}

func (cF CompositePotential[T]) EvaluateOnGridInPlace(fn, x []T) {
	compositeOnGridInPlace(cF.func1.EvaluateAt, cF.func2.EvaluateAt, fn, x)
}

func (cF CompositePotential[T]) ForceOnGridInPlace(fn, x []T) {
	compositeOnGridInPlace(cF.func1.ForceAt, cF.func2.ForceAt, fn, x)
}
