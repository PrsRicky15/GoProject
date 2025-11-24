package gridData

type TDfunc interface {
	EvaluateAtTime(x float64, t float64) float64
}

// TDPotentialOp General interface for the evaluating the potential on a grid
type TDPotentialOp interface {
	TDfunc
	EvaluateOnRGridTime(x []float64, t float64) []float64
	EvaluateOnRGridTimeInPlace(x, res []float64, t float64)
}

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
func (PF *ProductFunc) EvaluateAtTime(x, t float64) float64 {
	return PF.func1.EvaluateAt(x) * PF.func2.EvaluateAt(t)
}

type SumFunc struct {
	func1 Rfunc
	func2 Rfunc
}

func NewSumFunc(func1 Rfunc, func2 Rfunc) *ProductFunc {
	return &ProductFunc{func1, func2}
}

func (SF *SumFunc) Redefine(func1 Rfunc, func2 Rfunc) {
	SF.func1 = func1
	SF.func2 = func2
}

func (SF *SumFunc) EvaluateAt(x float64) float64 {
	return SF.func1.EvaluateAt(x) + SF.func2.EvaluateAt(x)
}

func (SF *SumFunc) EvaluateAtTime(x, t float64) float64 {
	return SF.func1.EvaluateAt(x) + SF.func2.EvaluateAt(t)
}

type FuncWithTDPert struct {
	tdFunc ProductFunc
	tiFunc Rfunc
}

func NewFuncWithTDPert(func1 ProductFunc, func2 Rfunc) *FuncWithTDPert {
	return &FuncWithTDPert{func1, func2}
}

func (SF *FuncWithTDPert) Redefine(func1 ProductFunc, func2 Rfunc) {
	SF.tdFunc = func1
	SF.tiFunc = func2
}

func (SF *FuncWithTDPert) EvaluateAtTime(x, t float64) float64 {
	return SF.tiFunc.EvaluateAt(x) + SF.tdFunc.EvaluateAtTime(x, t)
}

// CompositePotentialOp multiplication of two function
type CompositePotentialOp[T VarType] struct {
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

func (cF CompositePotentialOp[T]) EvaluateAt(x T) T {
	return cF.func1.EvaluateAt(x) * cF.func2.EvaluateAt(x)
}

func (cF CompositePotentialOp[T]) ForceAt(x T) T {
	return cF.func1.ForceAt(x) * cF.func2.ForceAt(x)
}

func (cF CompositePotentialOp[T]) EvaluateOnGrid(x []T) []T {
	return compositeOnGrid(cF.func1.EvaluateAt, cF.func2.EvaluateAt, x)
}

func (cF CompositePotentialOp[T]) ForceOnGrid(x []T) []T {
	return compositeOnGrid(cF.func1.ForceAt, cF.func2.ForceAt, x)
}

func (cF CompositePotentialOp[T]) EvaluateOnGridInPlace(fn, x []T) {
	compositeOnGridInPlace(cF.func1.EvaluateAt, cF.func2.EvaluateAt, fn, x)
}

func (cF CompositePotentialOp[T]) ForceOnGridInPlace(fn, x []T) {
	compositeOnGridInPlace(cF.func1.ForceAt, cF.func2.ForceAt, fn, x)
}
