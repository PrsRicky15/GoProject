package gridData

// compositeFunc multiplication of two function
type compositeFunc[T VarType] struct {
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

// CompositeF64 For float64 specialization
type CompositeF64 compositeFunc[float64]

// CompositeZ64 For complex64 specialization
type CompositeZ64 compositeFunc[complex128]

func (cF compositeFunc[T]) EvaluateAt(x T) T {
	return cF.func1.EvaluateAt(x) * cF.func2.EvaluateAt(x)
}

func (cF compositeFunc[T]) ForceAt(x T) T {
	return cF.func1.ForceAt(x) * cF.func2.ForceAt(x)
}

func (cF compositeFunc[T]) EvaluateOnGrid(x []T) []T {
	return compositeOnGrid(cF.func1.EvaluateAt, cF.func2.EvaluateAt, x)
}

func (cF compositeFunc[T]) ForceOnGrid(x []T) []T {
	return compositeOnGrid(cF.func1.ForceAt, cF.func2.ForceAt, x)
}
