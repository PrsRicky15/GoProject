package gridData

// CompositeFunc multiplication of two function
type CompositeFunc[T VarType] struct {
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

func (cF CompositeFunc[T]) EvaluateAt(x T) T {
	return cF.func1.EvaluateAt(x) * cF.func2.EvaluateAt(x)
}

func (cF CompositeFunc[T]) ForceAt(x T) T {
	return cF.func1.ForceAt(x) * cF.func2.ForceAt(x)
}

func (cF CompositeFunc[T]) EvaluateOnGrid(x []T) []T {
	return compositeOnGrid(cF.func1.EvaluateAt, cF.func2.EvaluateAt, x)
}

func (cF CompositeFunc[T]) ForceOnGrid(x []T) []T {
	return compositeOnGrid(cF.func1.ForceAt, cF.func2.ForceAt, x)
}

func (cF CompositeFunc[T]) EvaluateOnGridInPlace(fn, x []T) {
	compositeOnGridInPlace(cF.func1.EvaluateAt, cF.func2.EvaluateAt, fn, x)
}

func (cF CompositeFunc[T]) ForceOnGridInPlace(fn, x []T) {
	compositeOnGridInPlace(cF.func1.ForceAt, cF.func2.ForceAt, fn, x)
}
