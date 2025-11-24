package gridData

type TDfunc interface {
	EvaluateAt(x float64, t float64) float64
}

// TDPotentialOp General interface for the evaluating the potential on a grid
type TDPotentialOp interface {
	TDfunc
	EvaluateOnRGrid(x []float64, t float64) []float64
	EvaluateOnRGridInPlace(x, res []float64, t float64)
}

// TDPotential Idea: compose - PotentialOp[float64]*f(x,t)
type TDPotential struct {
	Static PotentialOp[float64]
	TDFunc func(x, t float64) float64
}

// NewTimeDependentPotential Constructor
func NewTimeDependentPotential(static PotentialOp[float64],
	tdTerm func(float64, float64) float64) *TDPotential {
	return &TDPotential{
		Static: static,
		TDFunc: tdTerm,
	}
}

// EvaluateAt returns V_static(x) + V_time(x, t)
func (tp *TDPotential) EvaluateAt(x float64, t float64) float64 {
	var vstat float64
	if tp.Static != nil {
		vstat = tp.Static.EvaluateAt(x)
	}
	if tp.TDFunc == nil {
		return vstat
	}
	return vstat + tp.TDFunc(x, t)
}

func (tp *TDPotential) EvaluateOnRGrid(x []float64, t float64) []float64 {
	res := make([]float64, len(x))
	tp.EvaluateOnRGridInPlace(x, res, t)
	return res
}

// EvaluateOnRGridInPlace res[i] = V_static(x[i]) + V_td(x[i], t)
func (tp *TDPotential) EvaluateOnRGridInPlace(x, res []float64, t float64) {
	if len(x) != len(res) {
		panic("EvaluateOnRGridInPlace: x and res must have equal length")
	}

	if tp.Static != nil {
		for i, xi := range x {
			res[i] = tp.Static.EvaluateAt(xi)
		}
	} else {
		for i := range res {
			res[i] = 0
		}
	}

	if tp.TDFunc != nil {
		for i, xi := range x {
			res[i] += tp.TDFunc(xi, t)
		}
	}
}
