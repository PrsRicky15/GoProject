package gridData

type Integrator interface {
	Step(s *ModelSystem) SystemInfo
}
