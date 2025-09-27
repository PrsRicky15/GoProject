package gridData

type InitialValueIntegrator interface {
	Step(s *ModelSystem) SystemInfo
}
