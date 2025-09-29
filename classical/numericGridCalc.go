package classical

type InitialValueIntegrator interface {
	Step(s *ModelSystem) SystemInfo
}
