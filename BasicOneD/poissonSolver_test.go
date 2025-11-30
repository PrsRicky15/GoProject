package BasicOneD

import (
	"GoProject/gridData"
	"testing"
)

func TestNewOneDimPoissonSolver(t *testing.T) {
	tests := []struct {
		name    string
		grid    *gridData.RadGrid
		fx      gridData.Rfunc
		wantErr bool
		errMsg  string
	}{
		{
			name:    "valid initialization",
			grid:    &gridData.RadGrid{},
			fx:      gridData.Gaussian[float64]{Cen: 1, Sigma: 2},
			wantErr: false,
		},
		{
			name:    "nil grid",
			grid:    nil,
			fx:      gridData.Gaussian[float64]{Cen: 1, Sigma: 1},
			wantErr: true,
			errMsg:  "grid cannot be nil",
		},
		{
			name:    "nil function",
			grid:    &gridData.RadGrid{},
			fx:      nil,
			wantErr: true,
			errMsg:  "function cannot be nil",
		},
		{
			name:    "both nil",
			grid:    nil,
			fx:      nil,
			wantErr: true,
			errMsg:  "grid cannot be nil",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			solver, err := NewOneDimPoissonSolver(tt.grid, tt.fx)

			if tt.wantErr {
				if err == nil {
					t.Errorf("expected error but got none")
				}
				if err != nil && err.Error() != tt.errMsg {
					t.Errorf("expected error '%s', got '%s'", tt.errMsg, err.Error())
				}
				if solver != nil {
					t.Errorf("expected nil solver on error")
				}
			} else {
				if err != nil {
					t.Errorf("unexpected error: %v", err)
				}
				if solver == nil {
					t.Errorf("expected valid solver")
				}
				if solver.grid != tt.grid {
					t.Errorf("grid not set correctly")
				}
				if solver.fx == nil {
					t.Errorf("function not set correctly")
				}
			}
		})
	}
}
