package Plots

import (
	"GoProject/gridData"
	"log"
	"slices"
)

// ==================== PLOT GENERATORS ====================

// generateMorsePotentialSurfaceData generates Morse potential plot data
func generateMorsePotentialSurfaceData(gridParams map[string]interface{}, params map[string]interface{}) PlotDataResponse {
	// Extract parameters
	gp := extractGridParams(gridParams)
	mp := extractMorseParams(params)

	// Create grid
	grid, err := gridData.NewRGrid(gp.RMin, gp.RMax, gp.NGrid)
	if err != nil {
		log.Printf("Error creating grid: %v", err)
		return PlotDataResponse{Status: "error"}
	}

	// Create Morse potential
	morse := gridData.Morse[float64]{
		De:    mp.D,
		Alpha: mp.A,
		Cen:   mp.R0,
	}

	// Calculate potential on grid
	x := grid.RValues()
	y := grid.PotentialOnGrid(morse)

	// Create plot elements
	trace := createTrace(x, y, "Potential Energy")
	layout := createBaseLayout("Morse Potential Energy Surface")
	setYAxisRange(layout, -10, mp.D+100)

	return PlotDataResponse{
		Data:   []map[string]interface{}{trace},
		Layout: layout,
		Status: statusCompleted,
	}
}

// generateSoftcorePotentialSurfaceData generates Softcore potential plot data
func generateSoftcorePotentialSurfaceData(gridParams map[string]interface{}, params map[string]interface{}) PlotDataResponse {
	// Extract parameters
	gp := extractGridParams(gridParams)
	sp := extractSoftcoreParams(params)

	// Create grid
	grid, err := gridData.NewRGrid(gp.RMin, gp.RMax, gp.NGrid)
	if err != nil {
		log.Printf("Error creating grid: %v", err)
		return PlotDataResponse{Status: "error"}
	}

	// Create Softcore potential
	sc := gridData.SoftCore[float64]{
		Charge:    sp.Charge,
		Centre:    sp.R0,
		SoftParam: sp.A,
	}

	// Calculate potential on grid
	x := grid.RValues()
	y := grid.PotentialOnGrid(sc)

	// Calculate y-axis range
	ymin := slices.Min(y)
	ymax := slices.Max(y)

	// Create plot elements
	trace := createTrace(x, y, "Potential Energy")
	layout := createBaseLayout("Softcore Potential Energy Surface")
	setYAxisRange(layout, ymin-1, ymax+1)

	return PlotDataResponse{
		Data:   []map[string]interface{}{trace},
		Layout: layout,
		Status: statusCompleted,
	}
}
