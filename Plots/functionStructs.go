package Plots

// ==================== TYPE DEFINITIONS ====================

type PlotRequest struct {
	Grid       map[string]interface{} `json:"grid"`
	PlotType   string                 `json:"plot_type"`
	Parameters map[string]interface{} `json:"parameters"`
}

type PlotDataResponse struct {
	Data   []map[string]interface{} `json:"data"`
	Layout map[string]interface{}   `json:"layout"`
	Status string                   `json:"status"`
}

// GridParams holds the grid configuration
type GridParams struct {
	RMin  float64
	RMax  float64
	NGrid uint32
}

// MorseParams holds Morse potential parameters
type MorseParams struct {
	D  float64 // Dissociation energy
	A  float64 // Width parameter
	R0 float64 // Equilibrium distance
}

// SoftcoreParams holds Softcore potential parameters
type SoftcoreParams struct {
	Charge float64
	A      float64
	R0     float64
}

// ==================== CONSTANTS ====================

const (
	defaultRMin  = 0.0
	defaultRMax  = 15.0
	defaultNGrid = 20

	defaultMorseD  = 100.0
	defaultMorseA  = 1.5
	defaultMorseR0 = 2.0

	defaultSoftcoreCharge = 1.0
	defaultSoftcoreA      = 1.5
	defaultSoftcoreR0     = 2.0

	plotTypeScatter = "scatter"
	plotModeLines   = "lines"
	lineColor       = "rgb(147, 51, 234)"
	lineWidth       = 3

	statusCompleted = "completed"
)

// ==================== PARAMETER EXTRACTION ====================

// extractGridParams extracts and validates grid parameters from the request
func extractGridParams(gridData map[string]interface{}) GridParams {
	params := GridParams{
		RMin:  defaultRMin,
		RMax:  defaultRMax,
		NGrid: defaultNGrid,
	}

	if val, ok := gridData["rMin"].(float64); ok {
		params.RMin = val
	}
	if val, ok := gridData["rMax"].(float64); ok {
		params.RMax = val
	}
	if val, ok := gridData["nGrid"].(float64); ok {
		params.NGrid = uint32(val)
	}

	return params
}

// ==================== LAYOUT BUILDERS ====================

// createBaseLayout creates the common layout configuration for all plots
func createBaseLayout(title string) map[string]interface{} {
	return map[string]interface{}{
		"title": map[string]interface{}{
			"text": title,
			"font": map[string]interface{}{
				"size":  24,
				"color": "white",
			},
		},
		"xaxis": map[string]interface{}{
			"title": map[string]interface{}{
				"text": "Distance (a.u.)",
				"font": map[string]interface{}{"color": "white"},
			},
			"gridcolor": "rgba(255,255,255,0.1)",
			"color":     "white",
		},
		"yaxis": map[string]interface{}{
			"title": map[string]interface{}{
				"text": "Energy (a.u.)",
				"font": map[string]interface{}{"color": "white"},
			},
			"gridcolor": "rgba(255,255,255,0.1)",
			"color":     "white",
		},
		"plot_bgcolor":  "rgba(0,0,0,0.1)",
		"paper_bgcolor": "rgba(0, 0, 0, 0.8)",
		"font": map[string]interface{}{
			"color": "white",
		},
	}
}

// setYAxisRange sets the y-axis range in the layout
func setYAxisRange(layout map[string]interface{}, min, max float64) {
	if yaxis, ok := layout["yaxis"].(map[string]interface{}); ok {
		yaxis["range"] = []float64{min, max}
	}
}

// ==================== TRACE BUILDERS ====================

// createTrace creates a plot trace with the given data
func createTrace(x, y []float64, name string) map[string]interface{} {
	return map[string]interface{}{
		"x":    x,
		"y":    y,
		"type": plotTypeScatter,
		"mode": plotModeLines,
		"name": name,
		"line": map[string]interface{}{
			"color": lineColor,
			"width": lineWidth,
		},
	}
}
