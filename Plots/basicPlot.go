package Plots

import (
	"encoding/json"
	"log"
	"net/http"
	"slices"

	"GoProject/gridData"
)

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

// Generate 1D potential energy surface data
func generateMorsePotentialSurfaceData(gridParams map[string]interface{}, params map[string]interface{}) PlotDataResponse {

	rMin := 0.
	rMax := 15.
	nGrid := uint32(20)
	if val, ok := gridParams["rMin"].(float64); ok {
		rMin = val
	}
	if val, ok := gridParams["a"].(float64); ok {
		rMax = val
	}
	if val, ok := gridParams["r0"].(float64); ok {
		nGrid = uint32(val)
	}

	grid, err := gridData.NewRGrid(rMin, rMax, nGrid)
	if err != nil {
		log.Fatal(err)
	}
	// Default parameters
	D := 100.0
	a := 1.5
	r0 := 2.0

	if val, ok := params["D"].(float64); ok {
		D = val
	}
	if val, ok := params["a"].(float64); ok {
		a = val
	}
	if val, ok := params["r0"].(float64); ok {
		r0 = val
	}

	morse := gridData.Morse[float64]{De: D, Alpha: a, Cen: r0}

	x := grid.RValues()
	y := grid.PotentialOnGrid(morse)

	trace := map[string]interface{}{
		"x":    x,
		"y":    y,
		"type": "scatter",
		"mode": "lines",
		"name": "Potential Energy",
		"line": map[string]interface{}{
			"color": "rgb(147, 51, 234)",
			"width": 3,
		},
	}

	layout := map[string]interface{}{
		"title": map[string]interface{}{
			"text": "Potential Energy Surface",
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
			"range":     []float64{-10, D + 100},
		},
		"plot_bgcolor":  "rgba(0,0,0,0)",
		"paper_bgcolor": "rgba(0,0,0,0)",
		"font": map[string]interface{}{
			"color": "white",
		},
	}

	return PlotDataResponse{
		Data:   []map[string]interface{}{trace},
		Layout: layout,
		Status: "completed",
	}
}

// Generate 1D potential energy surface data
func generateSoftcorePotentialSurfaceData(gridParams map[string]interface{}, params map[string]interface{}) PlotDataResponse {
	grid, _ := gridData.NewRGrid(gridParams["rMin"].(float64), gridParams["rMax"].(float64), gridParams["nGrid"].(uint32))
	// Default parameters
	charge := 1.0
	a := 1.5
	r0 := 2.0

	if val, ok := params["Charge"].(float64); ok {
		charge = val
	}
	if val, ok := params["a"].(float64); ok {
		a = val
	}
	if val, ok := params["r0"].(float64); ok {
		r0 = val
	}

	sc := gridData.SoftCore[float64]{Charge: charge, Centre: r0, SoftParam: a}

	x := grid.RValues()
	y := grid.PotentialOnGrid(sc)

	ymin := slices.Min(y)
	ymax := slices.Max(y)

	trace := map[string]interface{}{
		"x":    x,
		"y":    y,
		"type": "scatter",
		"mode": "lines",
		"name": "Potential Energy",
		"line": map[string]interface{}{
			"color": "rgb(147, 51, 234)",
			"width": 3,
		},
	}

	layout := map[string]interface{}{
		"title": map[string]interface{}{
			"text": "Potential Energy Surface",
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
			"range":     []float64{ymin - 1, ymax + 1},
		},
		"plot_bgcolor":  "rgba(0,0,0,0)",
		"paper_bgcolor": "rgba(0,0,0,0)",
		"font": map[string]interface{}{
			"color": "white",
		},
	}

	return PlotDataResponse{
		Data:   []map[string]interface{}{trace},
		Layout: layout,
		Status: "completed",
	}
}

// GeneratePlotData Main handler
func GeneratePlotData(w http.ResponseWriter, r *http.Request) {
	var req PlotRequest
	err := json.NewDecoder(r.Body).Decode(&req)
	if err != nil {
		http.Error(w, err.Error(), http.StatusBadRequest)
	}

	var response PlotDataResponse

	switch req.PlotType {
	case "Morse":
		response = generateMorsePotentialSurfaceData(req.Grid, req.Parameters)
	case "Softcore":
		response = generateSoftcorePotentialSurfaceData(req.Grid, req.Parameters)
	default:
		http.Error(w, "Unknown plot type", http.StatusBadRequest)
	}

	w.Header().Set("Content-Type", "application/json")
	err = json.NewEncoder(w).Encode(response)
	if err != nil {
		log.Println(err)
	}
}
