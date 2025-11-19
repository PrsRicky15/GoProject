package Plots

import (
    "encoding/json"
    "net/http"

    "GoProject/gridData"
)

type PlotRequest struct {
    PlotType   string                 `json:"plot_type"`
    Parameters map[string]interface{} `json:"parameters"`
}

type PlotDataResponse struct {
    Data   []map[string]interface{} `json:"data"`
    Layout map[string]interface{}   `json:"layout"`
    Status string                   `json:"status"`
}

// Generate 1D potential energy surface data
func generateMorsePotentialSurfaceData(params map[string]interface{}) PlotDataResponse {
    // Generate data points
    // numPoints := 200
    // x := make([]float64, numPoints)
    // y := make([]float64, numPoints)

    // Default grid for plots
    grid, _ := gridData.NewRGrid(0, 50, 200)

    // Get parameters with defaults
    D := 100.0  // Dissociation energy
    a := 1.5    // Width parameter
    r0 := 2.0   // Equilibrium distance

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

    // for i := 0; i < numPoints; i++ {
    //    x[i] = float64(i) * 0.05
    //    // Morse potential: D * (1 - exp(-a*(x-r0)))^2
    //    y[i] = D * math.Pow(1-math.Exp(-a*(x[i]-r0)), 2)
    // }

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

// Main handler
func GeneratePlotData(w http.ResponseWriter, r *http.Request) {
    var req PlotRequest
    err := json.NewDecoder(r.Body).Decode(&req)
    if err != nil {
        http.Error(w, err.Error(), http.StatusBadRequest)
        return
    }

    var response PlotDataResponse

    switch req.PlotType {
    case "Morse":
        response = generateMorsePotentialSurfaceData(req.Parameters)
    default:
        http.Error(w, "Unknown plot type", http.StatusBadRequest)
        return
    }

    w.Header().Set("Content-Type", "application/json")
    json.NewEncoder(w).Encode(response)
}