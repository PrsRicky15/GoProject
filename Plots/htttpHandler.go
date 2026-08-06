package Plots

import (
	"encoding/json"
	"log"
	"net/http"
)

// ==================== HTTP HANDLER ====================

// GeneratePlotData handles HTTP requests for plot generation
func GeneratePlotData(w http.ResponseWriter, r *http.Request) {
	// Parse request
	var req PlotRequest
	if err := json.NewDecoder(r.Body).Decode(&req); err != nil {
		http.Error(w, "Invalid request body: "+err.Error(), http.StatusBadRequest)
		return
	}

	// Route to appropriate plot generator
	var response PlotDataResponse
	switch req.PlotType {
	case "Morse":
		response = generateMorsePotentialSurfaceData(req.Grid, req.Parameters)
	case "Softcore":
		response = generateSoftcorePotentialSurfaceData(req.Grid, req.Parameters)
	default:
		http.Error(w, "Unknown plot type: "+req.PlotType, http.StatusBadRequest)
		return
	}

	// Send response
	w.Header().Set("Content-Type", "application/json")
	if err := json.NewEncoder(w).Encode(response); err != nil {
		log.Printf("Error encoding response: %v", err)
		http.Error(w, "Internal server error", http.StatusInternalServerError)
	}
}
