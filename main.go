package main

import (
	"GoProject/BasicOneD"
	"GoProject/Plots"
	"log"
	"net/http"

	"github.com/gorilla/mux"
	"github.com/rs/cors"
)

func webserver() {
	router := mux.NewRouter()

	// API routes
	router.HandleFunc(
		"/api/plots/data",
		Plots.GeneratePlotData,
	).Methods("POST")

	// Enable CORS
	c := cors.New(cors.Options{
		AllowedOrigins:   []string{"http://localhost:5173"},
		AllowedMethods:   []string{"GET", "POST", "OPTIONS"},
		AllowedHeaders:   []string{"Content-Type"},
		AllowCredentials: true,
	})
	handler := c.Handler(router)

	log.Println("Server starting on :8080")
	log.Fatal(http.ListenAndServe(":8080", handler))
}

func main() {
	webserver()
	err := BasicOneD.OneDimPoissonSolver()
	if err != nil {
		log.Fatal(err)
	}

	err = BasicOneD.BarrierPotential()
	if err != nil {
		log.Fatal(err)
	}
}
