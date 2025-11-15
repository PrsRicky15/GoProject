package main

import (
	"fmt"
    "log"
   	"net/http"

    "GoProject/Plots"
	"GoProject/gridData"
    "github.com/gorilla/mux"
    "github.com/rs/cors"
)

func main() {
	grid, _ := gridData.NewRGrid(0, 25, 100)
	fmt.Println(grid)

	router := mux.NewRouter()

    // API routes
    router.HandleFunc("/api/plots/data", Plots.GeneratePlotData).Methods("POST")

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
