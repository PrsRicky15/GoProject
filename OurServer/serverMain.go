package OurServer

import (
	"GoProject/Plots"
	"fmt"
	"log"
	"net/http"

	"github.com/gorilla/mux"
	"github.com/rs/cors"
)

const PortNumber = ":8080"

func WebServer() {
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

// Home is the about page handler
func Home(w http.ResponseWriter, r *http.Request) {
}

// About is the about page handler
func About(w http.ResponseWriter, r *http.Request) {
}

func BasicWebServer() {
	http.HandleFunc("/", Home)
	http.HandleFunc("/about", About)

	fmt.Println(fmt.Sprintf("Server starting on :%s", PortNumber))
	err := http.ListenAndServe(PortNumber, nil)
	if err != nil {
		fmt.Println(err)
	}
}
