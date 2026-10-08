package plotmath

import (
	"context"
	"fmt"
	"math"
)

type Request struct {
	Mode       string               `json:"mode"`
	Expression string               `json:"expression"`
	Resolution int                  `json:"resolution"`
	Bounds     map[string][]float64 `json:"bounds"`
	Level      float64              `json:"level"`
	Style      string               `json:"style"`
}
type Result struct {
	Settings Request              `json:"settings"`
	Coords   map[string][]float64 `json:"coords"`
	Values   []*float64           `json:"values"` // JSON null represents an undefined sample.
	Count    int                  `json:"count"`
	Invalid  int                  `json:"invalid"`
	Min      float64              `json:"min"`
	Max      float64              `json:"max"`
}
type limits struct{ dimensions, min, max int }

var modes = map[string]limits{"line": {1, 20, 10000}, "heatmap": {2, 10, 200}, "surface": {2, 10, 150}, "iso": {3, 10, 60}}
var axisNames = []string{"x", "y", "z"}

func Sample(ctx context.Context, request Request) (*Result, error) {
	limit, ok := modes[request.Mode]
	if !ok {
		return nil, fmt.Errorf("choose a supported plot type")
	}
	n := request.Resolution
	if n < limit.min || n > limit.max {
		return nil, fmt.Errorf("samples per axis must be an integer from %d to %d", limit.min, limit.max)
	}
	if request.Style != "heatmap" && request.Style != "contour" {
		return nil, fmt.Errorf("choose heatmap or contour display")
	}
	if !finite(request.Level) {
		return nil, fmt.Errorf("enter a finite isosurface level")
	}
	if err := ctx.Err(); err != nil {
		return nil, err
	}
	coords := make(map[string][]float64, limit.dimensions)
	count := 1
	for _, axis := range axisNames[:limit.dimensions] {
		pair := request.Bounds[axis]
		if len(pair) != 2 || !finite(pair[0]) || !finite(pair[1]) || pair[0] >= pair[1] || !finite(pair[1]-pair[0]) {
			return nil, fmt.Errorf("%s needs finite minimum and maximum bounds, with minimum < maximum", axis)
		}
		span := pair[1] - pair[0]
		if pair[0]+span/float64(n-1) == pair[0] {
			return nil, fmt.Errorf("%s range is too small to sample at this magnitude", axis)
		}
		coords[axis] = make([]float64, n)
		for i := range n {
			coords[axis][i] = pair[0] + float64(i)/float64(n-1)*span
		}
		coords[axis][n-1] = pair[1]
		count *= n
	}
	evaluate, err := Compile(request.Expression, limit.dimensions)
	if err != nil {
		return nil, err
	}
	result := &Result{Settings: request, Coords: coords, Values: make([]*float64, count), Count: count, Min: math.Inf(1), Max: math.Inf(-1)}
	// One contiguous allocation backs all non-null JSON values.
	values := make([]float64, count)
	for i := range count {
		if i%256 == 0 {
			if err := ctx.Err(); err != nil {
				return nil, err
			}
		}
		x := coords["x"][i%n]
		var y, z float64
		if limit.dimensions >= 2 {
			y = coords["y"][(i/n)%n]
		}
		if limit.dimensions == 3 {
			z = coords["z"][i/(n*n)]
		}
		v := evaluate(x, y, z)
		if !finite(v) {
			result.Invalid++
			continue
		}
		values[i] = v
		result.Values[i] = &values[i]
		result.Min = math.Min(result.Min, v)
		result.Max = math.Max(result.Max, v)
	}
	if result.Invalid == count {
		return nil, fmt.Errorf("no finite real values in this domain; change the expression or bounds")
	}
	if request.Mode == "iso" {
		if result.Invalid > 0 {
			return nil, fmt.Errorf("isosurfaces need a finite value at every grid point")
		}
		if request.Level <= result.Min || request.Level >= result.Max {
			return nil, fmt.Errorf("choose a level strictly between %.5g and %.5g (the sampled range)", result.Min, result.Max)
		}
	}
	return result, nil
}
