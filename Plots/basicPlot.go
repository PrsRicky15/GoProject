package Plots

// extractMorseParams extracts Morse potential parameters
func extractMorseParams(params map[string]interface{}) MorseParams {
	morse := MorseParams{
		D:  defaultMorseD,
		A:  defaultMorseA,
		R0: defaultMorseR0,
	}

	if val, ok := params["D"].(float64); ok {
		morse.D = val
	}
	if val, ok := params["a"].(float64); ok {
		morse.A = val
	}
	if val, ok := params["r0"].(float64); ok {
		morse.R0 = val
	}

	return morse
}

// extractSoftcoreParams extracts Softcore potential parameters
func extractSoftcoreParams(params map[string]interface{}) SoftcoreParams {
	sc := SoftcoreParams{
		Charge: defaultSoftcoreCharge,
		A:      defaultSoftcoreA,
		R0:     defaultSoftcoreR0,
	}

	if val, ok := params["Charge"].(float64); ok {
		sc.Charge = val
	}
	if val, ok := params["a"].(float64); ok {
		sc.A = val
	}
	if val, ok := params["r0"].(float64); ok {
		sc.R0 = val
	}

	return sc
}
