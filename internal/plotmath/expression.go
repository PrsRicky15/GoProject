// Package plotmath evaluates a bounded arithmetic grammar, never executable code.
package plotmath

import (
	"fmt"
	"math"
	"regexp"
	"strconv"
	"strings"
)

type Evaluator func(x, y, z float64) float64

var tokenPattern = regexp.MustCompile(`^(?:\s+|(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?|[a-zA-Z_][a-zA-Z_0-9]*|[+\-*/%^(),])`)

var unaryFunctions = map[string]func(float64) float64{
	"sin": math.Sin, "cos": math.Cos, "tan": math.Tan,
	"asin": math.Asin, "acos": math.Acos, "atan": math.Atan,
	"sinh": math.Sinh, "cosh": math.Cosh, "tanh": math.Tanh,
	"sqrt": math.Sqrt, "abs": math.Abs, "exp": math.Exp,
	"ln": math.Log, "log": math.Log10, "log10": math.Log10,
	"floor": math.Floor, "ceil": math.Ceil,
	// Match the calculator: ties round toward positive infinity.
	"round": func(v float64) float64 {
		f := math.Floor(v)
		if v-f >= .5 {
			return f + 1
		}
		return f
	},
}
var binaryFunctions = map[string]func(float64, float64) float64{
	"pow": math.Pow, "min": math.Min, "max": math.Max,
}

type operation struct {
	precedence int
	apply      func(float64, float64) float64
}

var operations = map[string]operation{
	"+": {1, func(a, b float64) float64 { return a + b }},
	"-": {1, func(a, b float64) float64 { return a - b }},
	"*": {2, func(a, b float64) float64 { return a * b }},
	"/": {2, func(a, b float64) float64 { return a / b }},
	"%": {2, math.Mod}, "^": {4, math.Pow},
}

type parser struct {
	tokens             []string
	cursor, dimensions int
}

// Compile parses once, then reuses the evaluation tree at every sample point.
func Compile(source string, dimensions int) (Evaluator, error) {
	if len(source) > 1024 {
		return nil, fmt.Errorf("keep expressions under 1,024 bytes")
	}
	input := strings.ReplaceAll(strings.ReplaceAll(source, "π", "pi"), "**", "^")
	p := parser{dimensions: dimensions}
	for len(input) > 0 {
		token := tokenPattern.FindString(input)
		if token == "" {
			return nil, fmt.Errorf("unexpected character; use * for multiplication")
		}
		input = input[len(token):]
		if strings.TrimSpace(token) != "" {
			p.tokens = append(p.tokens, token)
		}
		if len(p.tokens) > 256 {
			return nil, fmt.Errorf("expression is too complex (maximum 256 tokens)")
		}
	}
	if len(p.tokens) == 0 {
		return nil, fmt.Errorf("enter an expression")
	}
	result, err := p.parse(0)
	if err != nil {
		return nil, err
	}
	if p.cursor != len(p.tokens) {
		return nil, fmt.Errorf("unexpected %q; use * for multiplication", p.tokens[p.cursor])
	}
	return result, nil
}

func (p *parser) consume(want string) error {
	if p.cursor >= len(p.tokens) || p.tokens[p.cursor] != want {
		return fmt.Errorf("expected %q; check parentheses and arguments", want)
	}
	p.cursor++
	return nil
}

func (p *parser) parse(minimum int) (Evaluator, error) {
	if p.cursor >= len(p.tokens) {
		return nil, fmt.Errorf("expression is incomplete")
	}
	token := p.tokens[p.cursor]
	p.cursor++
	var left Evaluator
	switch {
	case token == "+" || token == "-":
		operand, err := p.parse(3)
		if err != nil {
			return nil, err
		}
		left = operand
		if token == "-" {
			left = func(x, y, z float64) float64 { return -operand(x, y, z) }
		}
	case token == "(":
		var err error
		left, err = p.parse(0)
		if err != nil {
			return nil, err
		}
		if err = p.consume(")"); err != nil {
			return nil, err
		}
	case token == "pi":
		left = func(_, _, _ float64) float64 { return math.Pi }
	case token == "e":
		left = func(_, _, _ float64) float64 { return math.E }
	case token == "x" && p.dimensions >= 1:
		left = func(x, _, _ float64) float64 { return x }
	case token == "y" && p.dimensions >= 2:
		left = func(_, y, _ float64) float64 { return y }
	case token == "z" && p.dimensions >= 3:
		left = func(_, _, z float64) float64 { return z }
	case unaryFunctions[token] != nil:
		fn := unaryFunctions[token]
		if err := p.consume("("); err != nil {
			return nil, err
		}
		operand, err := p.parse(0)
		if err != nil {
			return nil, err
		}
		if err = p.consume(")"); err != nil {
			return nil, err
		}
		left = func(x, y, z float64) float64 { return fn(operand(x, y, z)) }
	case binaryFunctions[token] != nil:
		fn := binaryFunctions[token]
		if err := p.consume("("); err != nil {
			return nil, err
		}
		a, err := p.parse(0)
		if err != nil {
			return nil, err
		}
		if err = p.consume(","); err != nil {
			return nil, err
		}
		b, err := p.parse(0)
		if err != nil {
			return nil, err
		}
		if err = p.consume(")"); err != nil {
			return nil, err
		}
		left = func(x, y, z float64) float64 { return fn(a(x, y, z), b(x, y, z)) }
	default:
		value, err := strconv.ParseFloat(token, 64)
		if err != nil || !finite(value) {
			return nil, fmt.Errorf("unknown name or invalid number: %q", token)
		}
		left = func(_, _, _ float64) float64 { return value }
	}
	for p.cursor < len(p.tokens) {
		token := p.tokens[p.cursor]
		op, ok := operations[token]
		if !ok || op.precedence < minimum {
			break
		}
		p.cursor++
		next := op.precedence + 1
		if token == "^" {
			next--
		}
		right, err := p.parse(next)
		if err != nil {
			return nil, err
		}
		previous := left
		left = func(x, y, z float64) float64 { return op.apply(previous(x, y, z), right(x, y, z)) }
	}
	return left, nil
}

func finite(v float64) bool { return !math.IsNaN(v) && !math.IsInf(v, 0) }
