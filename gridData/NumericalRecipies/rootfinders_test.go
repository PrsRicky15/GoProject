package NumericalRecipies

import (
    "math"
    "testing"

    "GoProject/gridData"
)

func approxEqual(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

func TestFindRoot_Sqrt2(t *testing.T) {
    f := gridData.RfuncFunc(func(x float64) float64 { return x*x - 2 })
    b := NewBisection(100)

    r, err := b.FindRoot(1, 2, 1e-12, f)
    if err != nil {
        t.Fatalf("unexpected error: %v", err)
    }
    if !approxEqual(r, math.Sqrt(2), 1e-9) {
        t.Fatalf("root mismatch: got %v want %v", r, math.Sqrt(2))
    }
}