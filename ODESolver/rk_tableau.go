package EquationSolver

import "math"

// RKMethod selects an explicit Runge-Kutta formula. Embedded pairs use their
// higher-order solution for both fixed and adaptive steps.
type RKMethod string

const (
	RKEuler             RKMethod = "euler"
	RKHeun2             RKMethod = "heun2"
	RKMidpoint2         RKMethod = "midpoint2"
	RKRalston2          RKMethod = "ralston2"
	RKRalston3          RKMethod = "ralston3"
	RKHeun3             RKMethod = "heun3"
	RKWray3             RKMethod = "wray3"
	RKSSP3              RKMethod = "ssp3"
	RKKutta3            RKMethod = "kutta3"
	RKClassic4          RKMethod = "rk4"
	RKThreeEighths4     RKMethod = "rk38"
	RKRalston4          RKMethod = "ralston4"
	RKNystrom5          RKMethod = "nystrom5"
	RKHeunEuler21       RKMethod = "heun-euler21"
	RKFehlberg21        RKMethod = "fehlberg21"
	RKBogackiShampine32 RKMethod = "bogacki-shampine32"
	RKFehlberg54        RKMethod = "fehlberg54"
	RKCashKarp54        RKMethod = "cash-karp54"
	RKDormandPrince54   RKMethod = "dormand-prince54"
)

type rkTableau struct {
	a               [][]float64
	b, c, low       []float64
	order, lowOrder int
	fsal            bool
}

func rkTable(m RKMethod) (q rkTableau, ok bool) {
	ok = true
	switch m {
	case RKEuler:
		q = rkTableau{a: [][]float64{{}}, b: []float64{1}, c: []float64{0}, order: 1}
	case RKHeun2, RKHeunEuler21:
		q = rkTableau{a: [][]float64{{}, {1}}, b: []float64{.5, .5}, c: []float64{0, 1}, order: 2}
		if m == RKHeunEuler21 {
			q.low = []float64{1, 0}
			q.lowOrder = 1
		}
	case RKMidpoint2:
		q = rkTableau{a: [][]float64{{}, {.5}}, b: []float64{0, 1}, c: []float64{0, .5}, order: 2}
	case RKRalston2:
		q = rkTableau{a: [][]float64{{}, {2. / 3}}, b: []float64{.25, .75}, c: []float64{0, 2. / 3}, order: 2}
	case RKRalston3:
		q = rkTableau{a: [][]float64{{}, {.5}, {0, .75}}, b: []float64{2. / 9, 1. / 3, 4. / 9}, c: []float64{0, .5, .75}, order: 3}
	case RKHeun3:
		q = rkTableau{a: [][]float64{{}, {1. / 3}, {0, 2. / 3}}, b: []float64{.25, 0, .75}, c: []float64{0, 1. / 3, 2. / 3}, order: 3}
	case RKWray3:
		q = rkTableau{a: [][]float64{{}, {8. / 15}, {.25, 5. / 12}}, b: []float64{.25, 0, .75}, c: []float64{0, 8. / 15, 2. / 3}, order: 3}
	case RKSSP3:
		q = rkTableau{a: [][]float64{{}, {1}, {.25, .25}}, b: []float64{1. / 6, 1. / 6, 2. / 3}, c: []float64{0, 1, .5}, order: 3}
	case RKKutta3:
		q = rkTableau{a: [][]float64{{}, {.5}, {-1, 2}}, b: []float64{1. / 6, 2. / 3, 1. / 6}, c: []float64{0, .5, 1}, order: 3}
	case RKClassic4:
		q = rkTableau{a: [][]float64{{}, {.5}, {0, .5}, {0, 0, 1}}, b: []float64{1. / 6, 1. / 3, 1. / 3, 1. / 6}, c: []float64{0, .5, .5, 1}, order: 4}
	case RKThreeEighths4:
		q = rkTableau{a: [][]float64{{}, {1. / 3}, {-1. / 3, 1}, {1, -1, 1}}, b: []float64{.125, .375, .375, .125}, c: []float64{0, 1. / 3, 2. / 3, 1}, order: 4}
	case RKRalston4:
		s := math.Sqrt(5)
		q = rkTableau{a: [][]float64{{}, {.4}, {(-2889 + 1428*s) / 1024, (3785 - 1620*s) / 1024}, {(-3365 + 2094*s) / 6040, (-975 - 3046*s) / 2552, (467040 + 203968*s) / 240845}}, b: []float64{(263 + 24*s) / 1812, (125 - 1000*s) / 3828, (3426304 + 1661952*s) / 5924787, (30 - 4*s) / 123}, c: []float64{0, .4, (14 - 3*s) / 16, 1}, order: 4}
	case RKNystrom5:
		q = rkTableau{a: [][]float64{{}, {1. / 3}, {4. / 25, 6. / 25}, {.25, -3, 15. / 4}, {2. / 27, 10. / 9, -50. / 81, 8. / 81}, {2. / 25, 12. / 25, 2. / 15, 8. / 75}}, b: []float64{23. / 192, 0, 125. / 192, 0, -27. / 64, 125. / 192}, c: []float64{0, 1. / 3, 2. / 5, 1, 2. / 3, 4. / 5}, order: 5}
	case RKFehlberg21:
		q = rkTableau{a: [][]float64{{}, {.5}, {1. / 256, 255. / 256}}, b: []float64{1. / 512, 255. / 256, 1. / 512}, low: []float64{1. / 256, 255. / 256, 0}, c: []float64{0, .5, 1}, order: 2, lowOrder: 1}
	case RKBogackiShampine32:
		q = rkTableau{a: [][]float64{{}, {.5}, {0, .75}, {2. / 9, 1. / 3, 4. / 9}}, b: []float64{2. / 9, 1. / 3, 4. / 9, 0}, low: []float64{7. / 24, .25, 1. / 3, .125}, c: []float64{0, .5, .75, 1}, order: 3, lowOrder: 2, fsal: true}
	case RKFehlberg54:
		q = rkTableau{a: [][]float64{{}, {.25}, {3. / 32, 9. / 32}, {1932. / 2197, -7200. / 2197, 7296. / 2197}, {439. / 216, -8, 3680. / 513, -845. / 4104}, {-8. / 27, 2, -3544. / 2565, 1859. / 4104, -11. / 40}}, b: []float64{16. / 135, 0, 6656. / 12825, 28561. / 56430, -9. / 50, 2. / 55}, low: []float64{25. / 216, 0, 1408. / 2565, 2197. / 4104, -.2, 0}, c: []float64{0, .25, 3. / 8, 12. / 13, 1, .5}, order: 5, lowOrder: 4}
	case RKCashKarp54:
		q = rkTableau{a: [][]float64{{}, {.2}, {3. / 40, 9. / 40}, {.3, -.9, 1.2}, {-11. / 54, 2.5, -70. / 27, 35. / 27}, {1631. / 55296, 175. / 512, 575. / 13824, 44275. / 110592, 253. / 4096}}, b: []float64{37. / 378, 0, 250. / 621, 125. / 594, 0, 512. / 1771}, low: []float64{2825. / 27648, 0, 18575. / 48384, 13525. / 55296, 277. / 14336, .25}, c: []float64{0, .2, .3, .6, 1, .875}, order: 5, lowOrder: 4}
	case RKDormandPrince54:
		q = rkTableau{a: [][]float64{{}, {.2}, {3. / 40, 9. / 40}, {44. / 45, -56. / 15, 32. / 9}, {19372. / 6561, -25360. / 2187, 64448. / 6561, -212. / 729}, {9017. / 3168, -355. / 33, 46732. / 5247, 49. / 176, -5103. / 18656}, {35. / 384, 0, 500. / 1113, 125. / 192, -2187. / 6784, 11. / 84}}, b: []float64{35. / 384, 0, 500. / 1113, 125. / 192, -2187. / 6784, 11. / 84, 0}, low: []float64{5179. / 57600, 0, 7571. / 16695, 393. / 640, -92097. / 339200, 187. / 2100, 1. / 40}, c: []float64{0, .2, .3, .8, 8. / 9, 1, 1}, order: 5, lowOrder: 4, fsal: true}
	default:
		ok = false
	}
	return
}
