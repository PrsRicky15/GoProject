export class Calculator {
    protected result = 0;

    add(a: number, b: number): number {
        this.result = a + b;
        return this.result;
    }

    subtract(a: number, b: number): number {
        this.result = a - b;
        return this.result;
    }

    multiply(a: number, b: number): number {
        this.result = a * b;
        return this.result;
    }

    divide(a: number, b: number): number {
        if (b === 0) throw new Error("Division by zero");
        this.result = a / b;
        return this.result;
    }

    getResult(): number {
        return this.result;
    }

    clear(): void {
        this.result = 0;
    }
}

export class ScientificCalculator extends Calculator {
    power(base: number, exponent: number): number {
        this.result = Math.pow(base, exponent);
        return this.result;
    }

    squareRoot(n: number): number {
        if (n < 0) throw new Error("Cannot calculate square root of negative number");
        this.result = Math.sqrt(n);
        return this.result;
    }

    nthRoot(n: number, root: number): number {
        if (n < 0 && root % 2 === 0) {
            throw new Error("Cannot calculate even root of negative number");
        }
        this.result = Math.pow(Math.abs(n), 1 / root) * (n < 0 ? -1 : 1);
        return this.result;
    }

    // Trigonometric functions
    sine(angle: number): number {
        this.result = Math.sin(angle);
        return this.result;
    }

    cosine(angle: number): number {
        this.result = Math.cos(angle);
        return this.result;
    }

    tangent(angle: number): number {
        this.result = Math.tan(angle);
        return this.result;
    }

    // Inverse trigonometric functions
    arcSine(x: number): number {
        if (x < -1 || x > 1) throw new Error("arcSine domain error: x must be in [-1, 1]");
        this.result = Math.asin(x);
        return this.result;
    }

    arcCosine(x: number): number {
        if (x < -1 || x > 1) throw new Error("arcCosine domain error: x must be in [-1, 1]");
        this.result = Math.acos(x);
        return this.result;
    }

    arcTangent(x: number): number {
        this.result = Math.atan(x);
        return this.result;
    }

    arcTangent2(y: number, x: number): number {
        this.result = Math.atan2(y, x);
        return this.result;
    }

    // Hyperbolic functions
    sinh(x: number): number {
        this.result = Math.sinh(x);
        return this.result;
    }

    cosh(x: number): number {
        this.result = Math.cosh(x);
        return this.result;
    }

    tanh(x: number): number {
        this.result = Math.tanh(x);
        return this.result;
    }

    // Inverse hyperbolic functions
    arcSinh(x: number): number {
        this.result = Math.asinh(x);
        return this.result;
    }

    arcCosh(x: number): number {
        if (x < 1) throw new Error("arcCosh domain error: x must be >= 1");
        this.result = Math.acosh(x);
        return this.result;
    }

    arcTanh(x: number): number {
        if (x <= -1 || x >= 1) throw new Error("arcTanh domain error: x must be in (-1, 1)");
        this.result = Math.atanh(x);
        return this.result;
    }

    // Logarithmic and exponential
    logarithm(n: number, base = 10): number {
        if (n <= 0) throw new Error("Logarithm requires positive input");
        this.result = Math.log(n) / Math.log(base);
        return this.result;
    }

    naturalLog(n: number): number {
        if (n <= 0) throw new Error("Natural log requires positive input");
        this.result = Math.log(n);
        return this.result;
    }

    exp(x: number): number {
        this.result = Math.exp(x);
        return this.result;
    }

    // Special functions
    factorial(n: number): number {
        if (n < 0 || !Number.isInteger(n)) {
            throw new Error("Factorial requires non-negative integer");
        }

        let fact = 1;
        for (let i = 2; i <= n; i++) {
            fact *= i;
        }

        this.result = fact;
        return this.result;
    }

    abs(x: number): number {
        this.result = Math.abs(x);
        return this.result;
    }

    ceil(x: number): number {
        this.result = Math.ceil(x);
        return this.result;
    }

    floor(x: number): number {
        this.result = Math.floor(x);
        return this.result;
    }

    round(x: number, decimals = 0): number {
        const multiplier = Math.pow(10, decimals);
        this.result = Math.round(x * multiplier) / multiplier;
        return this.result;
    }

    modulo(a: number, b: number): number {
        this.result = a % b;
        return this.result;
    }

    sign(x: number): number {
        this.result = Math.sign(x);
        return this.result;
    }

    // Statistical functions
    combination(n: number, k: number): number {
        if (n < 0 || k < 0 || !Number.isInteger(n) || !Number.isInteger(k)) {
            throw new Error("Combination requires non-negative integers");
        }
        if (k > n) return 0;

        this.result = this.factorial(n) / (this.factorial(k) * this.factorial(n - k));
        return this.result;
    }

    permutation(n: number, k: number): number {
        if (n < 0 || k < 0 || !Number.isInteger(n) || !Number.isInteger(k)) {
            throw new Error("Permutation requires non-negative integers");
        }
        if (k > n) return 0;

        this.result = this.factorial(n) / this.factorial(n - k);
        return this.result;
    }
}

export class MathFunctions {
    // Basic trigonometric
    static readonly sine = (x: number): number => Math.sin(x);
    static readonly cosine = (x: number): number => Math.cos(x);
    static readonly tangent = (x: number): number => Math.tan(x);

    // Inverse trigonometric
    static readonly arcSine = (x: number): number => Math.asin(x);
    static readonly arcCosine = (x: number): number => Math.acos(x);
    static readonly arcTangent = (x: number): number => Math.atan(x);

    // Hyperbolic
    static readonly sinh = (x: number): number => Math.sinh(x);
    static readonly cosh = (x: number): number => Math.cosh(x);
    static readonly tanh = (x: number): number => Math.tanh(x);

    // Inverse hyperbolic
    static readonly arcSinh = (x: number): number => Math.asinh(x);
    static readonly arcCosh = (x: number): number => Math.acosh(x);
    static readonly arcTanh = (x: number): number => Math.atanh(x);

    // Exponential and logarithmic
    static readonly exp = (x: number): number => Math.exp(x);
    static readonly log = (x: number): number => (x > 0 ? Math.log(x) : NaN);
    static readonly log10 = (x: number): number => (x > 0 ? Math.log10(x) : NaN);
    static readonly log2 = (x: number): number => (x > 0 ? Math.log2(x) : NaN);

    // Basic operations
    static readonly abs = (x: number): number => Math.abs(x);
    static readonly sign = (x: number): number => Math.sign(x);
    static readonly ceil = (x: number): number => Math.ceil(x);
    static readonly floor = (x: number): number => Math.floor(x);
    static readonly round = (x: number): number => Math.round(x);
    static readonly trunc = (x: number): number => Math.trunc(x);

    // Root functions
    static readonly sqrt = (x: number): number => Math.sqrt(x);
    static readonly cbrt = (x: number): number => Math.cbrt(x);

    // Piecewise functions
    static readonly heaviside = (x: number): number => (x >= 0 ? 1 : 0);
    static readonly rectPulse = (x: number, width = 1): number =>
        (Math.abs(x) <= width / 2 ? 1 : 0);

    // Clamping and modulo
    static clamp(min: number, max: number): (x: number) => number {
        return (x: number) => Math.max(min, Math.min(max, x));
    }

    static modulo(divisor: number): (x: number) => number {
        return (x: number) => ((x % divisor) + divisor) % divisor;
    }

    // Parameterized function factories
    static gaussian(sigma = 1, mu = 0): (x: number) => number {
        return (x: number) =>
            Math.exp(-Math.pow(x - mu, 2) / (2 * sigma * sigma)) /
            (sigma * Math.sqrt(2 * Math.PI));
    }

    static polynomial(coefficients: number[]): (x: number) => number {
        return (x: number) =>
            coefficients.reduce((sum, coeff, i) => sum + coeff * Math.pow(x, i), 0);
    }

    static power(exponent: number): (x: number) => number {
        return (x: number) => Math.pow(x, exponent);
    }

    static linear(slope: number, intercept = 0): (x: number) => number {
        return (x: number) => slope * x + intercept;
    }

    static sigmoid(k = 1): (x: number) => number {
        return (x: number) => 1 / (1 + Math.exp(-k * x));
    }

    static step(threshold = 0): (x: number) => number {
        return (x: number) => (x >= threshold ? 1 : 0);
    }

    static sawtooth(period = 1): (x: number) => number {
        return (x: number) => 2 * (x / period - Math.floor(x / period + 0.5));
    }

    static triangle(period = 1): (x: number) => number {
        return (x: number) => {
            const normalized = x / period - Math.floor(x / period);
            return normalized < 0.5 ? 4 * normalized - 1 : 3 - 4 * normalized;
        };
    }

    static square(period = 1, dutyCycle = 0.5): (x: number) => number {
        return (x: number) => {
            const normalized = x / period - Math.floor(x / period);
            return normalized < dutyCycle ? 1 : -1;
        };
    }

    // Exponential decay/growth
    static exponentialDecay(lambda: number): (x: number) => number {
        return (x: number) => Math.exp(-lambda * x);
    }

    static exponentialGrowth(lambda: number): (x: number) => number {
        return (x: number) => Math.exp(lambda * x);
    }

    // Damped oscillations
    static dampedSine(omega: number, gamma: number): (x: number) => number {
        return (x: number) => Math.exp(-gamma * x) * Math.sin(omega * x);
    }

    static dampedCosine(omega: number, gamma: number): (x: number) => number {
        return (x: number) => Math.exp(-gamma * x) * Math.cos(omega * x);
    }

    // Sinc function
    static sinc(x: number): number {
        if (x === 0) return 1;
        return Math.sin(Math.PI * x) / (Math.PI * x);
    }

    // Logistic function
    static logistic(L: number, k: number, x0: number): (x: number) => number {
        return (x: number) => L / (1 + Math.exp(-k * (x - x0)));
    }
}