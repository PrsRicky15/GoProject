
// ============================================================================
// CALCULATOR Namespace
// ============================================================================

export namespace calc {
    export class Calculator {
        protected result: number = 0;

        add(a: number, b: number): number {
            return this.result = a + b;
        }

        subtract(a: number, b: number): number {
            return this.result = a - b;
        }

        multiply(a: number, b: number): number {
            return this.result = a * b;
        }

        divide(a: number, b: number): number {
            if (b === 0) throw new Error("Division by zero");
            return this.result = a / b;
        }

        mod(a: number, b: number): number {
            return this.result = a % b;
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
            return this.result = Math.pow(base, exponent);
        }

        squareRoot(n: number): number {
            if (n < 0) throw new Error("Cannot calculate square root of negative number");
            return this.result = Math.sqrt(n);
        }

        nthRoot(n: number, nRoot: number): number {
            if (n < 0) throw new Error("Cannot calculate square root of negative number");
            return this.result = Math.pow(n, 1/nRoot);
        }

        sine(angle: number): number {
            return this.result = Math.sin(angle);
        }

        sineh(angle: number): number {
            return this.result = Math.sinh(angle);
        }

        cosec(angle: number): number {
            return this.result = 1./Math.sin(angle);
        }

        Asine(angle: number): number {
            return this.result = Math.asin(angle);
        }

        Asineh(angle: number): number {
            return this.result = Math.asinh(angle);
        }

        cosine(angle: number): number {
            return this.result = Math.cos(angle);
        }

        cosh(angle: number): number {
            return this.result = Math.cosh(angle);
        }

        sec(angle: number): number {
            return this.result = 1./Math.cos(angle);
        }

        Acos(angle: number): number {
            return this.result = Math.acos(angle);
        }

        Acosh(angle: number): number {
            return this.result = Math.acosh(angle);
        }

        tangent(angle: number): number {
            return this.result = Math.tan(angle);
        }

        cot(angle: number): number {
            return this.result = 1./Math.tan(angle);
        }

        tangenth(angle: number): number {
            return this.result = Math.tanh(angle);
        }

        atan(angle: number): number {
            return this.result = Math.atan(angle);
        }

        atanh(angle: number): number {
            return this.result = Math.atanh(angle);
        }

        logarithm(n: number, base: number = 10): number {
            if (n <= 0) throw new Error("Logarithm requires positive input");
            return this.result = Math.log(n) / Math.log(base);
        }

        naturalLog(n: number): number {
            if (n <= 0) throw new Error("Natural log requires positive input");
            return this.result = Math.log(n);
        }

        factorial(n: number): number {
            if (n < 0 || !Number.isInteger(n)) {
                throw new Error("Factorial requires non-negative integer");
            }
            let fact = 1;
            for (let i = 2; i <= n; i++) fact *= i;
            return this.result = fact;
        }
    }
}