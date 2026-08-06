// ============================================================================
// PLOT Namespace - Complete plotting functionality
// ============================================================================

import {Grid} from "./gridInformation";

export namespace plot {

    export interface PlotData {
        x: number[];
        y: number[];
        label?: string;
        color?: string;
    }

    export interface PlotOptions {
        title?: string;
        xlabel?: string;
        ylabel?: string;
        grid?: boolean;
        legend?: boolean;
        lineWidth?: number;
        colors?: string[];
    }

    export class FunctionEvaluator {
        private cache: Map<string, number[]> = new Map();

        evaluateOnGrid(fn: (x: number) => number, grid: Grid): PlotData {
            const { xMin, xMax, nx } = grid;
            const dx = (xMax - xMin) / (nx - 1);
            const x: number[] = [];
            const y: number[] = [];

            for (let i = 0; i < nx; i++) {
                const xi = xMin + i * dx;
                x.push(xi);
                y.push(fn(xi));
            }

            return { x, y };
        }

        evaluateMultiple(functions: Array<(x: number) => number>, grid: Grid): PlotData[] {
            return functions.map(fn => this.evaluateOnGrid(fn, grid));
        }

        clearCache(): void {
            this.cache.clear();
        }
    }

    export class Functions {
        // Basic trigonometric functions
        static sine = (x: number) => Math.sin(x);
        static cosine = (x: number) => Math.cos(x);
        static tangent = (x: number) => Math.tan(x);

        // Exponential and logarithmic functions
        static exp = (x: number) => Math.exp(x);
        static log = (x: number) => x > 0 ? Math.log(x) : NaN;
        static log10 = (x: number) => x > 0 ? Math.log10(x) : NaN;

        // Parameterized functions
        static gaussian(sigma: number = 1): (x: number) => number {
            return (x: number) => Math.exp(-Math.pow(x, 2) / (2 * sigma * sigma));
        }

        static polynomial(coefficients: number[]): (x: number) => number {
            return (x: number) => coefficients.reduce((sum, coeff, i) => sum + coeff * Math.pow(x, i), 0);
        }

        static power(exponent: number): (x: number) => number {
            return (x: number) => Math.pow(x, exponent);
        }

        static linear(slope: number, intercept: number = 0): (x: number) => number {
            return (x: number) => slope * x + intercept;
        }

        static sigmoid(k: number = 1): (x: number) => number {
            return (x: number) => 1 / (1 + Math.exp(-k * x));
        }

        static step(threshold: number = 0): (x: number) => number {
            return (x: number) => x >= threshold ? 1 : 0;
        }

        static abs(x: number): number {
            return Math.abs(x);
        }

        static sign(x: number): number {
            return Math.sign(x);
        }

        // Hyperbolic functions
        static sinh = (x: number) => Math.sinh(x);
        static cosh = (x: number) => Math.cosh(x);
        static tanh = (x: number) => Math.tanh(x);
    }

    // ========================================================================
    // PARAMETRIC PLOTTER
    // ========================================================================

    export interface ParametricCurve {
        x: (t: number) => number;
        y: (t: number) => number;
        tMin: number;
        tMax: number;
        nt: number;
    }

    export class ParametricPlotter {
        evaluateCurve(curve: ParametricCurve): PlotData {
            const { x: xFunc, y: yFunc, tMin, tMax, nt } = curve;
            const dt = (tMax - tMin) / (nt - 1);
            const x: number[] = [];
            const y: number[] = [];

            for (let i = 0; i < nt; i++) {
                const t = tMin + i * dt;
                x.push(xFunc(t));
                y.push(yFunc(t));
            }

            return { x, y };
        }

        circle(radius: number = 1, nt: number = 100): PlotData {
            return this.evaluateCurve({
                x: (t) => radius * Math.cos(t),
                y: (t) => radius * Math.sin(t),
                tMin: 0,
                tMax: 2 * Math.PI,
                nt
            });
        }

        ellipse(a: number, b: number, nt: number = 100): PlotData {
            return this.evaluateCurve({
                x: (t) => a * Math.cos(t),
                y: (t) => b * Math.sin(t),
                tMin: 0,
                tMax: 2 * Math.PI,
                nt
            });
        }

        lissajous(A: number, B: number, a: number, b: number, delta: number = 0, nt: number = 1000): PlotData {
            return this.evaluateCurve({
                x: (t) => A * Math.sin(a * t + delta),
                y: (t) => B * Math.sin(b * t),
                tMin: 0,
                tMax: 2 * Math.PI,
                nt
            });
        }

        spiral(a: number, b: number, nt: number = 1000): PlotData {
            return this.evaluateCurve({
                x: (t) => (a + b * t) * Math.cos(t),
                y: (t) => (a + b * t) * Math.sin(t),
                tMin: 0,
                tMax: 10 * Math.PI,
                nt
            });
        }
    }

    // ========================================================================
    // DATA ANALYSIS
    // ========================================================================
    export class DataAnalysis {
        static min(data: number[]): number {
            return Math.min(...data);
        }

        static max(data: number[]): number {
            return Math.max(...data);
        }

        static mean(data: number[]): number {
            return data.reduce((sum, val) => sum + val, 0) / data.length;
        }

        static median(data: number[]): number {
            if (data.length === 0) throw new Error("Cannot calculate median of empty array");
            const sorted = [...data].sort((a, b) => a - b);
            const mid = Math.floor(sorted.length / 2);
            return sorted.length % 2 ? sorted[mid]! : (sorted[mid - 1]! + sorted[mid]!) / 2;
        }

        static stdDev(data: number[]): number {
            const m = this.mean(data);
            const variance = data.reduce((sum, val) => sum + Math.pow(val - m, 2), 0) / data.length;
            return Math.sqrt(variance);
        }

        static variance(data: number[]): number {
            const m = this.mean(data);
            return data.reduce((sum, val) => sum + Math.pow(val - m, 2), 0) / data.length;
        }

        static range(data: number[]): number {
            return this.max(data) - this.min(data);
        }

        static sum(data: number[]): number {
            return data.reduce((sum, val) => sum + val, 0);
        }
    }

    // ========================================================================
    // GRID UTILITIES
    // ========================================================================
    export class GridUtils {
        static createLinearGrid(xMin: number, xMax: number, nx: number): Grid {
            return { xMin, xMax, nx };
        }

        static createLogGrid(xMin: number, xMax: number, nx: number): number[] {
            const logMin = Math.log10(xMin);
            const logMax = Math.log10(xMax);
            const dx = (logMax - logMin) / (nx - 1);
            return Array.from({ length: nx }, (_, i) => Math.pow(10, logMin + i * dx));
        }

        static mergeGrids(grid1: Grid, grid2: Grid): Grid {
            return {
                xMin: Math.min(grid1.xMin, grid2.xMin),
                xMax: Math.max(grid1.xMax, grid2.xMax),
                nx: Math.max(grid1.nx, grid2.nx)
            };
        }

        static refineMesh(x: number[], y: number[], threshold: number): { x: number[], y: number[] } {
            const newX: number[] = [x[0]!];
            const newY: number[] = [y[0]!];

            for (let i = 1; i < x.length; i++) {
                const dx = Math.abs(x[i]! - x[i - 1]!);
                const dy = Math.abs(y[i]! - y[i - 1]!);

                if (dx > threshold || dy > threshold) {
                    const midX = (x[i]! + x[i - 1]!) / 2;
                    const midY = (y[i]! + y[i - 1]!) / 2;
                    newX.push(midX);
                    newY.push(midY);
                }

                newX.push(x[i]!);
                newY.push(y[i]!);
            }

            return { x: newX, y: newY };
        }
    }


    // ========================================================================
    // TRANSFORMATIONS
    // ========================================================================
    export class Transformations {
        static translate(data: PlotData, dx: number, dy: number): PlotData {
            return {
                x: data.x.map(x => x + dx),
                y: data.y.map(y => y + dy),
                label: data.label,
                color: data.color
            };
        }

        static scale(data: PlotData, sx: number, sy: number): PlotData {
            return {
                x: data.x.map(x => x * sx),
                y: data.y.map(y => y * sy),
                label: data.label,
                color: data.color
            };
        }

        static rotate(data: PlotData, angle: number): PlotData {
            const cos = Math.cos(angle);
            const sin = Math.sin(angle);
            const x: number[] = [];
            const y: number[] = [];

            for (let i = 0; i < data.x.length; i++) {
                x.push(data.x[i]! * cos - data.y[i]! * sin);
                y.push(data.x[i]! * sin + data.y[i]! * cos);
            }

            return { x, y, label: data.label, color: data.color };
        }

        static normalize(data: PlotData): PlotData {
            const minY = DataAnalysis.min(data.y);
            const maxY = DataAnalysis.max(data.y);
            const range = maxY - minY;

            return {
                x: [...data.x],
                y: data.y.map(y => (y - minY) / range),
                label: data.label,
                color: data.color
            };
        }
    }
}