import { MathFunctions } from "./calculator_and_MathFunc";

const FUNCTION_ALIASES: Record<string, string> = {
    // exponentials
    "exp": "exp",
    "Exp": "exp",
    "EXP": "exp",
    "e^": "exp",
    "e**": "exp",

    // logarithms
    "ln": "log",
    "Ln": "log",
    "LOG": "log",
    "Log": "log",

    // trigonometric
    "sin": "sin",
    "Sin": "sin",
    "cos": "cos",
    "Cos": "cos",
    "tan": "tan",
    "Tan": "tan",
    "tg": "tan",
    "Tg": "tan",

    // inverse trig (optional)
    "asin": "asin",
    "acos": "acos",
    "atan": "atan"
};

function normalizeFunctions(expr: string): string {
    let result = expr;

    for (const [alias, canonical] of Object.entries(FUNCTION_ALIASES)) {
        const regex = new RegExp(`\\b${alias}\\s*\\(`, "g");
        result = result.replace(regex, `${canonical}(`);
    }

    return result;
}


export class Point {
    constructor(
        public x: number,
        public y: number
    ) { }

    translate(dx: number, dy: number): Point {
        return new Point(this.x + dx, this.y + dy);
    }

    shiftX(dx: number): Point {
        return new Point(this.x + dx, this.y);
    }

    shiftY(dy: number): Point {
        return new Point(this.x, this.y + dy);
    }

    scale(sx: number, sy: number): Point {
        return new Point(this.x * sx, this.y * sy);
    }

    rotate(angle: number): Point {
        const cos = Math.cos(angle);
        const sin = Math.sin(angle);
        return new Point(
            this.x * cos - this.y * sin,
            this.x * sin + this.y * cos
        );
    }

    static distance(p1: Point, p2: Point): number {
        const dx = p2.x - p1.x;
        const dy = p2.y - p1.y;
        return Math.sqrt(dx * dx + dy * dy);
    }
}

export class Point3D {
    constructor(
        public x: number,
        public y: number,
        public z: number
    ) { }

    translate(dx: number, dy: number, dz: number): Point3D {
        return new Point3D(this.x + dx, this.y + dy, this.z + dz);
    }

    shiftX(dx: number): Point3D {
        return new Point3D(this.x + dx, this.y, this.z);
    }

    shiftY(dy: number): Point3D {
        return new Point3D(this.x, this.y + dy, this.z);
    }

    shiftZ(dz: number): Point3D {
        return new Point3D(this.x, this.y, this.z + dz);
    }

    scale(sx: number, sy: number, sz: number): Point3D {
        return new Point3D(this.x * sx, this.y * sy, this.z * sz);
    }

    rotateZ(angle: number): Point3D {
        const cos = Math.cos(angle);
        const sin = Math.sin(angle);
        return new Point3D(
            this.x * cos - this.y * sin,
            this.x * sin + this.y * cos,
            this.z
        );
    }

    rotateX(angle: number): Point3D {
        const cos = Math.cos(angle);
        const sin = Math.sin(angle);
        return new Point3D(
            this.x,
            this.y * cos - this.z * sin,
            this.y * sin + this.z * cos
        );
    }

    rotateY(angle: number): Point3D {
        const cos = Math.cos(angle);
        const sin = Math.sin(angle);
        return new Point3D(
            this.x * cos + this.z * sin,
            this.y,
            -this.x * sin + this.z * cos
        );
    }

    static distance(p1: Point3D, p2: Point3D): number {
        const dx = p2.x - p1.x;
        const dy = p2.y - p1.y;
        const dz = p2.z - p1.z;
        return Math.sqrt(dx * dx + dy * dy + dz * dz);
    }
}

export interface PlotData {
    x: number[];
    y: number[];
    label?: string;
    color?: string;
}

export interface PlotData3D {
    x: number[];
    y: number[];
    z: number[];
    label?: string;
    color?: string;
}

export interface Grid {
    xMin: number;
    xMax: number;
    nx: number;
}

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
            nx: Math.max(grid1.nx, grid2.nx),
        };
    }

    static refineMesh(
        x: number[],
        y: number[],
        threshold: number
    ): { x: number[]; y: number[] } {
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

export interface PlotOptions {
    title?: string;
    xlabel?: string;
    ylabel?: string;
    grid?: boolean;
    legend?: boolean;
    lineWidth?: number;
    colors?: string[];
}

export interface ParametricCurve {
    x: (t: number) => number;
    y: (t: number) => number;
    tMin: number;
    tMax: number;
    nt: number;
}

export class FunctionEvaluator extends MathFunctions {
    private cache = new Map<string, number[]>();

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

    evaluateMultiple(
        functions: Array<(x: number) => number>,
        grid: Grid
    ): PlotData[] {
        return functions.map(fn => this.evaluateOnGrid(fn, grid));
    }

    clearCache(): void {
        this.cache.clear();
    }
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

    circle(radius = 1, nt = 100): PlotData {
        return this.evaluateCurve({
            x: (t) => radius * Math.cos(t),
            y: (t) => radius * Math.sin(t),
            tMin: 0,
            tMax: 2 * Math.PI,
            nt,
        });
    }

    ellipse(a: number, b: number, nt = 100): PlotData {
        return this.evaluateCurve({
            x: (t) => a * Math.cos(t),
            y: (t) => b * Math.sin(t),
            tMin: 0,
            tMax: 2 * Math.PI,
            nt,
        });
    }

    lissajous(
        A: number,
        B: number,
        a: number,
        b: number,
        delta = 0,
        nt = 1000
    ): PlotData {
        return this.evaluateCurve({
            x: (t) => A * Math.sin(a * t + delta),
            y: (t) => B * Math.sin(b * t),
            tMin: 0,
            tMax: 2 * Math.PI,
            nt,
        });
    }

    spiral(a: number, b: number, nt = 1000): PlotData {
        return this.evaluateCurve({
            x: (t) => (a + b * t) * Math.cos(t),
            y: (t) => (a + b * t) * Math.sin(t),
            tMin: 0,
            tMax: 10 * Math.PI,
            nt,
        });
    }
}

export class DataAnalysis {
    static min(data: number[]): number {
        if (data.length === 0) throw new Error("Cannot find min of empty array");
        return Math.min(...data);
    }

    static max(data: number[]): number {
        if (data.length === 0) throw new Error("Cannot find max of empty array");
        return Math.max(...data);
    }

    static mean(data: number[]): number {
        if (data.length === 0) throw new Error("Cannot find mean of empty array");
        return data.reduce((sum, val) => sum + val, 0) / data.length;
    }

    static median(data: number[]): number {
        if (data.length === 0) throw new Error("Cannot calculate median of empty array");

        const sorted = [...data].sort((a, b) => a - b);
        const mid = Math.floor(sorted.length / 2);

        return sorted.length % 2 === 0
            ? (sorted[mid - 1]! + sorted[mid]!) / 2
            : sorted[mid]!;
    }

    static variance(data: number[]): number {
        if (data.length === 0) throw new Error("Cannot calculate variance of empty array");

        const m = this.mean(data);
        return data.reduce((sum, val) => sum + Math.pow(val - m, 2), 0) / data.length;
    }

    static stdDev(data: number[]): number {
        return Math.sqrt(this.variance(data));
    }

    static range(data: number[]): number {
        return this.max(data) - this.min(data);
    }

    static sum(data: number[]): number {
        return data.reduce((sum, val) => sum + val, 0);
    }

    static interpolate(x: number[], y: number[], xi: number): number {
        if (x.length === 0 || y.length === 0) {
            throw new Error("Cannot interpolate with empty arrays");
        }
        if (x.length !== y.length) {
            throw new Error("Arrays must have same length");
        }

        if (xi <= x[0]!) return y[0]!;
        if (xi >= x[x.length - 1]!) return y[y.length - 1]!;

        for (let i = 0; i < x.length - 1; i++) {
            if (xi >= x[i]! && xi <= x[i + 1]!) {
                const t = (xi - x[i]!) / (x[i + 1]! - x[i]!);
                return y[i]! + t * (y[i + 1]! - y[i]!);
            }
        }

        return NaN;
    }

    static correlationCoefficient(x: number[], y: number[]): number {
        if (x.length !== y.length) {
            throw new Error("Arrays must have same length");
        }
        if (x.length === 0) {
            throw new Error("Cannot calculate correlation of empty arrays");
        }

        const meanX = this.mean(x);
        const meanY = this.mean(y);

        let numerator = 0;
        let denomX = 0;
        let denomY = 0;

        for (let i = 0; i < x.length; i++) {
            const dx = x[i]! - meanX;
            const dy = y[i]! - meanY;
            numerator += dx * dy;
            denomX += dx * dx;
            denomY += dy * dy;
        }

        return numerator / Math.sqrt(denomX * denomY);
    }
}

export class DataExporter {
    static toCSV(data: PlotData): string {
        const rows = data.x.map((x, i) => `${x},${data.y[i]}`);
        return `x,y\n${rows.join("\n")}`;
    }

    static toJSON(data: PlotData | PlotData[]): string {
        return JSON.stringify(data, null, 2);
    }

    static fromCSV(csv: string): PlotData {
        const lines = csv.trim().split("\n").slice(1);
        const x: number[] = [];
        const y: number[] = [];

        for (const line of lines) {
            const [xi, yi] = line.split(",").map(Number);
            if (xi !== undefined && yi !== undefined) {
                x.push(xi);
                y.push(yi);
            }
        }

        return { x, y };
    }

    static fromJSON(json: string): PlotData | PlotData[] {
        return JSON.parse(json);
    }

    static toTable(data: PlotData, precision = 4): string {
        let table = "x".padEnd(15) + "y".padEnd(15) + "\n";
        table += "-".repeat(30) + "\n";

        for (let i = 0; i < data.x.length; i++) {
            const xStr = data.x[i]!.toFixed(precision).padEnd(15);
            const yStr = data.y[i]!.toFixed(precision).padEnd(15);
            table += xStr + yStr + "\n";
        }

        return table;
    }
}