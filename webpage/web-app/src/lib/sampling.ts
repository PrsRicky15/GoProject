import type { Axis, Mode, PlotDraft, PlotRequest } from './types';

interface ModeConfig { label: string; formula: string; axes: Axis[]; min: number; max: number; resolution: number; expression: string; bounds: [number, number] }
export const MODES: Record<Mode, ModeConfig> = {
  line: { label: '1D curve', formula: 'y = f(x)', axes: ['x'], min: 20, max: 10000, resolution: 801, expression: 'exp(-x^2 / 2)', bounds: [-5, 5] },
  heatmap: { label: '2D field', formula: 'color = f(x, y)', axes: ['x', 'y'], min: 10, max: 200, resolution: 81, expression: 'sin(x) * cos(y)', bounds: [-5, 5] },
  surface: { label: '3D surface', formula: 'z = f(x, y)', axes: ['x', 'y'], min: 10, max: 150, resolution: 81, expression: 'exp(-(x^2 + y^2) / 2)', bounds: [-4, 4] },
  iso: { label: 'Isosurface', formula: 'f(x, y, z) = level', axes: ['x', 'y', 'z'], min: 10, max: 60, resolution: 31, expression: 'x^2 + y^2 + z^2', bounds: [-2, 2] },
};
export const MODE_KEYS = Object.keys(MODES) as Mode[];
export function defaults(mode: Mode): PlotDraft {
  const config = MODES[mode];
  return { mode, expression: config.expression, resolution: config.resolution,
    bounds: { x: [...config.bounds], y: [...config.bounds], z: [...config.bounds] }, level: 1, style: 'heatmap' };
}

// These checks provide immediate feedback. Go independently validates everything.
export function validateRequest(request: PlotDraft): PlotRequest {
  if (!Object.hasOwn(MODES, request.mode)) throw new Error('Choose a supported plot type.');
  const config = MODES[request.mode];
  const resolution = Number(request.resolution);
  if (!Number.isInteger(resolution) || resolution < config.min || resolution > config.max) {
    throw new Error(`Samples per axis must be an integer from ${config.min} to ${config.max}.`);
  }
  const bounds: PlotRequest['bounds'] = {};
  for (const axis of config.axes) {
    const pair = request.bounds[axis];
    if (!pair || pair.some((value) => String(value).trim() === '')) throw new Error(`Enter both ${axis} bounds.`);
    const [min, max] = pair.map(Number);
    if (!Number.isFinite(min) || !Number.isFinite(max) || min >= max || !Number.isFinite(max - min)) throw new Error(`${axis} minimum must be finite and smaller than its maximum.`);
    bounds[axis] = [min, max];
  }
  const level = request.mode === 'iso' ? Number(request.level) : 0;
  if (request.mode === 'iso' && (String(request.level).trim() === '' || !Number.isFinite(level))) throw new Error('Enter a finite isosurface level.');
  if (!request.expression.trim() || new TextEncoder().encode(request.expression).length > 1024) throw new Error('Enter an expression of 1–1,024 bytes.');
  return { ...request, resolution, bounds, level };
}