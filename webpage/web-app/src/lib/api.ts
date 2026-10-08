import { MODES } from './sampling';
import type { PlotRequest, PlotResult } from './types';

function record(value: unknown): value is Record<string, unknown> {
  return typeof value === 'object' && value !== null && !Array.isArray(value);
}
function finite(value: unknown): value is number { return typeof value === 'number' && Number.isFinite(value); }

// TypeScript types disappear at runtime: verify the API payload before rendering.
export function decodePlotResult(payload: unknown, settings: PlotRequest): PlotResult {
  const axes = MODES[settings.mode].axes;
  const expectedCount = settings.resolution ** axes.length;
  if (!record(payload) || !record(payload.coords) || payload.count !== expectedCount || !Array.isArray(payload.values) || payload.values.length !== expectedCount || !finite(payload.min) || !finite(payload.max)) throw new Error('The plot engine returned invalid data.');
  const coords: PlotResult['coords'] = { x: new Float64Array() };
  for (const axis of axes) {
    const samples: unknown = payload.coords[axis];
    if (!Array.isArray(samples) || samples.length !== settings.resolution || !samples.every(finite)) throw new Error('The plot engine returned invalid coordinates.');
    coords[axis] = Float64Array.from(samples);
  }
  let invalid = 0;
  const values = Float64Array.from(payload.values, (value: unknown) => {
    if (value === null) { invalid++; return NaN; }
    if (!finite(value)) throw new Error('The plot engine returned an invalid sample.');
    return value;
  });
  if (invalid !== payload.invalid || (settings.mode === 'iso' && invalid > 0)) throw new Error('The plot engine returned inconsistent data.');
  return { settings, coords, values, count: expectedCount, invalid, min: payload.min, max: payload.max };
}

export async function requestPlot(settings: PlotRequest, signal: AbortSignal): Promise<PlotResult> {
  let response: Response;
  try {
    response = await fetch('/api/plot', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(settings), signal });
  } catch (error) {
    if (signal.aborted) throw error;
    throw new Error('Cannot reach the Go plot engine. Start the Go server and try again.');
  }
  let payload: unknown;
  try { payload = await response.json(); }
  catch { throw new Error('The Go plot engine is unavailable. Check that the Go server is running.'); }
  if (!response.ok) throw new Error(record(payload) && typeof payload.error === 'string' ? payload.error : 'The Go plot engine could not complete this request.');
  return decodePlotResult(payload, settings);
}