import type { Data, Layout } from 'plotly.js';
import type { PlotResult } from './types';

// @types/plotly.js omits the isosurface fields and types `value` as a scalar.
// Keep those fields checked locally and bridge that upstream gap only here.
interface IsosurfaceTrace {
  type: 'isosurface'; x: Float64Array; y: Float64Array; z: Float64Array;
  value: Float64Array; isomin: number; isomax: number;
  surface: { count: number }; caps: Record<'x' | 'y' | 'z', { show: boolean }>;
  colorscale: string; showscale: boolean;
}

export function makeFigure(result: PlotResult): { data: Data[]; layout: Partial<Layout> } {
  const { settings, coords, values, count } = result;
  const { mode, resolution: n } = settings;
  const common = { colorscale: 'Viridis', colorbar: { title: { text: 'f' }, thickness: 12 }, showscale: true };
  let trace: Data;
  if (mode === 'line') {
    // Explicit holes keep undefined samples from being joined across a pole.
    trace = { type: 'scatter', mode: 'lines', x: coords.x,
      y: values,
      connectgaps: false, line: { color: '#7457e8', width: 2.5 }, name: 'f(x)' };
  } else if (mode === 'iso') {
    const x = new Float64Array(count), y = new Float64Array(count), z = new Float64Array(count);
    for (let i = 0; i < count; i++) {
      x[i] = coords.x[i % n]; y[i] = coords.y![Math.floor(i / n) % n]; z[i] = coords.z![Math.floor(i / (n * n))];
    }
    const iso: IsosurfaceTrace = { ...common, type: 'isosurface', x, y, z, value: values,
      isomin: settings.level, isomax: settings.level, surface: { count: 1 },
      caps: { x: { show: false }, y: { show: false }, z: { show: false } }, showscale: false };
    trace = iso as unknown as Data;
  } else {
    const z = Array.from({ length: n }, (_, row) => Array.from(values.subarray(row * n, (row + 1) * n)));
    trace = { ...common, type: mode === 'surface' ? 'surface' : settings.style,
      x: coords.x, y: coords.y, z, connectgaps: false };
  }
  const axis = (title: string) => ({ title: { text: title }, gridcolor: '#e8e8f0', zerolinecolor: '#b4b3c8', automargin: true });
  const spatial = mode === 'surface' || mode === 'iso';
  return { data: [trace], layout: {
    autosize: true, paper_bgcolor: '#ffffff', plot_bgcolor: '#ffffff',
    font: { family: 'system-ui, sans-serif', color: '#57566e', size: 12 },
    margin: spatial ? { l: 8, r: 24, t: 30, b: 8 } : { l: 56, r: 36, t: 36, b: 50 },
    xaxis: axis('x'), yaxis: { ...axis(mode === 'line' ? 'f(x)' : 'y'), ...(mode === 'heatmap' ? { scaleanchor: 'x', scaleratio: 1 } : {}) },
    scene: { xaxis: axis('x'), yaxis: axis('y'), zaxis: axis(mode === 'surface' ? 'f(x, y)' : 'z'), aspectmode: mode === 'iso' ? 'data' : 'auto' },
    uirevision: `${mode}:${JSON.stringify(settings.bounds)}`, showlegend: false,
  } };
}
