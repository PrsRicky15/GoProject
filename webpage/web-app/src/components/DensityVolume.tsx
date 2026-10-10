import { useId, useMemo, useState } from 'react';
import type { Data, Layout } from 'plotly.js';
import type { DensityVolume as Volume } from '../lib/simulation';
import type { IsosurfaceTrace } from '../lib/figure';
import PlotCanvas from './PlotCanvas';

export default function DensityVolume({ volume }: { volume: Volume }) {
  const [level, setLevel] = useState(20);
  const [status, setStatus] = useState({ busy: false, error: '' });
  const id = useId();
  const peak = useMemo(() => Math.max(...volume.values), [volume]);
  const minimum = useMemo(() => Math.min(...volume.values), [volume]);
  const uniform = peak - minimum < Math.max(1e-14, peak * 1e-8);
  const figure = useMemo(() => {
    const trace: IsosurfaceTrace = { type: 'isosurface', x: volume.x, y: volume.y, z: volume.z, value: volume.values, isomin: peak * level / 100, isomax: peak * .95, surface: { count: 3 }, caps: { x: { show: false }, y: { show: false }, z: { show: false } }, colorscale: 'Viridis', opacity: .55, colorbar: { thickness: 12, title: { text: '|ψ|²' } } };
    // Bridge the upstream Plotly isosurface typing gap after checking our schema.
    const data: Data[] = [trace as unknown as Data];
    const layout: Partial<Layout> = { autosize: true, margin: { l: 10, r: 30, t: 15, b: 10 }, paper_bgcolor: '#fff', scene: { xaxis: { title: { text: volume.labels[0] } }, yaxis: { title: { text: volume.labels[1] } }, zaxis: { title: { text: volume.labels[2] } }, aspectmode: 'data' }, uirevision: volume.title };
    return { data, layout };
  }, [volume, peak, level]);
  return <section className="panel density-volume"><h3>{volume.title}</h3><p className="field-note">Rotate to inspect equal-density surfaces. The threshold is a fraction of peak density, not enclosed probability.</p>{uniform ? <p className="field-note">This state has nearly uniform density, so it has no distinct isosurface. Use the density projections above.</p> : <><label htmlFor={id}>Lowest surface: {level}% of peak</label><input id={id} aria-label="Isosurface threshold" type="range" min={5} max={80} step={5} value={level} onChange={e => setLevel(Number(e.target.value))} /><div className="simulation-plot"><PlotCanvas figure={figure} label={volume.title} onStatus={setStatus} />{status.busy && <span className="chart-status" role="status">Drawing…</span>}{status.error && <p className="error" role="alert">{status.error}</p>}</div></>}</section>;
}
