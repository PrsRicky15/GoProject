import { useMemo, useState } from 'react';
import type { Data, Layout } from 'plotly.js';
import type { SimulationChart as Chart } from '../lib/simulation';
import PlotCanvas from './PlotCanvas';

export default function SimulationChart({ chart }: { chart: Chart }) {
  const [status, setStatus] = useState({ busy: false, error: '' });
  const figure = useMemo(() => {
    const data: Data[] = chart.z
      ? [{ type: 'heatmap', x: chart.x, y: chart.y, z: chart.z, colorscale: 'Viridis', colorbar: { thickness: 10, title: { text: '|ψ|²' } } }]
      : (chart.series ?? []).map(s => ({ type: 'scatter', mode: 'lines', name: s.name, x: s.x, y: s.y, line: { width: 2 } }));
    const layout: Partial<Layout> = { autosize: true, margin: { l: 60, r: chart.z ? 50 : 20, t: 20, b: 55 }, paper_bgcolor: '#fff', plot_bgcolor: '#fff', colorway: ['#7457e8', '#21a599', '#e59437', '#d35a87'], font: { family: 'system-ui', size: 11, color: '#676078' }, xaxis: { title: { text: chart.xLabel }, gridcolor: '#eeebf3', automargin: true }, yaxis: { title: { text: chart.yLabel }, gridcolor: '#eeebf3', automargin: true }, legend: { orientation: 'h', y: 1.12 }, uirevision: chart.title };
    return { data, layout };
  }, [chart]);
  return <section className="panel simulation-chart"><h3>{chart.title}</h3><div className="simulation-plot"><PlotCanvas figure={figure} label={chart.title} onStatus={setStatus} />{status.busy && <span className="chart-status" role="status">Drawing…</span>}{status.error && <p className="error" role="alert">{status.error}</p>}</div></section>;
}
