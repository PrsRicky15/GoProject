import { memo, useEffect, useRef } from 'react';
import type * as Plotly from 'plotly.js';
import type { PlotResult, RenderStatus } from '../lib/types';
import { makeFigure } from '../lib/figure';

const loadPlotly = (spatial: boolean) => spatial
  ? import('plotly.js-gl3d-dist-min').then((module) => module.default)
  : import('plotly.js-cartesian-dist-min').then((module) => module.default);

interface Controller { node: HTMLDivElement; alive: boolean; ready: boolean; queue: Promise<void>; plotly: typeof Plotly | null }
type Figure = { data: Plotly.Data[]; layout: Partial<Plotly.Layout> };
type Props = ({ result: PlotResult; figure?: never; label?: never } | { result?: never; figure: Figure; label: string }) & { onStatus: (status: RenderStatus) => void };
export default memo(function PlotCanvas({ result, figure: suppliedFigure, label, onStatus }: Props) {
  const host = useRef<HTMLDivElement>(null);
  const controller = useRef<Controller | null>(null);
  useEffect(() => {
    const node = document.createElement('div');
    node.className = 'plot-canvas';
    if (!host.current) return;
    host.current.appendChild(node);
    const state: Controller = { node, alive: true, ready: false, queue: Promise.resolve(), plotly: null };
    controller.current = state;
    let frame = 0;
    const observer = new ResizeObserver(() => {
      cancelAnimationFrame(frame);
      frame = requestAnimationFrame(() => {
        if (state.alive && state.plotly && state.ready) void Promise.resolve(state.plotly.Plots.resize(node)).catch(() => {});
      });
    });
    observer.observe(host.current);
    return () => {
      state.alive = false;
      observer.disconnect();
      cancelAnimationFrame(frame);
      node.remove();
      state.queue.finally(() => state.plotly?.purge(node));
    };
  }, []);

  useEffect(() => {
    const state = controller.current;
    if (!state) return;
    let current = true;
    // Serialize Plotly updates; cancelled renders cannot overwrite a newer result.
    state.queue = state.queue.then(async () => {
      if (!current || !state.alive) return;
      onStatus({ busy: true, error: '' });
      const spatial = result ? ['surface', 'iso'].includes(result.settings.mode) : false;
      const plotly = await loadPlotly(spatial);
      if (!current || !state.alive) return;
      state.plotly = plotly;
      if (spatial) {
        const canvas = document.createElement('canvas');
        const context = canvas.getContext('webgl2') || canvas.getContext('webgl');
        if (!context) throw new Error('3D plots need WebGL. Enable browser graphics acceleration, or use the 1D and 2D modes.');
        context.getExtension('WEBGL_lose_context')?.loseContext();
      }
      const figure = result ? makeFigure(result) : suppliedFigure!;
      await plotly.react(state.node, figure.data, figure.layout, {
        responsive: false, displaylogo: false, scrollZoom: false,
        toImageButtonOptions: { format: 'png', filename: `quantummld-${result?.settings.mode ?? 'calculation'}`, scale: 2 },
        modeBarButtonsToRemove: ['sendDataToCloud'],
      });
      state.ready = true;
      if (current && state.alive) onStatus({ busy: false, error: '' });
    }).catch((error) => {
      if (current && state.alive) onStatus({ busy: false, error: error.message || 'Unable to render this plot.' });
    });
    return () => { current = false; };
  }, [result, suppliedFigure, onStatus]);

  return <div className="plot-host" ref={host} role="img" aria-label={result ? `${result.settings.mode} plot of ${result.settings.expression}` : label} />;
});
