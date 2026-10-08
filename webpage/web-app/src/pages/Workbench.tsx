import { useEffect, useRef, useState } from 'react';
import type { FormEvent } from 'react';
import type { Mode, PlotDraft, PlotResult, RenderStatus, FieldStyle } from '../lib/types';
import { requestPlot } from '../lib/api';
import { Atom, ArrowLeft, ChartNoAxesCombined, Play, RotateCcw } from 'lucide-react';
import PlotCanvas from '../components/PlotCanvas';
import Calculator from '../components/Calculator';
import { defaults, MODES, MODE_KEYS, validateRequest } from '../lib/sampling';
import './workbench.css';

const PRESETS: Record<Mode, [string, string][]> = {
  line: [['Gaussian', 'exp(-x^2 / 2)'], ['Harmonic potential', '0.5 * x^2'], ['Oscillation', 'sin(2*x) * exp(-x^2/8)']],
  heatmap: [['Interference', 'sin(x) * cos(y)'], ['Gaussian', 'exp(-(x^2+y^2)/2)'], ['Saddle', 'x^2-y^2']],
  surface: [['Gaussian', 'exp(-(x^2+y^2)/2)'], ['Harmonic potential', '0.5*(x^2+y^2)'], ['Ripples', 'cos(sqrt(x^2+y^2))']],
  iso: [['Sphere', 'x^2+y^2+z^2'], ['Ellipsoid', 'x^2/2+y^2+2*z^2'], ['Orbital-like lobes', 'z^2*exp(-(x^2+y^2+z^2))']],
};

export default function Workbench() {
  const [drafts, setDrafts] = useState<Record<Mode, PlotDraft>>(() => ({ line: defaults('line'), heatmap: defaults('heatmap'), surface: defaults('surface'), iso: defaults('iso') }));
  const [mode, setMode] = useState<Mode>('line');
  const [request, setRequest] = useState(() => defaults('line'));
  const [result, setResult] = useState<PlotResult | null>(null);
  const [sampling, setSampling] = useState(true);
  const [error, setError] = useState('');
  const [render, setRender] = useState<RenderStatus>({ busy: false, error: '' });
  const activeRequest = useRef<AbortController | null>(null);
  const draft = drafts[mode];
  const config = MODES[mode];
  const dirty = JSON.stringify(draft) !== JSON.stringify(request);
  const update = (patch: Partial<PlotDraft>) => setDrafts((previous) => ({ ...previous, [mode]: { ...previous[mode], ...patch } }));

  useEffect(() => {
    const controller = new AbortController();
    activeRequest.current = controller;
    const timeout = window.setTimeout(() => {
      controller.abort(); setSampling(false); setError('The Go plot engine took too long. Reduce the resolution and try again.');
    }, 15000);
    requestPlot(validateRequest(request), controller.signal).then((data) => {
      if (controller.signal.aborted) return;
      setResult(data); setSampling(false); setError(''); setRender({ busy: true, error: '' });
    }).catch((problem: unknown) => {
      if (controller.signal.aborted) return;
      setSampling(false); setError(problem instanceof Error ? problem.message : 'Unable to generate the plot.');
    }).finally(() => window.clearTimeout(timeout));
    return () => { window.clearTimeout(timeout); controller.abort(); };
  }, [request]);

  function plot(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    try {
      validateRequest(draft);
      setError(''); setSampling(true); setRender({ busy: false, error: '' });
      setRequest(structuredClone(draft));
    } catch (problem) { setError(problem instanceof Error ? problem.message : 'Invalid plot settings.'); }
  }
  function cancel() {
    activeRequest.current?.abort();
    setSampling(false); setError('Plot cancelled. Change the settings or plot again.');
  }
  function chooseMode(next: Mode) {
    setMode(next); setError('');
  }
  const plotted = result?.settings;
  return <div className="workbench">
    <header className="workspace-header"><a className="brand" href="#home"><Atom size={27} /><span>Quantum<span>MLD</span></span></a><nav aria-label="Workspace"><a href="#home"><ArrowLeft size={15} /> Overview</a><a className="selected" href="#tools">Plotter & calculator</a><a href="#calculations">Calculations</a></nav><span className="local-badge"><i /> Go plot engine</span></header>
    <main className="workspace-main">
      <div className="workspace-intro"><div><span className="eyebrow">THE SCIENTIFIC WORKSPACE</span><h1>From expression to insight.</h1><p>Explore curves, fields, and surfaces. Make the mathematics visible.</p></div><span className="workspace-symbol"><ChartNoAxesCombined size={30} /></span></div>
      <div className="workspace-grid">
        <section className="plot-controls panel" aria-labelledby="plot-settings-title">
          <div className="section-heading"><h2 id="plot-settings-title">Build a plot</h2><span className="step-label">01 / DEFINE</span></div>
          <div className="mode-tabs" role="group" aria-label="Plot type">{MODE_KEYS.map((key) => <button key={key} type="button" aria-pressed={key === mode} className={key === mode ? 'active' : ''} onClick={() => chooseMode(key)}>{MODES[key].label}</button>)}</div>
          <form onSubmit={plot}>
            <label htmlFor="plot-expression">{config.formula}</label>
            <textarea id="plot-expression" value={draft.expression} onChange={(event) => update({ expression: event.target.value })} spellCheck="false" maxLength={1024} rows={2} />
            <div className="presets" aria-label="Example expressions">{PRESETS[mode].map(([label, expression]) => <button type="button" key={label} onClick={() => update({ expression, ...(mode === 'iso' ? { level: label === 'Orbital-like lobes' ? 0.1 : 1 } : {}) })}>{label}</button>)}</div>
            <fieldset className="domain"><legend>Domain</legend><div className="domain-head"><span>Axis</span><span>Minimum</span><span>Maximum</span></div>{config.axes.map((axis) => <div className="domain-row" key={axis}><span>{axis}</span>{[0, 1].map((side) => <input key={side} aria-label={`${axis} ${side ? 'maximum' : 'minimum'}`} type="number" step="any" required value={draft.bounds[axis][side]} onChange={(event) => update({ bounds: { ...draft.bounds, [axis]: [side === 0 ? event.target.value : draft.bounds[axis][0], side === 1 ? event.target.value : draft.bounds[axis][1]] } })} />)}</div>)}</fieldset>
            {mode === 'iso' && <div className="field"><label htmlFor="iso-level">Isosurface level</label><input id="iso-level" type="number" step="any" required value={draft.level} onChange={(event) => update({ level: event.target.value })} /></div>}
            {mode === 'heatmap' && <div className="field"><label htmlFor="field-style">Display</label><select id="field-style" value={draft.style} onChange={(event) => update({ style: event.target.value as FieldStyle })}><option value="heatmap">Heatmap</option><option value="contour">Contours</option></select></div>}
            <div className="field"><label htmlFor="resolution">Samples {mode !== 'line' && 'per axis'}</label><input id="resolution" type="number" min={config.min} max={config.max} step="1" required value={draft.resolution} onChange={(event) => update({ resolution: event.target.value })} /><small>{config.min}–{config.max} per axis · more samples add detail and take longer.</small></div>
            <div className="plot-actions"><button type="submit" className="primary-button"><Play size={15} /> {sampling ? 'Restart plot' : 'Plot expression'}</button><button type="button" className="reset-button" aria-label="Reset plot settings" title="Reset plot settings" onClick={() => { setDrafts((previous) => ({ ...previous, [mode]: defaults(mode) })); setError(''); }}><RotateCcw size={17} /></button></div>
            {sampling && <button type="button" className="cancel-button" onClick={cancel}>Cancel sampling</button>}
            {error && <p className="error" role="alert">{error}</p>}
          </form>
          <details className="math-help"><summary>Expression guide</summary><p>Use explicit multiplication: <code>2*x</code>. Powers: <code>x^2</code>. Constants: <code>pi</code>, <code>e</code>. Functions include sin, cos, tan, exp, sqrt, abs, ln and log (base 10). All angles are radians; values are real numbers.</p><p>These are uniformly sampled plots. Narrow features and discontinuities between samples can be missed. Increase resolution or narrow the domain to investigate.</p></details>
        </section>
        <section className="plot-panel panel" aria-labelledby="plot-title" aria-busy={sampling || render.busy}>
          <div className="plot-heading"><div><span className="eyebrow">02 / EXPLORE</span><h2 id="plot-title">{plotted ? MODES[plotted.mode].label : config.label}</h2></div><span className="tag">{dirty ? 'Unplotted changes' : 'Interactive view'}</span></div>
          <div className="plot-expression-label">{plotted ? `${MODES[plotted.mode].formula.split(' = ')[0]} = ${plotted.expression}${plotted.mode === 'iso' ? ` = ${plotted.level}` : ''}` : 'Preparing your visualization'}</div>
          <div className="plot-stage">
            {result && <PlotCanvas key={['surface', 'iso'].includes(result.settings.mode) ? '3d' : '2d'} result={result} onStatus={setRender} />}
            {(sampling || render.busy) && <div className="plot-message" role="status"><span className="loading-dot" />{sampling ? 'Sampling expression…' : 'Drawing plot…'}</div>}
            {render.error && <div className="plot-message error" role="alert">{render.error}</div>}
            {!result && !sampling && <div className="plot-message">Adjust the settings, then choose Plot expression.</div>}
          </div>
          <div className="plot-footer"><span>{result ? `${result.count.toLocaleString()} samples${result.invalid ? ` · ${result.invalid.toLocaleString()} undefined (shown as gaps)` : ''}` : 'Ready for your next idea'}</span><span>Drag to explore · toolbar to zoom or save PNG</span></div>
        </section>
        <aside className="panel calculator-container"><Calculator /></aside>
      </div>
      <footer className="workspace-footer"><span>QuantumMLD / Math tools</span><span>Plot sampling in Go · Calculator in TypeScript</span></footer>
    </main>
  </div>;
}
