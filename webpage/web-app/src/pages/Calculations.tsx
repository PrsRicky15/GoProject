import { useEffect, useRef, useState } from 'react';
import type { FormEvent } from 'react';
import { Atom, FlaskConical, Play, RotateCcw } from 'lucide-react';
import { CALCULATIONS, CALCULATION_KEYS, runSimulation, simulationDefaults, RK_METHODS, ADAPTIVE_METHODS, IMPLICIT_METHODS, isAdaptiveMethod, isImplicitMethod } from '../lib/simulation';
import type { CalculationKind, SimulationRequest, SimulationResult } from '../lib/simulation';
import SimulationChart from '../components/SimulationChart';
import ParticlePlayback from '../components/ParticlePlayback';
import NDControls from '../components/NDControls';
import DensityVolume from '../components/DensityVolume';
import './workbench.css';
import './calculations.css';

type NumberKey = { [K in keyof SimulationRequest]-?: SimulationRequest[K] extends number ? K : never }[keyof SimulationRequest];
const PROPAGATOR_GROUPS = [['Fixed-step Runge–Kutta', RK_METHODS], ['Adaptive Runge–Kutta', ADAPTIVE_METHODS], ['Implicit methods', IMPLICIT_METHODS]] as const;

export default function Calculations() {
  const [draft, setDraft] = useState(() => simulationDefaults('spectrum'));
  const [request, setRequest] = useState<SimulationRequest | null>(null);
  const [result, setResult] = useState<SimulationResult | null>(null);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState('');
  const controller = useRef<AbortController | null>(null);
  const quantum = draft.kind === 'quantum' || draft.kind === 'spectrum';
  const md = draft.kind === 'md';
  const nd = draft.kind === 'nd-spectrum';
  const evolving = draft.kind !== 'spectrum' && !nd;
  const adaptive = evolving && isAdaptiveMethod(draft.method);
  const implicit = evolving && isImplicitMethod(draft.method);
  const update = (patch: Partial<SimulationRequest>) => setDraft(d => ({ ...d, ...patch }));
  useEffect(() => {
    if (!request) return;
    const current = new AbortController(); controller.current = current;
    const timer = window.setTimeout(() => { current.abort(); setBusy(false); setError('Calculation timed out. Reduce the grid or steps and run again.'); }, 20000);
    runSimulation(request, current.signal).then(data => {
      if (!current.signal.aborted) { setResult(data); setBusy(false); setError(''); }
    }).catch((problem: unknown) => {
      if (!current.signal.aborted) { setBusy(false); setError(problem instanceof Error ? problem.message : 'Cannot reach the calculation service.'); }
    }).finally(() => window.clearTimeout(timer));
    return () => { current.abort(); window.clearTimeout(timer); };
  }, [request]);
  function submit(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    if (Object.values(draft).some(v => typeof v === 'number' && !Number.isFinite(v))) { setError('Complete all numerical fields.'); return; }
    if (draft.nd && (draft.nd.axes.some(axis => Object.values(axis).some(v => !Number.isFinite(v))) || Object.values(draft.nd).some(v => typeof v === 'number' && !Number.isFinite(v)))) { setError('Complete all coordinate and potential fields.'); return; }
    if (draft.nd && draft.nd.axes.reduce((n, axis) => n * axis.points, 1) > 8192) { setError('Reduce the grid to at most 8,192 total points.'); return; }
    setError(''); setBusy(true); setRequest({ ...draft });
  }
  function switchKind(kind: CalculationKind) { setDraft(simulationDefaults(kind)); setError(''); }
  function number(key: NumberKey, label: string, min: number, max: number, step: number | 'any' = 'any') {
    return <div className="field" key={key}><label htmlFor={`sim-${key}`}>{label}</label><input id={`sim-${key}`} type="number" required min={min} max={max} step={step} value={Number.isFinite(draft[key]) ? draft[key] : ''} onChange={e => update({ [key]: e.target.valueAsNumber })} /></div>;
  }
  return <div className="workbench calculation-workspace">
    <header className="workspace-header"><a className="brand" href="#home"><Atom size={27} /><span>Quantum<span>MLD</span></span></a><nav aria-label="Workspace"><a href="#tools">Plotter & calculator</a><a href="#calculations" className="selected">Calculations</a></nav><span className="local-badge"><i /> Go calculation engine</span></header>
    <main className="workspace-main">
      <div className="workspace-intro"><div><span className="eyebrow">QUANTUM & CLASSICAL LAB</span><h1>Define a model. Follow its dynamics.</h1><p>Energy levels, wave packets, and particle motion — computed by your Go solvers.</p></div><span className="workspace-symbol"><FlaskConical size={30} /></span></div>
      <div className="calculation-tabs" role="group" aria-label="Calculation type">{CALCULATION_KEYS.map(kind => <button type="button" key={kind} aria-pressed={kind === draft.kind} className={kind === draft.kind ? 'active' : ''} onClick={() => switchKind(kind)}>{CALCULATIONS[kind].name}</button>)}</div>
      <div className="calculation-grid">
        <section className="panel calculation-controls"><div className="section-heading"><h2>{CALCULATIONS[draft.kind].name}</h2><span className="step-label">01 / MODEL</span></div><p className="muted">{CALCULATIONS[draft.kind].description}</p>
          <form onSubmit={submit}>
            {nd && <NDControls request={draft} update={update} />}
            {!md && !nd && <>
              <div className="field"><label htmlFor="sim-potential">Potential</label><select id="sim-potential" value={draft.potential} disabled={draft.kind === 'molecular'} onChange={e => update({ potential: e.target.value })}><option value="harmonic">Harmonic oscillator</option><option value="morse">Morse bond</option><option value="double-well">Symmetric double well</option><option value="free">Free particle</option></select></div>
              <p className="model-formula">{draft.potential === 'harmonic' ? 'V(x) = ½ k (x − center)²' : draft.potential === 'morse' ? 'V(x) = Dₑ [1 − exp(−α(x − center))]²' : draft.potential === 'double-well' ? 'V(x) = strength (½ x⁴ − x²)' : 'V(x) = 0'}</p>
              <div className="parameter-pair">{number('mass', draft.kind === 'molecular' ? 'Reduced mass' : 'Mass', .1, 100)}{draft.potential !== 'free' && number('strength', draft.potential === 'harmonic' ? 'Spring constant k' : draft.potential === 'morse' ? 'Well depth Dₑ' : 'Strength', .01, 100)}</div>
              <div className="parameter-pair">{['harmonic', 'morse'].includes(draft.potential) && number('center', 'Potential center', -10, 10)}{draft.potential === 'morse' && number('alpha', 'Morse α', .05, 2)}</div>
            </>}
            {quantum && <>
              <div className="field"><label htmlFor="sim-basis">Hamiltonian basis</label><select id="sim-basis" value={draft.basis} onChange={e => update({ basis: e.target.value })}><option value="dvr">Sinc DVR · symmetric domain</option><option value="fourier">Fourier · periodic</option></select></div>
              <div className="parameter-pair">{number('halfWidth', 'Domain half-width L', 2, 20)}{number('points', 'Grid points N', 32, 192, 1)}</div>
              <p className="field-note">Domain [−L, L). Increase the grid and domain to check convergence.</p>
              {draft.kind === 'spectrum' && number('states', 'Eigenstates to show', 1, 8, 1)}
            </>}
            {evolving && <>
              <div className="field"><label htmlFor="sim-method">Propagator</label><select id="sim-method" disabled={md} value={draft.method} onChange={e => update({ method: e.target.value, ...(draft.kind === 'quantum' && RK_METHODS.some(([id]) => id === e.target.value) ? { dt: Math.min(draft.dt, e.target.value === 'euler' ? .0002 : .002) } : {}) })}>
                {draft.kind === 'quantum' ? <option value="split">Midpoint split-operator · 2nd order</option> : <><option value="verlet">Velocity-Verlet · 2nd order</option>{!md && <><option value="leapfrog">Leapfrog · 2nd order</option><option value="yoshida">Yoshida · 4th order</option></>}</>}
                {!md && PROPAGATOR_GROUPS.map(([label, methods]) => <optgroup key={label} label={label}>{methods.map(([id, name]) => <option key={id} value={id}>{name}</option>)}</optgroup>)}
              </select></div>
              <div className="parameter-pair">{number('dt', adaptive ? 'Output interval / maximum step' : 'Time step dt', .00001, md ? .01 : .05)}{number('steps', adaptive ? 'Output intervals' : 'Integration steps', 1, 2000, 1)}</div>
              <p className="field-note">Final time: {Number.isFinite(draft.dt * draft.steps) ? (draft.dt * draft.steps).toPrecision(6) : '—'}. {adaptive && 'Go may take smaller internal steps between outputs.'}</p>
              {(adaptive || implicit) && <fieldset className="drive-fields"><legend>{adaptive ? 'Adaptive error control' : 'Implicit solve accuracy'}</legend><div className="parameter-pair">{number('absoluteTolerance', 'Absolute tolerance', 1e-12, .001)}{number('relativeTolerance', 'Relative tolerance', 1e-12, .001)}</div><p className="field-note">{adaptive ? 'Tighter tolerances usually require more work. Local error control does not guarantee a global error bound.' : 'These tolerances control each equation solve. Reduce dt separately to check time accuracy.'}</p></fieldset>}
              {draft.kind === 'quantum' && draft.method !== 'split' && <p className="field-note">Norm is measured without renormalization. Explicit methods can drift; backward Euler and SDIRK2 can damp the state. Compare smaller time steps.</p>}
              {!md && <div className="parameter-pair">{number('x0', draft.kind === 'molecular' ? 'Initial bond coordinate' : 'Initial position x₀', -10, 10)}{number('p0', 'Initial momentum p₀', -20, 20)}</div>}
            </>}
            {draft.kind === 'quantum' && <>
              {number('sigma', 'Wave-packet width σ', .1, 4)}
              <fieldset className="drive-fields"><legend>Time-dependent field</legend><p className="field-note">Adds A x sin(ωt) to the potential. Set A = 0 for a static Hamiltonian.</p><div className="parameter-pair">{number('driveAmplitude', 'Amplitude A', -5, 5)}{number('driveFrequency', 'Frequency ω', 0, 10)}</div></fieldset>
            </>}
            {md && <><div className="parameter-pair">{number('particles', 'Particle count', 4, 64, 1)}{number('box', 'Box side length', 5, 30)}</div><div className="parameter-pair">{number('temperature', 'Initial temperature', 0, 2)}{number('seed', 'Random seed', 0, 1000000, 1)}</div><p className="field-note">Reduced units. Particles start on a lattice, with zero net momentum. No thermostat is used.</p></>}
            <div className="plot-actions"><button type="submit" disabled={busy} className="primary-button"><Play size={15} />{busy ? 'Calculating…' : 'Run calculation'}</button><button type="button" className="reset-button" title="Reset model" aria-label="Reset model" onClick={() => { setDraft(simulationDefaults(draft.kind)); setError(''); }}><RotateCcw size={17} /></button></div>
            {busy && <button className="cancel-button" type="button" onClick={() => { controller.current?.abort(); setBusy(false); setError('Calculation cancelled.'); }}>Cancel calculation</button>}
            {error && <p className="error" role="alert">{error}</p>}
          </form>
          <details className="math-help"><summary>Methods & scope</summary><p>Models use atomic units, except the Lennard–Jones fluid, which uses reduced units. All propagation runs in Go on the complete coupled state.</p><p>Quantum, classical and bond dynamics offer 13 fixed Runge–Kutta methods, six adaptive pairs and four implicit methods. Adaptive results report accepted/rejected steps and evaluation counts. Implicit quantum solves use matrix-free GMRES. The original split and symplectic methods remain available.</p><p>Wave-packet dynamics is one-dimensional in space; the n-dimensional page calculates stationary states. Bond vibration is a one-dimensional reduced-mass model. The particle mode is a two-dimensional model fluid using velocity-Verlet.</p></details>
        </section>
        <div className="calculation-results" aria-live="polite" aria-busy={busy}>
          {busy && <div className="panel run-status" role="status"><span className="loading-dot" /> Go is computing the requested model…</div>}
          {!result && <section className="panel calculation-empty"><Atom size={46} /><h2>Your model, made visible.</h2><p>Choose a calculation and its parameters, then run it. Results include plots and numerical diagnostics.</p><div><span>Quantum states</span><span>Time evolution</span><span>Particle trajectories</span></div></section>}
          {result && <>
            <section className="panel result-heading"><div><span className="eyebrow">02 / RESULTS</span><h2>{result.title}</h2><p>{result.units}</p></div><span className="tag">{JSON.stringify(result.request) === JSON.stringify(draft) ? 'Completed' : 'Previous run · settings changed'}</span></section>
            <div className="metrics-grid">{result.metrics.map(metric => <div className="panel metric" key={metric.label}><span>{metric.label}</span><strong>{metric.value.toPrecision(7)}</strong></div>)}</div>
            {result.frames && <ParticlePlayback key={JSON.stringify(result.request)} frames={result.frames} box={result.request.box} />}
            <div className="result-charts">{result.charts.map((chart, index) => <SimulationChart key={index} chart={chart} />)}</div>
            {result.volume && <DensityVolume key={JSON.stringify(result.request)} volume={result.volume} />}
            <section className="panel result-notes"><h3>How to interpret this run</h3>{result.notes.map(note => <p key={note}>{note}</p>)}<details><summary>Parameters used for this result</summary><pre>{JSON.stringify(result.request, null, 2)}</pre></details></section>
          </>}
        </div>
      </div>
      <footer className="workspace-footer"><span>QuantumMLD / Model calculations</span><span>Go solvers · TypeScript interface · Explicit units and diagnostics</span></footer>
    </main>
  </div>;
}
