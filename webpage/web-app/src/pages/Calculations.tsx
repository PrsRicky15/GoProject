import { useEffect, useRef, useState } from 'react';
import type { FormEvent } from 'react';
import { Atom, FlaskConical, Play, RotateCcw } from 'lucide-react';
import { CALCULATIONS, CALCULATION_KEYS, runSimulation, simulationDefaults } from '../lib/simulation';
import type { CalculationKind, SimulationRequest, SimulationResult } from '../lib/simulation';
import SimulationChart from '../components/SimulationChart';
import ParticlePlayback from '../components/ParticlePlayback';
import './workbench.css';
import './calculations.css';

type NumberKey = { [K in keyof SimulationRequest]: SimulationRequest[K] extends number ? K : never }[keyof SimulationRequest];

export default function Calculations() {
  const [draft, setDraft] = useState(() => simulationDefaults('spectrum'));
  const [request, setRequest] = useState<SimulationRequest | null>(null);
  const [result, setResult] = useState<SimulationResult | null>(null);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState('');
  const controller = useRef<AbortController | null>(null);
  const quantum = draft.kind === 'quantum' || draft.kind === 'spectrum';
  const md = draft.kind === 'md';
  const evolving = draft.kind !== 'spectrum';
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
            {!md && <>
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
              <div className="field"><label htmlFor="sim-method">Propagator</label><select id="sim-method" disabled={md} value={draft.method} onChange={e => update({ method: e.target.value, ...(e.target.value === 'ralston3' || e.target.value === 'nystrom5' ? { dt: .002 } : {}) })}>{draft.kind === 'quantum' ? <><option value="split">Midpoint split-operator (2nd order)</option><option value="rk4">Runge–Kutta (4th order)</option><option value="ralston3">Ralston (3rd order)</option><option value="nystrom5">Nyström (5th order)</option></> : <><option value="verlet">Velocity-Verlet (2nd order)</option>{!md && <><option value="leapfrog">Leapfrog (2nd order)</option><option value="yoshida">Yoshida (4th order)</option></>}</>}</select></div>
              <div className="parameter-pair">{number('dt', 'Time step dt', .00001, md ? .01 : .05)}{number('steps', 'Integration steps', 1, 2000, 1)}</div>
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
          <details className="math-help"><summary>Methods & scope</summary><p>Models use atomic units, except the Lennard–Jones fluid, which uses reduced units. The Hamiltonian and classical integrators are computed in Go.</p><p>The verified choices above include completed Ralston/Nyström grid steps and a new midpoint split propagator. Legacy adaptive and implicit solvers are not exposed here until their convergence behavior is verified.</p><p>Bond vibration is a one-dimensional reduced-mass model. The particle mode is a two-dimensional monatomic fluid, not a chemical force field.</p></details>
        </section>
        <div className="calculation-results" aria-live="polite" aria-busy={busy}>
          {busy && <div className="panel run-status" role="status"><span className="loading-dot" /> Go is computing the requested model…</div>}
          {!result && <section className="panel calculation-empty"><Atom size={46} /><h2>Your model, made visible.</h2><p>Choose a calculation and its parameters, then run it. Results include plots and numerical diagnostics.</p><div><span>Quantum states</span><span>Time evolution</span><span>Particle trajectories</span></div></section>}
          {result && <>
            <section className="panel result-heading"><div><span className="eyebrow">02 / RESULTS</span><h2>{result.title}</h2><p>{result.units}</p></div><span className="tag">{JSON.stringify(result.request) === JSON.stringify(draft) ? 'Completed' : 'Previous run · settings changed'}</span></section>
            <div className="metrics-grid">{result.metrics.map(metric => <div className="panel metric" key={metric.label}><span>{metric.label}</span><strong>{metric.value.toPrecision(7)}</strong></div>)}</div>
            {result.frames && <ParticlePlayback key={JSON.stringify(result.request)} frames={result.frames} box={result.request.box} />}
            <div className="result-charts">{result.charts.map((chart, index) => <SimulationChart key={index} chart={chart} />)}</div>
            <section className="panel result-notes"><h3>How to interpret this run</h3>{result.notes.map(note => <p key={note}>{note}</p>)}<details><summary>Parameters used for this result</summary><pre>{JSON.stringify(result.request, null, 2)}</pre></details></section>
          </>}
        </div>
      </div>
      <footer className="workspace-footer"><span>QuantumMLD / Model calculations</span><span>Go solvers · TypeScript interface · Explicit units and diagnostics</span></footer>
    </main>
  </div>;
}
