import { ndAxes } from '../lib/simulation';
import type { NDAxis, NDRequest, SimulationRequest } from '../lib/simulation';

const ND_POTENTIALS = {
  harmonic: { name: 'Harmonic oscillator', formula: 'V = ½ k Σ qᵢ²', strength: 1 },
  'coupled-harmonic': { name: 'Coupled harmonic oscillator', formula: 'V = k (½ Σ qᵢ² + c Σ qᵢ qᵢ₊₁)', strength: 1 },
  morse: { name: 'Morse · sum over coordinates', formula: 'V = Dₑ Σ [1 − exp(−α qᵢ)]²', strength: 2 },
  'double-well': { name: 'Quartic double wells', formula: 'V = ¼ k Σ (qᵢ² − 1)²', strength: .2 },
  gaussian: { name: 'Gaussian well', formula: 'V = −D exp(−r² / 2σ²)', strength: 5 },
  'super-gaussian': { name: 'Super-Gaussian well', formula: 'V = −D exp[−(r / σ)ᵖ]', strength: 5 },
  'soft-core': { name: 'Soft-core attraction', formula: 'V = −Z / √(r² + a²)', strength: 3 },
  'double-gaussian': { name: 'Two Gaussian wells', formula: 'V = −D Σ± exp[−|q ± (s/2)e₁|² / 2σ²]', strength: 4 },
  lattice: { name: 'Sinusoidal lattice', formula: 'V = D Σ sin²(κ qᵢ)', strength: 2 },
  free: { name: 'Free particle', formula: 'V = 0', strength: 1 },
  custom: { name: 'Custom formula', formula: 'Enter V(x1, …, xN)', strength: 1 },
} as const;

export default function NDControls({ request: r, update }: { request: SimulationRequest; update: (patch: Partial<SimulationRequest>) => void }) {
  const n = r.nd!;
  const model = ND_POTENTIALS[r.potential as keyof typeof ND_POTENTIALS];
  const custom = r.potential === 'custom';
  const radiusSquared = n.axes.map((_, d) => `x${d + 1}^2`).join(' + ');
  const examples = [
    { name: 'Coupled oscillator', expression: `0.5*(${radiusSquared}) + 0.2*x1*x2` },
    { name: 'Gaussian well', expression: `-5*exp(-(${radiusSquared})/(2*1.2^2))` },
    { name: 'Anharmonic', expression: `0.5*(${radiusSquared}) + 0.1*x1^4` },
  ];
  const set = (patch: Partial<NDRequest>) => update({ nd: { ...n, ...patch } });
  const total = n.axes.reduce((size, axis) => size * axis.points, 1);
  function field(id: string, label: string, value: number, min: number, max: number, change: (value: number) => void, step: number | 'any' = 'any') {
    return <div className="field"><label htmlFor={id}>{label}</label><input id={id} type="number" required min={min} max={max} step={step} value={Number.isFinite(value) ? value : ''} onChange={e => change(e.target.valueAsNumber)} /></div>;
  }
  function axisField(d: number, key: keyof NDAxis, label: string, min: number, max: number, step: number | 'any' = 'any') {
    return field(`nd-axis-${d}-${key}`, `${label} · x${d + 1}`, n.axes[d][key], min, max, value => set({ axes: n.axes.map((axis, i) => i === d ? { ...axis, [key]: value } : axis) }), step);
  }
  return <>
    <div className="parameter-pair">
      <div className="field"><label htmlFor="nd-dimensions">Coordinates</label><select id="nd-dimensions" value={n.axes.length} onChange={e => set({ axes: ndAxes(Number(e.target.value)), axisX: 0, axisY: 1 })}>{[2, 3, 4].map(d => <option key={d} value={d}>{d} dimensions</option>)}</select></div>
      <div className="field"><label htmlFor="nd-basis">Hamiltonian basis</label><select id="nd-basis" value={r.basis} onChange={e => update({ basis: e.target.value })}><option value="fourier">Fourier · periodic</option><option value="dvr">Sinc DVR</option></select></div>
    </div>
    <div className="field"><label htmlFor="nd-potential">N-dimensional potential</label><select id="nd-potential" value={r.potential} onChange={e => { const key = e.target.value as keyof typeof ND_POTENTIALS; update({ potential: key, strength: ND_POTENTIALS[key].strength, sigma: key === 'super-gaussian' ? 2 : 1.2, ...(key === 'custom' && n.expression === undefined ? { nd: { ...n, expression: examples[0].expression } } : {}) }); }}>{Object.entries(ND_POTENTIALS).map(([key, value]) => <option key={key} value={key}>{value.name}</option>)}</select></div>
    {custom ? <>
      <div className="field"><label htmlFor="nd-expression">Potential formula V</label><textarea id="nd-expression" rows={4} required maxLength={1024} spellCheck={false} aria-describedby="nd-expression-help" value={n.expression ?? ''} onChange={e => set({ expression: e.target.value })} /></div>
      <p className="field-note" id="nd-expression-help">Use x1 through x{n.axes.length}, the actual grid coordinates. Enter only the formula, without “V =”. Use * to multiply and ^ for powers. Go evaluates the formula and solves the Hamiltonian.</p>
      <div className="formula-examples" aria-label="Example potentials">{examples.map(example => <button type="button" key={example.name} onClick={() => set({ expression: example.expression })}>{example.name}</button>)}</div>
      <details className="formula-help"><summary>Functions and coordinate settings</summary><p>Functions: sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, sqrt, abs, exp, ln, log, log10, floor, ceil, round, pow, min, max. Constants: pi and e. ln is natural logarithm; log and log10 use base 10. pow, min and max take two arguments.</p><p>Write coefficients and coordinate shifts directly in the formula. Preset strength and coordinate scales do not apply. Slice positions fix hidden coordinates only in the potential plot. Changing dimensions keeps your formula; update its coordinates as needed.</p><p>Maximum 1,024 bytes and 256 tokens. Every sampled potential value must be finite and satisfy |V| ≤ 100,000, in atomic units.</p></details>
    </> : <p className="model-formula">{model.formula}<br />qᵢ = scaleᵢ (xᵢ − centerᵢ), r² = Σ qᵢ²</p>}
    {!custom && r.potential !== 'free' && field('nd-strength', r.potential.includes('harmonic') ? 'Spring constant k' : 'Strength / depth', r.strength, .01, 30, strength => update({ strength }))}
    {r.potential === 'coupled-harmonic' && field('nd-coupling', 'Neighbor coupling c', n.coupling, -.45, .45, coupling => set({ coupling }))}
    {r.potential === 'morse' && field('nd-alpha', 'Morse α', r.alpha, .05, 1, alpha => update({ alpha }))}
    {['gaussian', 'super-gaussian', 'double-gaussian'].includes(r.potential) && field('nd-sigma', 'Well width σ', r.sigma, .3, 4, sigma => update({ sigma }))}
    {r.potential === 'super-gaussian' && <div className="field"><label htmlFor="nd-order">Super-Gaussian order p</label><select id="nd-order" value={n.order} onChange={e => set({ order: Number(e.target.value) })}>{[2, 4, 6, 8].map(p => <option key={p}>{p}</option>)}</select></div>}
    {r.potential === 'soft-core' && field('nd-softening', 'Softening a', n.softening, .2, 3, softening => set({ softening }))}
    {r.potential === 'double-gaussian' && field('nd-separation', 'Well separation s (q₁)', n.separation, .2, 6, separation => set({ separation }))}
    {r.potential === 'lattice' && field('nd-wavenumber', 'Lattice wavenumber κ', n.wavenumber, .1, 3, wavenumber => set({ wavenumber }))}
    <p className={`nd-grid-size ${total > 8192 ? 'error' : ''}`}>{Number.isFinite(total) ? total.toLocaleString() : '—'} tensor points / 8,192 maximum</p>
    <p className="field-note">Changing dimensions resets the grid. Each axis spans [−L, L). Increase resolution to check convergence.</p>
    {n.axes.map((_, d) => <fieldset className="nd-axis" key={d}><legend>Coordinate x{d + 1}</legend><div className="parameter-pair">{axisField(d, 'points', 'Points', 4, 64, 1)}{axisField(d, 'halfWidth', 'Half-width L', 2, 12)}</div><div className="parameter-pair">{axisField(d, 'mass', 'Mass', .1, 10)}{!custom && axisField(d, 'scale', 'Coordinate scale', .25, 2)}</div>{axisField(d, 'center', custom ? 'Hidden slice position' : 'Center', -n.axes[d].halfWidth / 2, n.axes[d].halfWidth / 2)}</fieldset>)}
    <div className="parameter-pair">{field('nd-states', 'Eigenstates', r.states, 1, 4, states => update({ states, nd: { ...n, state: Math.min(n.state, states - 1) } }), 1)}<div className="field"><label htmlFor="nd-state">Density state</label><select id="nd-state" value={n.state} onChange={e => set({ state: Number(e.target.value) })}>{Array.from({ length: Number.isFinite(r.states) ? Math.max(0, Math.min(4, r.states)) : 0 }, (_, s) => <option key={s} value={s}>{s === 0 ? '0 · Ground state' : `${s} · Excited state`}</option>)}</select></div></div>
    <div className="parameter-pair">{(['axisX', 'axisY'] as const).map((key, index) => <div className="field" key={key}><label htmlFor={`nd-${key}`}>{index ? 'Vertical' : 'Horizontal'} plot axis</label><select id={`nd-${key}`} value={n[key]} onChange={e => { const value = Number(e.target.value); const other = key === 'axisX' ? 'axisY' : 'axisX'; set({ [key]: value, ...(value === n[other] ? { [other]: n[key] } : {}) }); }}>{n.axes.map((_, d) => <option key={d} value={d}>x{d + 1}</option>)}</select></div>)}</div>
    <p className="field-note">Densities integrate over hidden coordinates. For 3D/4D, the isosurface shows the density in x1, x2, x3.</p>
  </>;
}
