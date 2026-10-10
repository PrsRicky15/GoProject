export type CalculationKind = 'spectrum' | 'nd-spectrum' | 'quantum' | 'classical' | 'molecular' | 'md';
export interface NDAxis { points: number; halfWidth: number; mass: number; center: number; scale: number }
export interface NDRequest { expression?: string; axes: NDAxis[]; coupling: number; softening: number; separation: number; wavenumber: number; order: number; state: number; axisX: number; axisY: number }
export interface DensityVolume { title: string; labels: string[]; x: number[]; y: number[]; z: number[]; values: number[] }
export function ndAxes(dimensions: number): NDAxis[] {
  return Array.from({ length: dimensions }, () => ({ points: dimensions === 2 ? 24 : dimensions === 3 ? 16 : 8, halfWidth: 6, mass: 1, center: 0, scale: 1 }));
}
export interface SimulationRequest {
  nd?: NDRequest;
  kind: CalculationKind; basis: string; method: string; potential: string;
  mass: number; strength: number; alpha: number; center: number; halfWidth: number;
  points: number; states: number; dt: number; steps: number; x0: number; p0: number;
  sigma: number; driveAmplitude: number; driveFrequency: number;
  particles: number; box: number; temperature: number; seed: number;
}
export interface SimulationChart {
  valueLabel?: string;
  title: string; xLabel: string; yLabel: string;
  series?: { name: string; x: number[]; y: number[] }[];
  x?: number[]; y?: number[]; z?: number[][];
}
export interface ParticleFrame { time: number; x: number[]; y: number[] }
export interface SimulationResult {
  volume?: DensityVolume;
  request: SimulationRequest; title: string; units: string;
  metrics: { label: string; value: number }[]; charts: SimulationChart[]; notes: string[];
  frames?: ParticleFrame[];
}
export const CALCULATIONS: Record<CalculationKind, { name: string; description: string }> = {
  spectrum: { name: 'Energy levels', description: 'Diagonalize a 1D Hamiltonian and inspect its eigenstate densities.' },
  'nd-spectrum': { name: 'N-dimensional states', description: 'Solve a coupled 2D–4D Hamiltonian and inspect its energies, density projections and potential.' },
  quantum: { name: 'Quantum dynamics', description: 'Propagate a Gaussian wave packet, with an optional oscillating field.' },
  classical: { name: 'Classical dynamics', description: 'Explore a particle trajectory, phase space, and conserved energy.' },
  molecular: { name: 'Bond vibration', description: 'Follow a diatomic Morse bond using a reduced mass.' },
  md: { name: 'Molecular dynamics', description: 'Simulate a 2D Lennard–Jones model fluid with periodic boundaries.' },
};
export const CALCULATION_KEYS = Object.keys(CALCULATIONS) as CalculationKind[];
export function simulationDefaults(kind: CalculationKind): SimulationRequest {
  if (kind === 'nd-spectrum') return { ...simulationDefaults('spectrum'), kind, basis: 'fourier', sigma: 1.2, alpha: .25, nd: { axes: ndAxes(2), coupling: .2, softening: .7, separation: 3, wavenumber: 1, order: 4, state: 0, axisX: 0, axisY: 1 } };
  return { kind, basis: 'dvr', method: kind === 'quantum' ? 'split' : 'verlet', potential: kind === 'molecular' ? 'morse' : 'harmonic', mass: 1, strength: kind === 'molecular' ? 5 : 1, alpha: .5, center: 0, halfWidth: 8, points: 64, states: 4, dt: kind === 'md' ? .002 : .01, steps: 400, x0: 1, p0: 0, sigma: .7, driveAmplitude: 0, driveFrequency: 1, particles: 16, box: 6, temperature: .2, seed: 42 };
}

function object(value: unknown): value is Record<string, unknown> { return typeof value === 'object' && value !== null && !Array.isArray(value); }
function numbers(value: unknown, limit: number): value is number[] { return Array.isArray(value) && value.length <= limit && value.every(v => typeof v === 'number' && Number.isFinite(v)); }
function text(value: unknown): value is string { return typeof value === 'string' && value.length < 3000; }

export function decodeSimulation(value: unknown, request: SimulationRequest): SimulationResult {
  if (!object(value) || !text(value.title) || !text(value.units) || !Array.isArray(value.metrics) || value.metrics.length > 16 || !Array.isArray(value.charts) || value.charts.length > 6 || !Array.isArray(value.notes) || value.notes.length > 12 || !value.notes.every(text)) throw new Error('Invalid calculation response.');
  const metrics = value.metrics.map(m => { if (!object(m) || !text(m.label) || typeof m.value !== 'number' || !Number.isFinite(m.value)) throw new Error('Invalid calculation diagnostics.'); return { label: m.label, value: m.value }; });
  const charts: SimulationChart[] = value.charts.map(c => {
    if (!object(c) || !text(c.title) || !text(c.xLabel) || !text(c.yLabel)) throw new Error('Invalid chart response.');
    const chart: SimulationChart = { title: c.title, xLabel: c.xLabel, yLabel: c.yLabel };
    if (c.valueLabel !== undefined) { if (!text(c.valueLabel)) throw new Error('Invalid chart units.'); chart.valueLabel = c.valueLabel; }
    if (c.z !== undefined) {
      if (!numbers(c.x, 192) || !numbers(c.y, 152) || !Array.isArray(c.z) || c.z.length !== c.y.length || !c.z.every(row => numbers(row, 192) && row.length === (c.x as number[]).length)) throw new Error('Invalid density grid.');
      chart.x = c.x; chart.y = c.y; chart.z = c.z as number[][];
    } else {
      if (!Array.isArray(c.series) || c.series.length > 8) throw new Error('Invalid result series.');
      chart.series = c.series.map(s => { if (!object(s) || !text(s.name) || !numbers(s.x, 402) || !numbers(s.y, 402) || s.x.length !== s.y.length) throw new Error('Invalid result samples.'); return { name: s.name, x: s.x, y: s.y }; });
    }
    return chart;
  });
  let frames: ParticleFrame[] | undefined;
  if (value.frames !== undefined) {
    if (!Array.isArray(value.frames) || value.frames.length > 152) throw new Error('Invalid particle frames.');
    frames = value.frames.map(f => { if (!object(f) || typeof f.time !== 'number' || !Number.isFinite(f.time) || !numbers(f.x, 64) || !numbers(f.y, 64) || f.x.length !== request.particles || f.y.length !== request.particles) throw new Error('Invalid particle positions.'); return { time: f.time, x: f.x, y: f.y }; });
  }
  let volume: DensityVolume | undefined;
  if (value.volume !== undefined) {
    const v = value.volume;
    const n = request.nd?.axes.slice(0, 3).reduce((product, axis) => product * axis.points, 1);
    if (request.kind !== 'nd-spectrum' || !request.nd || request.nd.axes.length < 3 || !object(v) || !text(v.title) || !Array.isArray(v.labels) || v.labels.length !== 3 || !v.labels.every(text) || !numbers(v.x, 8192) || !numbers(v.y, 8192) || !numbers(v.z, 8192) || !numbers(v.values, 8192) || v.values.some(x => x < 0) || v.x.length !== n || v.y.length !== n || v.z.length !== n || v.values.length !== n) throw new Error('Invalid 3D density data.');
    volume = { title: v.title, labels: v.labels, x: v.x, y: v.y, z: v.z, values: v.values };
  }
  return { request, title: value.title, units: value.units, metrics, charts, notes: value.notes, frames, volume };
}

export async function runSimulation(request: SimulationRequest, signal: AbortSignal): Promise<SimulationResult> {
  const response = await fetch('/api/simulations', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(request), signal });
  let payload: unknown;
  try { payload = await response.json(); } catch { throw new Error('Calculation service is unavailable. Restart the Go server with the latest code.'); }
  if (!response.ok) throw new Error(object(payload) && typeof payload.error === 'string' ? payload.error : 'Calculation failed.');
  return decodeSimulation(payload, request);
}
