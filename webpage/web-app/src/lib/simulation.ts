export type CalculationKind = 'spectrum' | 'quantum' | 'classical' | 'molecular' | 'md';
export interface SimulationRequest {
  kind: CalculationKind; basis: string; method: string; potential: string;
  mass: number; strength: number; alpha: number; center: number; halfWidth: number;
  points: number; states: number; dt: number; steps: number; x0: number; p0: number;
  sigma: number; driveAmplitude: number; driveFrequency: number;
  particles: number; box: number; temperature: number; seed: number;
}
export interface SimulationChart {
  title: string; xLabel: string; yLabel: string;
  series?: { name: string; x: number[]; y: number[] }[];
  x?: number[]; y?: number[]; z?: number[][];
}
export interface ParticleFrame { time: number; x: number[]; y: number[] }
export interface SimulationResult {
  request: SimulationRequest; title: string; units: string;
  metrics: { label: string; value: number }[]; charts: SimulationChart[]; notes: string[];
  frames?: ParticleFrame[];
}
export const CALCULATIONS: Record<CalculationKind, { name: string; description: string }> = {
  spectrum: { name: 'Energy levels', description: 'Diagonalize a 1D Hamiltonian and inspect its eigenstate densities.' },
  quantum: { name: 'Quantum dynamics', description: 'Propagate a Gaussian wave packet, with an optional oscillating field.' },
  classical: { name: 'Classical dynamics', description: 'Explore a particle trajectory, phase space, and conserved energy.' },
  molecular: { name: 'Bond vibration', description: 'Follow a diatomic Morse bond using a reduced mass.' },
  md: { name: 'Molecular dynamics', description: 'Simulate a 2D Lennard–Jones model fluid with periodic boundaries.' },
};
export const CALCULATION_KEYS = Object.keys(CALCULATIONS) as CalculationKind[];
export function simulationDefaults(kind: CalculationKind): SimulationRequest {
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
  return { request, title: value.title, units: value.units, metrics, charts, notes: value.notes, frames };
}

export async function runSimulation(request: SimulationRequest, signal: AbortSignal): Promise<SimulationResult> {
  const response = await fetch('/api/simulations', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(request), signal });
  let payload: unknown;
  try { payload = await response.json(); } catch { throw new Error('Calculation service is unavailable. Restart the Go server with the latest code.'); }
  if (!response.ok) throw new Error(object(payload) && typeof payload.error === 'string' ? payload.error : 'Calculation failed.');
  return decodeSimulation(payload, request);
}
