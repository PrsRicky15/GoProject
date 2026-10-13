import { describe, expect, test } from 'bun:test';
import { decodeSimulation, simulationDefaults, ndAxes, RK_METHODS, ADAPTIVE_METHODS, IMPLICIT_METHODS, isAdaptiveMethod, isImplicitMethod } from './simulation';

describe('calculation response validation', () => {
 test('offers the complete vector method catalog with correct controls', () => {
  expect(RK_METHODS).toHaveLength(13); expect(ADAPTIVE_METHODS).toHaveLength(6); expect(IMPLICIT_METHODS).toHaveLength(4);
  expect(new Set([...RK_METHODS, ...ADAPTIVE_METHODS, ...IMPLICIT_METHODS].map(([id]) => id)).size).toBe(23);
  for (const [id] of ADAPTIVE_METHODS) { expect(isAdaptiveMethod(id)).toBe(true); expect(isImplicitMethod(id)).toBe(false); }
  for (const [id] of IMPLICIT_METHODS) { expect(isImplicitMethod(id)).toBe(true); expect(isAdaptiveMethod(id)).toBe(false); }
  for (const [id] of RK_METHODS) { expect(isImplicitMethod(id)).toBe(false); expect(isAdaptiveMethod(id)).toBe(false); }
  const defaults = simulationDefaults('quantum'); expect(defaults.absoluteTolerance).toBe(1e-9); expect(defaults.relativeTolerance).toBe(1e-7);
 });
 const request = simulationDefaults('quantum');
 const response = { title: 'Quantum dynamics', units: 'Atomic units', metrics: [{ label: 'Norm', value: 1 }], notes: [], charts: [{ title: 'Density', xLabel: 'x', yLabel: 't', x: [0, 1], y: [0], z: [[1, 0]] }] };
 test('retains density orientation and the submitted settings', () => {
  const result = decodeSimulation(response, request);
  expect(result.charts[0].z).toEqual([[1, 0]]);
  expect(result.request).toEqual(request);
 });
 test('rejects non-finite diagnostics and inconsistent density rows', () => {
  expect(() => decodeSimulation({ ...response, metrics: [{ label: 'Norm', value: Infinity }] }, request)).toThrow();
  expect(() => decodeSimulation({ ...response, charts: [{ ...response.charts[0], z: [[1]] }] }, request)).toThrow();
  expect(() => decodeSimulation({ ...response, charts: [{ ...response.charts[0], z: [[1, null]] }] }, request)).toThrow();
 });
 test('rejects particle frames that do not match the requested system', () => {
  expect(() => decodeSimulation({ ...response, frames: [{ time: 0, x: [1], y: [1] }] }, simulationDefaults('md'))).toThrow();
 });
 test('validates multidimensional volume lengths and keeps potential units', () => {
  const ndRequest = simulationDefaults('nd-spectrum');
  ndRequest.nd = { ...ndRequest.nd!, axes: ndAxes(3).map(axis => ({ ...axis, points: 4 })) };
  const volume = { title: '3D density', labels: ['x1', 'x2', 'x3'], x: Array(64).fill(0), y: Array(64).fill(0), z: Array(64).fill(0), values: Array(64).fill(.1) };
  const ndResponse = { ...response, volume, charts: [{ ...response.charts[0], valueLabel: 'V' }] };
  const result = decodeSimulation(ndResponse, ndRequest);
  expect(result.volume?.values.length).toBe(64);
  expect(result.charts[0].valueLabel).toBe('V');
  expect(() => decodeSimulation({ ...ndResponse, volume: { ...volume, values: [1] } }, ndRequest)).toThrow();
  expect(() => decodeSimulation({ ...ndResponse, volume: { ...volume, values: Array(64).fill(-1) } }, ndRequest)).toThrow();
  expect(() => decodeSimulation({ ...ndResponse, volume: { ...volume, x: Array(64).fill(NaN) } }, ndRequest)).toThrow();
  expect(() => decodeSimulation(ndResponse, request)).toThrow();
 });
 test('dimension presets fit the web tensor budget and own their arrays', () => {
  for (const d of [2, 3, 4]) expect(ndAxes(d).reduce((n, axis) => n * axis.points, 1)).toBeLessThanOrEqual(8192);
  const first = simulationDefaults('nd-spectrum'); const second = simulationDefaults('nd-spectrum');
  first.nd!.axes[0].points = 4;
  expect(second.nd!.axes[0].points).toBe(24);
 });
});
