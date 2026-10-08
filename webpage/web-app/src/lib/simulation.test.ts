import { describe, expect, test } from 'bun:test';
import { decodeSimulation, simulationDefaults } from './simulation';

describe('calculation response validation', () => {
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
});
