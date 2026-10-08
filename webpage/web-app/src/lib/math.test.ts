import { describe, expect, test } from 'bun:test';
import fixtures from '../../../../testdata/math-expressions.json';
import { calculateExpression, compileExpression } from './expression';
import { defaults, validateRequest } from './sampling';
import { decodePlotResult } from './api';
import { makeFigure } from './figure';
import type { PlotResult } from './types';

describe('TypeScript calculator', () => {
  test('matches Go arithmetic fixtures', () => {
    for (const { expression, expected } of fixtures) expect(calculateExpression(expression)).toBeCloseTo(expected, 12);
    expect(calculateExpression('ans+2', 5)).toBe(7);
  });
  test('rejects executable, malformed and non-finite expressions', () => {
    for (const source of ['globalThis', 'constructor(1)', 'Math.sin(1)', 'alert(1)', 'x=1', '2x', 'sin()', '1+', 'pow(2)', 'min(1,2,3)', '(2', '2)', '1;2', '2[0]', '', '1e999', '1+'.repeat(200)+'1', ' '.repeat(1025)]) expect(() => compileExpression(source)).toThrow();
    for (const source of ['1/0', 'sqrt(-1)', 'exp(1000)']) expect(() => calculateExpression(source)).toThrow();
  });
});

describe('Go API boundary', () => {
  const settings = validateRequest({ ...defaults('line'), resolution: 21 });
  const payload = { coords: { x: Array.from({length:21}, (_,i) => i) }, values: Array.from({length:21}, (_,i) => i === 10 ? null : i), count:21, invalid:1, min:0, max:20 };
  test('converts JSON null to a gap without coercing it to zero', () => {
    const result = decodePlotResult(payload, settings);
    expect(result.values[10]).toBeNaN();
    expect(result.values[11]).toBe(11);
    expect(result.coords.x.length).toBe(21);
  });
  test('rejects malformed data rather than trusting TypeScript annotations', () => {
    for (const bad of [null, {}, {...payload, count:22}, {...payload, values:[1]}, {...payload, invalid:0}, {...payload, coords:{x:['1']}}, {...payload, values:payload.values.map(() => '1')}]) expect(() => decodePlotResult(bad,settings)).toThrow();
  });
  test('normalizes form numbers and rejects invalid ranges', () => {
    expect(validateRequest({...defaults('line'), resolution:'801'}).resolution).toBe(801);
    expect(() => validateRequest({...defaults('iso'), resolution:1000000})).toThrow();
    expect(() => validateRequest({...defaults('line'), bounds:{x:['',2],y:[0,1],z:[0,1]}})).toThrow();
  });
});

test('renderer retains the Go grid orientation', () => {
  const settings = validateRequest({...defaults('surface'), resolution:10});
  const result: PlotResult = { settings, coords:{x:Float64Array.from({length:10},(_,i)=>i),y:Float64Array.from({length:10},(_,i)=>i)}, values:Float64Array.from({length:100},(_,i)=>i), count:100,invalid:0,min:0,max:99 };
  const trace = makeFigure(result).data[0];
  expect(trace.type).toBe('surface');
  expect('z' in trace && (trace.z as number[][])[2][3]).toBe(23);
});