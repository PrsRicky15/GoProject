import type { FormEvent } from 'react';
import { useRef, useState } from 'react';
import { calculateExpression } from '../lib/expression';

const keys = ['sin(', 'cos(', 'tan(', 'ln(', 'log(', 'sqrt(', '(', ')', '^', 'pi', '7', '8', '9', '/', 'e', '4', '5', '6', '*', 'ans', '1', '2', '3', '-', 'x²', '0', '.', '%', '+', '1/x'];

export default function Calculator() {
  const [expression, setExpression] = useState('');
  const [answer, setAnswer] = useState(0);
  const [history, setHistory] = useState<{ expression: string; value: number }[]>([]);
  const [error, setError] = useState('');
  const input = useRef<HTMLInputElement>(null);
  function insert(key: string) {
    setError('');
    if (key === 'x²' || key === '1/x') {
      setExpression(key === 'x²' ? `(${expression || 'ans'})^2` : `1/(${expression || 'ans'})`);
      input.current?.focus();
      return;
    }
    const start = input.current?.selectionStart ?? expression.length;
    const end = input.current?.selectionEnd ?? start;
    setExpression(expression.slice(0, start) + key + expression.slice(end));
    requestAnimationFrame(() => { input.current?.focus(); input.current?.setSelectionRange(start + key.length, start + key.length); });
  }
  function calculate(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    try {
      const value = calculateExpression(expression, answer);
      setAnswer(value);
      setHistory((previous) => [{ expression, value }, ...previous].slice(0, 8));
      setError('');
    } catch (problem) { setError(problem instanceof Error ? problem.message : 'Invalid expression.'); }
  }
  return <section className="calculator-panel" aria-labelledby="calculator-title">
    <div className="section-heading"><div><span className="eyebrow">QUICK COMPUTATION</span><h2 id="calculator-title">Scientific calculator</h2></div><span className="tag">Radians</span></div>
    <p className="muted">Evaluate expressions alongside your plot. Use <code>ans</code> for your last answer.</p>
    <form onSubmit={calculate}>
      <label htmlFor="calculation">Expression</label>
      <input ref={input} id="calculation" autoComplete="off" spellCheck="false" placeholder="sqrt(2) / 2" value={expression} onChange={(event) => setExpression(event.target.value)} maxLength={1024} />
      <div className="answer" aria-live="polite"><span>Answer</span><output>{String(answer)}</output></div>
      {error && <p className="error" role="alert">{error}</p>}
      <div className="calculator-keys">{keys.map((key) => <button type="button" key={key} onClick={() => insert(key)} aria-label={key === 'x²' ? 'Square expression' : key === '1/x' ? 'Reciprocal expression' : undefined}>{key === 'pi' ? 'π' : key.replace('(', key === '(' ? '(' : '')}</button>)}</div>
      <div className="calculator-actions"><button type="button" onClick={() => { setExpression(''); setError(''); input.current?.focus(); }}>Clear</button><button className="primary-button" type="submit">Calculate <span>=</span></button></div>
    </form>
    <details className="math-help"><summary>Functions & notation</summary><p>Use * to multiply, ^ for powers, % for remainder. log is base 10; ln is natural log. Trigonometric functions use radians.</p><p>Also available: exp, abs, asin, acos, atan, sinh, cosh, tanh, floor, ceil, round, pow(a,b), min(a,b), max(a,b).</p></details>
    {history.length > 0 && <div className="calculation-history"><h3>Recent calculations</h3>{history.map((item, index) => <button type="button" key={index} onClick={() => { setExpression(item.expression); setError(''); input.current?.focus(); }}><span>{item.expression}</span><strong>{String(item.value)}</strong></button>)}</div>}
  </section>;
}
