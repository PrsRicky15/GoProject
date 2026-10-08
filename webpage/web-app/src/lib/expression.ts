// A small arithmetic grammar: user input never becomes JavaScript source.
type Scope = Record<string, number>;
type Evaluator = (scope: Scope) => number;
type MathFunction = [((a: number) => number), 1] | [((a: number, b: number) => number), 2];
const FUNCTIONS: Readonly<Record<string, MathFunction>> = Object.freeze({
  sin: [Math.sin, 1], cos: [Math.cos, 1], tan: [Math.tan, 1],
  asin: [Math.asin, 1], acos: [Math.acos, 1], atan: [Math.atan, 1],
  sinh: [Math.sinh, 1], cosh: [Math.cosh, 1], tanh: [Math.tanh, 1],
  sqrt: [Math.sqrt, 1], abs: [Math.abs, 1], exp: [Math.exp, 1],
  ln: [Math.log, 1], log: [Math.log10, 1], log10: [Math.log10, 1],
  floor: [Math.floor, 1], ceil: [Math.ceil, 1], round: [Math.round, 1],
  pow: [Math.pow, 2], min: [Math.min, 2], max: [Math.max, 2],
});
const CONSTANTS: Readonly<Record<string, number>> = Object.freeze({ pi: Math.PI, e: Math.E });
const OPS: Record<string, [number, (a: number, b: number) => number]> = {
  '+': [1, (a, b) => a + b], '-': [1, (a, b) => a - b],
  '*': [2, (a, b) => a * b], '/': [2, (a, b) => a / b],
  '%': [2, (a, b) => a % b], '^': [4, (a, b) => a ** b],
};

export function compileExpression(source: string, variables: string[] = []): Evaluator {
  if (typeof source !== 'string' || !source.trim()) throw new Error('Enter an expression.');
  if (source.length > 1024) throw new Error('Keep expressions under 1,024 characters.');
  const input = source.replaceAll('π', 'pi').replaceAll('**', '^');
  const tokens: string[] = [];
  let offset = 0;
  while (offset < input.length) {
    const match = /^(?:\s+|(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?|[a-zA-Z_][a-zA-Z_0-9]*|[+\-*/%^(),])/.exec(input.slice(offset));
    if (!match) throw new Error(`Unexpected character at position ${offset + 1}. Use * for multiplication.`);
    if (match[0].trim()) tokens.push(match[0]);
    offset += match[0].length;
  }
  if (tokens.length > 256) throw new Error('Expression is too complex (maximum 256 tokens).');
  const allowed = new Set(variables);
  let cursor = 0;
  const consume = (token: string) => {
    if (tokens[cursor++] !== token) throw new Error(`Expected "${token}". Check parentheses and arguments.`);
  };
  function parse(minimum = 0): Evaluator {
    const token = tokens[cursor++];
    let left: Evaluator;
    if (token === '+' || token === '-') {
      const operand = parse(3);
      left = token === '-' ? (scope) => -operand(scope) : operand;
    } else if (token === '(') {
      left = parse();
      consume(')');
    } else if (token && /^(?:\d|\.)/.test(token)) {
      const number = Number(token);
      if (!Number.isFinite(number)) throw new Error('Number is too large.');
      left = () => number;
    } else if (Object.hasOwn(CONSTANTS, token)) {
      left = () => CONSTANTS[token];
    } else if (allowed.has(token)) {
      left = (scope) => scope[token];
    } else if (Object.hasOwn(FUNCTIONS, token)) {
      const [fn, count] = FUNCTIONS[token];
      consume('(');
      const first = parse();
      if (count === 2) {
        consume(',');
        const second = parse();
        left = (scope) => fn(first(scope), second(scope));
      } else left = (scope) => fn(first(scope));
      consume(')');
    } else throw new Error(token ? `Unknown name or unexpected token: "${token}".` : 'Expression is incomplete.');
    while (Object.hasOwn(OPS, tokens[cursor])) {
      const operator = tokens[cursor];
      const [precedence, fn] = OPS[operator];
      if (precedence < minimum) break;
      cursor++;
      const right = parse(precedence + (operator === '^' ? 0 : 1));
      const previous = left;
      left = (scope) => fn(previous(scope), right(scope));
    }
    return left;
  }
  const evaluate = parse();
  if (cursor !== tokens.length) throw new Error(`Unexpected "${tokens[cursor]}". Use * for multiplication.`);
  return evaluate;
}

export function calculateExpression(source: string, ans = 0) {
  const result = compileExpression(source, ['ans'])({ ans });
  if (!Number.isFinite(result)) throw new Error('Result is undefined or outside the real-number range.');
  return result;
}
