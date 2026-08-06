import React, { useState, useMemo } from 'react';
import { ScientificCalculator } from './basicPlotter/calculator_and_MathFunc.ts';
import './basicPlotter/Calculator.css';

export default function CalculatorApp() {
    // ===== STATE MANAGEMENT =====
    const calculator = useMemo(() => new ScientificCalculator(), []);
    const [display, setDisplay] = useState('0');
    const [currentInput, setCurrentInput] = useState('');
    const [previousValue, setPreviousValue] = useState<number | null>(null);
    const [operation, setOperation] = useState<string | null>(null);
    const [showScientific, setShowScientific] = useState(false);
    const [error, setError] = useState<string | null>(null);

    const resetCalculator = () => {
        calculator.clear();
        setDisplay('0');
        setCurrentInput('');
        setPreviousValue(null);
        setOperation(null);
        setError(null);
    };

    const updateDisplay = (value: string) => {
        setDisplay(value);
        setCurrentInput(value);
    };

    const canAddDecimal = (): boolean => {
        return !currentInput.includes('.');
    };

    const shouldReplaceDisplay = (num: string): boolean => {
        return display === '0' && num !== '.';
    };

    const formatNumber = (num: number): string => {
        // If it's a whole number, show as is
        if (Number.isInteger(num)) return String(num);

        // If it has many decimals, round to 8 places
        return num.toFixed(8).replace(/\.?0+$/, '');
    };

    const showError = (message: string) => {
        setError(message);
        setDisplay('Error');
    };

    const clearError = () => {
        setError(null);
    };

    const handleNumberClick = (num: string) => {
        clearError();

        // Prevent multiple decimal points
        if (num === '.' && !canAddDecimal()) return;

        const newValue = shouldReplaceDisplay(num) ? num : display + num;
        updateDisplay(newValue);
    };

    const handleBasicOperation = (op: string) => {
        if (currentInput === '') return;

        clearError();
        const current = parseFloat(currentInput);

        // If we already have a pending operation, calculate it first
        if (previousValue !== null && operation) {
            executeOperation();
        } else {
            setPreviousValue(current);
        }

        setOperation(op);
        setCurrentInput('');
        setDisplay('0');
    };

    const executeOperation = () => {
        if (previousValue === null || currentInput === '' || !operation) return;

        const current = parseFloat(currentInput);

        try {
            let result: number;

            switch (operation) {
                case '+':
                    result = calculator.add(previousValue, current);
                    break;
                case '-':
                    result = calculator.subtract(previousValue, current);
                    break;
                case '*':
                    result = calculator.multiply(previousValue, current);
                    break;
                case '/':
                    result = calculator.divide(previousValue, current);
                    break;
                default:
                    return;
            }

            const formattedResult = formatNumber(result);
            setDisplay(formattedResult);
            setCurrentInput(formattedResult);
            setPreviousValue(null);
            setOperation(null);

        } catch (err) {
            showError(err instanceof Error ? err.message : 'Error');
        }
    };

    const handleScientificFunction = (func: string) => {
        if (currentInput === '') return;

        clearError();
        const value = parseFloat(currentInput);

        try {
            let result: number;

            switch (func) {
                case 'sqrt':
                    result = calculator.squareRoot(value);
                    break;
                case 'square':
                    result = calculator.power(value, 2);
                    break;
                case 'sin':
                    result = calculator.sine(value);
                    break;
                case 'cos':
                    result = calculator.cosine(value);
                    break;
                case 'tan':
                    result = calculator.tangent(value);
                    break;
                case 'ln':
                    result = calculator.naturalLog(value);
                    break;
                case 'log':
                    result = calculator.logarithm(value, 10);
                    break;
                case 'factorial':
                    result = calculator.factorial(Math.floor(value));
                    break;
                case 'abs':
                    result = calculator.abs(value);
                    break;
                default:
                    return;
            }

            const formattedResult = formatNumber(result);
            updateDisplay(formattedResult);

        } catch (err) {
            showError(err instanceof Error ? err.message : 'Error');
        }
    };

    const handleBackspace = () => {
        if (currentInput.length <= 1) {
            updateDisplay('0');
            setCurrentInput('');
        } else {
            const newValue = currentInput.slice(0, -1);
            updateDisplay(newValue);
        }
    };

    const toggleMode = () => {
        setShowScientific(!showScientific);
    };

    const renderButton = (
        label: string,
        onClick: () => void,
        className: string = 'btn-number',
        spanCols: number = 1
    ) => {
        const spanClass = spanCols > 1 ? 'btn-span-2' : '';
        return (
            <button
                onClick={onClick}
                className={`calc-btn ${className} ${spanClass}`}
            >
                {label}
            </button>
        );
    };

    const renderScientificPanel = () => {
        if (!showScientific) return null;

        const scientificButtons = [
            { label: 'sin', func: 'sin' },
            { label: 'cos', func: 'cos' },
            { label: 'tan', func: 'tan' },
            { label: 'sqrt', func: 'sqrt' },
            { label: 'x^2', func: 'square' },
            { label: 'ln', func: 'ln' },
            { label: 'log', func: 'log' },
            { label: 'n!', func: 'factorial' },
        ];

        return (
            <div className="scientific-panel">
                {scientificButtons.map((btn) => (
                    <button
                        key={btn.func}
                        onClick={() => handleScientificFunction(btn.func)}
                        className="scientific-btn"
                    >
                        {btn.label}
                    </button>
                ))}
            </div>
        );
    };

    return (
        <div className="calculator-container">
            <div className="calculator-wrapper">
                {/* Header */}
                <div className="calculator-header">
                    <h1 className="calculator-title">Calculator</h1>
                    <button onClick={toggleMode} className="mode-toggle-btn">
                        {showScientific ? 'Basic' : 'Scientific'}
                    </button>
                </div>

                {/* Display */}
                <div className="calculator-display">
                    <div className="display-operation">
                        {operation && previousValue !== null && `${previousValue} ${operation}`}
                    </div>
                    <div className={`display-value ${error ? 'error' : ''}`}>
                        {error || display}
                    </div>
                </div>

                {/* Scientific Panel (conditional) */}
                {renderScientificPanel()}

                {/* Main Button Grid */}
                <div className="button-grid">
                    {/* Row 1 */}
                    {renderButton('C', resetCalculator, 'btn-clear')}
                    {renderButton('del', handleBackspace, 'btn-special')}
                    {renderButton('|x|', () => handleScientificFunction('abs'), 'btn-special')}
                    {renderButton('/', () => handleBasicOperation('/'), 'btn-operator')}

                    {/* Row 2 */}
                    {renderButton('7', () => handleNumberClick('7'))}
                    {renderButton('8', () => handleNumberClick('8'))}
                    {renderButton('9', () => handleNumberClick('9'))}
                    {renderButton('*', () => handleBasicOperation('*'), 'btn-operator')}

                    {/* Row 3 */}
                    {renderButton('4', () => handleNumberClick('4'))}
                    {renderButton('5', () => handleNumberClick('5'))}
                    {renderButton('6', () => handleNumberClick('6'))}
                    {renderButton('−', () => handleBasicOperation('-'), 'btn-operator')}

                    {/* Row 4 */}
                    {renderButton('1', () => handleNumberClick('1'))}
                    {renderButton('2', () => handleNumberClick('2'))}
                    {renderButton('3', () => handleNumberClick('3'))}
                    {renderButton('+', () => handleBasicOperation('+'), 'btn-operator')}

                    {/* Row 5 */}
                    {renderButton('0', () => handleNumberClick('0'), 'btn-number', 2)}
                    {renderButton('.', () => handleNumberClick('.'))}
                    {renderButton('=', executeOperation, 'btn-equals')}
                </div>
            </div>
        </div>
    );
}