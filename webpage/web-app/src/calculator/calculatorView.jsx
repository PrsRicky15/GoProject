import React, { useState } from 'react';
import '../components/plotDark.css';

export default function CalculatorView() {
    const [display, setDisplay] = useState('0');

    const appendValue = (value) => {
        if (display === '0') {
            setDisplay(value);
        } else {
            setDisplay(display + value);
        }
    };

    const clearDisplay = () => {
        setDisplay('0');
    };

    const calculate = () => {
        try {
            // Replace Math functions and constants for evaluation
            let expression = display
                .replace(/Math\.PI/g, Math.PI)
                .replace(/Math\.E/g, Math.E);

            // Evaluate the expression safely
            const result = eval(expression);
            setDisplay(result.toString());
        } catch (error) {
            setDisplay('Error');
        }
    };

    return (
        <div className="Calculator">
            <div className="CalcContainer">
                <aside className="sidebar-calc">
                    <div className="Scientific-Calculator">
                        {/* Display */}
                        <label htmlFor="display"></label>
                        <input
                            type="text"
                            placeholder="0"
                            id="display"
                            value={display}
                            readOnly
                        />

                        {/* Row 1 */}
                        <div className="row">
                            <button onClick={() => appendValue('Math.sin(')}>sin</button>
                            <button onClick={() => appendValue('Math.cos(')}>cos</button>
                            <button onClick={() => appendValue('Math.tan(')}>tan</button>
                            <button onClick={() => appendValue('Math.log(')}>log</button>
                            <button onClick={() => appendValue('Math.exp(')}>exp</button>
                        </div>

                        {/* Row 2 */}
                        <div className="row">
                            <button onClick={() => appendValue('7')}>7</button>
                            <button onClick={() => appendValue('8')}>8</button>
                            <button onClick={() => appendValue('9')}>9</button>
                            <button onClick={() => appendValue('/')}>/</button>
                            <button onClick={() => appendValue('(')}>(</button>
                        </div>

                        {/* Row 3 */}
                        <div className="row">
                            <button onClick={() => appendValue('4')}>4</button>
                            <button onClick={() => appendValue('5')}>5</button>
                            <button onClick={() => appendValue('6')}>6</button>
                            <button onClick={() => appendValue('*')}>*</button>
                            <button onClick={() => appendValue(')')}>)</button>
                        </div>

                        {/* Row 4 */}
                        <div className="row">
                            <button onClick={() => appendValue('1')}>1</button>
                            <button onClick={() => appendValue('2')}>2</button>
                            <button onClick={() => appendValue('3')}>3</button>
                            <button onClick={() => appendValue('-')}>-</button>
                            <button onClick={clearDisplay}>C</button>
                        </div>

                        {/* Row 5 */}
                        <div className="row">
                            <button onClick={() => appendValue('0')}>0</button>
                            <button onClick={() => appendValue('.')}>.</button>
                            <button onClick={() => appendValue('+')}>+</button>
                            <button onClick={() => appendValue('Math.sqrt(')}>√</button>
                            <button onClick={() => appendValue('Math.pow(')}>^</button>
                        </div>

                        {/* Row 6 */}
                        <div className="row">
                            <button onClick={() => appendValue('Math.PI')}>π</button>
                            <button onClick={() => appendValue('Math.E')}>e</button>
                            <button onClick={() => appendValue('%')}>%</button>
                            <button onClick={() => appendValue('1/')}>1/x</button>
                            <button onClick={() => appendValue('**2')}>x²</button>
                        </div>

                        {/* Row 7 */}
                        <div className="row">
                            <button onClick={() => appendValue('Math.floor(')}>floor</button>
                            <button onClick={() => appendValue('Math.ceil(')}>ceil</button>
                            <button onClick={() => appendValue('Math.abs(')}>abs</button>
                            <button onClick={() => appendValue('Math.round(')}>round</button>
                            <button onClick={() => appendValue('Math.random()')}>rand</button>
                        </div>

                        {/* Row 8 */}
                        <div className="row">
                            <button onClick={() => appendValue('Math.sin(Math.PI/2)')}>sin90°</button>
                            <button onClick={() => appendValue('Math.cos(Math.PI)')}>cos180°</button>
                            <button onClick={() => appendValue('Math.tan(Math.PI/4)')}>tan45°</button>
                            <button onClick={() => appendValue('Math.log10(')}>log10</button>
                            <button onClick={() => appendValue('Math.exp(1)')}>e¹</button>
                        </div>

                        {/* Row 9 */}
                        <div className="row">
                            <button onClick={calculate} style={{backgroundColor:'#4CAF50', color:'white'}}>=</button>
                            <button onClick={clearDisplay} style={{backgroundColor:'#f44336', color:'white'}}>AC</button>
                            <button onClick={() => appendValue('**0.5')}>√x</button>
                            <button onClick={() => appendValue('Math.pow(')}>xʸ</button>
                            <button onClick={() => appendValue('Math.log(')}>ln</button>
                        </div>

                        <div className="row">
                            <button>Ac</button>
                            <button>Ce</button>
                            <button onClick={() => appendValue('0')}>0</button>
                            <button onClick={() => appendValue('00')}>00</button>
                        </div>
                    </div>
                </aside>

                <main className="plot-area">
                    <div id="plot-container">
                        <h2>Graph Plotter</h2>
                        <p>Plot will be displayed here...</p>
                    </div>
                </main>
            </div>
        </div>
    );
}