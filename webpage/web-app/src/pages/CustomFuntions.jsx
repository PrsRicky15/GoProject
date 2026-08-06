import React, { useState, useMemo, useCallback } from 'react';
import * as math from 'mathjs';
import Plot from 'react-plotly.js';

const presetFunctions = [
    { label: 'Quadratic', expr: 'x^2' },
    { label: 'Sine Wave', expr: 'sin(x)' },
    { label: 'Cosine Wave', expr: 'cos(x)' },
    { label: 'Exponential', expr: 'e^x' },
    { label: 'Logarithm', expr: 'log(x)' },
    { label: 'Cubic', expr: 'x^3 - 3x' },
    { label: 'Tangent', expr: 'tan(x)' },
    { label: 'Square Root', expr: 'sqrt(x)' },
    { label: 'Gaussian', expr: 'e^(-x^2)' },
    { label: 'Sinc', expr: 'sin(x)/x' },
];

const mathHelp = [
    { category: 'Basic', ops: ['+', '-', '*', '/', '^'] },
    { category: 'Trig', ops: ['sin', 'cos', 'tan', 'asin', 'acos', 'atan'] },
    { category: 'Other', ops: ['sqrt', 'log', 'exp', 'abs', 'floor', 'ceil'] },
    { category: 'Constants', ops: ['pi', 'e'] },
];

export default function FunctionVisualizer() {
    const [expression, setExpression] = useState('sin(x) * x');
    const [xMin, setXMin] = useState(-10);
    const [xMax, setXMax] = useState(10);
    const [nPoints, setNpoints] = useState(200);
    const [error, setError] = useState(null);



    const parseExpression = useCallback((expr) => {
        try {
            // Pre-process the expression for common mathematical conventions
            let processed = expr
                .replace(/(\d)([a-zA-Z])/g, '$1*$2')  // 2x -> 2*x
                .replace(/([a-zA-Z])(\d)/g, '$1*$2')  // x2 -> x*2 (less common but handle it)
                .replace(/\)\(/g, ')*(')              // )(  -> )*(
                .replace(/(\d)\(/g, '$1*(')           // 2(  -> 2*(
                .replace(/\)(\d)/g, ')*$1')           // )2  -> )*2
                .replace(/([a-zA-Z])\(/g, (match, p1) => {
                    // Don't add * for known functions
                    const funcs = ['sin', 'cos', 'tan', 'asin', 'acos', 'atan', 'sinh', 'cosh', 'tanh', 
                                'sqrt', 'log', 'ln', 'exp', 'abs', 'floor', 'ceil', 'round', 'sign'];
                    if (funcs.some(f => expr.includes(f + '('))) {
                        return match;
                    }
                    return p1 + '*(';
                });
      
            const compiled = math.compile(processed);
            return { compiled, processed, error: null };
        } catch (e) {
            return { compiled: null, processed: expr, error: e.message };
        }
    }, []);

    const { xData, yData, parsedExpr } = useMemo(() => {
        const { compiled, processed, error } = parseExpression(expression);
    
        if (error) {
            setError(error);
            return { xData: [], yData: [], parsedExpr: processed };
        }
    
        setError(null);
        const xPoints = [];
        const yPoints = [];
        const step = (xMax - xMin) / nPoints;
    
        for (let i = 0; i <= nPoints; i++) {
            const x = xMin + i * step;
            try {
                let y = compiled.evaluate({ x });
        
                // Handle special cases
                if (!isFinite(y) || isNaN(y)) {
                    y = null;
                } else if (Math.abs(y) > 1e10) {
                    y = null; // Clip extremely large values
                }
        
                xPoints.push(parseFloat(x.toFixed(6)));
                yPoints.push(y);
            } catch {
                xPoints.push(parseFloat(x.toFixed(6)));
                yPoints.push(null);
            }
        }
    
        return { xData: xPoints, yData:yPoints, parsedExpr: processed, parseError: null };
    }, [expression, xMin, xMax, nPoints, parseExpression]);

    const plotData = [
        {
            x: xData,
            y: yData,
            type: 'scatter',
            mode: 'lines',
            name: `f(x) = ${parsedExpr}`,
            line: {
                color: '#dbe72ee0',
                width: 2.5,
                shape: 'spline',
            },
            connectgaps: false, // Important: don't connect across discontinuities
            hovertemplate: 'x: %{x:.4f}<br>y: %{y:.4f}<extra></extra>',
        },
    ];

    const plotLayout = {
        paper_bgcolor: 'rgba(26, 32, 44, 0)',
        plot_bgcolor: 'rgba(10, 10, 15, 0.8)',
        font: {
            family: '"JetBrains Mono", "Fira Code", monospace',
            color: '#000000ff',
        },
        title: {
            text: `y = ${parsedExpr}`,
            font: { size: 24, color: '#000000ff' },
            x: 0.5,
        },
        xaxis: {
            title: { text: 'x', font: { size: 20 } },
            gridcolor: 'rgba(99, 179, 237, 0.15)',
            zerolinecolor: 'rgba(99, 179, 237, 0.5)',
            zerolinewidth: 2,
            tickfont: { size: 16 },
            range: [xMin, xMax],
        },
        yaxis: {
            title: { text: 'f(x)', font: { size: 20 } },
            gridcolor: 'rgba(99, 179, 237, 0.15)',
            zerolinecolor: 'rgba(99, 179, 237, 0.5)',
            zerolinewidth: 2,
            tickfont: { size: 16 },
            autorange: true,
        },
        margin: { t: 50, r: 30, b: 50, l: 60 },
        hovermode: 'closest',
        dragmode: 'zoom',
        showlegend: false,
    };

    const plotConfig = {
        responsive: true,
        displaylogo: false,
        modeBarButtonsToRemove: ['lasso2d', 'select2d'],
        toImageButtonOptions: {
            format: 'svg',
            filename: `function_${expression.replace(/[^a-zA-Z0-9]/g, '_')}`,
        },
    };

    return (
        <div style={{
            padding: '24px',
            
        }}>
            <div style={{
                backgroundSize: '40px 40px',
                pointerEvents: 'auto',
                position: 'fixed',
                inset: 0,
                overflowX: 'auto',
            }}>
                <h1 style={{ textAlign: 'center', marginBottom: '32px', padding: '20px'}}>Custom Function Visualizer</h1>
                <label style={{ 
                    paddingLeft: '30px',
                    display: 'block', 
                    fontSize: '1.5rem', 
                    color: '#63b3ed',
                    letterSpacing: '0.15em',
                    textTransform: 'uppercase',
                    marginBottom: '8px',
                }}>
                    Function Expression
                </label>

                <div style={{ display: 'flex', gap: '12px', alignItems: 'center', padding: '30px' }}>
                    <span style={{ 
                        fontSize: '1.5rem', 
                        color: '#dbe72ee0',
                        fontStyle: 'italic',
                    }}>
                        f(x) =
                    </span>
                
                    <input
                        type="text"
                        value={expression}
                        onChange={(e) => setExpression(e.target.value)}
                        placeholder="e.g., sin(x) * x^2"
                        style={{
                            flex: 1,
                            background: 'rgba(0, 0, 0, 0.4)',
                            border: error ? '2px solid #fc8181' : '2px solid rgba(99, 179, 237, 0.3)',
                            borderRadius: '8px',
                            padding: '16px 20px',
                            fontSize: '1.25rem',
                            color: '#fff',
                            fontFamily: 'inherit',
                            outline: 'none',
                            transition: 'border-color 0.2s, box-shadow 0.2s',
                        }}
                        onFocus={(e) => {
                            e.target.style.borderColor = '#63b3ed';
                            e.target.style.boxShadow = '0 0 20px rgba(99, 179, 237, 0.3)';
                        }}
                        onBlur={(e) => {
                            e.target.style.borderColor = error ? '#fc8181' : 'rgba(99, 179, 237, 0.3)';
                            e.target.style.boxShadow = 'none';
                        }}
                    />
                </div>
                {error && (
                    <div style={{
                        marginTop: '12px',
                        padding: '12px 16px',
                        background: 'rgba(252, 129, 129, 0.1)',
                        border: '1px solid rgba(252, 129, 129, 0.3)',
                        borderRadius: '8px',
                        color: '#fc8181',
                        fontSize: '0.875rem',
                    }}>
                        ⚠ Parse Error: {error}
                    </div>
                )}

                {/* Preset Functions */}
                <div style={{ marginTop: '0px', paddingLeft: '30px' }}>
                    <span style={{ 
                        fontSize: '0.7rem', 
                        color: '#718096',
                        letterSpacing: '0.1em',
                        textTransform: 'uppercase',
                    }}>
                        Quick Presets:
                    </span>
                    <div style={{ 
                        display: 'flex',  
                        gap: '20px', 
                        marginTop: '8px' 
                    }}>
                        {presetFunctions.map((preset) => (
                            <button
                                key={preset.expr}
                                onClick={() => setExpression(preset.expr)}
                                style={{
                                    border: '1px solid rgba(99, 179, 237, 0.3)',
                                    borderRadius: '20px',
                                    padding: '6px 14px',
                                    color: expression === preset.expr ? '#fff' : '#a0aec0',
                                    fontSize: '0.75rem',
                                    cursor: 'pointer',
                                    transition: 'all 0.2s',
                                    fontFamily: 'inherit',
                                }}
                                onMouseEnter={(e) => {
                                    if (expression !== preset.expr) {
                                        e.target.style.background = 'rgba(99, 179, 237, 0.2)';
                                        e.target.style.color = '#fff';
                                    }
                                }}
                                onMouseLeave={(e) => {
                                    if (expression !== preset.expr) {
                                        e.target.style.background = 'rgba(255, 255, 255, 0.05)';
                                        e.target.style.color = '#a0aec0';
                                    }
                                }}
                            >
                                {preset.label}
                            </button>
                        ))}
                    </div>
                </div>
                
                {/*Plotting area*/}
                <div style={{
                    display: 'flex',
                    flexDirection: 'row',
                    paddingTop: '30px'
                }}>
                    {/* Controls Column */}
                    <div style={{
                        display: 'flex',
                        flexDirection: 'column',
                        gap: '16px',
                        marginBottom: '24px',
                        alignItems: 'flex-start',
                        padding: '30px'
                    }}>
                        {/* X Range Controls */}
                        <div style={{
                            background: 'rgba(26, 32, 44, 0.6)',
                            borderRadius: '12px',
                            border: '1px solid rgba(99, 179, 237, 0.15)',
                            padding: '16px',
                            alignContent: 'center',
                        }}>
                            <label style={{ 
                                fontSize: '0.7rem', 
                                color: '#63b3ed',
                                letterSpacing: '0.1em',
                                textTransform: 'uppercase',
                            }}>
                                X Range
                            </label>
                            <div style={{ display: 'flex', gap: '12px', marginTop: '8px', alignItems: 'center', paddingBottom: '8px'}}>
                                <input
                                    type="number"
                                    value={xMin}
                                    onChange={(e) => setXMin(parseFloat(e.target.value) || -10)}
                                    style={{
                                        width: '80px',
                                        background: 'rgba(0, 0, 0, 0.4)',
                                        border: '1px solid rgba(99, 179, 237, 0.3)',
                                        borderRadius: '6px',
                                        padding: '8px 12px',
                                        color: '#fff',
                                        fontSize: '0.875rem',
                                        fontFamily: 'inherit',
                                        textAlign: 'center',
                                    }}
                                />
                                <span style={{ color: '#fff' }}>to</span>
                                <input
                                    type="number"
                                    value={xMax}
                                    onChange={(e) => setXMax(parseFloat(e.target.value) || 10)}
                                    style={{
                                        width: '80px',
                                        background: 'rgba(0, 0, 0, 0.4)',
                                        border: '1px solid rgba(99, 179, 237, 0.3)',
                                        borderRadius: '6px',
                                        padding: '8px 12px',
                                        color: '#fff',
                                        fontSize: '0.875rem',
                                        fontFamily: 'inherit',
                                        textAlign: 'center',
                                    }}
                                />
                            </div>

                            <label style={{ 
                                paddingTop: '20px',
                                fontSize: '0.7rem', 
                                color: '#63b3ed',
                                letterSpacing: '0.1em',
                                textTransform: 'uppercase',
                            }}>
                                Resolution: {nPoints} points
                            </label>
                            <div style={{ display: 'flex', gap: '12px', marginTop: '8px', alignItems: 'center' }}>
                                <input
                                    type="number"
                                    value={nPoints}
                                    onChange={(e) => setNpoints(parseInt(e.target.value) || 200)}
                                    style={{
                                        width: '80px',
                                        background: 'rgba(0, 0, 0, 0.4)',
                                        border: '1px solid rgba(99, 179, 237, 0.3)',
                                        borderRadius: '6px',
                                        padding: '8px 12px',
                                        color: '#fff',
                                        fontSize: '0.875rem',
                                        fontFamily: 'inherit',
                                        textAlign: 'center',
                                    }}
                                />
                            </div>
                        </div>

                        {/* Math Reference */}
                        <div style={{
                            background: 'rgba(26, 32, 44, 0.6)',
                            borderRadius: '12px',
                            border: '1px solid rgba(99, 179, 237, 0.15)',
                            padding: '16px',
                            gridColumn: 'span 1',
                        }}>
                            <label style={{ 
                                fontSize: '1rem', 
                                color: '#63b3ed',
                                letterSpacing: '0.2em',
                                textTransform: 'uppercase',
                            }}>
                                Operators & Functions
                            </label>
                            <div style={{ marginTop: '8px', fontSize: '1rem' }}>
                                {mathHelp.map((cat) => (
                                    <div key={cat.category} style={{ marginBottom: '4px' }}>
                                        <span style={{ color: '#718096' }}>{cat.category}: </span>
                                        <span style={{ color: '#e2e8f0' }}>{cat.ops.join(', ')}</span>
                                    </div>
                                ))}
                            </div>
                        </div>
                    </div>

                    {/* Plot Container */}
                    <Plot
                        data={plotData}
                        layout={plotLayout}
                        config={plotConfig}
                        style={{
                            width: '100%',
                            height: 'auto', 
                        }}
                        useResizeHandler={true}
                    />
                </div>
            </div>    
        </div>
        
  );
}