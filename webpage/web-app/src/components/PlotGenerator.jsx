import React, { useState } from 'react';
import Plot from 'react-plotly.js';
import { BarChart3, Loader, Download, RefreshCw, Code } from 'lucide-react';
import './plotDark.css';

// ============================================================================
// CONSTANTS & CONFIGURATION
// ============================================================================

const PLOT_TYPES = [
    { value: 'Custom', label: 'Custom Function' },
    { value: 'Harmonic', label: 'Harmonic' },
    { value: 'Polynomial', label: 'Polynomial' },
    { value: 'Barrier', label: 'Barrier' },
    { value: 'Gaussian', label: 'Gaussian' },
    { value: 'MultiGaussian', label: 'Multiple Gaussian' },
    { value: 'SuperGaussian', label: 'Super Gaussian' },
    { value: 'Sinusoidal', label: 'Sinusoidal' },
    { value: 'Morse', label: 'Morse' },
    { value: 'Softcore', label: 'Softcore' },
    { value: 'Demo-plot', label: 'Demo Plot' },
    { value: 'surface_3d', label: '3D Surface Plot' },
    { value: 'energy_levels', label: 'Energy Level Diagram' },
  ];

const PARAM_CONFIGS = {
    RectangleBarrier: [
        { key: 'D', label: 'Charge (q)', step: '1' },
        { key: 'a', label: 'Softcore Parameter (a)', step: '0.1' },
        { key: 'r0', label: 'Center (r₀)', step: '0.1' },
    ],
    Sinusoidal: [
        { key: 'D', label: 'Charge (q)', step: '1' },
        { key: 'a', label: 'Softcore Parameter (a)', step: '0.1' },
        { key: 'r0', label: 'Center (r₀)', step: '0.1' },
    ],
    Harmonic: [
        { key: 'x0', label: 'Center', step: '0.1' },
        { key: 'k', label: 'Force Constant', step: '1' },
    ],
    Polynomial: [
        { key: 'x0', label: 'Center', step: '0.1' },
        { key: 'k', label: 'Force Constant', step: '1' },
    ],
    Morse: [
        { key: 'D0', label: 'Dissociation Energy (D)', step: '10' },
        { key: 'a', label: 'Width Parameter (a)', step: '0.1' },
        { key: 'rEq', label: 'Equilibrium Distance', step: '0.1' },
    ],
    Softcore: [
        { key: 'q', label: 'Charge (q)', step: '1' },
        { key: 'a', label: 'Softcore Parameter (a)', step: '0.1' },
        { key: 'r0', label: 'Center (r₀)', step: '0.1' },
    ],
    Gaussian: [
        { key: 'x0', label: 'Center', step: '1' },
        { key: 'a', label: 'Softcore Parameter (a)', step: '0.1' },
        { key: 'r0', label: 'Center (r₀)', step: '0.1' },
    ],
    MultiGaussian: [
        { key: 'x0', label: 'Center', step: '1' },
        { key: 'a', label: 'Softcore Parameter (a)', step: '0.1' },
        { key: 'r0', label: 'Center (r₀)', step: '0.1' },
    ],
    SuperGaussian: [
        { key: 'x0', label: 'Center', step: '1' },
        { key: 'a', label: 'Softcore Parameter (a)', step: '0.1' },
        { key: 'r0', label: 'Center (r₀)', step: '0.1' },
    ],
};

const EXAMPLE_FUNCTIONS = [
    { label: 'Sine', expr: 'sin(x)' },
    { label: 'Parabola', expr: 'x^2' },
    { label: 'Damped', expr: 'exp(-x/5)*cos(2*x)' },
    { label: 'Gaussian', expr: 'exp(-x^2/2)' },
    { label: 'Cubic', expr: 'x^3 - 3*x' },
];

const INITIAL_GRID_PARAMS = { rMin: 0.0, rMax: 10.0, nGrid: 100 };
const INITIAL_POT_PARAMS = { D: 100.0, a: 1.5, r0: 2.0, D0: 100.0, rEq: 2.0, x0: 5.0, k: 1.0, q: 1.0 };

const PLOT_CONFIG = {
    responsive: true,
    displayModeBar: true,
    displaylogo: false,
    modeBarButtonsToRemove: ['lasso2d', 'select2d'],
    toImageButtonOptions: {
        format: 'png',
        filename: 'quantum_plot',
        height: 800,
        width: 1200,
        scale: 2,
    },
};

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

const evaluateExpression = (expr, x) => {
    try {
        const mathFunctions = {
            sin: Math.sin,
            cos: Math.cos,
            tan: Math.tan,
            exp: Math.exp,
            log: Math.log,
            sqrt: Math.sqrt,
            abs: Math.abs,
            pow: Math.pow,
            PI: Math.PI,
            E: Math.E,
        };

        let processedExpr = expr
            .replace(/\^/g, '**')
            .replace(/(\d+)x/g, '$1*x')
            .replace(/x(\d+)/g, 'x*$1')
            .replace(/\)(\d)/g, ')*$1')
            .replace(/(\d)\(/g, '$1*(');

        const func = new Function('x', 'math', `
      with(math) {
        return ${processedExpr};
      }
    `);

        return func(x, mathFunctions);
    } catch (error) {
        throw new Error(`Invalid expression: ${error.message}`);
    }
};

// ============================================================================
// SUB-COMPONENTS
// ============================================================================

const NumberInput = ({ label, value, onChange, step = "0.1" }) => (
    <div className="input-group">
        <label className="input-label">{label}</label>
        <input
            type="number"
            step={step}
            value={value}
            onChange={(e) => onChange(parseFloat(e.target.value))}
            className="input-field"
        />
    </div>
);

const GridParametersSection = ({ gridParams, updateGridParam }) => (
    <div>
        <NumberInput
            label="rMin"
            value={gridParams.rMin}
            onChange={(val) => updateGridParam('rMin', val)}
        />
        <NumberInput
            label="rMax"
            value={gridParams.rMax}
            onChange={(val) => updateGridParam('rMax', val)}
        />
        <NumberInput
            label="nGrid"
            value={gridParams.nGrid}
            onChange={(val) => updateGridParam('nGrid', val)}
            step="1"
        />
    </div>
);

const PlotTypeSelector = ({ plotType, setPlotType }) => (
    <div className="input-group">
        <label className="input-label">Plot Type</label>
        <select
            value={plotType}
            onChange={(e) => setPlotType(e.target.value)}
            className="select-field"
        >
            {PLOT_TYPES.map((type) => (
                <option key={type.value} value={type.value}>
                    {type.label}
                </option>
            ))}
        </select>
    </div>
);

const CustomFunctionInput = ({ customFunction, setCustomFunction, loadExample, error }) => (
    <div className="input-group">
        <label className="input-label">
            <Code className="icon" style={{ display: 'inline', verticalAlign: 'middle', marginRight: '6px' }} />
            Function f(x)
        </label>
        <textarea
            value={customFunction}
            onChange={(e) => setCustomFunction(e.target.value)}
            placeholder="e.g., sin(x) * exp(-x/5)"
            className="input-field"
            style={{
                fontFamily: 'monospace',
                minHeight: '80px',
                resize: 'vertical'
            }}
        />

        <div style={{ marginTop: '8px' }}>
            <p style={{ fontSize: '12px', opacity: 0.7, marginBottom: '6px' }}>
                Quick Examples:
            </p>
            <div style={{ display: 'flex', flexWrap: 'wrap', gap: '6px' }}>
                {EXAMPLE_FUNCTIONS.map((ex) => (
                    <button
                        key={ex.expr}
                        onClick={() => loadExample(ex.expr)}
                        className="input-field"
                        style={{
                            padding: '4px 8px',
                            fontSize: '11px',
                            cursor: 'pointer',
                            width: 'auto',
                            minHeight: 'auto'
                        }}
                    >
                        {ex.label}
                    </button>
                ))}
            </div>
        </div>

        <div className="info-box" style={{ marginTop: '12px', fontSize: '11px' }}>
            <p style={{ fontWeight: 'bold', marginBottom: '4px' }}>Supported:</p>
            <p>sin, cos, tan, exp, log, sqrt, abs, pow, ^, +, -, *, /, (, )</p>
        </div>
    </div>
);

const DynamicParametersSection = ({ plotType, potParams, updatePotParam }) => {
    const paramConfig = PARAM_CONFIGS[plotType];

    if (!paramConfig) return null;

    return (
        <div className="params-section">
            {paramConfig.map(({ key, label, step }) => (
                <NumberInput
                    key={key}
                    label={label}
                    value={potParams[key] || 0}
                    step={step}
                    onChange={(val) => updatePotParam(key, val)}
                />
            ))}
        </div>
    );
};

const GenerateButton = ({ isGenerating, onClick }) => (
    <button
        onClick={onClick}
        disabled={isGenerating}
        className="generate-btn"
    >
        {isGenerating ? (
            <>
                <Loader className="icon animate-spin" />
                <span>Generating...</span>
            </>
        ) : (
            <>
                <RefreshCw className="icon" />
                <span>Generate Plot</span>
            </>
        )}
    </button>
);

const InfoBox = () => (
    <div className="info-box">
        <p className="info-title">✨ Interactive Features:</p>
        <ul className="info-list">
            <li>• Zoom: Drag to select area</li>
            <li>• Pan: Shift + drag</li>
            <li>• Reset: Double click</li>
            <li>• Hover for values</li>
            <li>• Download via toolbar</li>
        </ul>
    </div>
);

const PlotDisplay = ({ plotData, isGenerating, plotLayout, downloadPlot }) => (
    <div>
        <div className="plot-header">
            <h2 className="plot-title">Visualization</h2>
            {plotData && (
                <button onClick={downloadPlot} className="download-btn">
                    <Download className="icon" />
                    <span>Download PNG</span>
                </button>
            )}
        </div>

        <div className="plot-container">
            {isGenerating ? (
                <div className="plot-placeholder">
                    <div className="placeholder-content">
                        <Loader className="placeholder-icon animate-spin" />
                        <p className="placeholder-text">Generating plot...</p>
                    </div>
                </div>
            ) : plotData ? (
                <Plot
                    data={plotData.data}
                    layout={plotLayout}
                    config={PLOT_CONFIG}
                    style={{ width: '100%', height: '600px' }}
                    useResizeHandler={true}
                />
            ) : (
                <div className="plot-placeholder">
                    <div className="placeholder-content">
                        <BarChart3 className="placeholder-icon-large" />
                        <p className="placeholder-text">
                            Configure parameters and click "Generate Plot" to visualize
                        </p>
                    </div>
                </div>
            )}
        </div>

        {plotData && (
            <div className="alert alert-success">
                ✓ Plot generated successfully - Interact with the plot using your mouse!
            </div>
        )}
    </div>
);

// ============================================================================
// MAIN COMPONENT
// ============================================================================

export default function PlotGenerator() {
    const [gridParams, setGridParams] = useState(INITIAL_GRID_PARAMS);
    const [potParams, setPotParams] = useState(INITIAL_POT_PARAMS);
    const [plotType, setPlotType] = useState('Morse');
    const [customFunction, setCustomFunction] = useState('sin(x) + cos(2*x)');
    const [isGenerating, setIsGenerating] = useState(false);
    const [plotData, setPlotData] = useState(null);
    const [error, setError] = useState(null);

    const updateGridParam = useCallback((key, value) => {
        setGridParams(prev => ({ ...prev, [key]: value }));
    }, []);

    const updatePotParam = useCallback((key, value) => {
        setPotParams(prev => ({ ...prev, [key]: value }));
    }, []);

    const generateCustomPlot = useCallback(() => {
        const { rMin, rMax, nGrid } = gridParams;
        const step = (rMax - rMin) / (nGrid - 1);
        const xValues = [];
        const yValues = [];

        for (let i = 0; i < nGrid; i++) {
            const x = rMin + i * step;
            try {
                const y = evaluateExpression(customFunction, x);
                if (isFinite(y)) {
                    xValues.push(x);
                    yValues.push(y);
                }
            } catch (err) {
                throw new Error(`Error at x=${x.toFixed(2)}: ${err.message}`);
            }
        }

        return {
            data: [{
                x: xValues,
                y: yValues,
                type: 'scatter',
                mode: 'lines',
                line: { color: '#3b82f6', width: 2 },
                name: `f(x) = ${customFunction}`
            }],
            layout: {
                title: `f(x) = ${customFunction}`,
                xaxis: { title: 'x' },
                yaxis: { title: 'f(x)' },
                paper_bgcolor: 'rgba(0,0,0,0)',
                plot_bgcolor: 'rgba(0,0,0,0.02)'
            }
        };
    }, [gridParams, customFunction]);

    const generatePlot = useCallback(async () => {
        setIsGenerating(true);
        setError(null);

        try {
            if (plotType === 'Custom') {
                const result = generateCustomPlot();
                setPlotData(result);
            } else {
                // For demo purposes - replace with actual API call
                await new Promise(resolve => setTimeout(resolve, 1000));
                const x = Array.from({ length: 100 }, (_, i) => i / 10);
                const y = x.map(val => Math.sin(val) * Math.exp(-val / 10));
                setPlotData({
                    data: [{ x, y, type: 'scatter', mode: 'lines', line: { color: '#3b82f6' } }],
                    layout: {
                        title: `${plotType} Potential`,
                        xaxis: { title: 'r' },
                        yaxis: { title: 'V(r)' }
                    }
                });
            }
        } catch (err) {
            setError(err.message);
        } finally {
            setIsGenerating(false);
        }
    }, [gridParams, plotType, potParams, generateCustomPlot]);

    const downloadPlot = useCallback(() => {
        const plotElement = document.querySelector('.js-plotly-plot');
        if (plotElement && window.Plotly) {
            window.Plotly.downloadImage(plotElement, {
                format: 'png',
                width: 1200,
                height: 800,
                filename: `plot_${Date.now()}`,
            });
        }
    }, []);

    const loadExample = useCallback((expr) => {
        setCustomFunction(expr);
        setError(null);
    }, []);

    const plotLayout = useMemo(() => ({
        ...plotData?.layout,
        autosize: true,
        height: 600,
    }), [plotData]);

    return (
        <div className="plot-generator">
            {/* Header */}
            <header className="header">
                <h1 className="title">Interactive Plot Generator</h1>
                <p className="subtitle">Generate plots for functions</p>
            </header>

            <div className="CalcContainer">
                <div className="content-grid">

                    {/* Control Panel */}
                    <aside className="sidebar-calc">
                        <h2 className="panel-title">
                            <BarChart3 className="icon" />
                            <span>Plot Settings</span>
                        </h2>

                        <GridParametersSection
                            gridParams={gridParams}
                            updateGridParam={updateGridParam}
                        />

                        <PlotTypeSelector
                            plotType={plotType}
                            setPlotType={setPlotType}
                        />

                        {plotType === 'Custom' && (
                            <CustomFunctionInput
                                customFunction={customFunction}
                                setCustomFunction={setCustomFunction}
                                loadExample={loadExample}
                                error={error}
                            />
                        )}

                        {plotType !== 'Custom' && (
                            <DynamicParametersSection
                                plotType={plotType}
                                potParams={potParams}
                                updatePotParam={updatePotParam}
                            />
                        )}

                        <GenerateButton
                            isGenerating={isGenerating}
                            onClick={generatePlot}
                        />

                        {error && (
                            <div className="alert alert-error">
                                {error}
                            </div>
                        )}

                        <InfoBox />
                    </aside>

                    {/* Plot Display */}
                    <main className="plot-area">
                        <PlotDisplay
                            plotData={plotData}
                            isGenerating={isGenerating}
                            plotLayout={plotLayout}
                            downloadPlot={downloadPlot}
                        />
                    </main>
                </div>
            </div>
        </div>
    );
}