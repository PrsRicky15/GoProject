import React, { useState } from 'react';
import Plot from 'react-plotly.js';
import { BarChart3, Loader, Download, RefreshCw } from 'lucide-react';
import plotApi from '../api/plotApi';

export default function PlotGenerator() {
    const [gridParams, setGridParams] = useState({
        rMin: -0.,
        rMax: 10.,
        nGrid: 100,
    });

    const [plotType, setPlotType] = useState('Morse');
    const [isGenerating, setIsGenerating] = useState(false);
    const [plotData, setPlotData] = useState(null);
    const [error, setError] = useState(null);
    const [potParams, setPotParams] = useState({
        D: 100.0,
        a: 1.5,
        r0: 2.0,
    });


  const plotTypes = [
    { value: 'Morse', label: 'Morse Potential' },
    { value: 'Softcore', label: 'Softcore Potential' },
    { value: 'surface_3d', label: '3D Surface Plot' },
    { value: 'energy_levels', label: 'Energy Level Diagram' },
  ];

  const generatePlot = async () => {
    setIsGenerating(true);
    setError(null);

    try {
      const response = await plotApi.generatePlotData(plotType, potParams);
      setPlotData(response);
    } catch (err) {
      setError(err.message);
    } finally {
      setIsGenerating(false);
    }
  };

  const downloadPlot = () => {
    if (plotData) {
      // Plotly has built-in download functionality via the modebar
      // Or you can trigger it programmatically
      const plotElement = document.querySelector('.js-plotly-plot');
      if (plotElement) {
        window.Plotly.downloadImage(plotElement, {
          format: 'png',
          width: 1200,
          height: 800,
          filename: `plot_${Date.now()}`,
        });
      }
    }
  };

  return (
    <div className="min-h-screen min-w-screen bg-gradient-to-br from-slate-900 via-purple-900 to-slate-900 text-white p-6">
      <div className="max-w-7xl mx-auto">
        <div className="text-center mb-8">
          <h1 className="text-4xl font-bold mb-2">Interactive Plot Generator</h1>
          <p className="text-gray-400">Generate and explore quantum chemistry visualizations</p>
        </div>

        <div className="grid lg:grid-cols-4 gap-6">
          {/* Control Panel */}
          <div className="lg:col-span-1 bg-white/5 border border-white/10 rounded-2xl p-6 space-y-6 h-fit">
            <div>
              <h2 className="text-xl font-semibold mb-4 flex items-center space-x-2">
                <BarChart3 className="w-5 h-5" />
                <span>Plot Settings</span>
              </h2>
            </div>

              {/* Plot Type Selector */}
              <div>
                  <label className="block text-sm text-gray-400 mb-2">rMin</label>
                  <input
                      type="number"
                      step="0.1"
                      value={gridParams.rMin}
                      onChange={(e) =>
                          setGridParams({ ...gridParams, rMin: parseFloat(e.target.value) })
                      }
                      className="w-full px-4 py-2 bg-white/5 border border-white/10 rounded-lg text-white focus:border-purple-500 outline-none"
                  />
              </div>

              <div>
                  <label className="block text-sm text-gray-400 mb-2">rMax</label>
                  <input
                      type="number"
                      step="0.1"
                      value={gridParams.rMax}
                      onChange={(e) =>
                          setGridParams({ ...gridParams, rMax: parseFloat(e.target.value) })
                      }
                      className="w-full px-4 py-2 bg-white/5 border border-white/10 rounded-lg text-white focus:border-purple-500 outline-none"
                  />
              </div>

              <div>
                  <label className="block text-sm text-gray-400 mb-2">nGrid</label>
                  <input
                      type="number"
                      step="0.1"
                      value={gridParams.nGrid}
                      onChange={(e) =>
                          setGridParams({ ...gridParams, nGrid: parseFloat(e.target.value) })
                      }
                      className="w-full px-4 py-2 bg-white/5 border border-white/10 rounded-lg text-white focus:border-purple-500 outline-none"
                  />
              </div>

            {/* Plot Type Selector */}
            <div>
              <label className="block text-sm text-gray-400 mb-2">Plot Type</label>
              <select
                value={plotType}
                onChange={(e) => setPlotType(e.target.value)}
                className="w-full px-4 py-2 bg-white/5 border border-white/10 rounded-lg text-white focus:border-purple-500 outline-none"
              >
                {plotTypes.map((type) => (
                  <option key={type.value} value={type.value}>
                    {type.label}
                  </option>
                ))}
              </select>
            </div>

            {/* Parameters - Show only for relevant plot types */}
            {plotType === 'Morse' && (
              <div className="space-y-4">
                <div>
                  <label className="block text-sm text-gray-400 mb-2">
                    Dissociation Energy (D)
                  </label>
                  <input
                    type="number"
                    step="10"
                    value={potParams.D}
                    onChange={(e) =>
                      setPotParams({ ...potParams, D: parseFloat(e.target.value) })
                    }
                    className="w-full px-4 py-2 bg-white/5 border border-white/10 rounded-lg text-white focus:border-purple-500 outline-none"
                  />
                </div>

                <div>
                  <label className="block text-sm text-gray-400 mb-2">
                    Width Parameter (a)
                  </label>
                  <input
                    type="number"
                    step="0.1"
                    value={potParams.a}
                    onChange={(e) =>
                      setPotParams({ ...potParams, a: parseFloat(e.target.value) })
                    }
                    className="w-full px-4 py-2 bg-white/5 border border-white/10 rounded-lg text-white focus:border-purple-500 outline-none"
                  />
                </div>

                <div>
                  <label className="block text-sm text-gray-400 mb-2">
                    Equilibrium Distance (r₀)
                  </label>
                  <input
                    type="number"
                    step="0.1"
                    value={potParams.r0}
                    onChange={(e) =>
                      setPotParams({ ...potParams, r0: parseFloat(e.target.value) })
                    }
                    className="w-full px-4 py-2 bg-white/5 border border-white/10 rounded-lg text-white focus:border-purple-500 outline-none"
                  />
                </div>
              </div>
            )}

            {plotType === 'Softcore' && (
              <div className="space-y-4">
                <div>
                  <label className="block text-sm text-gray-400 mb-2">
                    Charge (q)
                  </label>
                  <input
                    type="number"
                    step="1"
                    value={potParams.D}
                    onChange={(e) =>
                      setPotParams({ ...potParams, D: parseFloat(e.target.value) })
                    }
                    className="w-full px-4 py-2 bg-white/5 border border-white/10 rounded-lg text-white focus:border-purple-500 outline-none"
                  />
                </div>

                <div>
                  <label className="block text-sm text-gray-400 mb-2">
                    Softcore Parameter (a)
                  </label>
                  <input
                    type="number"
                    step="0.1"
                    value={potParams.a}
                    onChange={(e) =>
                      setPotParams({ ...potParams, a: parseFloat(e.target.value) })
                    }
                    className="w-full px-4 py-2 bg-white/5 border border-white/10 rounded-lg text-white focus:border-purple-500 outline-none"
                  />
                </div>

                <div>
                  <label className="block text-sm text-gray-400 mb-2">
                    Center (r₀)
                  </label>
                  <input
                    type="number"
                    step="0.1"
                    value={potParams.r0}
                    onChange={(e) =>
                      setPotParams({ ...potParams, r0: parseFloat(e.target.value) })
                    }
                    className="w-full px-4 py-2 bg-white/5 border border-white/10 rounded-lg text-white focus:border-purple-500 outline-none"
                  />
                </div>
              </div>
            )}

            {/* Generate Button */}
            <button
              onClick={generatePlot}
              disabled={isGenerating}
              className="w-full px-6 py-3 bg-gradient-to-r from-purple-600 to-blue-600 rounded-lg font-medium flex items-center justify-center space-x-2 hover:shadow-lg disabled:opacity-50 disabled:cursor-not-allowed transition"
            >
              {isGenerating ? (
                <>
                  <Loader className="w-5 h-5 animate-spin" />
                  <span>Generating...</span>
                </>
              ) : (
                <>
                  <RefreshCw className="w-5 h-5" />
                  <span>Generate Plot</span>
                </>
              )}
            </button>

            {/* Error Display */}
            {error && (
              <div className="p-3 bg-red-500/10 border border-red-500/30 rounded-lg text-red-400 text-sm">
                {error}
              </div>
            )}

            {/* Info Box */}
            <div className="p-4 bg-blue-500/10 border border-blue-500/30 rounded-lg text-sm space-y-2">
              <p className="text-blue-400 font-medium">✨ Interactive Features:</p>
              <ul className="text-gray-300 space-y-1 text-xs">
                <li>• Zoom: Drag to select area</li>
                <li>• Pan: Shift + drag</li>
                <li>• Reset: Double click</li>
                <li>• Hover for values</li>
                <li>• Download via toolbar</li>
              </ul>
            </div>
          </div>

          {/* Plot Display Area */}
          <div className="lg:col-span-3 bg-white/10 border border-white/10 rounded-2xl p-6">
            <div className="flex items-center justify-between mb-4">
              <h2 className="text-xl font-semibold">Visualization</h2>
              {plotData && (
                <button
                  onClick={downloadPlot}
                  className="px-4 py-2 bg-white/10 hover:bg-white/20 rounded-lg flex items-center space-x-2 transition"
                >
                  <Download className="w-4 h-4" />
                  <span className="text-sm">Download PNG</span>
                </button>
              )}
            </div>

            <div className="bg-white/5 rounded-xl overflow-hidden">
              {isGenerating ? (
                <div className="min-h-[600px] flex items-center justify-center">
                  <div className="text-center space-y-4">
                    <Loader className="w-12 h-12 animate-spin mx-auto text-purple-400" />
                    <p className="text-gray-400">Generating plot...</p>
                  </div>
                </div>
              ) : plotData ? (
                <Plot
                  data={plotData.data}
                  layout={{
                    ...plotData.layout,
                    autosize: true,
                    height: 600,
                  }}
                  config={{
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
                  }}
                  style={{ width: '100%', height: '600px' }}
                  useResizeHandler={true}
                />
              ) : (
                <div className="min-h-[600px] flex items-center justify-center">
                  <div className="text-center space-y-4">
                    <BarChart3 className="w-16 h-16 mx-auto text-gray-600" />
                    <p className="text-gray-400">
                      Configure parameters and click "Generate Plot" to visualize
                    </p>
                  </div>
                </div>
              )}
            </div>

            {plotData && (
              <div className="mt-4 p-4 bg-green-500/10 border border-green-500/30 rounded-lg">
                <p className="text-sm text-green-400">
                  ✓ Plot generated successfully - Interact with the plot using your mouse!
                </p>
              </div>
            )}
          </div>
        </div>
      </div>
    </div>
  );
}