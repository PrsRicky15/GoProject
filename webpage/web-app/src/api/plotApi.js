const API_BASE_URL = 'http://localhost:8080/api';

class PlotAPI {
  async generatePlotData(setGrid, plotType, potParams) {
    const response = await fetch(`${API_BASE_URL}/plots/data`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
          grid: setGrid,
        plot_type: plotType,
        parameters: potParams,
      }),
    });

    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }

    return await response.json();
  }
}

export default new PlotAPI();