import {calc} from "./calculator";
import {plot} from "./plotter";

export const examples = {
    basicCalculator: () => {
        const calculator = new calc.ScientificCalculator();
        calculator.add(5, 3);
        calculator.power(2, 8);
        return calculator.getResult();
    },

    plotSine: () => {
        const evaluator = new plot.FunctionEvaluator();
        const grid = plot.GridUtils.createLinearGrid(0, 2 * Math.PI, 100);
        return evaluator.evaluateOnGrid(plot.Functions.sine, grid);
    },

    plotMultiple: () => {
        const evaluator = new plot.FunctionEvaluator();
        const grid = plot.GridUtils.createLinearGrid(-5, 5, 200);
        return evaluator.evaluateMultiple([
            plot.Functions.sine,
            plot.Functions.cosine,
            plot.Functions.gaussian(1)
        ], grid);
    },

    parametricCircle: () => {
        const plotter = new plot.ParametricPlotter();
        return plotter.circle(5, 100);
    },

    dataAnalysis: () => {
        const data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        return {
            mean: plot.DataAnalysis.mean(data),
            median: plot.DataAnalysis.median(data),
            stdDev: plot.DataAnalysis.stdDev(data),
            range: plot.DataAnalysis.range(data)
        };
    }
};