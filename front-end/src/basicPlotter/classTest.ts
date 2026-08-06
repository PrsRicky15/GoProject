import {
    FunctionEvaluator, ParametricPlotter,
    DataAnalysis, GridUtils, type PlotData,
} from "./plotter";

import { ScientificCalculator, MathFunctions } from "./calculator_and_MathFunc";

class Hello {
    constructor(private greeting: string) { }
    greet() {
        console.log(this.greeting);
    }
}

class User extends Hello {
    public hobbies: string[] = [];
    public homeTown: string = "";
    private name: string;
    private age: number;

    constructor(name: string, age: number) {
        super("Hello User!");
        this.name = name;
        this.age = age;
    }

    override greet() {
        super.greet();
        console.log(`Hello ${this.name}, you are ${this.age} years old`);
    }

    getName() {
        return this.name;
    }
    getAge() {
        return this.age;
    }
    getHomeTown() {
        return this.homeTown;
    }
    getHobbies() {
        return this.hobbies;
    }
    setName(name: string) {
        this.name = name;
    }
    setAge(age: number) {
        this.age = age;
    }

    print() {
        console.log(`Name: ${this.name}, Age: ${this.age}, HomeTown: ${this.homeTown}`);
    }

    printHobbies() {
        console.log(`Hobbies: ${this.hobbies}`);
    }

    static count: number = 0;
    static incrementCount() {
        User.count++;
    }
    static decrementCount() {
        User.count--;
    }
    static getCount() {
        return User.count;
    }

    static printCount() {
        console.log(`Count: ${User.count}`);
    }
}

export const examples = {
    basicCalculator: (): number => {
        const calculator = new ScientificCalculator();
        calculator.add(5, 3);
        calculator.power(2, 8);
        return calculator.getResult();
    },

    advancedCalculator: () => {
        const calculator = new ScientificCalculator();
        return {
            sinh: calculator.sinh(1),
            arcSine: calculator.arcSine(0.5),
            arcTanh: calculator.arcTanh(0.5),
            combination: calculator.combination(10, 3),
            abs: calculator.abs(-42),
        };
    },

    plotSine: (): PlotData => {
        const evaluator = new FunctionEvaluator();
        const grid = GridUtils.createLinearGrid(0, 2 * Math.PI, 100);
        return evaluator.evaluateOnGrid(MathFunctions.sine, grid);
    },

    plotMultiple: (): PlotData[] => {
        const evaluator = new FunctionEvaluator();
        const grid = GridUtils.createLinearGrid(-5, 5, 200);
        return evaluator.evaluateMultiple(
            [
                MathFunctions.sine,
                MathFunctions.cosine,
                MathFunctions.gaussian(1),
            ],
            grid
        );
    },

    plotAdvancedFunctions: (): PlotData[] => {
        const evaluator = new FunctionEvaluator();
        const grid = GridUtils.createLinearGrid(-5, 5, 200);
        return evaluator.evaluateMultiple(
            [
                MathFunctions.abs,
                MathFunctions.sinh,
                MathFunctions.arcTanh,
                MathFunctions.dampedSine(2, 0.5),
                MathFunctions.sigmoid(1),
            ],
            grid
        );
    },

    parametricCircle: (): PlotData => {
        const plotter = new ParametricPlotter();
        return plotter.circle(5, 100);
    },

    dataAnalysis: () => {
        const data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        return {
            mean: DataAnalysis.mean(data),
            median: DataAnalysis.median(data),
            stdDev: DataAnalysis.stdDev(data),
            range: DataAnalysis.range(data),
        };
    },
};