export const display = document.getElementById('display');

// Public Calculator API
export const Calculator = {
    appendValue(values) {
        display.value += values;
    },

    clearDisplay() {
        display.value = '';
    },

    calculate() {
        try {
            display.value = eval(display.value);
        } catch (error) {
            display.value = 'Error';
        }
    }
};