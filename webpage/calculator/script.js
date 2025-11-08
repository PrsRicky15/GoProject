// Get the calculator display
const display = document.getElementById('display');

// Add value to the display when a button is clicked
function appendValue(value) {
    display.value += value;
}

// Clear the display (for 'C' button)
function clearDisplay() {
    display.value = '';
}

// Evaluate the expression on the display
function calculate() {
    try {
        // Use eval to compute the result
        // Math functions like Math.sin(), Math.log() are supported
        display.value = eval(display.value);
    } catch (error) {
        display.value = 'Error';
    }
}
