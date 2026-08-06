function greeting() {
    const greetEl = document.getElementById("greet");
    greetEl.style.transform = "scale(1.1)";
    greetEl.innerHTML = "GoodBye";
    setTimeout(() => {
        greetEl.style.transform = "scale(1)";
    }, 300);
}

document.addEventListener('mousemove', (e) => {
    const items = document.querySelectorAll('.grid-item');
    const x = e.clientX / window.innerWidth;
    const y = e.clientY / window.innerHeight;

    items.forEach((item, index) => {
        const speed = (index % 3 + 1) * 2;
        const xOffset = (x - 0.5) * speed;
        const yOffset = (y - 0.5) * speed;
        item.style.transform = `translate(${xOffset}px, ${yOffset}px)`;
    });
});