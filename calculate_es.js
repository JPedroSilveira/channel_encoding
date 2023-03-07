const coordenadas = [
    [1,0], 
    [1/Math.sqrt(2),  1/Math.sqrt(2)], 
    [0, 1],
    [-1/Math.sqrt(2), 1/Math.sqrt(2)],
    [-1, 0],
    [-1/Math.sqrt(2), -1/Math.sqrt(2)],
    [0, -1],
    [1/Math.sqrt(2), -1/Math.sqrt(2)],
    [Math.sqrt(2 + Math.sqrt(2))/2*Math.sqrt(2), Math.sqrt(2 - Math.sqrt(2))/2*Math.sqrt(2)],
    [Math.sqrt(2 - Math.sqrt(2))/2*Math.sqrt(2), Math.sqrt(2 + Math.sqrt(2))/2*Math.sqrt(2)],
    [-Math.sqrt(2 - Math.sqrt(2))/2*Math.sqrt(2), Math.sqrt(2 + Math.sqrt(2))/2*Math.sqrt(2)],
    [-Math.sqrt(2 + Math.sqrt(2))/2*Math.sqrt(2), Math.sqrt(2 - Math.sqrt(2))/2*Math.sqrt(2)],
    [-Math.sqrt(2 + Math.sqrt(2))/2*Math.sqrt(2), -Math.sqrt(2 - Math.sqrt(2))/2*Math.sqrt(2)],
    [-Math.sqrt(2 - Math.sqrt(2))/2*Math.sqrt(2), -Math.sqrt(2 + Math.sqrt(2))/2*Math.sqrt(2)],
    [Math.sqrt(2 - Math.sqrt(2))/2*Math.sqrt(2), -Math.sqrt(2 + Math.sqrt(2))/2*Math.sqrt(2)],
    [Math.sqrt(2 + Math.sqrt(2))/2*Math.sqrt(2), -Math.sqrt(2 - Math.sqrt(2))/2*Math.sqrt(2)]
]

let total = 0;
coordenadas.forEach(coordenada => {
    total += Math.pow(coordenada[0], 2) + Math.pow(coordenada[1], 2);
});

console.log(total);
