const coordenadas = [
    [-3, -3],
    [-3, -1],
    [-3, 3],
    [-3, 1],
    [-1, -3],
    [-1, -1],
    [-1, 3],
    [-1, 1],
    [3, -3],
    [3, -1],
    [3, 3],
    [3, 1],
    [1, -3],
    [1, -1],
    [1, 3],
    [1, 1],
]

let total = 0;
coordenadas.forEach(coordenada => {
    total += Math.pow(coordenada[0], 2) + Math.pow(coordenada[1], 2);
});

console.log(total/coordenadas.length);
