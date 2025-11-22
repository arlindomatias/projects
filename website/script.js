document.getElementById("calcular").addEventListener("click", () => {
    const texto = document.getElementById("dados").value;

    const valores = texto
        .split(",")
        .map(v => parseFloat(v.trim()))
        .filter(v => !isNaN(v));

    if (valores.length === 0) {
        document.getElementById("resultado").innerHTML =
            "Nenhum valor numérico válido foi encontrado.";
        return;
    }

    const media = valores.reduce((a, b) => a + b, 0) / valores.length;

    const variancia =
        valores.reduce((acc, v) => acc + Math.pow(v - media, 2), 0) /
        (valores.length - 1);

    const desvio = Math.sqrt(variancia);

    document.getElementById("resultado").innerHTML = `
        <p><strong>N =</strong> ${valores.length}</p>
        <p><strong>Média =</strong> ${media.toFixed(4)}</p>
        <p><strong>Desvio-padrão =</strong> ${desvio.toFixed(4)}</p>
    `;
});
