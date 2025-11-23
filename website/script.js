// ==========================
// Funções auxiliares
// ==========================

// converte texto → números aceitando várias formas de separação
function parseNumbers(text) {
    let cleaned = text
        .replace(/,/g, '.')      // vírgula → ponto decimal
        .replace(/\t/g, ' ')     // tab → espaço
        .replace(/;/g, ' ')      // ponto e vírgula → espaço
        .replace(/\n/g, ' ')     // quebra de linha → espaço
        .replace(/\s+/g, ' ');   // múltiplos espaços → um só

    return cleaned
        .trim()
        .split(' ')
        .map(Number)
        .filter(x => !isNaN(x));
}

// média
function mean(arr) {
    return arr.reduce((a,b) => a+b, 0) / arr.length;
}

// mediana
function median(arr) {
    let a = [...arr].sort((a,b)=>a-b);
    let mid = Math.floor(a.length/2);
    return a.length % 2 === 0 ? (a[mid-1]+a[mid])/2 : a[mid];
}

// desvio padrão populacional
function sd(arr) {
    let m = mean(arr);
    return Math.sqrt(arr.reduce((s,x) => s + (x-m)**2, 0) / arr.length);
}

// erro padrão da média
function sem(arr) {
    return sd(arr) / Math.sqrt(arr.length);
}

// moda
function mode(arr) {
    let counts = {};
    arr.forEach(x => counts[x] = (counts[x]||0)+1);
    let maxCount = Math.max(...Object.values(counts));
    let modes = Object.keys(counts).filter(x => counts[x] === maxCount);
    return modes.join(', ');
}

// intervalo interquartil (Q3 - Q1)
function iqr(arr) {
    let a = [...arr].sort((a,b)=>a-b);
    let q1 = median(a.slice(0, Math.floor(a.length/2)));
    let q3 = median(a.slice(Math.ceil(a.length/2)));
    return q3 - q1;
}

// ==========================
// Controle dos grupos
// ==========================

let groupCount = 0;

function addGroup() {
    groupCount++;
    const div = document.createElement("div");
    div.className = "group";
    div.id = "group" + groupCount;

    div.innerHTML = `
        <h4>Grupo ${groupCount}</h4>
        <textarea placeholder="Insira os dados"></textarea>
        <br>
        <button onclick="removeGroup(${groupCount})">Remover grupo</button>
    `;

    document.getElementById("groups").appendChild(div);
}

function removeGroup(n) {
    document.getElementById("group" + n).remove();
}

// ==========================
// Análise dos dados
// ==========================

function analyze() {
    let groups = [];
    let labels = [];

    document.querySelectorAll(".group").forEach((group, index) => {
        let text = group.querySelector("textarea").value;
        let values = parseNumbers(text);

        if (values.length > 0) {
            labels.push("Grupo " + (index + 1));

            groups.push({
                mean: mean(values),
                median: median(values),
                sd: sd(values),
                sem: sem(values),
                mode: mode(values),
                iqr: iqr(values),
                min: Math.min(...values),
                max: Math.max(...values),
                sum: values.reduce((a,b)=>a+b,0),
                n: values.length
            });
        }
    });

    if (groups.length === 0) {
        document.getElementById("results").innerHTML = 
            "<p>Nenhum grupo válido encontrado.</p>";
        return;
    }

    let rows = [
        ["Média",            g => g.mean.toFixed(4)],
        ["Mediana",          g => g.median.toFixed(4)],
        ["Desvio Padrão",    g => g.sd.toFixed(4)],
        ["Erro-Padrão",      g => g.sem.toFixed(4)],
        ["Moda",             g => g.mode],
        ["Intervalo IQ",     g => g.iqr.toFixed(4)],
        ["Mínimo",           g => g.min],
        ["Máximo",           g => g.max],
        ["Soma",             g => g.sum.toFixed(4)],
        ["N",                g => g.n]
    ];

    let table = `
        <table>
            <tr>
                <th>Estatística</th>
                ${labels.map(l => `<th>${l}</th>`).join("")}
            </tr>
    `;

    rows.forEach(([label, getter]) => {
        table += `
            <tr>
                <th>${label}</th>
                ${groups.map(g => `<td>${getter(g)}</td>`).join("")}
            </tr>
        `;
    });

    table += "</table>";

    document.getElementById("results").innerHTML = table;
}

addGroup();  // adiciona um grupo inicial

// ==========================
// Requisição ao backend R para gerar gráficos
// ==========================
function collectGroups() {
    const result = [];

    document.querySelectorAll(".group textarea").forEach(t => {
        const values = parseNumbers(t.value);
        if (values.length > 0) result.push(values);
    });

    return result;
}

async function requestPlots() {
    const groups = collectGroups();  // você já possui isso
    if (groups.length === 0) {
        alert("Insira ao menos um grupo.");
        return;
    }

    try {
        const resp = await fetch("http://127.0.0.1:8000/plots", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ groups })
        });

        if (!resp.ok) {
            alert("Erro do servidor: " + resp.status);
            return;
        }

        // backend retorna PNG -> precisamos ler como blob
        const blob = await resp.blob();
        const url = URL.createObjectURL(blob);

        // abrir em nova aba
        window.open(url, "_blank");

    } catch (err) {
        console.error(err);
        alert("Falha ao comunicar com o backend: " + err.message);
    }
}
