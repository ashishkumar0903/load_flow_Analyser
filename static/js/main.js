function addBusRow() {
    const table = document.getElementById("busTable").getElementsByTagName('tbody')[0];
    const newRow = table.insertRow();
    newRow.innerHTML = `
        <td><input type="number" class="b-id" value="${table.rows.length + 1}"></td>
        <td>
            <select class="b-type">
                <option value="slack">Slack Bus</option>
                <option value="pv">PV Bus</option>
                <option value="pq" selected>PQ Bus</option>
            </select>
        </td>
        <td><input type="number" step="0.01" value="1.0" class="b-v"></td>
        <td><input type="number" step="0.1" value="0.0" class="b-ang"></td>
        <td><input type="number" step="0.1" value="0" class="b-pg"></td>
        <td><input type="number" step="0.1" value="0" class="b-qg"></td>
        <td><input type="number" step="0.1" value="0" class="b-pl"></td>
        <td><input type="number" step="0.1" value="0" class="b-ql"></td>
    `;
    updateKPICounts();
}

function addLineRow() {
    const table = document.getElementById("lineTable").getElementsByTagName('tbody')[0];
    const newRow = table.insertRow();
    newRow.innerHTML = `
        <td><input type="number" class="l-from" value="1"></td>
        <td><input type="number" class="l-to" value="2"></td>
        <td><input type="number" step="0.001" value="0.00" class="l-r"></td>
        <td><input type="number" step="0.001" value="0.00" class="l-x"></td>
        <td><input type="number" step="0.001" value="0.00" class="l-b"></td>
    `;
    updateKPICounts();
}

function updateKPICounts() {
    const busCount = document.querySelectorAll("#busTable tbody tr").length;
    const lineCount = document.querySelectorAll("#lineTable tbody tr").length;
    
    let totalPl = 0;
    document.querySelectorAll(".b-pl").forEach(input => {
        totalPl += parseFloat(input.value) || 0;
    });

    document.getElementById("kpiBuses").innerText = busCount;
    document.getElementById("kpiLines").innerText = lineCount;
    document.getElementById("kpiMwLoad").innerText = totalPl.toFixed(2) + " MW";
}

function runCalculation() {
    const method = document.getElementById("method").value;
    
    const busRows = document.querySelectorAll("#busTable tbody tr");
    const buses = [];
    busRows.forEach(row => {
        buses.push({
            id: row.querySelector(".b-id").value,
            type: row.querySelector(".b-type").value,
            v: row.querySelector(".b-v").value,
            ang: row.querySelector(".b-ang").value,
            pg: row.querySelector(".b-pg").value,
            qg: row.querySelector(".b-qg").value,
            pl: row.querySelector(".b-pl").value,
            ql: row.querySelector(".b-ql").value
        });
    });

    const lineRows = document.querySelectorAll("#lineTable tbody tr");
    const lines = [];
    lineRows.forEach(row => {
        lines.push({
            from: row.querySelector(".l-from").value,
            to: row.querySelector(".l-to").value,
            r: row.querySelector(".l-r").value,
            x: row.querySelector(".l-x").value,
            b: row.querySelector(".l-b").value
        });
    });

    fetch('https://load-flow-analysis.onrender.com  /calculate', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ method, buses, lines })
    })
    .then(res => res.json())
    .then(data => {
        document.getElementById("resultsCard").style.display = "block";
        document.getElementById("resStatus").innerText = data.status;
        document.getElementById("resIter").innerText = data.iterations;

        const tbody = document.querySelector("#resultsTable tbody");
        tbody.innerHTML = "";
        data.bus_results.forEach(res => {
            const row = tbody.insertRow();
            row.innerHTML = `
                <td><strong>${res.bus_id}</strong></td>
                <td>${res.voltage.toFixed(4)}</td>
                <td>${res.angle_deg.toFixed(2)}°</td>
                <td>${res.p_mw.toFixed(2)}</td>
                <td>${res.q_mvar.toFixed(2)}</td>
            `;
        });
    });
}

// Initialize KPI metrics on load
document.addEventListener("DOMContentLoaded", updateKPICounts);