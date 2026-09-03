# ⚡ Load flow Analyser

An interactive, full-stack web application designed for real-time power system load flow calculations, grid stability diagnostics, and numerical solver benchmarking.

---

## 📌 Introduction

**Load flow Analyser** is a browser-accessible power system analysis platform that eliminates the need for heavy desktop software installations. Built on a high-performance Python backend and an interactive JavaScript frontend, the tool solves non-linear power flow equations for multi-bus networks using standard numerical methods. 

Engineers, researchers, and students can model network topologies, evaluate bus voltage profiles, monitor active and reactive transmission losses, perform $V$-$P$ stability stress testing, and instantly detect line thermal overloads.

---

## ✨ Key Features

* **⚡ Multi-Algorithm Solvers:** Executes **Newton-Raphson (NR)**, **Gauss-Seidel (GS)**, and **Fast Decoupled Load Flow (FDLF)** algorithms with configurable convergence tolerances and iteration limits.
* **🕸 Interactive Network Graph:** Renders real-time node-based network topology diagrams using `vis-network` to visualize Slack, PV, and PQ buses along with connecting transmission lines.
* **📈 Voltage Stability Stress Testing:** Automatically scales system load factors ($50\%$ to $150\%$) to calculate system stress and plot $V$-$P$ nose curves for voltage collapse threshold identification.
* **🚨 Automated Health Diagnostics:** Instantly flags operational limit violations, including bus over/under-voltages ($V_{\text{min}} / V_{\text{max}}$) and line thermal capacities ($>100\%$ MVA limits).
* **📊 Visual Analytics:** Interactive Chart.js charts providing iteration-by-iteration convergence tracking, voltage magnitude profiles, and active/reactive power loss distributions.

---

## 🛠 Tech Stack

* **Backend:** Python 3, Flask REST API, NumPy (Complex admittance matrix computations & linear algebra)
* **Frontend:** HTML5, CSS3 (Modern Flexbox/Grid), JavaScript (ES6+)
* **Visualization Libraries:** Chart.js (Data analytics), Vis.js (Interactive graph layout)
* **Data Format:** JSON-based network import/export schemas

---

## 📐 System Architecture
┌─────────────────────────────────────────────────────────────┐
│                     Client Browser (UI/UX)                   │
│  • HTML5 / Responsive CSS                                   │
│  • Chart.js (Convergence, Voltage, Power & Loss Graphs)     │
│  • vis-network (Interactive Bus/Line Topology Diagram)       │
└──────────────────────────────┬──────────────────────────────┘
│
REST API / JSON
│
┌──────────────────────────────▼──────────────────────────────┐
│                      Flask Python Backend                   │
│  • Python Engine & NumPy Linear Algebra Engine               │
│  • Y-Bus Admittance Matrix Construction                     │
│  • Newton-Raphson, Gauss-Seidel & Fast Decoupled Solvers     │
│  • Voltage Stability & Thermal Limit Diagnostic Engines     │
└─────────────────────────────────────────────────────────────┘

---
