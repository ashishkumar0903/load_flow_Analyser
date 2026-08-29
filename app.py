from flask import Flask, render_template, request, jsonify
import numpy as np
import time

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

def build_ybus(n, lines):
    Ybus = np.zeros((n, n), dtype=complex)
    line_params = []
    for l in lines:
        try:
            fb = int(l["from"]) - 1
            tb = int(l["to"]) - 1
            r = float(l.get("r", 0.0))
            x = float(l.get("x", 0.0001))
            b_half = float(l.get("b", 0.0))
            rate_a = float(l.get("rate_a", 100.0) or 100.0)
            
            z = complex(r, x)
            y = 1.0 / z
            y_sh = complex(0, b_half)
            
            if 0 <= fb < n and 0 <= tb < n:
                Ybus[fb, fb] += y + y_sh
                Ybus[tb, tb] += y + y_sh
                Ybus[fb, tb] -= y
                Ybus[tb, fb] -= y
            line_params.append({'from': fb, 'to': tb, 'y': y, 'b_half': y_sh, 'rate_a': rate_a})
        except (ValueError, IndexError):
            continue
    return Ybus, line_params

def validate_inputs(buses, lines, base_mva, max_iter, tol):
    errors = []
    if base_mva <= 0: errors.append("Base MVA must be > 0.")
    if max_iter <= 0: errors.append("Max Iterations must be > 0.")
    if tol <= 0: errors.append("Tolerance must be > 0.")
    if not buses: return ["At least one bus must exist."]

    bus_ids = set()
    slack_count = 0

    for idx, b in enumerate(buses):
        b_id = b.get("id")
        b_type = (b.get("type") or "").lower()
        try:
            b_id_int = int(b_id)
            if b_id_int in bus_ids: errors.append(f"Duplicate Bus ID: {b_id}.")
            bus_ids.add(b_id_int)
        except (ValueError, TypeError):
            errors.append(f"Bus row {idx+1} has an invalid Bus ID.")

        if b_type not in ["slack", "pv", "pq"]:
            errors.append(f"Bus {b_id} has an invalid type '{b_type}'.")
        if b_type == "slack": slack_count += 1

    if slack_count != 1:
        errors.append(f"Exactly one Slack Bus is required (Found: {slack_count}).")

    if not lines: errors.append("At least one line is required.")
    seen = set()
    for idx, l in enumerate(lines):
        try:
            fb, tb = int(l.get("from")), int(l.get("to"))
            if fb not in bus_ids: errors.append(f"Line {idx+1}: From-Bus {fb} does not exist.")
            if tb not in bus_ids: errors.append(f"Line {idx+1}: To-Bus {tb} does not exist.")
            if fb == tb: errors.append(f"Line {idx+1}: Self-loop detected at Bus {fb}.")
            pair = tuple(sorted([fb, tb]))
            if pair in seen: errors.append(f"Line {idx+1}: Duplicate connection between Bus {fb} and {tb}.")
            seen.add(pair)
        except (ValueError, TypeError):
            errors.append(f"Line {idx+1} contains non-numeric node IDs.")

    return errors

def solve_newton_raphson(n, Ybus, V, delta, P_spec, Q_spec, bus_types, max_iter, tol):
    G, B = Ybus.real, Ybus.imag
    pv = [i for i, t in enumerate(bus_types) if t == "pv"]
    pq = [i for i, t in enumerate(bus_types) if t == "pq"]
    non_slack = pv + pq
    history, converged = [], False

    for it in range(1, max_iter + 1):
        P_calc, Q_calc = np.zeros(n), np.zeros(n)
        for i in range(n):
            for j in range(n):
                t_ij = delta[i] - delta[j]
                P_calc[i] += V[i] * V[j] * (G[i, j] * np.cos(t_ij) + B[i, j] * np.sin(t_ij))
                Q_calc[i] += V[i] * V[j] * (G[i, j] * np.sin(t_ij) - B[i, j] * np.cos(t_ij))

        dP, dQ = P_spec - P_calc, Q_spec - Q_calc
        mismatch = np.concatenate([dP[non_slack], dQ[pq]])
        max_m = float(np.max(np.abs(mismatch))) if len(mismatch) > 0 else 0.0
        history.append({"iteration": it, "mismatch": round(max_m, 6)})

        if max_m < tol:
            converged = True
            break

        n_ns, n_pq = len(non_slack), len(pq)
        J11, J12 = np.zeros((n_ns, n_ns)), np.zeros((n_ns, n_pq))
        J21, J22 = np.zeros((n_pq, n_ns)), np.zeros((n_pq, n_pq))

        for i_idx, i in enumerate(non_slack):
            for j_idx, j in enumerate(non_slack):
                if i == j: J11[i_idx, j_idx] = -Q_calc[i] - (V[i]**2) * B[i, i]
                else:
                    t_ij = delta[i] - delta[j]
                    J11[i_idx, j_idx] = V[i] * V[j] * (G[i, j] * np.sin(t_ij) - B[i, j] * np.cos(t_ij))

        for i_idx, i in enumerate(non_slack):
            for j_idx, j in enumerate(pq):
                if i == j: J12[i_idx, j_idx] = (P_calc[i] / V[i]) + G[i, i] * V[i]
                else:
                    t_ij = delta[i] - delta[j]
                    J12[i_idx, j_idx] = V[i] * (G[i, j] * np.cos(t_ij) + B[i, j] * np.sin(t_ij))

        for i_idx, i in enumerate(pq):
            for j_idx, j in enumerate(non_slack):
                if i == j: J21[i_idx, j_idx] = P_calc[i] - (V[i]**2) * G[i, i]
                else:
                    t_ij = delta[i] - delta[j]
                    J21[i_idx, j_idx] = -V[i] * V[j] * (G[i, j] * np.cos(t_ij) + B[i, j] * np.sin(t_ij))

        for i_idx, i in enumerate(pq):
            for j_idx, j in enumerate(pq):
                if i == j: J22[i_idx, j_idx] = (Q_calc[i] / V[i]) - B[i, i] * V[i]
                else:
                    t_ij = delta[i] - delta[j]
                    J22[i_idx, j_idx] = V[i] * (G[i, j] * np.sin(t_ij) - B[i, j] * np.cos(t_ij))

        J = np.block([[J11, J12], [J21, J22]])
        try: dx = np.linalg.solve(J, mismatch)
        except np.linalg.LinAlgError: break

        delta[non_slack] += dx[:n_ns]
        V[pq] += dx[n_ns:]

    return V, delta, converged, history

def solve_gauss_seidel(n, Ybus, V, delta, P_spec, Q_spec, bus_types, max_iter, tol):
    V_complex = V * np.exp(1j * delta)
    history, converged = [], False
    for it in range(1, max_iter + 1):
        max_diff = 0.0
        for i in range(n):
            if bus_types[i] == "slack": continue
            if bus_types[i] == "pv":
                Q_calc = -np.imag(np.conj(V_complex[i]) * sum(Ybus[i, j] * V_complex[j] for j in range(n)))
                S_spec = complex(P_spec[i], Q_calc)
            else: S_spec = complex(P_spec[i], Q_spec[i])

            sum_YV = sum(Ybus[i, j] * V_complex[j] for j in range(n) if j != i)
            V_new = (1.0 / Ybus[i, i]) * ((np.conj(S_spec) / np.conj(V_complex[i])) - sum_YV)
            if bus_types[i] == "pv": V_new = V[i] * (V_new / abs(V_new))

            diff = abs(V_new - V_complex[i])
            if diff > max_diff: max_diff = diff
            V_complex[i] = V_new

        history.append({"iteration": it, "mismatch": round(float(max_diff), 6)})
        if max_diff < tol:
            converged = True
            break
    return np.abs(V_complex), np.angle(V_complex), converged, history

def solve_fast_decoupled(n, Ybus, V, delta, P_spec, Q_spec, bus_types, max_iter, tol):
    B_mat = Ybus.imag
    pv = [i for i, t in enumerate(bus_types) if t == "pv"]
    pq = [i for i, t in enumerate(bus_types) if t == "pq"]
    non_slack = pv + pq

    B_prime = -B_mat[np.ix_(non_slack, non_slack)]
    B_dprime = -B_mat[np.ix_(pq, pq)]
    history, converged = [], False

    for it in range(1, max_iter + 1):
        G, B = Ybus.real, Ybus.imag
        P_calc, Q_calc = np.zeros(n), np.zeros(n)
        for i in range(n):
            for j in range(n):
                t_ij = delta[i] - delta[j]
                P_calc[i] += V[i] * V[j] * (G[i, j] * np.cos(t_ij) + B[i, j] * np.sin(t_ij))
                Q_calc[i] += V[i] * V[j] * (G[i, j] * np.sin(t_ij) - B[i, j] * np.cos(t_ij))

        dP, dQ = (P_spec - P_calc) / V, (Q_spec - Q_calc) / V
        mismatch = np.concatenate([dP[non_slack], dQ[pq]])
        max_m = float(np.max(np.abs(mismatch))) if len(mismatch) > 0 else 0.0
        history.append({"iteration": it, "mismatch": round(max_m, 6)})

        if max_m < tol:
            converged = True
            break

        try:
            delta[non_slack] += np.linalg.solve(B_prime, dP[non_slack])
            V[pq] += np.linalg.solve(B_dprime, dQ[pq])
        except np.linalg.LinAlgError: break

    return V, delta, converged, history

def execute_solver_core(method, buses, lines, base_mva, max_iter, tol, load_factor=1.0):
    n = len(buses)
    Ybus, line_params = build_ybus(n, lines)

    V, delta = np.zeros(n), np.zeros(n)
    P_spec, Q_spec = np.zeros(n), np.zeros(n)
    bus_types = []

    for i, b in enumerate(buses):
        b_type = b.get("type", "pq").lower()
        bus_types.append(b_type)
        v_val = float(b.get("v", 1.0)) if b.get("v") not in [None, "", "—"] else 1.0
        ang_val = float(b.get("ang", 0.0)) if b.get("ang") not in [None, "", "—"] else 0.0

        if b_type == "slack": V[i], delta[i] = v_val, np.radians(ang_val)
        elif b_type == "pv": V[i], delta[i] = v_val, 0.0
        else: V[i], delta[i] = (1.0 if v_val == 0 else v_val), 0.0

        pl = (float(b.get("pl", 0) or 0) * load_factor)
        ql = (float(b.get("ql", 0) or 0) * load_factor)
        P_spec[i] = (float(b.get("pg", 0) or 0) - pl) / base_mva
        Q_spec[i] = (float(b.get("qg", 0) or 0) - ql) / base_mva

    t0 = time.time()
    if method == "gauss_seidel":
        V, delta, converged, history = solve_gauss_seidel(n, Ybus, V, delta, P_spec, Q_spec, bus_types, max_iter, tol)
    elif method == "fast_decoupled":
        V, delta, converged, history = solve_fast_decoupled(n, Ybus, V, delta, P_spec, Q_spec, bus_types, max_iter, tol)
    else:
        V, delta, converged, history = solve_newton_raphson(n, Ybus, V, delta, P_spec, Q_spec, bus_types, max_iter, tol)
    exec_time = round(time.time() - t0, 4)

    return V, delta, converged, history, exec_time, line_params

@app.route('/calculate', methods=['POST'])
def calculate():
    data = request.json
    method = data.get('method', 'newton_raphson')
    buses, lines = data.get('buses', []), data.get('lines', [])
    
    try:
        base_mva = float(data.get('baseMva', 100))
        max_iter = int(data.get('maxIter', 20))
        tol = float(data.get('tolerance', 0.001))
    except (ValueError, TypeError):
        return jsonify({"errors": ["Parameters must be numeric."]}), 400

    errors = validate_inputs(buses, lines, base_mva, max_iter, tol)
    if errors: return jsonify({"errors": errors}), 400

    n = len(buses)
    V, delta, converged, history, exec_time, line_params = execute_solver_core(method, buses, lines, base_mva, max_iter, tol)

    Ybus, _ = build_ybus(n, lines)
    G, B = Ybus.real, Ybus.imag
    bus_results = []
    total_pg, total_pl, total_qg, total_ql = 0.0, 0.0, 0.0, 0.0

    low_v_count, high_v_count, normal_v_count = 0, 0, 0

    for i in range(n):
        pl_mw = float(buses[i].get("pl", 0) or 0)
        ql_mvar = float(buses[i].get("ql", 0) or 0)
        v_min = float(buses[i].get("v_min", 0.94) or 0.94)
        v_max = float(buses[i].get("v_max", 1.06) or 1.06)

        P_calc_i, Q_calc_i = 0.0, 0.0
        for j in range(n):
            t_ij = delta[i] - delta[j]
            P_calc_i += V[i] * V[j] * (G[i, j] * np.cos(t_ij) + B[i, j] * np.sin(t_ij))
            Q_calc_i += V[i] * V[j] * (G[i, j] * np.sin(t_ij) - B[i, j] * np.cos(t_ij))

        b_type = (buses[i].get("type") or "pq").lower()
        if b_type == "slack":
            pg_mw = round((P_calc_i + pl_mw / base_mva) * base_mva, 2)
            qg_mvar = round((Q_calc_i + ql_mvar / base_mva) * base_mva, 2)
        elif b_type == "pv":
            pg_mw = float(buses[i].get("pg", 0) or 0)
            qg_mvar = round((Q_calc_i + ql_mvar / base_mva) * base_mva, 2)
        else:
            pg_mw = float(buses[i].get("pg", 0) or 0)
            qg_mvar = float(buses[i].get("qg", 0) or 0)

        total_pg += pg_mw; total_pl += pl_mw
        total_qg += qg_mvar; total_ql += ql_mvar

        v_mag = round(float(V[i]), 4)
        if v_mag < v_min:
            v_status = "LOW"; low_v_count += 1
        elif v_mag > v_max:
            v_status = "HIGH"; high_v_count += 1
        else:
            v_status = "NORMAL"; normal_v_count += 1

        bus_results.append({
            "bus_id": buses[i]["id"], "type": buses[i]["type"],
            "voltage": v_mag, "angle_deg": round(float(np.degrees(delta[i])), 3),
            "v_min": v_min, "v_max": v_max, "v_status": v_status,
            "p_gen": pg_mw, "q_gen": qg_mvar, "p_load": pl_mw, "q_load": ql_mvar
        })

    line_results = []
    total_p_loss, total_q_loss = 0.0, 0.0
    overloaded_lines_count = 0

    for lp in line_params:
        fb, tb = lp['from'], lp['to']
        V_f = V[fb] * np.exp(1j * delta[fb])
        V_t = V[tb] * np.exp(1j * delta[tb])

        I_ft = (V_f - V_t) * lp['y'] + V_f * lp['b_half']
        S_ft = V_f * np.conj(I_ft) * base_mva
        I_tf = (V_t - V_f) * lp['y'] + V_t * lp['b_half']
        S_tf = V_t * np.conj(I_tf) * base_mva

        S_loss = S_ft + S_tf
        p_loss, q_loss = round(abs(S_loss.real), 3), round(abs(S_loss.imag), 3)
        total_p_loss += p_loss; total_q_loss += q_loss

        mva_flow = abs(S_ft)
        loading_pct = round((mva_flow / lp['rate_a']) * 100, 1)

        if loading_pct > 100:
            load_status = "OVERLOADED"; overloaded_lines_count += 1
        elif loading_pct >= 90: load_status = "CRITICAL"
        elif loading_pct >= 70: load_status = "WARNING"
        else: load_status = "NORMAL"

        line_results.append({
            "from_bus": buses[fb]["id"], "to_bus": buses[tb]["id"],
            "p_flow": round(S_ft.real, 2), "q_flow": round(S_ft.imag, 2),
            "line_current": round(abs(I_ft), 4), "p_loss": p_loss, "q_loss": q_loss,
            "rate_a": lp['rate_a'], "loading_pct": loading_pct, "load_status": load_status
        })

    v_array = [b["voltage"] for b in bus_results]
    min_v = min(v_array) if v_array else 0
    max_v = max(v_array) if v_array else 0

    return jsonify({
        "status": "CONVERGED" if converged else "FAILED",
        "iterations": len(history), "exec_time": exec_time, "convergence_history": history,
        "total_pg": round(total_pg, 2), "total_pl": round(total_pl, 2),
        "total_qg": round(total_qg, 2), "total_ql": round(total_ql, 2),
        "total_p_loss": round(total_p_loss, 3), "total_q_loss": round(total_q_loss, 3),
        "min_voltage": min_v, "max_voltage": max_v,
        "voltage_summary": { "normal": normal_v_count, "low": low_v_count, "high": high_v_count, "violations": low_v_count + high_v_count },
        "overloaded_lines_count": overloaded_lines_count,
        "bus_results": bus_results, "line_results": line_results
    })

@app.route('/stress_analysis', methods=['POST'])
def stress_analysis():
    data = request.json
    method = data.get('method', 'newton_raphson')
    buses, lines = data.get('buses', []), data.get('lines', [])
    start_pct = float(data.get('start_pct', 50))
    end_pct = float(data.get('end_pct', 150))
    step_pct = float(data.get('step_pct', 5))

    base_mva = float(data.get('baseMva', 100))
    max_iter = int(data.get('maxIter', 20))
    tol = float(data.get('tolerance', 0.001))

    factors = np.arange(start_pct, end_pct + step_pct / 2.0, step_pct)
    curve_data = []

    for pct in factors:
        lf = pct / 100.0
        V, _, converged, history, _, _ = execute_solver_core(method, buses, lines, base_mva, max_iter, tol, load_factor=lf)
        point = {"load_pct": round(pct, 1), "converged": converged, "iterations": len(history), "buses": {}}
        for i, b in enumerate(buses):
            point["buses"][b["id"]] = round(float(V[i]), 4)
        curve_data.append(point)

    return jsonify({"stress_curve": curve_data})

@app.route('/compare_methods', methods=['POST'])
def compare_methods():
    data = request.json
    buses, lines = data.get('buses', []), data.get('lines', [])
    base_mva = float(data.get('baseMva', 100))
    max_iter = int(data.get('maxIter', 100))
    tol = float(data.get('tolerance', 0.001))

    methods = [
        ("Newton-Raphson", "newton_raphson"),
        ("Gauss-Seidel", "gauss_seidel"),
        ("Fast Decoupled", "fast_decoupled")
    ]
    comparison = []

    for name, key in methods:
        V, _, converged, history, exec_time, _ = execute_solver_core(key, buses, lines, base_mva, max_iter, tol)
        last_mismatch = history[-1]["mismatch"] if history else 0.0
        comparison.append({
            "method": name, "iterations": len(history),
            "final_error": last_mismatch, "status": "CONVERGED" if converged else "FAILED",
            "exec_time": exec_time
        })

    return jsonify({"comparison": comparison})

if __name__ == '__main__':
    app.run(debug=True)