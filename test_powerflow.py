import unittest
import numpy as np
from app import validate_inputs, build_ybus, execute_solver_core

class TestPowerFlowSuite(unittest.TestCase):

    # -------------------------------------------------------------
    # 1. INPUT VALIDATION TESTS
    # -------------------------------------------------------------
    def test_missing_slack_bus(self):
        buses = [
            {"id": 1, "type": "pq", "v": 1.0, "ang": 0, "pl": 10, "ql": 5},
            {"id": 2, "type": "pq", "v": 1.0, "ang": 0, "pl": 20, "ql": 10}
        ]
        lines = [{"from": 1, "to": 2, "r": 0.01, "x": 0.05, "b": 0.0, "rate_a": 100}]
        errors = validate_inputs(buses, lines, base_mva=100, max_iter=20, tol=0.001)
        self.assertIn("Exactly one Slack Bus is required (Found: 0).", errors)

    def test_duplicate_bus_id(self):
        buses = [
            {"id": 1, "type": "slack", "v": 1.0, "ang": 0},
            {"id": 1, "type": "pq", "v": 1.0, "ang": 0}
        ]
        lines = [{"from": 1, "to": 1, "r": 0.01, "x": 0.05, "b": 0.0, "rate_a": 100}]
        errors = validate_inputs(buses, lines, base_mva=100, max_iter=20, tol=0.001)
        self.assertTrue(any("Duplicate Bus ID" in err for err in errors))

    def test_invalid_line_connection(self):
        buses = [{"id": 1, "type": "slack", "v": 1.0, "ang": 0}]
        lines = [{"from": 1, "to": 99, "r": 0.01, "x": 0.05, "b": 0.0, "rate_a": 100}]
        errors = validate_inputs(buses, lines, base_mva=100, max_iter=20, tol=0.001)
        self.assertTrue(any("To-Bus 99 does not exist" in err for err in errors))

    # -------------------------------------------------------------
    # 2. Y-BUS ADMITTANCE MATRIX VERIFICATION
    # -------------------------------------------------------------
    def test_ybus_calculation(self):
        # 2-bus test system
        lines = [{"from": 1, "to": 2, "r": 0.0, "x": 0.1, "b": 0.0, "rate_a": 100}]
        Ybus, _ = build_ybus(2, lines)
        
        # Expected: Y11 = -10j, Y12 = 10j, Y21 = 10j, Y22 = -10j
        expected_Y11 = complex(0, -10.0)
        expected_Y12 = complex(0, 10.0)
        
        np.testing.assert_almost_equal(Ybus[0, 0], expected_Y11, decimal=4)
        np.testing.assert_almost_equal(Ybus[0, 1], expected_Y12, decimal=4)

    # -------------------------------------------------------------
    # 3. POWER FLOW SOLVER INTEGRITY (3-BUS SYSTEM)
    # -------------------------------------------------------------
    def test_solver_convergence(self):
        buses = [
            {"id": 1, "type": "slack", "v": 1.06, "ang": 0.0, "pg": 0, "qg": 0, "pl": 0, "ql": 0},
            {"id": 2, "type": "pv", "v": 1.045, "ang": 0.0, "pg": 40, "qg": 0, "pl": 20, "ql": 10},
            {"id": 3, "type": "pq", "v": 1.0, "ang": 0.0, "pg": 0, "qg": 0, "pl": 45, "ql": 15}
        ]
        lines = [
            {"from": 1, "to": 2, "r": 0.02, "x": 0.06, "b": 0.03, "rate_a": 100},
            {"from": 1, "to": 3, "r": 0.08, "x": 0.24, "b": 0.025, "rate_a": 100},
            {"from": 2, "to": 3, "r": 0.06, "x": 0.18, "b": 0.02, "rate_a": 100}
        ]
        
        V, delta, converged, history, exec_time, _ = execute_solver_core(
            "newton_raphson", buses, lines, base_mva=100, max_iter=20, tol=0.001
        )
        
        self.assertTrue(converged)
        self.assertLess(len(history), 10)
        self.assertAlmostEqual(V[0], 1.06, places=2)
        self.assertAlmostEqual(V[1], 1.045, places=2)

    # -------------------------------------------------------------
    # 4. IEEE 14-BUS REFERENCE BENCHMARK
    # -------------------------------------------------------------
    def test_ieee_14_bus_system(self):
        # standard IEEE 14-bus test case setup
        buses = [
            {"id": 1, "type": "slack", "v": 1.06, "ang": 0, "pg": 0, "qg": 0, "pl": 0, "ql": 0},
            {"id": 2, "type": "pv", "v": 1.045, "ang": 0, "pg": 40, "qg": 0, "pl": 21.7, "ql": 12.7},
            {"id": 3, "type": "pv", "v": 1.01, "ang": 0, "pg": 0, "qg": 0, "pl": 94.2, "ql": 19.0},
            {"id": 4, "type": "pq", "v": 1.0, "ang": 0, "pg": 0, "qg": 0, "pl": 47.8, "ql": -3.9},
            {"id": 5, "type": "pq", "v": 1.0, "ang": 0, "pg": 0, "qg": 0, "pl": 7.6, "ql": 1.6},
            {"id": 6, "type": "pv", "v": 1.07, "ang": 0, "pg": 0, "qg": 0, "pl": 11.2, "ql": 7.5},
            {"id": 7, "type": "pq", "v": 1.0, "ang": 0, "pg": 0, "qg": 0, "pl": 0, "ql": 0},
            {"id": 8, "type": "pv", "v": 1.09, "ang": 0, "pg": 0, "qg": 0, "pl": 0, "ql": 0},
            {"id": 9, "type": "pq", "v": 1.0, "ang": 0, "pg": 0, "qg": 0, "pl": 29.5, "ql": 16.6},
            {"id": 10, "type": "pq", "v": 1.0, "ang": 0, "pg": 0, "qg": 0, "pl": 9.0, "ql": 5.8},
            {"id": 11, "type": "pq", "v": 1.0, "ang": 0, "pg": 0, "qg": 0, "pl": 3.5, "ql": 1.8},
            {"id": 12, "type": "pq", "v": 1.0, "ang": 0, "pg": 0, "qg": 0, "pl": 6.1, "ql": 1.6},
            {"id": 13, "type": "pq", "v": 1.0, "ang": 0, "pg": 0, "qg": 0, "pl": 13.5, "ql": 5.8},
            {"id": 14, "type": "pq", "v": 1.0, "ang": 0, "pg": 0, "qg": 0, "pl": 14.9, "ql": 5.0}
        ]
        
        lines = [
            {"from": 1, "to": 2, "r": 0.01938, "x": 0.05917, "b": 0.0264, "rate_a": 100},
            {"from": 1, "to": 5, "r": 0.05403, "x": 0.22304, "b": 0.0246, "rate_a": 100},
            {"from": 2, "to": 3, "r": 0.04699, "x": 0.19797, "b": 0.0219, "rate_a": 100},
            {"from": 2, "to": 4, "r": 0.05811, "x": 0.17632, "b": 0.0187, "rate_a": 100},
            {"from": 2, "to": 5, "r": 0.05695, "x": 0.17388, "b": 0.0170, "rate_a": 100},
            {"from": 3, "to": 4, "r": 0.06701, "x": 0.17103, "b": 0.0173, "rate_a": 100},
            {"from": 4, "to": 5, "r": 0.01335, "x": 0.04211, "b": 0.0064, "rate_a": 100},
            {"from": 4, "to": 7, "r": 0.0, "x": 0.20912, "b": 0.0, "rate_a": 100},
            {"from": 4, "to": 9, "r": 0.0, "x": 0.55618, "b": 0.0, "rate_a": 100},
            {"from": 5, "to": 6, "r": 0.0, "x": 0.25202, "b": 0.0, "rate_a": 100},
            {"from": 6, "to": 11, "r": 0.09498, "x": 0.19890, "b": 0.0, "rate_a": 100},
            {"from": 6, "to": 12, "r": 0.12291, "x": 0.25581, "b": 0.0, "rate_a": 100},
            {"from": 6, "to": 13, "r": 0.06615, "x": 0.13027, "b": 0.0, "rate_a": 100},
            {"from": 7, "to": 8, "r": 0.0, "x": 0.17615, "b": 0.0, "rate_a": 100},
            {"from": 7, "to": 9, "r": 0.0, "x": 0.11001, "b": 0.0, "rate_a": 100},
            {"from": 9, "to": 10, "r": 0.03181, "x": 0.08450, "b": 0.0, "rate_a": 100},
            {"from": 9, "to": 14, "r": 0.12711, "x": 0.27038, "b": 0.0, "rate_a": 100},
            {"from": 10, "to": 11, "r": 0.08205, "x": 0.19207, "b": 0.0, "rate_a": 100},
            {"from": 12, "to": 13, "r": 0.22092, "x": 0.19988, "b": 0.0, "rate_a": 100},
            {"from": 13, "to": 14, "r": 0.17093, "x": 0.34802, "b": 0.0, "rate_a": 100}
        ]

        V, delta, converged, history, exec_time, _ = execute_solver_core(
            "newton_raphson", buses, lines, base_mva=100, max_iter=20, tol=0.001
        )

        self.assertTrue(converged, "IEEE 14-bus failed to converge.")
        self.assertLess(len(history), 6, "IEEE 14-bus took too many iterations.")
        # Verify Bus 14 solved voltage magnitude stays within realistic physical bounds (0.9 to 1.1 p.u.)
        self.assertTrue(0.9 <= V[13] <= 1.1)

if __name__ == '__main__':
    unittest.main()