#!/usr/bin/env python3
"""Reproduce the independent 3x3 dt/order scan quoted in the technical report."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

try:
    from .exact_reference import (
        ExactConfig, continuous_evolve, evolve, initial_states, phase_aligned_error,
        total_flux, with_dt,
    )
except ImportError:  # Direct execution: python validation/trotter_scan.py
    from exact_reference import (
        ExactConfig, continuous_evolve, evolve, initial_states, phase_aligned_error,
        total_flux, with_dt,
    )


def scan() -> pd.DataFrame:
    rows = []
    for K in (0.5, 1.5, 3.0, 3.04433, 4.0, 6.0):
        base = ExactConfig(nx=3, ny=3, K=K, Gamma=1.0, defect=(2, 2))
        scv0, loop0 = initial_states(base)
        scv_continuous = continuous_evolve(scv0, base, 1.0)
        loop_continuous = continuous_evolve(loop0, base, 1.0)
        delta_continuous = total_flux(loop_continuous, base) - total_flux(scv_continuous, base)
        for order in ("lie_x_zz", "strang"):
            for dt in (0.1, 0.05, 0.025):
                cfg = with_dt(base, dt, order)
                steps = round(1.0 / dt)
                scv = evolve(scv0, cfg, steps)
                loop = evolve(loop0, cfg, steps)
                delta = total_flux(loop, cfg) - total_flux(scv, cfg)
                rows.append({
                    "nx": 3,
                    "ny": 3,
                    "K": K,
                    "Gamma": 1.0,
                    "order": order,
                    "dt": dt,
                    "steps": steps,
                    "time": 1.0,
                    "delta_phi_continuous": delta_continuous,
                    "delta_phi_trotter": delta,
                    "delta_phi_abs_error": abs(delta - delta_continuous),
                    "loop_state_phase_aligned_l2_error": phase_aligned_error(loop_continuous, loop),
                })
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("trotter_scan_3x3.csv"))
    args = parser.parse_args()
    frame = scan()
    frame.to_csv(args.output, index=False)
    print(args.output)


if __name__ == "__main__":
    main()
