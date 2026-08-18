#!/usr/bin/env python3
"""Materialize the documented staged matrix as one validated TOML per run."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Run:
    nx: int
    ny: int
    K: float
    dt: float
    order: str = "lie_x_zz"
    backend: str = "bp"
    maxdim: int = 24
    tolerance: float = 1e-8
    maxiter: int = 50
    tmax: float = 1.0
    boundary: str = "open_dual"


def runs_for(stage: str, selected_dimension: int) -> list[Run]:
    critical = 3.04433
    if stage == "stage0_exact":
        return [Run(3, 3, K, dt, order, "exact", 16, 1e-10, 100)
                for K in (0.5, 1.5, 3.0, critical, 4.0, 6.0)
                for order in ("lie_x_zz", "strang") for dt in (0.1, 0.05, 0.025)]
    if stage == "stage0_bp_exact":
        return [Run(4, 4, K, 0.05, "lie_x_zz", backend, D, 1e-10, 100, 0.5)
                for K in (0.5, critical, 6.0) for backend in ("exact", "bp")
                for D in ((16,) if backend == "exact" else (4, 8, 16))]
    if stage == "stage1_coarse":
        return [Run(8, 6, K, 0.05, maxdim=16)
                for K in (0.5, 2.0, 3.0, critical, 4.0, 6.0)]
    if stage == "stage1_D":
        return [Run(8, 6, K, 0.05, maxdim=D)
                for K in (0.5, critical, 4.0, 6.0) for D in (8, 12, 16, 24, 32)]
    if stage == "stage1_BP":
        return [Run(8, 6, K, 0.05, maxdim=24, tolerance=tol, maxiter=maxiter)
                for K in (0.5, critical, 4.0)
                for tol, maxiter in ((1e-6, 25), (1e-8, 50), (1e-10, 100))]
    if stage == "stage1_trotter":
        return [Run(8, 6, K, dt, order, maxdim=24)
                for K in (0.5, critical, 4.0, 6.0)
                for order, dt in (("lie_x_zz", 0.1), ("lie_x_zz", 0.05),
                                  ("lie_x_zz", 0.025), ("strang", 0.05), ("strang", 0.025))]
    if stage == "stage2_size":
        return [Run(nx, ny, K, 0.05, maxdim=selected_dimension)
                for K in (critical, 4.0) for nx, ny in ((6, 6), (8, 6), (10, 8), (12, 10))]
    if stage == "stage3":
        main = [Run(12, 10, K, 0.05, maxdim=selected_dimension)
                for K in (3.0, critical, 4.0, 0.5, 2.0, 6.0)]
        reference = [Run(12, 10, K, 0.025, "strang", maxdim=selected_dimension)
                     for K in (critical, 4.0)]
        return main + reference
    raise ValueError(stage)


def slug(value: float) -> str:
    return f"{value:.8g}".replace(".", "p").replace("-", "m")


def render(run: Run, stage: str) -> str:
    steps = round(run.tmax / run.dt)
    if abs(steps * run.dt - run.tmax) > 1e-12:
        raise ValueError(f"tmax={run.tmax} is not commensurate with dt={run.dt}")
    defect_x, defect_y = (run.nx + 1) // 2, (run.ny + 1) // 2
    label = (f"{stage}_{run.nx}x{run.ny}_K{slug(run.K)}_{run.order}_dt{slug(run.dt)}_"
             f"D{run.maxdim}_tol{slug(run.tolerance)}")
    return f'''[simulation]
nx = {run.nx}
ny = {run.ny}
K = {run.K}
Gamma = 1.0
dt = {run.dt}
steps = {steps}
trotter = "{run.order}"
boundary = "{run.boundary}"
backend = "{run.backend}"
defect_x = {defect_x}
defect_y = {defect_y}

[peps]
maxdim = {run.maxdim}
cutoff = 1.0e-10
normalize_tensors = true

[bp]
tolerance = {run.tolerance}
maxiter = {run.maxiter}
fail_on_nonconvergence = true

[observables]
measure_every = 1
radial_bin_width = 0.5
angular_bins = 16
front_quantile = 0.9
correlation_max_distance = 4
wilson_max_area = 4

[output]
root = "results"
label = "{label}"
run_id = ""
overwrite = false
'''


def main() -> None:
    stages = ("stage0_exact", "stage0_bp_exact", "stage1_coarse", "stage1_D",
              "stage1_BP", "stage1_trotter", "stage2_size", "stage3")
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=stages, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--selected-dimension", type=int, default=24,
                        help="Accepted D cap from Stage 1; used only for Stages 2/3.")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    runs = runs_for(args.stage, args.selected_dimension)
    paths = []
    for index, run in enumerate(runs, start=1):
        filename = (f"{index:03d}_{run.backend}_{run.nx}x{run.ny}_K{slug(run.K)}_"
                    f"{run.order}_dt{slug(run.dt)}_D{run.maxdim}_tol{slug(run.tolerance)}.toml")
        path = args.output / filename
        if path.exists() and not args.overwrite:
            raise SystemExit(f"refusing to overwrite {path}; pass --overwrite deliberately")
        path.write_text(render(run, args.stage), encoding="utf-8")
        paths.append(path)
    print(f"materialized {len(paths)} configurations in {args.output}")


if __name__ == "__main__":
    main()
