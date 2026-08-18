#!/usr/bin/env python3
"""Fail-fast validation of one standardized BP/exact/PP/QPU run directory."""

from __future__ import annotations

import argparse
import sys
import tomllib
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_FILES = (
    "metadata.toml",
    "global.csv",
    "links.csv",
    "sites.csv",
    "profiles.csv",
    "wilson.csv",
    "correlations.csv",
    "diagnostics.csv",
    "bp_updates.csv",
    "pauli.csv",
)

REQUIRED_COLUMNS = {
    "global.csv": {"run_id", "backend", "step", "time", "phi_scv", "phi_loop", "delta_phi",
                   "loop_survival_scv", "loop_survival_loop", "front_rms", "front_quantile"},
    "links.csv": {"run_id", "backend", "step", "time", "link_id", "n_scv", "n_loop", "delta_n"},
    "sites.csv": {"run_id", "backend", "step", "time", "x", "y", "delta_X",
                  "local_flux_scv", "local_flux_loop", "delta_local_flux"},
    "profiles.csv": {"run_id", "backend", "step", "time", "profile_kind", "count", "delta_n_mean"},
    "wilson.csv": {"run_id", "backend", "step", "time", "observable_id", "delta_value"},
    "correlations.csv": {"run_id", "backend", "step", "time", "kind", "observable_id", "delta_value"},
    "diagnostics.csv": {"run_id", "backend", "step", "time", "state", "truncation_max", "truncation_sum",
                        "bp_update_count", "bp_final_residual_max", "bp_all_converged",
                        "physical_link_violations"},
    "bp_updates.csv": {"run_id", "backend", "step", "time", "state", "update_id", "iterations",
                       "final_residual", "converged"},
    "pauli.csv": {"run_id", "backend", "step", "time", "observable_id", "basis", "qpu_setting",
                  "pauli", "support", "value_scv", "value_loop", "delta_value"},
}


def validate_run(path: Path, strict_bp: bool = True) -> list[str]:
    errors: list[str] = []
    for filename in REQUIRED_FILES:
        if not (path / filename).is_file():
            errors.append(f"missing {filename}")
    if errors:
        return errors

    with (path / "metadata.toml").open("rb") as handle:
        metadata = tomllib.load(handle)
    tables = {filename: pd.read_csv(path / filename) for filename in REQUIRED_COLUMNS}
    for filename, required in REQUIRED_COLUMNS.items():
        missing = required - set(tables[filename].columns)
        if missing:
            errors.append(f"{filename} missing columns: {sorted(missing)}")
    if errors:
        return errors

    global_df = tables["global.csv"]
    links_df = tables["links.csv"]
    sites_df = tables["sites.csv"]
    profiles = tables["profiles.csv"]
    diagnostics = tables["diagnostics.csv"]
    bp_updates = tables["bp_updates.csv"]
    pauli = tables["pauli.csv"]

    run_id = metadata["schema"]["run_id"]
    backend = metadata["simulation"]["backend"]
    dt = float(metadata["simulation"]["dt"])
    for filename, frame in tables.items():
        if len(frame) and set(frame.run_id.astype(str)) != {run_id}:
            errors.append(f"{filename} contains a run_id different from metadata")
        if len(frame) and set(frame.backend.astype(str)) != {backend}:
            errors.append(f"{filename} contains a backend different from metadata")
        if len(frame) and not np.allclose(frame.time, frame.step * dt, atol=1e-12, rtol=1e-12):
            errors.append(f"{filename} has time != step*dt")

    expected = float(metadata["conventions"]["expected_delta_phi_t0"])
    at_zero = global_df.loc[global_df["step"] == 0]
    if len(at_zero) != 1:
        errors.append(f"global.csv must contain one t=0 row, found {len(at_zero)}")
    else:
        row = at_zero.iloc[0]
        if not np.isclose(row.phi_scv, 0.0, atol=1e-8):
            errors.append(f"phi_scv(0)={row.phi_scv}, expected 0")
        if not np.isclose(row.delta_phi, expected, atol=1e-8):
            errors.append(f"delta_phi(0)={row.delta_phi}, expected {expected}")

    if global_df.duplicated(["step"]).any():
        errors.append("global.csv contains duplicate steps")
    if links_df.duplicated(["step", "link_id"]).any():
        errors.append("links.csv contains duplicate step/link_id keys")
    if sites_df.duplicated(["step", "x", "y"]).any():
        errors.append("sites.csv contains duplicate step/site keys")
    if pauli.duplicated(["step", "observable_id"]).any():
        errors.append("pauli.csv contains duplicate step/observable_id keys")

    bad_links = links_df[(links_df.n_scv < -1e-8) | (links_df.n_scv > 1 + 1e-8) |
                         (links_df.n_loop < -1e-8) | (links_df.n_loop > 1 + 1e-8)]
    if len(bad_links):
        errors.append(f"{len(bad_links)} link expectations lie outside [0,1]")

    grouped = links_df.groupby("step", sort=True)[["n_scv", "n_loop"]].sum()
    merged = global_df.set_index("step").join(grouped, how="left")
    if not np.allclose(merged.phi_scv, merged.n_scv, atol=1e-8, rtol=1e-8):
        errors.append("phi_scv does not equal the link sum")
    if not np.allclose(merged.phi_loop, merged.n_loop, atol=1e-8, rtol=1e-8):
        errors.append("phi_loop does not equal the link sum")

    if not np.allclose(global_df.delta_phi, global_df.phi_loop - global_df.phi_scv,
                       atol=1e-10, rtol=1e-10):
        errors.append("global delta_phi is not phi_loop-phi_scv")
    if not np.allclose(links_df.delta_n, links_df.n_loop - links_df.n_scv,
                       atol=1e-10, rtol=1e-10):
        errors.append("link delta_n is not n_loop-n_scv")

    expected_links = int(metadata["conventions"]["measured_links"])
    link_counts = links_df.groupby("step").size()
    if not (link_counts == expected_links).all():
        errors.append("one or more measured steps does not contain every configured link")
    expected_sites = int(metadata["simulation"]["nx"]) * int(metadata["simulation"]["ny"])
    site_counts = sites_df.groupby("step").size()
    if not (site_counts == expected_sites).all():
        errors.append("one or more measured steps does not contain every dual site")
    expected_paulis = int(metadata["conventions"]["raw_pauli_observables"])
    pauli_counts = pauli.groupby("step").size()
    if not (pauli_counts == expected_paulis).all():
        errors.append("one or more measured steps does not contain every raw Pauli observable")
    if set(pauli.qpu_setting) - {"global_X", "global_Z"}:
        errors.append("pauli.csv contains a non-principal QPU measurement setting")
    wrong_setting = pauli[((pauli.basis == "X") & (pauli.qpu_setting != "global_X")) |
                          ((pauli.basis == "Z") & (pauli.qpu_setting != "global_Z"))]
    if len(wrong_setting):
        errors.append("pauli.csv basis and QPU setting disagree")
    manifest_variation = pauli.groupby("observable_id")[["basis", "qpu_setting", "pauli", "support"]].nunique()
    if (manifest_variation.to_numpy() > 1).any():
        errors.append("the raw Pauli manifest changes between measured times")
    bad_paulis = pauli[(pauli.value_scv < -1 - 1e-8) | (pauli.value_scv > 1 + 1e-8) |
                       (pauli.value_loop < -1 - 1e-8) | (pauli.value_loop > 1 + 1e-8)]
    if len(bad_paulis):
        errors.append(f"{len(bad_paulis)} raw Pauli expectations lie outside [-1,1]")
    if not np.allclose(pauli.delta_value, pauli.value_loop - pauli.value_scv,
                       atol=1e-10, rtol=1e-10):
        errors.append("Pauli delta_value is not value_loop-value_scv")
    for kind in ("radial", "angular"):
        counts = profiles[profiles.profile_kind == kind].groupby("step")["count"].sum()
        if len(counts) != len(global_df) or not (counts == expected_links).all():
            errors.append(f"{kind} profile bins do not partition all measured links")

    local_sums = sites_df.groupby("step")[["local_flux_scv", "local_flux_loop"]].sum()
    merged_local = global_df.set_index("step").join(local_sums, how="left")
    if not np.allclose(merged_local.phi_scv, merged_local.local_flux_scv, atol=1e-8, rtol=1e-8):
        errors.append("site-centered SCV flux does not sum to phi_scv")
    if not np.allclose(merged_local.phi_loop, merged_local.local_flux_loop, atol=1e-8, rtol=1e-8):
        errors.append("site-centered loop flux does not sum to phi_loop")

    if (diagnostics.truncation_max < -1e-15).any() or (diagnostics.truncation_sum < -1e-15).any():
        errors.append("negative truncation diagnostics recorded")
    if diagnostics.duplicated(["step", "state"]).any():
        errors.append("diagnostics.csv contains duplicate step/state keys")
    expected_diagnostic_rows = 2 * (int(metadata["simulation"]["steps"]) + 1)
    if len(diagnostics) != expected_diagnostic_rows:
        errors.append(f"diagnostics.csv has {len(diagnostics)} rows, expected {expected_diagnostic_rows}")

    if strict_bp and backend == "bp":
        unconverged = diagnostics[diagnostics.bp_all_converged.astype(str).str.lower() != "true"]
        if len(unconverged):
            errors.append(f"{len(unconverged)} state/step rows contain unconverged BP updates")
        tolerance = float(metadata["bp"]["tolerance"])
        finite = diagnostics.bp_final_residual_max.dropna()
        if len(finite) and (finite > tolerance * (1 + 1e-10)).any():
            errors.append("a converged BP summary residual exceeds the configured tolerance")
        recorded = bp_updates.groupby(["step", "state"]).size().rename("recorded")
        declared = diagnostics.set_index(["step", "state"])["bp_update_count"]
        joined = declared.to_frame().join(recorded, how="left").fillna({"recorded": 0})
        if not np.array_equal(joined.bp_update_count.to_numpy(int), joined.recorded.to_numpy(int)):
            errors.append("bp_updates.csv counts disagree with diagnostics.csv")
    if (diagnostics.physical_link_violations > 0).any():
        errors.append("diagnostics records physical link-bound violations")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_directory", type=Path)
    parser.add_argument("--allow-unconverged-bp", action="store_true")
    args = parser.parse_args()
    errors = validate_run(args.run_directory, strict_bp=not args.allow_unconverged_bp)
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(f"validated {args.run_directory}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
