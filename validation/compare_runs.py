#!/usr/bin/env python3
"""Compare commensurate runs while keeping approximation channels separate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import tomllib
from typing import Any

import numpy as np
import pandas as pd


CHANNELS = ("bond_dimension", "bp_environment", "trotter", "backend", "finite_size")


def _load_metadata(path: Path) -> dict[str, Any]:
    with (path / "metadata.toml").open("rb") as handle:
        return tomllib.load(handle)


def _flatten_configuration(metadata: dict[str, Any]) -> dict[str, Any]:
    flat: dict[str, Any] = {}
    for section in ("simulation", "peps", "bp", "observables", "conventions"):
        for key, value in metadata.get(section, {}).items():
            flat[f"{section}.{key}"] = value
    return flat


def metadata_mismatches(reference: dict[str, Any], candidate: dict[str, Any], channel: str) -> list[str]:
    """Return configuration differences that are not permitted for ``channel``."""
    if channel not in CHANNELS:
        raise ValueError(f"unknown comparison channel {channel}")
    allowed = {
        "bond_dimension": {"peps.maxdim"},
        "bp_environment": {"bp.tolerance", "bp.maxiter"},
        "trotter": {"simulation.dt", "simulation.steps", "simulation.trotter",
                    "observables.measure_every"},
        # PEPS/BP controls do not apply to exact/PP/QPU backends and are therefore
        # permitted to differ in a backend comparison.
        "backend": {"simulation.backend", "peps.maxdim", "peps.cutoff",
                    "peps.normalize_tensors", "bp.tolerance", "bp.maxiter"},
        "finite_size": {"simulation.nx", "simulation.ny", "simulation.defect_x",
                        "simulation.defect_y", "conventions.internal_bonds",
                        "conventions.measured_links", "conventions.raw_pauli_observables",
                        "conventions.boundary_distance_sites"},
    }[channel]
    ref = _flatten_configuration(reference)
    can = _flatten_configuration(candidate)
    mismatches = []
    for key in sorted(set(ref) | set(can)):
        if key in allowed:
            continue
        if ref.get(key) != can.get(key):
            mismatches.append(f"{key}: {ref.get(key)!r} != {can.get(key)!r}")

    ref_major = str(reference.get("schema", {}).get("version", "0")).split(".")[0]
    can_major = str(candidate.get("schema", {}).get("version", "0")).split(".")[0]
    if ref_major != can_major:
        mismatches.append(f"schema major: {ref_major} != {can_major}")
    if channel == "trotter":
        ref_sim, can_sim = reference["simulation"], candidate["simulation"]
        ref_tmax = float(ref_sim["dt"]) * int(ref_sim["steps"])
        can_tmax = float(can_sim["dt"]) * int(can_sim["steps"])
        if not np.isclose(ref_tmax, can_tmax, rtol=1e-12, atol=1e-12):
            mismatches.append(f"simulation.tmax: {ref_tmax!r} != {can_tmax!r}")
    return mismatches


def _time_key(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["_time_key"] = np.round(out["time"].to_numpy(float), 12)
    return out


def _metrics(values: pd.Series | np.ndarray, prefix: str) -> dict[str, float | None]:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if not len(array):
        return {f"{prefix}_max_abs": None, f"{prefix}_rms": None}
    return {
        f"{prefix}_max_abs": float(np.max(np.abs(array))),
        f"{prefix}_rms": float(np.sqrt(np.mean(array**2))),
    }


def _optional_long_table(reference: Path, candidate: Path, filename: str,
                         keys: list[str], value: str, prefix: str) -> dict[str, float | None]:
    ref = _time_key(pd.read_csv(reference / filename))
    can = _time_key(pd.read_csv(candidate / filename))
    merged = ref.merge(can, on=["_time_key", *keys], suffixes=("_ref", "_candidate"),
                       validate="one_to_one")
    return _metrics(merged[f"{value}_candidate"] - merged[f"{value}_ref"], prefix)


def _finite_size_compare(reference: Path, candidate: Path, reference_metadata: dict[str, Any],
                         candidate_metadata: dict[str, Any]) -> dict[str, float | int | None]:
    ref_g = _time_key(pd.read_csv(reference / "global.csv"))
    can_g = _time_key(pd.read_csv(candidate / "global.csv"))
    merged_g = ref_g.merge(can_g, on="_time_key", suffixes=("_ref", "_candidate"),
                           validate="one_to_one").sort_values("_time_key")
    if len(merged_g) != len(ref_g) or len(merged_g) != len(can_g):
        raise ValueError("finite-size comparison requires identical measured physical times")
    delta_error = merged_g.delta_phi_candidate - merged_g.delta_phi_ref
    report: dict[str, float | int | None] = {
        "matched_time_points": int(len(merged_g)),
        **_metrics(delta_error, "delta_phi"),
        "delta_phi_time_integrated_abs": float(np.trapezoid(
            np.abs(delta_error.to_numpy(float)), merged_g._time_key.to_numpy(float))),
    }
    for column in ("front_rms", "front_quantile", "delta_loop_survival"):
        report.update(_metrics(merged_g[f"{column}_candidate"] - merged_g[f"{column}_ref"], column))

    def relative_links(path: Path, metadata: dict[str, Any]) -> pd.DataFrame:
        frame = _time_key(pd.read_csv(path / "links.csv"))
        defect = metadata["simulation"]
        frame["rel_mid_x"] = np.round(frame.mid_x - float(defect["defect_x"]), 12)
        frame["rel_mid_y"] = np.round(frame.mid_y - float(defect["defect_y"]), 12)
        return frame

    ref_l = relative_links(reference, reference_metadata)
    can_l = relative_links(candidate, candidate_metadata)
    link_keys = ["_time_key", "orientation", "rel_mid_x", "rel_mid_y"]
    merged_l = ref_l.merge(can_l, on=link_keys, suffixes=("_ref", "_candidate"), validate="one_to_one")
    report.update(_metrics(merged_l.delta_n_candidate - merged_l.delta_n_ref, "common_local_delta_n"))
    common_per_time = merged_l.groupby("_time_key").size()
    smaller_per_time = min(ref_l.groupby("_time_key").size().min(), can_l.groupby("_time_key").size().min())
    report["common_links_per_time"] = int(common_per_time.min())
    report["common_link_fraction_of_smaller"] = float(common_per_time.min() / smaller_per_time)

    def relative_sites(path: Path, metadata: dict[str, Any]) -> pd.DataFrame:
        frame = _time_key(pd.read_csv(path / "sites.csv"))
        defect = metadata["simulation"]
        frame["rel_x"] = frame.x - int(defect["defect_x"])
        frame["rel_y"] = frame.y - int(defect["defect_y"])
        return frame

    ref_s = relative_sites(reference, reference_metadata)
    can_s = relative_sites(candidate, candidate_metadata)
    merged_s = ref_s.merge(can_s, on=["_time_key", "rel_x", "rel_y"],
                           suffixes=("_ref", "_candidate"), validate="one_to_one")
    report.update(_metrics(merged_s.delta_X_candidate - merged_s.delta_X_ref, "common_local_delta_X"))
    report.update(_metrics(merged_s.delta_local_flux_candidate - merged_s.delta_local_flux_ref,
                           "common_site_centered_delta_flux"))
    report["common_sites_per_time"] = int(merged_s.groupby("_time_key").size().min())

    def relative_wilson(path: Path, metadata: dict[str, Any]) -> pd.DataFrame:
        frame = _time_key(pd.read_csv(path / "wilson.csv"))
        defect = metadata["simulation"]
        frame["rel_x0"] = frame.x0 - int(defect["defect_x"])
        frame["rel_y0"] = frame.y0 - int(defect["defect_y"])
        return frame

    ref_w = relative_wilson(reference, reference_metadata)
    can_w = relative_wilson(candidate, candidate_metadata)
    merged_w = ref_w.merge(can_w, on=["_time_key", "width", "height", "rel_x0", "rel_y0"],
                           suffixes=("_ref", "_candidate"), validate="one_to_one")
    report.update(_metrics(merged_w.delta_value_candidate - merged_w.delta_value_ref,
                           "common_wilson_delta"))

    principal_kinds = {"electric_string", "magnetic_connected"}
    ref_c = _time_key(pd.read_csv(reference / "correlations.csv"))
    can_c = _time_key(pd.read_csv(candidate / "correlations.csv"))
    ref_c = ref_c[ref_c.kind.isin(principal_kinds)]
    can_c = can_c[can_c.kind.isin(principal_kinds)]
    merged_c = ref_c.merge(can_c, on=["_time_key", "kind", "observable_id"],
                           suffixes=("_ref", "_candidate"), validate="one_to_one")
    report.update(_metrics(merged_c.delta_value_candidate - merged_c.delta_value_ref,
                           "common_string_magnetic_delta"))
    return report


def compare(reference: Path, candidate: Path, channel: str | None = None,
            reference_metadata: dict[str, Any] | None = None,
            candidate_metadata: dict[str, Any] | None = None) -> dict[str, float | int | None]:
    if channel == "finite_size":
        reference_metadata = reference_metadata or _load_metadata(reference)
        candidate_metadata = candidate_metadata or _load_metadata(candidate)
        return _finite_size_compare(reference, candidate, reference_metadata, candidate_metadata)
    ref_g = _time_key(pd.read_csv(reference / "global.csv"))
    can_g = _time_key(pd.read_csv(candidate / "global.csv"))
    merged_g = ref_g.merge(can_g, on="_time_key", suffixes=("_ref", "_candidate"),
                           validate="one_to_one").sort_values("_time_key")
    if merged_g.empty:
        raise ValueError("the runs have no common measured physical times")
    if channel == "trotter":
        if len(merged_g) != min(len(ref_g), len(can_g)):
            raise ValueError("the coarser Trotter measurement grid is not a subset of the finer grid")
    elif channel is not None and (len(merged_g) != len(ref_g) or len(merged_g) != len(can_g)):
        raise ValueError(f"{channel} comparison requires identical measured physical times")
    delta_error = merged_g.delta_phi_candidate - merged_g.delta_phi_ref
    report: dict[str, float | int | None] = {
        "matched_time_points": int(len(merged_g)),
        "reference_time_points": int(len(ref_g)),
        "candidate_time_points": int(len(can_g)),
        **_metrics(delta_error, "delta_phi"),
        "delta_phi_time_integrated_abs": float(np.trapezoid(
            np.abs(delta_error.to_numpy(float)), merged_g._time_key.to_numpy(float))),
    }
    for column in ("front_rms", "front_quantile", "delta_loop_survival"):
        if f"{column}_ref" in merged_g and f"{column}_candidate" in merged_g:
            report.update(_metrics(merged_g[f"{column}_candidate"] - merged_g[f"{column}_ref"], column))

    ref_l = _time_key(pd.read_csv(reference / "links.csv"))
    can_l = _time_key(pd.read_csv(candidate / "links.csv"))
    merged_l = ref_l.merge(can_l, on=["_time_key", "link_id"], suffixes=("_ref", "_candidate"),
                           validate="one_to_one")
    report.update(_metrics(merged_l.delta_n_candidate - merged_l.delta_n_ref, "local_delta_n"))

    ref_s = _time_key(pd.read_csv(reference / "sites.csv"))
    can_s = _time_key(pd.read_csv(candidate / "sites.csv"))
    merged_s = ref_s.merge(can_s, on=["_time_key", "x", "y"], suffixes=("_ref", "_candidate"),
                           validate="one_to_one")
    report.update(_metrics(merged_s.delta_X_candidate - merged_s.delta_X_ref, "local_delta_X"))
    report.update(_metrics(merged_s.delta_local_flux_candidate - merged_s.delta_local_flux_ref,
                           "site_centered_delta_flux"))

    report.update(_optional_long_table(reference, candidate, "wilson.csv", ["observable_id"],
                                       "delta_value", "wilson_delta"))
    report.update(_optional_long_table(reference, candidate, "correlations.csv",
                                       ["kind", "observable_id"], "delta_value", "correlation_delta"))
    report.update(_optional_long_table(reference, candidate, "pauli.csv", ["observable_id"],
                                       "delta_value", "raw_pauli_delta"))
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--channel", required=True, choices=CHANNELS,
                        help="The one approximation channel allowed to differ; never a generic error bar.")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    ref_metadata = _load_metadata(args.reference)
    can_metadata = _load_metadata(args.candidate)
    mismatches = metadata_mismatches(ref_metadata, can_metadata, args.channel)
    if mismatches:
        details = "\n  - ".join(mismatches)
        raise SystemExit(f"runs are not commensurate for {args.channel}:\n  - {details}")

    report = {"comparison_channel": args.channel, **compare(
        args.reference, args.candidate, args.channel, ref_metadata, can_metadata)}
    text = json.dumps(report, indent=2, sort_keys=True, allow_nan=False)
    print(text)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
