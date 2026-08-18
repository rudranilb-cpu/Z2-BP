#!/usr/bin/env python3
"""Publication-oriented figures for standardized BP, exact, PP, and QPU runs."""

from __future__ import annotations

import argparse
from pathlib import Path
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import pandas as pd

try:
    from .common import RunData, channel_signature, fit_velocity, load_run, nearest_steps
except ImportError:  # Direct execution: python analysis/plot_results.py ...
    from common import RunData, channel_signature, fit_velocity, load_run, nearest_steps


def set_style() -> None:
    mpl.rcParams.update({
        "font.family": "serif",
        "font.size": 10,
        "axes.labelsize": 10,
        "axes.titlesize": 10,
        "legend.fontsize": 8,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def save(fig: plt.Figure, out: Path, stem: str) -> None:
    fig.savefig(out / f"{stem}.pdf")
    fig.savefig(out / f"{stem}.png")
    plt.close(fig)


def plot_delta_flux(runs: list[RunData], out: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    for run in runs:
        df = run.global_df
        errors = df.delta_phi_stderr.to_numpy(float)
        if np.isfinite(errors).any():
            ax.errorbar(df.time, df.delta_phi, yerr=errors, marker="o", ms=3, lw=1.2, label=run.label)
        else:
            ax.plot(df.time, df.delta_phi, marker="o", ms=3, lw=1.2, label=run.label)
    ax.axhline(0, color="0.5", lw=0.7)
    ax.set(xlabel=r"$t$", ylabel=r"$\Delta\Phi=\Phi_{\rm loop}-\Phi_{\rm SCV}$")
    ax.legend(frameon=False)
    save(fig, out, "delta_flux_vs_time")


def _link_segments(frame: pd.DataFrame) -> list[list[tuple[float, float]]]:
    segments = []
    for row in frame.itertuples():
        x2 = row.v2_x if np.isfinite(row.v2_x) else 2 * row.mid_x - row.v1_x
        y2 = row.v2_y if np.isfinite(row.v2_y) else 2 * row.mid_y - row.v1_y
        segments.append([(row.v1_x, row.v1_y), (x2, y2)])
    return segments


def plot_link_heatmaps(run: RunData, out: Path, requested_times: list[float]) -> None:
    steps = nearest_steps(run.global_df, requested_times)
    vmax = max(abs(run.links.loc[run.links.step.isin(steps), "delta_n"]).max(), 1e-12)
    fig, axes = plt.subplots(1, len(steps), figsize=(3.3 * len(steps), 3.0), squeeze=False)
    for ax, step in zip(axes[0], steps):
        frame = run.links[run.links.step == step]
        collection = LineCollection(_link_segments(frame), array=frame.delta_n.to_numpy(), cmap="RdBu_r",
                                    norm=mpl.colors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax), linewidths=6)
        ax.add_collection(collection)
        vertices = set(zip(frame.v1_x, frame.v1_y))
        vertices.update((x, y) for x, y in zip(frame.v2_x, frame.v2_y) if np.isfinite(x) and np.isfinite(y))
        vx, vy = zip(*sorted(vertices))
        ax.scatter(vx, vy, s=6, color="0.15", zorder=3)
        ax.autoscale()
        ax.set_aspect("equal")
        time = float(frame.time.iloc[0])
        ax.set(title=fr"$t={time:.3g}$", xlabel="dual x", ylabel="dual y")
    fig.colorbar(collection, ax=axes.ravel().tolist(), label=r"$\Delta n_b$", shrink=0.82)
    save(fig, out, "local_flux_heatmaps")


def _dimension_errors(run: RunData, reference: RunData) -> pd.DataFrame:
    frame = run.global_df[["time", "delta_phi", "front_rms"]].merge(
        reference.global_df[["time", "delta_phi", "front_rms"]], on="time",
        suffixes=("", "_ref"), validate="one_to_one")
    frame["delta_phi_error"] = (frame.delta_phi - frame.delta_phi_ref).abs()
    frame["front_rms_error"] = (frame.front_rms - frame.front_rms_ref).abs()

    links = run.links[["time", "link_id", "delta_n"]].merge(
        reference.links[["time", "link_id", "delta_n"]], on=["time", "link_id"],
        suffixes=("", "_ref"), validate="one_to_one")
    links["error"] = links.delta_n - links.delta_n_ref
    link_metrics = links.groupby("time", as_index=False).agg(
        local_link_max_error=("error", lambda x: np.max(np.abs(x))),
        local_link_rms_error=("error", lambda x: np.sqrt(np.mean(x**2))))
    frame = frame.merge(link_metrics, on="time", validate="one_to_one")

    for table_name, keys, prefix in (("wilson", ["observable_id"], "wilson"),
                                     ("correlations", ["kind", "observable_id"], "correlation")):
        left = getattr(run, table_name)
        right = getattr(reference, table_name)
        merged = left[["time", *keys, "delta_value"]].merge(
            right[["time", *keys, "delta_value"]], on=["time", *keys],
            suffixes=("", "_ref"), validate="one_to_one")
        merged["error"] = (merged.delta_value - merged.delta_value_ref).abs()
        metric = merged.groupby("time", as_index=False).agg(**{f"{prefix}_max_error": ("error", "max")})
        frame = frame.merge(metric, on="time", how="left", validate="one_to_one")
    frame[["wilson_max_error", "correlation_max_error"]] = frame[
        ["wilson_max_error", "correlation_max_error"]].fillna(0.0)
    return frame


def plot_radial_front(run: RunData, out: Path, requested_times: list[float], tmin: float | None,
                      tmax: float | None, front_column: str) -> None:
    steps = nearest_steps(run.global_df, requested_times)
    radial = run.profiles[(run.profiles.profile_kind == "radial") & run.profiles.step.isin(steps)]
    angular = run.profiles[(run.profiles.profile_kind == "angular") & run.profiles.step.isin(steps)]
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.2))
    for step, frame in radial.groupby("step", sort=True):
        axes[0].plot(frame.bin_center, frame.abs_delta_n_mean, marker="o", ms=3,
                     label=fr"$t={frame.time.iloc[0]:.3g}$")
    axes[0].set(xlabel="bond-midpoint radius", ylabel=r"radial mean $|\Delta n_b|$")
    axes[0].legend(frameon=False)

    for step, frame in angular.groupby("step", sort=True):
        axes[1].plot(frame.bin_center, frame.abs_delta_n_mean, marker="o", ms=3,
                     label=fr"$t={frame.time.iloc[0]:.3g}$")
    axes[1].set(xlabel=r"bond angle $\theta$", ylabel=r"angular mean $|\Delta n_b|$")
    axes[1].legend(frameon=False)

    df = run.global_df
    fits = []
    columns = [front_column, "front_quantile" if front_column == "front_rms" else "front_rms"]
    labels = {"front_rms": "RMS", "front_quantile": "weighted quantile"}
    for index, column in enumerate(columns):
        try:
            fit = fit_velocity(run.global_df, tmin=tmin, tmax=tmax, column=column)
        except ValueError:
            continue
        boundary_distance = float(run.metadata["conventions"]["boundary_distance_sites"])
        fit_window = df[(df.time >= fit["tmin"]) & (df.time <= fit["tmax"])]
        fit["boundary_distance_sites"] = boundary_distance
        fit["minimum_boundary_margin"] = boundary_distance - float(fit_window[column].max())
        if fit["minimum_boundary_margin"] <= 0:
            warnings.warn(f"{column} velocity window reaches the nearest open boundary")
        fits.append(fit)
        axes[2].plot(df.time, df[column], "o-", ms=3, label=f"{labels[column]} front")
        grid = np.linspace(fit["tmin"], fit["tmax"], 100)
        axes[2].plot(grid, fit["velocity"] * grid + fit["intercept"], "--",
                     label=fr"{labels[column]} fit $v={fit['velocity']:.3g}$")
        if index == 0:
            axes[2].axvspan(fit["tmin"], fit["tmax"], color="0.8", alpha=0.25)
    if not fits:
        raise ValueError("neither front definition has enough finite points for a velocity fit")
    axes[2].set(xlabel=r"$t$", ylabel="front radius")
    axes[2].legend(frameon=False)
    pd.DataFrame(fits).to_csv(out / "velocity_fit.csv", index=False)
    save(fig, out, "radial_front_and_velocity")


def plot_bond_convergence(runs: list[RunData], out: Path) -> None:
    bp_runs = [run for run in runs if run.metadata["simulation"]["backend"] == "bp"]
    groups: dict[tuple, list[RunData]] = {}
    for run in bp_runs:
        groups.setdefault(channel_signature(run, "bond_dimension"), []).append(run)
    groups = {key: group for key, group in groups.items()
              if len({run.metadata["peps"]["maxdim"] for run in group}) >= 2}
    if not groups:
        warnings.warn("bond-dimension convergence plot skipped: supply at least two BP dimensions")
        return
    summary = []
    for group_id, group in enumerate(groups.values(), start=1):
        by_dimension = {run.metadata["peps"]["maxdim"]: run for run in group}
        group = [by_dimension[dimension] for dimension in sorted(by_dimension)]
        reference = group[-1].global_df[["time", "delta_phi"]].rename(columns={"delta_phi": "reference"})
        reference_run = group[-1]
        fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.2))
        per_dimension = {}
        for run in group:
            dimension = run.metadata["peps"]["maxdim"]
            axes[0].plot(run.global_df.time, run.global_df.delta_phi, marker="o", ms=2.5, label=f"D={dimension}")
            merged = run.global_df.merge(reference, on="time")
            error = merged.delta_phi - merged.reference
            axes[1].plot(merged.time, np.abs(error), marker="o", ms=2.5, label=f"D={dimension}")
            metrics = _dimension_errors(run, reference_run)
            per_dimension[dimension] = metrics
            summary.append({
                "group_id": group_id, "maxdim": dimension,
                "max_abs_delta_phi_vs_largest_D": metrics.delta_phi_error.max(),
                "rms_delta_phi_vs_largest_D": np.sqrt(np.mean(metrics.delta_phi_error**2)),
                "max_local_delta_n_vs_largest_D": metrics.local_link_max_error.max(),
                "rms_local_delta_n_vs_largest_D": np.sqrt(np.mean(metrics.local_link_rms_error**2)),
                "max_front_rms_vs_largest_D": metrics.front_rms_error.max(),
                "max_wilson_vs_largest_D": metrics.wilson_max_error.max(),
                "max_correlation_vs_largest_D": metrics.correlation_max_error.max(),
            })
        required_rows = []
        for time in reference_run.global_df.time:
            chosen = None
            for dimension in sorted(per_dimension):
                row = per_dimension[dimension].loc[np.isclose(per_dimension[dimension].time, time)].iloc[0]
                passes = (row.delta_phi_error <= 0.02 and row.local_link_rms_error <= 0.005 and
                          row.local_link_max_error <= 0.02 and row.front_rms_error <= 0.05 and
                          row.wilson_max_error <= 0.02 and row.correlation_max_error <= 0.02)
                if passes:
                    chosen = (dimension, row)
                    break
            if chosen is None:
                dimension = max(per_dimension)
                row = per_dimension[dimension].loc[np.isclose(per_dimension[dimension].time, time)].iloc[0]
                chosen = (dimension, row)
            dimension, row = chosen
            required_rows.append({
                "group_id": group_id, "time": time, "minimum_tested_D_passing": dimension,
                "largest_D_reference": max(per_dimension), "delta_phi_error": row.delta_phi_error,
                "local_link_rms_error": row.local_link_rms_error,
                "local_link_max_error": row.local_link_max_error, "front_rms_error": row.front_rms_error,
                "wilson_max_error": row.wilson_max_error, "correlation_max_error": row.correlation_max_error,
            })
        required_frame = pd.DataFrame(required_rows)
        axes[2].step(required_frame.time, required_frame.minimum_tested_D_passing, where="mid")
        axes[0].set(xlabel=r"$t$", ylabel=r"$\Delta\Phi$")
        axes[1].set(xlabel=r"$t$", ylabel="absolute difference from largest D")
        if any((frame.delta_phi_error > 0).any() for frame in per_dimension.values()):
            axes[1].set_yscale("log")
        else:
            axes[1].text(0.5, 0.5, "identical in supplied data", ha="center", va="center",
                         transform=axes[1].transAxes)
        axes[2].set(xlabel=r"$t$", ylabel="minimum tested D passing screen")
        for ax in axes[:2]:
            ax.legend(frameon=False)
        if group_id == 1:
            all_required = required_frame
        else:
            all_required = pd.concat([all_required, required_frame], ignore_index=True)
        stem = "bond_dimension_convergence" if len(groups) == 1 else f"bond_dimension_convergence_g{group_id}"
        save(fig, out, stem)
    pd.DataFrame(summary).to_csv(out / "bond_dimension_convergence_metrics.csv", index=False)
    all_required.to_csv(out / "minimum_tested_D_vs_time.csv", index=False)


def plot_trotter_convergence(runs: list[RunData], out: Path) -> None:
    groups: dict[tuple, list[RunData]] = {}
    for run in runs:
        groups.setdefault(channel_signature(run, "trotter"), []).append(run)
    groups = {key: group for key, group in groups.items() if len({
        (run.metadata["simulation"]["dt"], run.metadata["simulation"]["trotter"])
        for run in group}) >= 2}
    if not groups:
        warnings.warn("dt/Trotter convergence plot skipped: supply multiple dt/order runs")
        return
    for group_id, group in enumerate(groups.values(), start=1):
        fig, ax = plt.subplots(figsize=(5.2, 3.4))
        for run in sorted(group, key=lambda r: (r.metadata["simulation"]["trotter"], r.metadata["simulation"]["dt"])):
            sim = run.metadata["simulation"]
            ax.plot(run.global_df.time, run.global_df.delta_phi, marker="o", ms=2.5,
                    label=f"{sim['trotter']}, dt={sim['dt']}")
        ax.set(xlabel=r"$t$", ylabel=r"$\Delta\Phi$")
        ax.legend(frameon=False)
        stem = "dt_trotter_convergence" if len(groups) == 1 else f"dt_trotter_convergence_g{group_id}"
        save(fig, out, stem)


def plot_loop_spreading(run: RunData, out: Path) -> None:
    df = run.global_df
    fig, axes = plt.subplots(2, 2, figsize=(8.0, 6.0), sharex=True)
    axes[0, 0].plot(df.time, df.loop_survival_loop, "o-", ms=3, label="loop initial state")
    axes[0, 0].plot(df.time, df.loop_survival_scv, "o-", ms=3, label="SCV")
    axes[0, 0].set(ylabel=r"$\langle S_{\rm loop}\rangle$")
    axes[0, 0].legend(frameon=False)
    axes[0, 1].plot(df.time, df.central_star_delta_flux, "o-", ms=3)
    axes[0, 1].set(ylabel="central-star differential flux")
    axes[1, 0].plot(df.time, df.front_participation, "o-", ms=3)
    axes[1, 0].set(xlabel=r"$t$", ylabel="effective bond participation")
    axes[1, 1].plot(df.time, df.front_anisotropy, "o-", ms=3)
    axes[1, 1].set(xlabel=r"$t$", ylabel="second-moment anisotropy")
    save(fig, out, "loop_survival_and_spreading")


def plot_runtime_diagnostics(runs: list[RunData], out: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.1))
    for run in runs:
        grouped = run.diagnostics.groupby("step", as_index=False).agg(
            time=("time", "first"), evolution_seconds=("evolution_seconds", "sum"),
            max_bond_dimension=("max_bond_dimension", "max"), live_memory_mb=("live_memory_mb", "max"))
        axes[0].plot(grouped.time, grouped.evolution_seconds, marker="o", ms=2.5, label=run.label)
        axes[1].plot(grouped.time, grouped.max_bond_dimension, marker="o", ms=2.5)
        axes[2].plot(grouped.time, grouped.live_memory_mb, marker="o", ms=2.5)
    axes[0].set(xlabel=r"$t$", ylabel="evolution seconds / step")
    axes[1].set(xlabel=r"$t$", ylabel="maximum bond dimension")
    axes[2].set(xlabel=r"$t$", ylabel="Julia live heap [MiB]")
    axes[0].legend(frameon=False)
    save(fig, out, "runtime_and_resource_diagnostics")


def plot_classical_scaling(runs: list[RunData], out: Path) -> None:
    rows = []
    for run in runs:
        sim = run.metadata["simulation"]
        evolved = run.diagnostics[run.diagnostics.step > 0]
        per_step = evolved.groupby("step").agg(
            evolution_seconds=("evolution_seconds", "sum"),
            measurement_seconds=("measurement_seconds", "sum"))
        rows.append({
            "label": run.label, "backend": sim["backend"], "nx": sim["nx"], "ny": sim["ny"],
            "nsites": sim["nx"] * sim["ny"], "maxdim_cap": run.metadata["peps"]["maxdim"],
            "median_evolution_seconds_per_pair_step": per_step.evolution_seconds.median(),
            "p90_evolution_seconds_per_pair_step": per_step.evolution_seconds.quantile(0.9),
            "median_measurement_seconds": per_step.measurement_seconds.median(),
            "peak_achieved_bond_dimension": run.diagnostics.max_bond_dimension.max(),
            "peak_live_memory_mb": run.diagnostics.live_memory_mb.max(),
            "peak_maxrss_mb": run.diagnostics.maxrss_mb.max(),
            "max_bp_iterations_per_state_step": run.diagnostics.bp_iterations.max(),
            "max_bp_final_residual": run.diagnostics.bp_final_residual_max.max(),
            "max_local_truncation": run.diagnostics.truncation_max.max(),
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(out / "classical_difficulty_summary.csv", index=False)
    if summary.nsites.nunique() < 2:
        return
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.2))
    for row in summary.itertuples():
        axes[0].scatter(row.nsites, row.median_evolution_seconds_per_pair_step,
                        label=f"{row.nx}x{row.ny}, D={row.maxdim_cap}")
        axes[1].scatter(row.nsites, row.peak_live_memory_mb)
    axes[0].set(xlabel="dual-spin sites", ylabel="median evolution seconds / pair step",
                xscale="log", yscale="log")
    axes[1].set(xlabel="dual-spin sites", ylabel="peak Julia live heap [MiB]", xscale="log", yscale="log")
    axes[0].legend(frameon=False)
    save(fig, out, "lattice_size_runtime_memory_scaling")


def plot_numerical_diagnostics(runs: list[RunData], out: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(8.0, 6.0), sharex=True)
    aggregated = []
    for run in runs:
        grouped = run.diagnostics.groupby("step", as_index=False).agg(
            time=("time", "first"), truncation_max=("truncation_max", "max"),
            bp_final_residual_max=("bp_final_residual_max", "max"), norm_drift=("norm_drift", "max"),
            physical_link_violations=("physical_link_violations", "sum"))
        aggregated.append(grouped)
        axes[0, 0].plot(grouped.time, grouped.truncation_max, marker="o", ms=2.5, label=run.label)
        axes[0, 1].plot(grouped.time, grouped.bp_final_residual_max, marker="o", ms=2.5)
        axes[1, 0].plot(grouped.time, grouped.norm_drift, marker="o", ms=2.5)
        axes[1, 1].plot(grouped.time, grouped.physical_link_violations, marker="o", ms=2.5)
    axes[0, 0].set(ylabel="max local discarded weight")
    if any((frame.truncation_max.abs() > 0).any() for frame in aggregated):
        axes[0, 0].set_yscale("symlog", linthresh=1e-16)
    axes[0, 1].set(ylabel="max final BP residual")
    residuals = pd.concat([frame.bp_final_residual_max for frame in aggregated], ignore_index=True)
    if np.isfinite(residuals).any() and (residuals.dropna() > 0).any():
        axes[0, 1].set_yscale("log")
    else:
        axes[0, 1].text(0.5, 0.5, "not applicable", ha="center", va="center", transform=axes[0, 1].transAxes)
    axes[1, 0].set(xlabel=r"$t$", ylabel="norm drift")
    if any((frame.norm_drift.abs() > 0).any() for frame in aggregated):
        axes[1, 0].set_yscale("symlog", linthresh=1e-16)
    axes[1, 1].set(xlabel=r"$t$", ylabel="link-bound violations")
    axes[0, 0].legend(frameon=False)
    save(fig, out, "numerical_error_diagnostics")


def plot_selected_observables(run: RunData, out: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.1))
    for observable_id, frame in run.wilson.groupby("observable_id"):
        axes[0].plot(frame.time, frame.delta_value, marker="o", ms=2.5, label=observable_id)
    axes[0].set(xlabel=r"$t$", ylabel=r"$\Delta W(C)$")
    axes[0].legend(frameon=False, ncol=2)

    strings = run.correlations[run.correlations.kind == "electric_string"]
    for observable_id, frame in strings.groupby("observable_id"):
        axes[1].plot(frame.time, frame.delta_value, marker="o", ms=2.5, label=observable_id)
    axes[1].set(xlabel=r"$t$", ylabel="differential electric string")
    axes[1].legend(frameon=False, ncol=2)

    magnetic = run.correlations[run.correlations.kind == "magnetic_connected"]
    for observable_id, frame in magnetic.groupby("observable_id"):
        axes[2].plot(frame.time, frame.delta_value, marker="o", ms=2.5, label=observable_id)
    axes[2].set(xlabel=r"$t$", ylabel="differential connected magnetic correlation")
    axes[2].legend(frameon=False, ncol=2)
    save(fig, out, "wilson_string_correlations")

    flux = run.correlations[run.correlations.kind == "flux_connected"]
    if len(flux):
        envelope = flux.assign(abs_delta=flux.delta_value.abs()).groupby("time", as_index=False).agg(
            max_abs_connected_delta=("abs_delta", "max"), rms_connected_delta=("delta_value", lambda x: np.sqrt(np.mean(x**2))))
        fig, ax = plt.subplots(figsize=(5.2, 3.4))
        ax.plot(envelope.time, envelope.max_abs_connected_delta, "o-", ms=3, label="maximum")
        ax.plot(envelope.time, envelope.rms_connected_delta, "o-", ms=3, label="RMS over central anchors/links")
        ax.set(xlabel=r"$t$", ylabel="differential connected flux correlation")
        ax.legend(frameon=False)
        save(fig, out, "connected_flux_correlations")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("runs", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, default=Path("figures"))
    parser.add_argument("--times", default="", help="comma-separated heatmap/profile times; default 0, midpoint, final")
    parser.add_argument("--velocity-tmin", type=float)
    parser.add_argument("--velocity-tmax", type=float)
    parser.add_argument("--front-column", choices=("front_rms", "front_quantile"), default="front_rms")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    set_style()
    args.output.mkdir(parents=True, exist_ok=True)
    runs = [load_run(path) for path in args.runs]
    primary = runs[0]
    if args.times:
        requested_times = [float(value) for value in args.times.split(",")]
    else:
        final = float(primary.global_df.time.max())
        requested_times = [0.0, final / 2, final]

    plot_delta_flux(runs, args.output)
    plot_link_heatmaps(primary, args.output, requested_times)
    try:
        plot_radial_front(primary, args.output, requested_times, args.velocity_tmin, args.velocity_tmax, args.front_column)
    except ValueError as error:
        warnings.warn(f"front velocity plot skipped: {error}")
    plot_bond_convergence(runs, args.output)
    plot_trotter_convergence(runs, args.output)
    plot_runtime_diagnostics(runs, args.output)
    plot_classical_scaling(runs, args.output)
    plot_numerical_diagnostics(runs, args.output)
    plot_loop_spreading(primary, args.output)
    plot_selected_observables(primary, args.output)
    print(f"figures written to {args.output.resolve()}")


if __name__ == "__main__":
    main()
