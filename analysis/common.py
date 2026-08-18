"""Shared loader, consistency checks, and front-velocity fits."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import tomllib

import numpy as np
import pandas as pd


@dataclass
class RunData:
    path: Path
    metadata: dict
    global_df: pd.DataFrame
    links: pd.DataFrame
    sites: pd.DataFrame
    profiles: pd.DataFrame
    wilson: pd.DataFrame
    correlations: pd.DataFrame
    diagnostics: pd.DataFrame
    pauli: pd.DataFrame

    @property
    def label(self) -> str:
        sim = self.metadata["simulation"]
        peps = self.metadata["peps"]
        ratio = np.inf if sim["Gamma"] == 0 else sim["K"] / sim["Gamma"]
        return (
            f"{sim['backend']}: {sim['nx']}x{sim['ny']}, "
            f"K/G={ratio:.4g}, dt={sim['dt']:.4g}, "
            f"{sim['trotter']}, D={peps['maxdim']}"
        )


def load_run(path: str | Path) -> RunData:
    path = Path(path)
    with (path / "metadata.toml").open("rb") as handle:
        metadata = tomllib.load(handle)
    version = metadata["schema"]["version"]
    if version.split(".")[0] != "1":
        raise ValueError(f"unsupported schema version {version}")
    return RunData(
        path=path,
        metadata=metadata,
        global_df=pd.read_csv(path / "global.csv"),
        links=pd.read_csv(path / "links.csv"),
        sites=pd.read_csv(path / "sites.csv"),
        profiles=pd.read_csv(path / "profiles.csv"),
        wilson=pd.read_csv(path / "wilson.csv"),
        correlations=pd.read_csv(path / "correlations.csv"),
        diagnostics=pd.read_csv(path / "diagnostics.csv"),
        pauli=pd.read_csv(path / "pauli.csv"),
    )


def fit_velocity(global_df: pd.DataFrame, tmin: float | None = None, tmax: float | None = None,
                 column: str = "front_quantile") -> dict[str, float]:
    clean = global_df[["time", column]].replace([np.inf, -np.inf], np.nan).dropna()
    if tmin is None:
        positive = clean.loc[clean.time > 0, "time"]
        tmin = float(positive.min()) if len(positive) else float(clean.time.min())
    if tmax is None:
        tmax = float(clean.time.quantile(0.65))
    window = clean[(clean.time >= tmin) & (clean.time <= tmax)]
    if len(window) < 3:
        eligible = clean[clean.time >= tmin]
        if len(eligible) >= 3:
            tmax = float(eligible.time.iloc[2])
            window = clean[(clean.time >= tmin) & (clean.time <= tmax)]
    if len(window) < 3:
        raise ValueError("velocity fit needs at least three finite time points")
    x = window.time.to_numpy(float)
    y = window[column].to_numpy(float)
    design = np.column_stack([x, np.ones_like(x)])
    slope, intercept = np.linalg.lstsq(design, y, rcond=None)[0]
    residual = y - (slope * x + intercept)
    dof = len(x) - 2
    sigma2 = float(residual @ residual / dof) if dof > 0 else np.nan
    covariance = sigma2 * np.linalg.inv(design.T @ design) if dof > 0 else np.full((2, 2), np.nan)
    total = float(np.sum((y - y.mean()) ** 2))
    r2 = np.nan if total == 0 else 1.0 - float(residual @ residual) / total
    return {
        "velocity": float(slope),
        "velocity_fit_stderr": float(np.sqrt(covariance[0, 0])),
        "intercept": float(intercept),
        "r_squared": float(r2),
        "tmin": float(tmin),
        "tmax": float(tmax),
        "npoints": int(len(window)),
        "front_column": column,
    }


def nearest_steps(global_df: pd.DataFrame, requested_times: list[float]) -> list[int]:
    available = global_df[["step", "time"]].drop_duplicates().sort_values("time")
    steps: list[int] = []
    for time in requested_times:
        index = (available.time - time).abs().idxmin()
        step = int(available.loc[index, "step"])
        if step not in steps:
            steps.append(step)
    return steps


def convergence_key(run: RunData, omit: tuple[str, ...] = ()) -> tuple:
    sim = run.metadata["simulation"]
    fields = ("nx", "ny", "K", "Gamma", "dt", "steps", "trotter", "boundary", "defect_x", "defect_y")
    return tuple((field, sim[field]) for field in fields if field not in omit)


def channel_signature(run: RunData, channel: str) -> tuple:
    """Immutable metadata signature after removing only the named error channel."""
    allowed = {
        "bond_dimension": {"peps.maxdim"},
        "trotter": {"simulation.dt", "simulation.steps", "simulation.trotter",
                    "observables.measure_every"},
    }[channel]
    entries = []
    for section in ("simulation", "peps", "bp", "observables", "conventions"):
        for key, value in sorted(run.metadata.get(section, {}).items()):
            path = f"{section}.{key}"
            if path not in allowed:
                entries.append((path, repr(value)))
    return tuple(entries)
