"""Independent NumPy/SciPy state-vector reference for the frozen dual-TFIM convention."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Iterable

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply


@dataclass(frozen=True)
class ExactConfig:
    nx: int = 3
    ny: int = 3
    K: float = 3.0
    Gamma: float = 1.0
    dt: float = 0.05
    trotter: str = "lie_x_zz"
    boundary: str = "open_dual"
    defect: tuple[int, int] | None = None

    @property
    def nsites(self) -> int:
        return self.nx * self.ny

    @property
    def resolved_defect(self) -> tuple[int, int]:
        return self.defect or ((self.nx + 1) // 2, (self.ny + 1) // 2)


@dataclass(frozen=True)
class Link:
    v1: tuple[int, int]
    v2: tuple[int, int] | None
    orientation: str


def site_index(v: tuple[int, int], nx: int) -> int:
    """Zero-based qubit index; x is the fast coordinate."""
    x, y = v
    return (x - 1) + (y - 1) * nx


def links(cfg: ExactConfig) -> list[Link]:
    out: list[Link] = []
    for y in range(1, cfg.ny + 1):
        for x in range(1, cfg.nx):
            out.append(Link((x, y), (x + 1, y), "x"))
    for y in range(1, cfg.ny):
        for x in range(1, cfg.nx + 1):
            out.append(Link((x, y), (x, y + 1), "y"))
    if cfg.boundary == "fixed_exterior":
        for y in range(1, cfg.ny + 1):
            out.append(Link((1, y), None, "left"))
            out.append(Link((cfg.nx, y), None, "right"))
        for x in range(1, cfg.nx + 1):
            out.append(Link((x, 1), None, "bottom"))
            out.append(Link((x, cfg.ny), None, "top"))
    return out


def initial_states(cfg: ExactConfig) -> tuple[np.ndarray, np.ndarray]:
    dim = 1 << cfg.nsites
    scv = np.zeros(dim, dtype=np.complex128)
    loop = np.zeros(dim, dtype=np.complex128)
    scv[0] = 1.0
    q = site_index(cfg.resolved_defect, cfg.nx)
    loop[1 << q] = 1.0
    return scv, loop


def _zsign(bits: int, q: int) -> float:
    return 1.0 if bits & (1 << q) == 0 else -1.0


def apply_x_layer(state: np.ndarray, cfg: ExactConfig, scale: float = 1.0) -> None:
    alpha = cfg.K * cfg.dt * scale
    c, s = np.cos(alpha), 1j * np.sin(alpha)
    for q in range(cfg.nsites):
        stride = 1 << q
        block = stride << 1
        for base in range(0, state.size, block):
            a0 = state[base : base + stride].copy()
            a1 = state[base + stride : base + block].copy()
            state[base : base + stride] = c * a0 + s * a1
            state[base + stride : base + block] = s * a0 + c * a1


def apply_diagonal_layer(state: np.ndarray, cfg: ExactConfig, scale: float = 1.0) -> None:
    alpha = cfg.Gamma * cfg.dt * scale
    all_links = links(cfg)
    for bits in range(state.size):
        eigenvalue = 0.0
        for link in all_links:
            z1 = _zsign(bits, site_index(link.v1, cfg.nx))
            if link.v2 is None:
                eigenvalue += z1
            else:
                z2 = _zsign(bits, site_index(link.v2, cfg.nx))
                eigenvalue += z1 * z2
        state[bits] *= np.exp(1j * alpha * eigenvalue)


def trotter_step(state: np.ndarray, cfg: ExactConfig) -> None:
    if cfg.trotter == "lie_x_zz":
        apply_x_layer(state, cfg)
        apply_diagonal_layer(state, cfg)
    elif cfg.trotter == "lie_zz_x":
        apply_diagonal_layer(state, cfg)
        apply_x_layer(state, cfg)
    elif cfg.trotter == "strang":
        apply_x_layer(state, cfg, 0.5)
        apply_diagonal_layer(state, cfg)
        apply_x_layer(state, cfg, 0.5)
    else:
        raise ValueError(f"unknown Trotter convention {cfg.trotter}")


def evolve(state: np.ndarray, cfg: ExactConfig, steps: int) -> np.ndarray:
    out = state.copy()
    for _ in range(steps):
        trotter_step(out, cfg)
    return out


def z_product(state: np.ndarray, cfg: ExactConfig, vertices: Iterable[tuple[int, int]]) -> float:
    qubits = [site_index(v, cfg.nx) for v in vertices]
    probabilities = np.abs(state) ** 2
    value = 0.0
    for bits, probability in enumerate(probabilities):
        sign = np.prod([_zsign(bits, q) for q in qubits], dtype=float)
        value += probability * sign
    return float(value)


def link_fluxes(state: np.ndarray, cfg: ExactConfig) -> np.ndarray:
    values = []
    for link in links(cfg):
        vertices = [link.v1] if link.v2 is None else [link.v1, link.v2]
        values.append((1.0 - z_product(state, cfg, vertices)) / 2.0)
    return np.asarray(values)


def total_flux(state: np.ndarray, cfg: ExactConfig) -> float:
    return float(link_fluxes(state, cfg).sum())


def hamiltonian(cfg: ExactConfig) -> sparse.csr_matrix:
    dim = 1 << cfg.nsites
    rows: list[int] = []
    cols: list[int] = []
    data: list[complex] = []
    all_links = links(cfg)
    for bits in range(dim):
        diagonal = 0.0
        for link in all_links:
            z1 = _zsign(bits, site_index(link.v1, cfg.nx))
            z2 = 1.0 if link.v2 is None else _zsign(bits, site_index(link.v2, cfg.nx))
            diagonal += -cfg.Gamma * z1 * z2
        rows.append(bits)
        cols.append(bits)
        data.append(diagonal)
        for q in range(cfg.nsites):
            rows.append(bits ^ (1 << q))
            cols.append(bits)
            data.append(-cfg.K)
    return sparse.coo_matrix((data, (rows, cols)), shape=(dim, dim)).tocsr()


def continuous_evolve(state: np.ndarray, cfg: ExactConfig, final_time: float) -> np.ndarray:
    return expm_multiply((-1j * final_time) * hamiltonian(cfg), state)


def phase_aligned_error(reference: np.ndarray, approximation: np.ndarray) -> float:
    overlap = np.vdot(reference, approximation)
    if abs(overlap) == 0:
        return float(np.linalg.norm(reference - approximation))
    aligned = approximation * np.exp(-1j * np.angle(overlap))
    return float(np.linalg.norm(reference - aligned))


def with_dt(cfg: ExactConfig, dt: float, trotter: str | None = None) -> ExactConfig:
    return replace(cfg, dt=dt, trotter=trotter or cfg.trotter)
