# Z2-BP: controlled BP/PEPS benchmarks for 2+1D Z2 gauge dynamics

This repository contains a reproducible classical benchmark pipeline for real-time
\(2+1\)D pure \(\mathbb Z_2\) lattice gauge dynamics in the dual square-lattice
transverse-field Ising representation. The main use case is to prepare controlled
Belief-Propagation/PEPS reference data before Pauli-propagation and quantum-hardware
runs, using observables that can be matched directly across all three methods.

## Physics convention

The finite dual model is

\[
H=-K\sum_r X_r-\Gamma\sum_{\langle r r'\rangle} Z_r Z_{r'},
\]

with electric-link occupation

\[
n_\ell=\frac{1-Z_rZ_{r'}}{2},
\qquad
\Delta\Phi=\Phi_{\rm loop}-\Phi_{\rm SCV}.
\]

A single flipped bulk dual spin represents an elementary four-link electric loop, so
\(\Delta\Phi(0)=4\) exactly. The reference first-order circuit applies the \(X\) layer
followed by the \(ZZ\) layer (`lie_x_zz`), i.e.
\(U_{\rm step}=U_{ZZ}(\delta t)U_X(\delta t)\).
See `docs/CONVENTIONS.md` for boundary and duality details.

## Frozen BP benchmark Set 1

The first validated working point is

```text
dual lattice: 8 x 6, open_dual
K/Gamma:      2.2
Gamma:        1
D:            16
dt:           0.05
Trotter:      lie_x_zz
window:       0 <= t <= 0.5
```

The strongest quantitative validation is through \(t\lesssim0.3\). At \(t=0.4\),
exact-vs-BP errors on the 4x4 validation system remain small for total flux, local
\(X\), and endpoint \(ZZ\) strings, while connected \(XX\) and Wilson observables
reach roughly \(3\times10^{-2}\). At \(t=0.5\), those two richer families reach
roughly \(7\times10^{-2}\), so the last part of the window is retained as a
semi-quantitative extension rather than treated as uniformly controlled.

For the 8x6 baseline, the recorded flux/spatial results are identical to CSV precision
for \(D=8,16,24\), and for BP tolerances \(10^{-4},10^{-6},10^{-8}\). Relative to
Lie `dt=0.05`, halving the Lie step to `dt=0.025` changes \(\Delta\Phi\) by at most
about 0.023 over \(t\le0.5\).

The committed benchmark data are under `benchmarks/set1/`, and the corresponding
paper-ready plots are under `benchmark_figures/set1/`. A numerical summary and the
exact acceptance interpretation are in `docs/BENCHMARK_SET1.md`.

## QPU-compatible observable set

The benchmark archive contains:

- total and link-resolved differential electric flux;
- radial flux profiles;
- local dual \(X_r\) (direct-lattice plaquette/magnetic operator);
- short electric strings, represented by endpoint \(Z_{r_0}Z_r\);
- connected \(XX\) correlators;
- direct-lattice Wilson loops represented as products of dual \(X\)'s over the
  enclosed plaquette area.

These are chosen so that the same Pauli observables can be evaluated with BP, Pauli
propagation, and QPU measurements. See `docs/OBSERVABLES.md`.

## Installation

The tested Julia environment is Julia 1.10.10.

```bash
juliaup add 1.10.10
julia +1.10.10 --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

For validation and plotting:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements-analysis.txt
```

## Tests

```bash
julia +1.10.10 --project=. -e 'using Pkg; Pkg.test()'
Z2_RUN_BP_TESTS=1 julia +1.10.10 --project=. -e 'using Pkg; Pkg.test()'

source .venv/bin/activate
python -m unittest -v \
  validation.test_exact_reference \
  validation.test_compare_runs \
  analysis.test_common
```

The BP integration test checks the actual TNQS gate path, the \(t=0\) identities,
first-step sign convention, BP convergence flags, and truncation diagnostics.

## Reproducing Set 1

The exact scripts used for the lightweight production families are retained in
`scripts/`. The recommended sequence is:

```bash
bash scripts/run_set1_controls.sh
bash scripts/run_set1_qpu_observables.sh
bash scripts/run_set1_validation.sh
```

To regenerate the figures from the **committed** benchmark CSV files:

```bash
source .venv/bin/activate
python scripts/make_benchmark_figures_set1.py
```

Detailed commands and expected outputs are in `docs/REPRODUCTION.md`.

## General simulation package

The reusable package lives in `src/`. A configurable run can be launched with

```bash
julia +1.10.10 --project=. bin/run.jl configs/smoke_exact.toml
```

or, for several materialized TOML configurations,

```bash
julia +1.10.10 --project=. bin/run_sweep.jl config1.toml config2.toml
```

The general full-schema runner writes metadata, global/link/site observables, profiles,
Wilson/string/correlation data, BP updates, truncation diagnostics, and the raw Pauli
manifest. It is intentionally more expensive than the lightweight Set-1 scripts.

## Repository layout

```text
src/                  reusable Julia package
bin/                  configurable run entry points
configs/              smoke/reproduction/production templates
scripts/              lightweight Set-1 and validation drivers
test/                  Julia tests
validation/            independent exact/reference checks
analysis/              general analysis utilities
docs/                  conventions, observables, schemas, benchmark notes
benchmarks/set1/       curated accepted CSV benchmark archive
benchmark_figures/     curated Set-1 figures
legacy/                original scripts retained for provenance
```

Raw local runs, logs, virtual environments, and intermediate figures are deliberately
ignored by Git. Only the curated benchmark archive is versioned.

## Interpretation

BP-assisted simple update is an approximation. A tiny BP residual proves convergence
of the message iteration, not exactness of the PEPS environment. Bond dimension,
BP-environment convergence, Trotter error, finite-size effects, and exact-vs-BP
discrepancy are therefore kept as separate validation channels.

The current repository establishes a controlled classical baseline; it does **not**
make a quantum-advantage claim. The next use is same-observable comparison with
Pauli propagation and improved-noise square-lattice QPU data.
