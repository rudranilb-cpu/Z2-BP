# Z2-BP: controlled PEPS/BP benchmarks for 2+1D Z2 gauge dynamics

This repository contains the Julia/Python benchmark pipeline for real-time dynamics of the (2+1)-dimensional pure Z2 lattice gauge theory in its dual transverse-field Ising representation. The numerical program starts from a localized elementary electric-flux loop and produces observables designed for direct comparison with Pauli-propagation and quantum-hardware calculations.

## Model and convention

The dual Hamiltonian is

```text
H = -K sum_r X_r - Gamma sum_<rr'> Z_r Z_r'
```

The differential electric flux is defined as `loop - SCV`. For a single flipped bulk dual spin, `DeltaPhi(0) = 4`. The principal finite-lattice calculations use the open dual-lattice convention.

## Frozen BP Set 1

```text
Nx x Ny = 8 x 6
K/Gamma = 2.2
Gamma = 1
D = 16
dt = 0.05
Trotter = Lie (X then ZZ)
```

The full observable set is treated as quantitatively controlled through approximately `t <= 0.4`; data in `0.4 < t <= 0.5` are retained as a semi-quantitative extension with explicit BP/simple-update systematic uncertainty.

Controls include `D = 8,16,24`, BP tolerances `1e-4,1e-6,1e-8`, Lie timesteps `0.10,0.05,0.025`, Strang references, and exact-vs-BP comparisons on `4 x 4`.

## QPU-compatible observables

- local electric-link flux and total differential flux;
- spatial flux maps and radial profiles;
- local dual `X` (plaquette/magnetic) observables;
- short electric strings represented by endpoint `ZZ`;
- connected `XX` correlations;
- direct-lattice Wilson loops represented as products of dual `X` over enclosed plaquettes.

## Layout

```text
src/                  Julia package code
bin/                  command-line drivers
scripts/              production, validation, and plotting scripts
configs/              reproducible run configurations
test/                  Julia tests
validation/            exact-reference and comparison tests
analysis/              analysis utilities
benchmarks/set1/       curated validated benchmark outputs
benchmark_figures/     frozen figure set generated from curated data
```

Raw/intermediate `results/`, `figures/`, and `logs/` are intentionally ignored by Git.

## Environment

The benchmarked environment uses Julia 1.10.10.

```bash
juliaup add 1.10.10
julia +1.10.10 --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

For the Python analysis utilities:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements-analysis.txt
```

## Tests

```bash
julia +1.10.10 --project=. -e 'using Pkg; Pkg.test()'
Z2_RUN_BP_TESTS=1 julia +1.10.10 --project=. -e 'using Pkg; Pkg.test()'
python -m unittest -v validation.test_exact_reference validation.test_compare_runs analysis.test_common
```

See [BENCHMARKS.md](BENCHMARKS.md) for numerical controls, observable definitions, and Set-1 provenance.
