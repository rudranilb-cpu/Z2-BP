#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-$HOME/Desktop/Z2_BP}"
REPO_URL="https://github.com/rudranilb-cpu/Z2-BP.git"
BRANCH="import-benchmark-pipeline"

cd "$ROOT"

if [[ ! -f Project.toml || ! -d src || ! -d scripts ]]; then
  echo "ERROR: $ROOT does not look like the Z2_BP project root." >&2
  exit 1
fi

mkdir -p benchmarks/set1 benchmark_figures/set1

# Curated benchmark directories: compact CSV/TOML/JSON/TXT outputs only.
bench_dirs=(
  paper_K2p2_set1_light
  light_qpu_X_ZZ_K2p2_8x6
  light_qpu_XX_K2p2_8x6
  light_wilson_K2p2_8x6
  light_wilson_K2p2_8x6_D24_Lie_dt005
  light_wilson_K2p2_8x6_D16_Strang_dt0025
  exact_vs_bp_wilson_4x4
  fast_exact_vs_bp_4x4_K2p2_dt005
  fast_exact_vs_bp_4x4_K2p2_dt0025
  fast_flux_8x6_D8 fast_flux_8x6_D16 fast_flux_8x6_D24
  fast_flux_8x6_dt010 fast_flux_8x6_dt005 fast_flux_8x6_dt0025
  fast_flux_8x6_bp1e-4 fast_flux_8x6_bp1e-6 fast_flux_8x6_bp1e-8
  fast_spatial_8x6_K0p5_G1 fast_spatial_8x6_K1_G1 fast_spatial_8x6_K2_G1
  fast_spatial_8x6_K3_G1 fast_spatial_8x6_K4_G1 fast_spatial_8x6_K6_G1
  fast_spatial_8x6_K2p2_G1 fast_spatial_8x6_K2p4_G1
  fast_spatial_8x6_K2p6_G1 fast_spatial_8x6_K2p8_G1
)

rm -rf benchmarks/set1/*
for d in "${bench_dirs[@]}"; do
  if [[ -d "results/$d" ]]; then
    mkdir -p "benchmarks/set1/$d"
    while IFS= read -r -d '' f; do
      rel="${f#results/$d/}"
      mkdir -p "benchmarks/set1/$d/$(dirname "$rel")"
      cp "$f" "benchmarks/set1/$d/$rel"
    done < <(find "results/$d" -type f \( -name '*.csv' -o -name '*.toml' -o -name '*.json' -o -name '*.txt' \) -print0)
  fi
done

if [[ -f results/exact_vs_bp_qpuobs_4x4.csv ]]; then
  cp results/exact_vs_bp_qpuobs_4x4.csv benchmarks/set1/
fi

rm -rf benchmark_figures/set1/*
if [[ -d figures/paper_BP_set1 ]]; then
  find figures/paper_BP_set1 -maxdepth 1 -type f \( -name '*.png' -o -name '*.pdf' -o -name '*.txt' \) -exec cp {} benchmark_figures/set1/ \;
fi

cat > .gitignore <<'EOF'
.DS_Store
*.swp
*.swo
.venv/
__pycache__/
*.pyc
logs/
results/
figures/
*.log
*.tmp
*.aux
*.out
*.toc
scripts/_paper_case.jl
EOF

cat > README.md <<'EOF'
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
EOF

cat > BENCHMARKS.md <<'EOF'
# Benchmark Set 1

Set 1 is the controlled classical BP/PEPS baseline frozen before the main Pauli-propagation and QPU campaign. It is restricted to observables that translate directly into Pauli measurements on the dual-spin register.

## Principal point

- dual lattice: 8 x 6
- `H = -K sum X - Gamma sum ZZ`
- `K/Gamma = 2.2`, `Gamma = 1`
- SCV and one flipped bulk dual spin (elementary electric loop)
- differential convention: `loop - SCV`
- `D = 16`
- `dt = 0.05`
- Lie ordering: X layer followed by ZZ layer
- conservative controlled window: `t <= 0.4`
- `0.4 < t <= 0.5` retained with explicit BP/simple-update systematics

## Numerical controls

The curated outputs contain checks of

1. `D = 8, 16, 24`;
2. BP tolerance `1e-4, 1e-6, 1e-8`;
3. `dt = 0.10, 0.05, 0.025`;
4. Lie versus Strang splitting;
5. exact-versus-BP dynamics on `4 x 4` for `DeltaPhi`, local `X`, endpoint `ZZ` strings, connected `XX`, and Wilson loops.

The bond-dimension and BP-tolerance channels are substantially tighter than the late-time BP/simple-update systematic. Exact-vs-BP comparisons therefore set the conservative time window.

## Gauge/dual observable dictionary

```text
B_p                         <-> X_r
electric link               <-> (1 - Z_r Z_r')/2
short electric-string ends  <-> Z_r Z_r'
W(partial R)                <-> product_{r in R} X_r
```

The Wilson-loop CSVs are therefore area products of dual `X` operators, not products around a dual-lattice boundary.

## Parameter-selection scans

The included `fast_spatial_8x6_K*_G1` directories document the broad and refined scans used to select `K/Gamma = 2.2` as a working point with visible outward propagation before strong local revival dominates.

## Data policy

Only compact, decision-relevant CSV/TOML/JSON/text outputs are versioned. Full raw run directories and logs remain local and can be regenerated from the committed scripts and configurations.
EOF

# Initialize from the GitHub bootstrap commit so histories are related.
if [[ ! -d .git ]]; then
  git init
fi
if git remote get-url origin >/dev/null 2>&1; then
  git remote set-url origin "$REPO_URL"
else
  git remote add origin "$REPO_URL"
fi

git fetch origin main
git checkout -B "$BRANCH" origin/main

# Remove any abandoned import scratch material if present.
if [[ -e .import ]]; then
  git rm -rf .import || rm -rf .import
fi

# Local identity only if one has not already been configured.
if ! git config user.name >/dev/null 2>&1; then
  git config user.name "Rudranil Basu"
fi
if ! git config user.email >/dev/null 2>&1; then
  git config user.email "rudranilb-cpu@users.noreply.github.com"
fi

git add -- .gitignore README.md BENCHMARKS.md Project.toml Manifest.toml
[[ -f requirements-analysis.txt ]] && git add -- requirements-analysis.txt
for p in src bin scripts configs test validation analysis docs benchmarks benchmark_figures; do
  [[ -e "$p" ]] && git add -- "$p"
done

echo "--- Files to be committed ---"
git status --short

git commit -m "Establish validated Z2 BP benchmark pipeline"

# Fast-forward main from the feature branch; never force-push.
git push origin "$BRANCH"
git push origin "$BRANCH":main

echo
echo "Published validated Z2 BP pipeline to $REPO_URL"
echo "Branch retained as $BRANCH for provenance."
