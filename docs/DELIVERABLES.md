# Deliverable index

| Requested item | Location |
|---|---|
| A. Modular Julia codebase | `Project.toml`, `Manifest.toml`, `src/`, `bin/`, `configs/` |
| B. Setup/run README | `README.md`, `requirements-analysis.txt`, `scripts/setup.jl` |
| C. Validation suite | `test/runtests.jl`, `validation/` |
| D. Analysis/plotting | `analysis/common.py`, `analysis/plot_results.py` |
| E. BP production matrix | `docs/PRODUCTION_MATRIX.md`, `scripts/materialize_production_matrix.py`, templates in `configs/` |
| F. Technical report | `docs/TECHNICAL_REPORT.md`, with detailed source audit in `docs/AUDIT.md` |
| G. QPU observable priorities | `docs/OBSERVABLES.md`; raw two-setting manifest in every `pauli.csv` |
| H. Executed tests | `docs/TEST_RESULTS.md` |

Supporting definitions are frozen in `docs/CONVENTIONS.md`; the cross-backend schema
and separate approximation channels are specified in `docs/OUTPUT_SCHEMA.md`. Original
public scripts are preserved unchanged under `legacy/`, with compatibility wrappers at
their former root paths.
