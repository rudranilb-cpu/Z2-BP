# Output schema v1.1.0

The unit of interchange is one run directory. PP and QPU producers should write the
same filenames and columns; unavailable deterministic/statistical fields are `nan`,
not silently omitted. `run_id`, `backend`, `step`, and `time` are the common keys.

## Metadata

`metadata.toml` records the full input, schema/git identity (including dirty state),
Julia/TNQS/runtime environment, link/boundary counts, t=0 expectation,
Hamiltonian/duality/gate/sign conventions, and the meanings of six distinct
uncertainty channels. QPU importers should add backend-specific job IDs,
shots, calibration timestamp, mitigation method, and physicality/conditioning flags
without changing the frozen convention fields.

## Tables

- `global.csv`: one row per measured step. `phi_*`, `delta_phi`, normalized densities,
  central-star flux, central-loop survival, absolute-signal weight, mean/RMS/weighted-
  quantile fronts, effective bond participation, and second-moment anisotropy.
  Statistical or bootstrap standard errors have dedicated columns.
- `links.csv`: one row per measured link and step. Stable `link_id`, orientation,
  endpoints, midpoint polar coordinates, `n_scv`, `n_loop`, `delta_n`, and standard
  errors. Blank second endpoints identify fixed-exterior boundary links.
- `sites.csv`: one row per dual site and step. Local X, diagnostic Z, and site-centered
  flux density for both states and their difference, with dedicated uncertainty columns.
- `profiles.csv`: long-format radial/angular bins with bounds, population, signed
  differential mean, absolute differential mean, and bootstrap-error fields.
- `wilson.csv`: rectangle geometry and X-product expectation for each state/difference.
- `correlations.csv`: long-format `electric_string`, `magnetic_connected`, and
  `flux_connected` rows. Anchor/target IDs retain the exact geometry.
- `diagnostics.csv`: one row per state and step. Evolution/measurement wall time,
  achieved D, local truncation statistics, BP summary, approximate norm drift,
  link-bound violations, live Julia heap, and maximum RSS.
- `bp_updates.csv`: one row per explicit message-convergence episode with iteration
  count, final average message change, and convergence boolean.
- `pauli.csv`: one row per raw Pauli word and measured step. Stable observable ID,
  basis, global QPU setting, support, both state values/difference, and uncertainty
  columns make this the direct PP/QPU interchange table; the other tables are derived
  physics views. Only `global_X` and `global_Z` occur, i.e. two measurement settings
  per prepared state and circuit time before any hardware-specific twirling/mitigation.

## Required invariants

`validation/check_run.py` enforces:

1. exactly one t=0 global row;
2. `phi_scv(0)=0` and the metadata-specified initial differential flux;
3. stable run/backend/time keys, no duplicate physical keys, and complete link/site/
   Pauli/profile row counts at every measured step;
4. every physical link occupation lies in `[0,1]` and every raw Pauli mean in
   `[-1,1]` within tolerance;
5. global totals equal both emitted link sums and site-centered sums, and every
   differential column equals `loop-SCV`;
6. the raw Pauli manifest is time-independent and uses only consistent `global_X`
   and `global_Z` settings;
7. nonnegative truncation diagnostics and exactly two state-diagnostic rows per step;
8. no unconverged BP update in strict mode, residuals below the configured tolerance,
   and exact agreement between BP summary/update row counts;
9. no recorded link-bound violation.

An ODR/mitigated QPU estimator can algebraically leave the physical interval. Such a
row must be flagged and must not replace the raw physical estimate in a quantitative
plot. It should fail strict validation by design.

## Comparison rule

Never combine the following into one error bar:

- QPU shot/statistical uncertainty;
- physical finite-size/boundary-window variation;
- PEPS bond-dimension variation;
- BP environment/contraction variation;
- local SVD truncation diagnostics;
- dt/Trotter variation.

`validation/compare_runs.py` requires the caller to name the channel under test and
reports max, RMS, and integrated differences. A convergence claim needs a stable
observable under each relevant channel, not a single visually overlapping plot.
When multiple D runs are plotted, `minimum_tested_D_vs_time.csv` applies the explicit
project thresholds in `PRODUCTION_MATRIX.md` against the largest supplied D. It is a
conditional D-screen only; it is not proof that the largest D, BP environment, or
Trotterization is converged.
The plotting driver also writes `classical_difficulty_summary.csv` for every supplied
run and adds a lattice-size runtime/memory plot when at least two sizes are present.
