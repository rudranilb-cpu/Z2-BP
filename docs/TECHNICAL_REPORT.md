# Technical report: controlled BP baseline status

## What changed

The original two ad-hoc scripts have been replaced by a Julia package with one frozen
Hamiltonian/duality/sign convention, explicit open-boundary choices, TOML inputs, an
exact small-lattice backend, controlled TNQS/BP simple update, a common observable
plan, versioned outputs, validation tools, and publication plotting. Untouched originals
remain in `legacy/`; compatibility wrappers point users to the new drivers.

The non-cached `dual_gauge.jl` was removed as a validation backend because source
inspection establishes that it used environment-free simple update for two-site gates
and BP by default for loopy contractions. It was neither exact PEPS nor independent of
BP. The controlled comparison is now exact state vector versus explicitly converged,
version-pinned BP.

## Physics corrections

- The reference circuit is X then ZZ, mathematically `U_ZZ U_X`; the original scripts
  used opposite Lie orderings.
- `Delta Phi` is always `loop-SCV`. A bulk one-spin loop has `Delta Phi(0)=4`.
- The public model is an open dual TFIM with 82 internal bonds on 8x6 and 218 on
  12x10. A full open direct lattice requires additional fixed-exterior boundary Z
  fields/fluxes, implemented as a separate configuration.
- `K=3, Gamma=1` is near the square-lattice TFIM critical ratio, approximately
  3.044330(6), rather than a representative deep strong-coupling post-quench point.
- “Curves overlap across D” has been downgraded to observed D stability. It is not a
  convergence conclusion without BP, dt/order, and exact/smaller-system checks.

## Validation completed in this environment

The independent NumPy/SciPy reference and analysis unit tests were executed successfully:

```text
validation/test_exact_reference.py: 5/5 passed
validation/test_compare_runs.py:    4/4 passed
analysis/test_common.py:             2/2 passed
Python byte-compilation:             passed
```

The five physics tests cover link/boundary counting, `Delta Phi(0)=4`, the analytic
one-step formula

\[
\Delta\Phi(dt)=4\cos^2(2Kdt)
\]

for the reference Lie layer, separation of the two Lie orderings after two steps, and
first-/second-order state convergence under dt halving. The reproducible 3x3 scan is
in `validation/trotter_scan_3x3.csv`.

As one useful warning from that small system, at `K/Gamma=3.04433`, `t=1`, and
`dt=.05`, the Lie differential total differs from continuous evolution by about
`8.24e-3`, while the loop-state phase-aligned L2 error is about `7.81e-2`. At
`K/Gamma=6`, the same total-flux error is accidentally only `8.54e-4` even though the
state error is about `8.21e-2`. Therefore total differential flux alone can seriously
underdiagnose Trotter error; local maps, strings, and X/Wilson observables must also be
checked.

Julia 1.10.10 was installed transiently for validation. The checked-in manifest was
instantiated and precompiled, the standard `Pkg.test()` suite passed (21/21
non-BP assertions), and the opt-in 3x3 TNQS/BP integration test passed (9/9). That
integration test covers explicit BP initialization,
the t=0 product-state identities, one controlled `lie_x_zz` step, convergence flags,
nonnegative truncation diagnostics, and agreement with
`4*cos(2*K*dt)^2`. It also compares the sign-sensitive local X result from the actual
TNQS `Rzz` gate path to the exact backend at `1e-10` tolerance. The final BP test took 144.8
seconds in a fresh Julia test process; compilation/setup timing is not reported as
steady-state PEPS timing. The 3x3/4x4 open-dual and 3x3 fixed-exterior exact CLI smoke
runs passed strict schema validation. The relative-coordinate finite-size comparator
and all applicable plotting paths also passed.

The commands remain the required first action on a new production host:

```bash
julia scripts/setup.jl
julia --project=. -e 'using Pkg; Pkg.test()'
Z2_RUN_BP_TESTS=1 julia --project=. -e 'using Pkg; Pkg.test()'
```

## What remains approximate

1. BP-assisted PEPS uses a product-message environment for simple update. A small
   residual establishes a fixed point of that approximation; it does not make the
   environment exact.
2. Finite D introduces local SVD truncation. Per-gate discarded weights are diagnostics,
   not a rigorous accumulated state or observable error.
3. BP contraction of observables is approximate even when D is sufficient. It must be
   compared with exact contraction/state vector where feasible and with tighter BP
   settings independently of D.
4. First-order Trotter dynamics differs from the target Hamiltonian at finite dt. The
   same circuit can still be a valid BP/PP/QPU comparison, but a Hamiltonian-dynamics
   claim requires dt/Strang control.
5. Finite open grids become boundary sensitive after the front reaches the edge. Front
   fits must use a declared pre-boundary window and be checked across sizes.
6. Runtime step 1 includes Julia compilation unless a host-level warm-up protocol is
   used. Steady-state runtime conclusions should exclude/label it.

These channels are emitted separately and are never collapsed into a generic error bar.

## Datasets required before PP/QPU begins

At minimum generate and validate:

1. exact 3x3 and exact/BP 4x4 same-circuit comparisons at confined, critical, and
   deconfined ratios;
2. corrected 8x6 K=3 reproduction with raw local/link data and full diagnostics;
3. 8x6 D sequences and BP-tolerance sequences at `K/Gamma=0.5,3.04433,4.0,6.0`;
4. Lie dt and Strang reference sequences at matched physical times;
5. 6x6/8x6/10x8/12x10 finite-size comparison at the near-critical and spreading
   points;
6. accepted-D 12x10 same-circuit data at K/Gamma `3.0`, `3.04433`, and `4.0`, plus
   the `0.5` easy control; add `2.0` and `6.0` after the first four pass;
7. for every accepted run, all schema tables, raw metadata, convergence comparison
   reports, and the generated figures—not only a plotted central curve.

The main PP/QPU shortlist should be the easy confined calibration (`0.5`), the existing
continuity point (`3.0`), the near-critical stress point (`3.04433`), and the moderately
deconfined spreading point (`4.0`). Select the final flagship only after the recorded
classical difficulty and the measured signal/noise are both favorable. Classical runtime
alone is not evidence of quantum advantage.

## Acceptance state

The package, exact backend, minimal controlled BP path, schema validator, and plotting
pipeline are smoke-tested and ready for staged data generation. It is **not yet valid
to call the old late-time BP data controlled**, and no new 8x6/12x10 production BP
dataset has been generated in this container. The production matrix and fail-fast
validator make the remaining decision points explicit.
