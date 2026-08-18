# Configuration files

- `smoke_exact.toml`: fast 3x3 exact state-vector sanity check.
- `smoke_exact_4x4.toml`: 4x4 exact integration/finite-size comparison smoke test.
- `smoke_fixed_exterior_exact.toml`: checks the optional full-open direct-lattice boundary convention.
- `reproduce_8x6_bp.toml`: corrected version of the public 8x6, K=3, Gamma=1 BP run.
- `production_12x10_nearcritical.toml`: a starting point, not a convergence claim. Run the validation matrix first.

The reference first-order ordering is `lie_x_zz`: sequential X gates followed by
diagonal ZZ gates, hence `U_step = U_ZZ U_X`. `strang` is the preferred higher-order
reference. `lie_zz_x` exists only to reproduce and quantify the old baseline mismatch.
