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
