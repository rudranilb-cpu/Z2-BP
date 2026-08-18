# Source audit and dispositions

Audited inputs:

- public repository `main` at `5f679859...` (three files);
- `legacy/dual_gauge_original.jl` (original non-cached script);
- `legacy/Belief_Propagation_Dual_gauge_original.jl` (original TNQS/BP script);
- the attached consolidated manuscript and `CONSISTENCY_AUDIT.md`;
- TensorNetworkQuantumSimulator.jl source at
  `b5d4089849de1cc23806aa8325e8db56a55f2e0b` (package version 0.4.4);
- ITensorNetworks.jl source at
  `69ef2f64f090ade119768b5c1f6b5eb6131514d1` for behavior verification only.

The tables below account for every executable block in the two original scripts.

## Original `dual_gauge.jl`

| Original lines | Audit | Disposition |
|---|---|---|
| 5-13 imports | Dependencies were neither pinned nor activated. `Observers` and `LinearAlgebra.norm` were unused directly. | Replaced by a package project and pinned TNQS setup. Original archived unchanged. |
| 15-23 globals | All physics, PEPS, and run parameters were global; functions silently captured `K_coup` and `Gamma`. | Replaced by validated `SimulationConfig` loaded from TOML. |
| 26-30 lattice | `named_grid((Nx,Ny))` is an open **dual** grid, not a full open direct lattice. Coordinate labels are `(x,y)`. | Retained as `open_dual`; the alternative `fixed_exterior` model is explicit. Stable link enumeration is tested. |
| 33-39 ZZ gates | `h_j=-Gamma ZZ` and `exp(-i dt h_j)=exp(+i Gamma dt ZZ)` have the correct sign. Gates were appended first. | Sign retained. Order is exposed as `lie_zz_x`, but is no longer the reference. |
| 40-44 X gates | `h_j=-K X` gives the correct `exp(+i K dt X)`. Gates were appended after ZZ. | Sign retained; reference `lie_x_zz` lists X first. |
| 33-45 layer API | `dt` was an argument but `K_coup` and `Gamma` were globals. No Strang layer existed. | All couplings/configuration are explicit; Lie in both orders and Strang are implemented. |
| 48-65 flux | `(1-<ZZ>)/2` is the correct internal-link flux. Division by an estimated norm is required because the PEPS tensors are not globally normalized. | Same operator, but all links are emitted individually before forming totals. Bounds and link sums are checked. |
| 52,60 contraction | `inner` did not specify an algorithm. Current ITensorNetworks defaults loopy contractions to BP; `maxiter=100` alone does not establish a residual tolerance. | This script is not treated as an independent/exact PEPS reference. Exact state vector and controlled cached BP replace it. |
| 55-60 measurement loop | One approximate contraction was performed per edge, with no shared convergence record. | TNQS cache is converged explicitly and all requested Pauli values use the same converged cache. |
| 75 lattice construction | Correct open dual geometry. | Retained. |
| 78 SCV | All `Up` is the dual representative of the no-electric-flux state. | Retained. |
| 79 MID | `(4,3)` on 8x6 is a bulk degree-four site, but not a unique geometric center on an even-by-even grid. | Default is `(cld(Nx,2),cld(Ny,2))`; exact coordinate stored in metadata. |
| 81 gate list | Constructed ZZ then X, opposite to the BP script. Current ITensorNetworks applies vectors sequentially. | Mismatch made selectable and tested to differ after two steps. Reference frozen to X then ZZ. |
| 89-99 t=0 | Correct sign `MID-SCV`; no hard assertion that it equals four. | Fail-fast t=0 assertion added, generalized to defect coordination/boundary choice. |
| 102-106 update | Vector application is sequential. Current two-site application without `envs` is environment-free simple update, not exact PEPS evolution. No truncation callback was used. | Removed from controlled backends rather than presenting it as independent validation. Original retained for provenance. |
| 108-118 measurement | Same BP-contraction caveat; no residual or truncation output. | Replaced by controlled diagnostics. |
| 120-133 export | CSV omitted physical time, metadata, D-independent IDs, local observables, errors, and conventions. | Replaced by schema v1 directory with metadata and eight tables. |

## Original `Belief_Propagation_Dual_gauge.jl`

| Original lines | Audit | Disposition |
|---|---|---|
| 3-4 activation | Hard-coded `C:\Users\dell\TNQS_Project`; non-portable and not the repository environment. | Removed. `scripts/setup.jl` pins the tested TNQS commit. |
| 6-11 imports | No project/manifest; package behavior could change under the same script. | Direct dependencies and compat are recorded; TNQS source revision pinned. |
| 13-20 globals | Fixed to one 8x6, K=3, Gamma=1, dt=.05 run. | TOML configuration covers geometry, physics, D, BP, defect, output, profiles, and method. |
| 22-27 X layer | For TNQS 0.4.4, `Rx(theta)=exp(-i theta X/2)`; `theta=-2Kdt` is correct. X is listed first. | Retained exactly. Unit test checks the analytic first-step flux. |
| 28-31 ZZ layer | TNQS 0.4.4 deliberately exposes Qiskit `Rzz(theta)` and rescales to the ITensor half-angle internally. `theta=-2 Gamma dt` is correct. Four colors suffice for a square grid; ideal ZZ gates commute. | Retained. Controlled driver records BP/truncation consequences of each scheduled layer. |
| 36-42 flux | Correct internal-link definition and batching. | Retained and expanded to local maps, profiles, strings, loops, and correlations. |
| 44-53 run setup | Geometry/state are correct. `(4,3)` has degree four but is only one of four central sites. | Configurable and metadata-backed. |
| 55-57 caches | Caches were created but not explicitly converged before t=0 measurement. Default messages happen to give exact product-state t=0 values, but this is not a general measurement protocol. | Both caches are explicitly warm-started to tolerance before any measurement. |
| 59 apply kwargs | D/cutoff are meaningful local SVD controls. Tensor normalization is a gauge/numerical operation, not proof of global normalization. | Retained; approximate BP norm drift and local discarded weights are recorded separately. |
| 64-74 t=0 | Correct differential sign and values in ideal arithmetic; no exact assertion. Timing measured only measurement. | Fail-fast check; measurement time and evolution time have separate columns. |
| 77-84 evolution | Upstream `apply_gates` performs BP-assisted simple update and refreshes messages between overlapping two-site updates, but does not return BP iterations/residuals. | Version-pinned controlled driver reproduces the scheduling while performing one BP iteration at a time and recording the final average projective message change. |
| 86-90 sign/storage | `MID-SCV` is correct and conflicts with legacy hardware plots labeled `SCV-Mid`. | Canonical sign frozen in metadata and column name `delta_*`; legacy opposite sign is not silently imported. |
| 91 truncation | Maximum local SVD discarded weight is useful but not a global state/observable error. | Max, sum, and count are saved as diagnostics only. |
| 92 max D | Actual achieved virtual dimension is useful. | Saved per state and step. |
| 94-95 log | Diagnostics were printed but mostly not archived. | All archived in machine-readable tables. |
| 98 return | Only totals/timing returned. | Complete standardized output. |
| 101-103 sweep/output | D list was hard-coded, excluded the D=34 curve used elsewhere, and output used a Windows desktop. | Sweep takes arbitrary TOML files; output root/run IDs portable and unique. |
| 108-115 CSV | No time, BP residuals, truncation, local data, or metadata. | Replaced by schema v1. |

## Analytical/manuscript discrepancies corrected

1. The two scripts implemented opposite finite-dt Lie products. The difference is not
   stylistic and is nonzero after two steps; the independent test checks this.
2. The canonical sign is `loop - SCV`; some legacy figures use its negative.
3. A bulk one-spin defect gives exactly four units of differential flux on every grid
   large enough to contain a degree-four site. Lattice size alone cannot change it.
4. The public finite model omits direct-lattice boundary electric links. Interpreting
   its total as the flux of a full open direct lattice is incorrect.
5. The old non-cached script is neither an exact PEPS evolution nor a BP-independent
   contraction reference.
6. Visual agreement across D is only **observed D stability**. It says nothing by
   itself about the BP environment, Trotter error, or late-time physical accuracy.
7. `K=3, Gamma=1` is close to the thermodynamic TFIM critical ratio, not an arbitrary
   generic point. This makes it a useful stress point but not a phase scan.
8. Absolute fluxes across different lattice sizes are extensive and not directly
   commensurate. Use same-size runs, flux density, or pre-boundary local/differential
   data.
9. The D=34 curve and its raw data are absent from the public repository; it cannot be
   treated as reproducible evidence until archived with configuration/metadata.
