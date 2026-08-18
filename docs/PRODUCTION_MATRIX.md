# Recommended staged BP production matrix

Set `Gamma=1` throughout so times and ratios are comparable. The phase labels below
refer to the post-quench Hamiltonian; every run still starts from the same Z-polarized
SCV/one-loop pair. The square-lattice critical ratio is
`K/Gamma = 3.044330(6)`.

Do not run the full 12x10 matrix first. Advance a point only after its smaller-stage
checks pass.

## Stage 0: exact convention and circuit validation

| Purpose | Lattice | K/Gamma | dt | tmax | Order |
|---|---:|---:|---:|---:|---|
| signs, t=0, one/two steps | 3x3 | 0.5, 1.5, 3.0, 3.04433, 4.0, 6.0 | 0.10, 0.05, 0.025 | 1.0 | both Lie + Strang |
| BP versus exact same circuit | 4x4 | 0.5, 3.04433, 6.0 | 0.05 | 0.5 | `lie_x_zz` |
| boundary interpretation | 3x3 | 3.04433 | 0.05 | 0.5 | Strang; both boundary modes |

For 4x4 BP use `D=4,8,16`, BP tolerance `1e-10`, maxiter 100. At a D with no
discarded weight, any remaining exact/BP observable difference is an environment/
contraction effect, not a bond-dimension error.

## Stage 1: 8x6 control and coarse physics scan

### 1A. Coarse regimes

Run `D=16`, `dt=0.05`, `tmax=1.0`, BP tolerance `1e-8`, maxiter 50:

| K/Gamma | Regime/purpose |
|---:|---|
| 0.5 | deep confined; near-eigenstate/easy control |
| 2.0 | interacting confined quench |
| 3.0 | exact continuity with existing figures; near-critical stress point |
| 3.04433 | thermodynamic critical benchmark |
| 4.0 | moderately deconfined, strong spreading candidate |
| 6.0 | deep transverse-field/deconfined control; fast local rotation |

The points likely to be hardest are the near-critical/intermediate quenches, but that
must be diagnosed from D growth, BP residuals, and operator/entanglement growth—not
assumed from equilibrium criticality alone.

### 1B. Bond-dimension channel

At `K/Gamma = 0.5, 3.04433, 4.0, 6.0`, run
`D = 8, 12, 16, 24, 32` at fixed `dt=0.05`, `tmax=1.0`, tolerance `1e-8`.
If D32 is not stable against D24 in the principal observables, add D40 only at the
failing points/times rather than everywhere.

At `K/Gamma=3.04433`, `D=24`, also repeat `cutoff=1e-8,1e-10,1e-12` through the
accepted time window. This is a truncation-control sensitivity test, not a substitute
for the D sequence. The plotting driver writes a threshold-based
`minimum_tested_D_vs_time.csv`; interpret it only relative to the largest supplied D.

### 1C. BP-environment channel

At `D=24` and `K/Gamma = 0.5, 3.04433, 4.0`, compare:

| tolerance | maxiter |
|---:|---:|
| 1e-6 | 25 |
| 1e-8 | 50 |
| 1e-10 | 100 |

A run that reaches maxiter without meeting tolerance is **not** an error-bar estimate;
it is an unconverged contraction/update and must be extended or excluded.

### 1D. Trotter channel

At `D=24`, compare Lie `dt=0.10,0.05,0.025` and Strang `dt=0.05,0.025` at fixed
physical times for `K/Gamma = 0.5,3.04433,4.0,6.0`. The exact 3x3 tests show why a
single total-flux curve is insufficient: at K=6 and t=1, the Lie `dt=.05` loop-state
vector error is appreciable even when the small-lattice differential total happens to
show a cancellation. Check local maps, Wilson terms, and strings as well.

## Stage 2: finite-size/boundary channel

At `K/Gamma = 3.04433` and `4.0`, use the smallest D already demonstrated stable and
run `6x6, 8x6, 10x8, 12x10` at `dt=.05`, `tmax=1`. Compare only common central regions
and differential/radial quantities using `compare_runs.py --channel finite_size`.
Keep D/BP/dt fixed within each size pair so the channel is not mixed. Extend t only while:

- the weighted front remains separated from the boundary;
- the 10x8 and 12x10 central/link profiles agree at matched coordinates;
- D and BP criteria remain satisfied.

## Stage 3: 12x10 data to generate before PP/QPU

The minimum same-circuit archive is:

| Priority | K/Gamma | dt/order | tmax | Role |
|---:|---:|---|---:|---|
| 1 | 3.0 | .05 Lie | 1.0 | reproduce/replace present BP curve with full metadata |
| 1 | 3.04433 | .05 Lie | 1.0 | near-critical stress/crossover target |
| 1 | 4.0 | .05 Lie | 1.0 | deconfined spreading target |
| 2 | 0.5 | .05 Lie | 1.0 | confined/easy hardware control |
| 2 | 2.0 | .05 Lie | 1.0 | intermediate confined comparison |
| 2 | 6.0 | .05 Lie | 1.0 | high-field control; scrutinize Trotter error |
| reference | 3.04433, 4.0 | .025 Strang | 1.0 | higher-order Hamiltonian reference |

For every row, choose D from Stage 1/2 rather than fixing D24 by fiat. Archive at
least the largest accepted D and its next-lower comparison. Only after t<=1 passes
should the hardest one or two points be extended to `t=1.5` and then `t=2.0`.

## Quantitative acceptance targets

These are project tolerances, not statistical confidence intervals:

- exact t=0 identities pass at machine/BP tolerance;
- no link occupation outside `[0,1]` beyond `1e-8` for deterministic data;
- every recorded BP update converges, and tightening tolerance changes principal
  observables less than the chosen physics tolerance;
- between the two largest D values over the accepted window:
  `max |delta Phi_D-delta Phi_ref| <= 0.02`, local-link RMS <= `5e-3`, local-link
  max <= `0.02`, front-radius difference <= `0.05`, and selected Wilson/string
  differences <= `0.02`;
- halving dt and comparing to Strang meets the same observable tolerances if the
  result is claimed as Hamiltonian dynamics. Otherwise label it explicitly as the
  implemented Trotter circuit benchmark;
- local discarded weights decrease with D and remain small, but are never promoted
  to a global error bar;
- velocity is stable under at least two front definitions and reasonable pre-boundary
  fit-window shifts.

If a scientifically meaningful signal is smaller than these targets, tighten the
relevant channel rather than normalizing the discrepancy away.
