# Test record for this delivery

Date: 2026-08-15

Executed in the supplied work container:

```text
$ python -m unittest -v validation/test_exact_reference.py
Ran 5 tests in 0.295s — OK

$ python -m unittest -v analysis/test_common.py
Ran 2 tests in 0.029s — OK

$ python -m unittest -v validation.test_exact_reference \
    validation.test_compare_runs analysis.test_common
Ran 11 tests — OK

$ python -m py_compile analysis/common.py analysis/plot_results.py \
    validation/exact_reference.py validation/check_run.py validation/compare_runs.py \
    scripts/materialize_production_matrix.py
exit code 0

$ python scripts/materialize_production_matrix.py --stage stage1_D \
    --output /tmp/z2-matrix-stage1D
20 TOMLs materialized; all parsed by Python tomllib and Julia load_config

$ python validation/trotter_scan.py --output /tmp/z2-trotter-scan-final.csv
byte-identical to validation/trotter_scan_3x3.csv

$ JULIA_DEPOT_PATH=/tmp/julia-depot julia --project=. \
    -e 'using Pkg; Pkg.test()'
lattice and t=0 identities: 11/11 passed
gate sign and first Lie step: 5/5 passed
fixed exterior counting: 5/5 passed
Z2GaugeBenchmarks tests passed

$ Z2_RUN_BP_TESTS=1 JULIA_DEPOT_PATH=/tmp/julia-depot julia --project=. \
    -e 'using Pkg; Pkg.test()'
all 21 non-BP assertions passed
BP minimal integration: 9/9 passed (144.8 s in the final fresh test process)

$ JULIA_DEPOT_PATH=/tmp/julia-depot julia --project=. \
    bin/run.jl configs/smoke_exact.toml
completed 3x3 open-dual exact run (steps 0 through 4)

$ JULIA_DEPOT_PATH=/tmp/julia-depot julia --project=. \
    bin/run.jl configs/smoke_fixed_exterior_exact.toml
completed 3x3 fixed-exterior exact run (steps 0 through 4)

$ JULIA_DEPOT_PATH=/tmp/julia-depot julia --project=. \
    bin/run.jl configs/smoke_exact_4x4.toml
completed 4x4 open-dual exact run (steps 0 through 4)

$ python validation/check_run.py results/<smoke-run>
validated

$ python validation/compare_runs.py results/<3x3> results/<4x4> \
    --channel finite_size
relative-coordinate finite-size report generated with all 3x3 links/sites matched

$ python analysis/plot_results.py results/<smoke-run> --output /tmp/<figures>
all single-run physics/diagnostic figures generated, including loop survival,
spreading, and connected flux; D and dt comparison panels correctly skipped because
no comparison runs were supplied
```

Julia was installed transiently as version 1.10.10, matching `Manifest.toml`. These
smoke tests certify the exercised API paths, not late-time or large-D accuracy. Repeat
the Julia tests on the production host before starting the BP matrix.
