using Pkg
Pkg.activate(".")
using Z2GaugeBenchmarks
using Printf

const Z2 = Z2GaugeBenchmarks

function main()
    common = (
        nx=4, ny=4,
        K=2.2, Gamma=1.0,
        dt=0.05, steps=20,
        trotter=:lie_x_zz,
        boundary=:open_dual,
        defect_x=2, defect_y=2,
        maxdim=16,
        cutoff=1e-10,
        normalize_tensors=true,
        bp_tolerance=1e-6,
        bp_maxiter=20,
        bp_fail_on_nonconvergence=false,
    )

    cfg_bp = SimulationConfig(; common..., backend=:bp)
    cfg_ex = SimulationConfig(; common..., backend=:exact)

    links = dual_links(cfg_bp)
    obs = [("ZZ", [l.v1, something(l.v2)]) for l in links]

    pair_bp, _ = Z2.initialize_bp(cfg_bp)
    pair_ex = Z2.initialize_exact(cfg_ex)
    circuit = build_trotter_layer(cfg_bp, pair_bp.graph)

    bp_flux(p) = [0.5*(1-real(v)) for v in Z2.TNQS.expect(p, obs)]

    function ex_flux(state)
        [0.5*(1-Z2._expect_z(state,
            (l.v1, something(l.v2)), cfg_ex)) for l in links]
    end

    outdir = "results/fast_exact_vs_bp_4x4_K2p2_dt005"
    rm(outdir; recursive=true, force=true)
    mkpath(outdir)

    open(joinpath(outdir, "global.csv"), "w") do f
        println(f, "step,time,delta_exact,delta_bp,abs_error,max_local_flux_error")

        for step in 0:cfg_bp.steps
            es = ex_flux(pair_ex.scv)
            el = ex_flux(pair_ex.loop)
            bs = bp_flux(pair_bp.scv)
            bl = bp_flux(pair_bp.loop)

            de = sum(el) - sum(es)
            db = sum(bl) - sum(bs)
            localerr = maximum(abs.((bl .- bs) .- (el .- es)))

            @printf("step %2d t=%.2f  exact=%.8f  BP=%.8f  |err|=%.3e  local=%.3e\n",
                    step, step*cfg_bp.dt, de, db, abs(db-de), localerr)

            @printf(f, "%d,%.8f,%.12f,%.12f,%.12e,%.12e\n",
                    step, step*cfg_bp.dt, de, db, abs(db-de), localerr)
            flush(f)

            step == cfg_bp.steps && break

            pair_ex, _ = Z2.evolve_exact!(pair_ex, cfg_ex)
            pair_bp, _ = Z2.evolve_bp!(pair_bp, circuit, cfg_bp)
        end
    end

    println("DONE: $outdir/global.csv")
end

main()
