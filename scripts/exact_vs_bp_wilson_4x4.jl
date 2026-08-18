using Pkg
Pkg.activate(".")
using Z2GaugeBenchmarks
using Printf

const Z2 = Z2GaugeBenchmarks

function main()

    common = (
        nx=4, ny=4,
        K=2.2, Gamma=1.0,
        dt=0.05, steps=10,
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

    x,y = resolved_defect(cfg_bp)

    loops = [
        ("W1x1","1x1",[(x,y)]),

        ("W2x1_R","2x1",[(x,y),(x+1,y)]),
        ("W2x1_L","2x1",[(x-1,y),(x,y)]),

        ("W1x2_U","1x2",[(x,y),(x,y+1)]),
        ("W1x2_D","1x2",[(x,y-1),(x,y)]),

        ("W2x2_NE","2x2",
            [(x,y),(x+1,y),(x,y+1),(x+1,y+1)]),

        ("W2x2_NW","2x2",
            [(x-1,y),(x,y),(x-1,y+1),(x,y+1)]),

        ("W2x2_SE","2x2",
            [(x,y-1),(x+1,y-1),(x,y),(x+1,y)]),

        ("W2x2_SW","2x2",
            [(x-1,y-1),(x,y-1),(x-1,y),(x,y)])
    ]

    obs = [
        (repeat("X",length(support)),support)
        for (_,_,support) in loops
    ]

    pair_ex = Z2.initialize_exact(cfg_ex)
    pair_bp, _ = Z2.initialize_bp(cfg_bp)
    circuit = build_trotter_layer(cfg_bp,pair_bp.graph)

    outdir = "results/exact_vs_bp_wilson_4x4"
    rm(outdir; recursive=true, force=true)
    mkpath(outdir)

    f = open(joinpath(outdir,"wilson_exact_vs_bp.csv"),"w")

    println(
        f,
        "step,time,id,shape,exact_scv,exact_loop,exact_delta," *
        "bp_scv,bp_loop,bp_delta,delta_abs_error"
    )

    for step in 0:cfg_bp.steps

        t = step*cfg_bp.dt

        if iseven(step)

            ex_scv = [
                Z2._expect_x(pair_ex.scv,Tuple(support),cfg_ex)
                for (_,_,support) in loops
            ]

            ex_loop = [
                Z2._expect_x(pair_ex.loop,Tuple(support),cfg_ex)
                for (_,_,support) in loops
            ]

            bp_scv = real.(Z2.TNQS.expect(pair_bp.scv,obs))
            bp_loop = real.(Z2.TNQS.expect(pair_bp.loop,obs))

            println()
            @printf("t = %.2f\n",t)

            for i in eachindex(loops)

                id,shape,_ = loops[i]

                de = ex_loop[i] - ex_scv[i]
                db = bp_loop[i] - bp_scv[i]
                err = abs(db-de)

                @printf(
                    "  %-8s exact=% .7f BP=% .7f err=%.3e\n",
                    id,de,db,err
                )

                @printf(
                    f,
                    "%d,%.8f,%s,%s,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12e\n",
                    step,t,id,shape,
                    ex_scv[i],ex_loop[i],de,
                    bp_scv[i],bp_loop[i],db,err
                )
            end

            flush(f)
        end

        step == cfg_bp.steps && break

        pair_ex, _ = Z2.evolve_exact!(pair_ex,cfg_ex)
        pair_bp, _ = Z2.evolve_bp!(pair_bp,circuit,cfg_bp)
    end

    close(f)

    println()
    println("DONE: $outdir/wilson_exact_vs_bp.csv")
end

main()
