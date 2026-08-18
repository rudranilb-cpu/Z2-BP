using Pkg
Pkg.activate(".")
using Z2GaugeBenchmarks
using Printf
using Statistics

const Z2 = Z2GaugeBenchmarks

function main()

    cfg = SimulationConfig(
        nx=8, ny=6,
        K=2.2, Gamma=1.0,
        dt=0.025, steps=20,
        trotter=:strang,
        boundary=:open_dual,
        backend=:bp,
        defect_x=4, defect_y=3,
        maxdim=16,
        cutoff=1e-10,
        normalize_tensors=true,
        bp_tolerance=1e-6,
        bp_maxiter=20,
        bp_fail_on_nonconvergence=false,
    )

    r = resolved_defect(cfg)
    x,y = r

    loops = [
        ("W1x1",       "1x1", "center",
            [(x,y)], 1, 4),

        ("W2x1_R",      "2x1", "right",
            [(x,y),(x+1,y)], 2, 6),

        ("W2x1_L",      "2x1", "left",
            [(x-1,y),(x,y)], 2, 6),

        ("W1x2_U",      "1x2", "up",
            [(x,y),(x,y+1)], 2, 6),

        ("W1x2_D",      "1x2", "down",
            [(x,y-1),(x,y)], 2, 6),

        ("W2x2_NE",     "2x2", "NE",
            [(x,y),(x+1,y),(x,y+1),(x+1,y+1)], 4, 8),

        ("W2x2_NW",     "2x2", "NW",
            [(x-1,y),(x,y),(x-1,y+1),(x,y+1)], 4, 8),

        ("W2x2_SE",     "2x2", "SE",
            [(x,y-1),(x+1,y-1),(x,y),(x+1,y)], 4, 8),

        ("W2x2_SW",     "2x2", "SW",
            [(x-1,y-1),(x,y-1),(x-1,y),(x,y)], 4, 8)
    ]

    obs = [
        (repeat("X",length(support)),support)
        for (_,_,_,support,_,_) in loops
    ]

    outdir = "results/light_wilson_K2p2_8x6_D16_Strang_dt0025"
    rm(outdir; recursive=true, force=true)
    mkpath(outdir)

    fw = open(joinpath(outdir,"wilson_loops.csv"),"w")
    fa = open(joinpath(outdir,"wilson_averages.csv"),"w")
    fd = open(joinpath(outdir,"diagnostics.csv"),"w")

    println(fw,
        "step,time,id,shape,placement,area,perimeter,support,scv,loop,delta")

    println(fa,
        "step,time,shape,nplacements,scv_mean,loop_mean,delta_mean,delta_std")

    println(fd,
        "step,time,state,evolution_seconds,trunc_max,trunc_sum,bp_updates,bp_iterations,residual_max,all_converged")

    println("Initializing BP states...")
    pair, _ = Z2.initialize_bp(cfg)
    circuit = build_trotter_layer(cfg,pair.graph)

    for step in 0:cfg.steps

        t = step*cfg.dt

        if iseven(step)

            println()
            @printf("Wilson measurement step %d/%d, t=%.2f\n",
                    step,cfg.steps,t)

            elapsed = @elapsed begin

                vscv  = real.(Z2.TNQS.expect(pair.scv,obs))
                vloop = real.(Z2.TNQS.expect(pair.loop,obs))
                delta = vloop .- vscv

                for i in eachindex(loops)

                    id,shape,placement,support,area,perimeter = loops[i]

                    supp = join(
                        ["$(v[1]):$(v[2])" for v in support],
                        "|"
                    )

                    @printf(
                        fw,
                        "%d,%.8f,%s,%s,%s,%d,%d,%s,%.12f,%.12f,%.12f\n",
                        step,t,id,shape,placement,
                        area,perimeter,supp,
                        vscv[i],vloop[i],delta[i]
                    )
                end

                for shape in ("1x1","2x1","1x2","2x2")

                    inds = findall(i -> loops[i][2] == shape,
                                   eachindex(loops))

                    ds = delta[inds]

                    @printf(
                        fa,
                        "%d,%.8f,%s,%d,%.12f,%.12f,%.12f,%.12f\n",
                        step,t,shape,length(inds),
                        mean(vscv[inds]),
                        mean(vloop[inds]),
                        mean(ds),
                        length(ds) > 1 ? std(ds) : 0.0
                    )
                end

                flush(fw)
                flush(fa)
            end

            @printf("  Wilson family finished in %.2f s\n",elapsed)
        end

        step == cfg.steps && break

        println("Evolving to step $(step+1)...")

        pair,evo = Z2.evolve_bp!(pair,circuit,cfg)

        for (name,e) in (("scv",evo.scv),("loop",evo.loop))

            truncmax = isempty(e.errors) ? 0.0 : maximum(e.errors)
            truncsum = sum(e.errors)

            residual = isempty(e.bp) ? NaN :
                maximum(s.final_residual for s in e.bp)

            iterations = sum(s.iterations for s in e.bp)
            converged = all(s.converged for s in e.bp)

            @printf(
                fd,
                "%d,%.8f,%s,%.8f,%.12e,%.12e,%d,%d,%.12e,%s\n",
                step+1,(step+1)*cfg.dt,name,e.runtime,
                truncmax,truncsum,length(e.bp),
                iterations,residual,string(converged)
            )
        end

        flush(fd)
    end

    close(fw)
    close(fa)
    close(fd)

    println()
    println("DONE: $outdir")
end

main()
