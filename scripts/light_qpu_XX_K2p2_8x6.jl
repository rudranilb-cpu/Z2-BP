using Pkg
Pkg.activate(".")
using Z2GaugeBenchmarks
using Printf

const Z2 = Z2GaugeBenchmarks

function main()

    cfg = SimulationConfig(
        nx=8, ny=6,
        K=2.2, Gamma=1.0,
        dt=0.05, steps=10,
        trotter=:lie_x_zz,
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

    defect = resolved_defect(cfg)

    targets = Tuple{String,Int,Tuple{Int,Int}}[]

    for (dir,dx,dy) in [
        ("right",  1, 0),
        ("left",  -1, 0),
        ("up",     0, 1),
        ("down",   0,-1)
    ]
        for d in 1:3
            v = (defect[1] + d*dx, defect[2] + d*dy)

            if 1 <= v[1] <= cfg.nx &&
               1 <= v[2] <= cfg.ny
                push!(targets,(dir,d,v))
            end
        end
    end

    # Only the local X values actually required.
    required_sites = unique(vcat([defect], [v for (_,_,v) in targets]))

    local_obs = [("X",[v]) for v in required_sites]
    xx_obs    = [("XX",[defect,v]) for (_,_,v) in targets]

    # Measure all X-family quantities in one batch.
    xobs = vcat(local_obs,xx_obs)
    nlocal = length(local_obs)

    outdir = "results/light_qpu_XX_K2p2_8x6"
    rm(outdir; recursive=true, force=true)
    mkpath(outdir)

    fc = open(joinpath(outdir,"connected_XX.csv"),"w")
    fd = open(joinpath(outdir,"diagnostics.csv"),"w")

    println(
        fc,
        "step,time,id,direction,distance,target_x,target_y," *
        "scv_xx,loop_xx,delta_xx," *
        "scv_connected,loop_connected,delta_connected"
    )

    println(
        fd,
        "step,time,state,evolution_seconds,trunc_max,trunc_sum," *
        "bp_updates,bp_iterations,residual_max,all_converged"
    )

    println("Initializing BP states...")
    pair, _ = Z2.initialize_bp(cfg)
    circuit = build_trotter_layer(cfg,pair.graph)

    for step in 0:cfg.steps

        t = step * cfg.dt

        if iseven(step)

            println()
            @printf(
                "XX measurement step %d/%d, t=%.2f\n",
                step,cfg.steps,t
            )

            elapsed = @elapsed begin

                vals_scv  = real.(Z2.TNQS.expect(pair.scv,xobs))
                vals_loop = real.(Z2.TNQS.expect(pair.loop,xobs))

                local_scv = Dict(
                    required_sites[i] => vals_scv[i]
                    for i in eachindex(required_sites)
                )

                local_loop = Dict(
                    required_sites[i] => vals_loop[i]
                    for i in eachindex(required_sites)
                )

                for i in eachindex(targets)

                    dir,d,v = targets[i]
                    j = nlocal + i

                    xx_scv  = vals_scv[j]
                    xx_loop = vals_loop[j]

                    c_scv =
                        xx_scv -
                        local_scv[defect] * local_scv[v]

                    c_loop =
                        xx_loop -
                        local_loop[defect] * local_loop[v]

                    @printf(
                        fc,
                        "%d,%.8f,CXX_%s_d%d,%s,%d,%d,%d,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f\n",
                        step,t,
                        dir,d,dir,d,v[1],v[2],
                        xx_scv,
                        xx_loop,
                        xx_loop-xx_scv,
                        c_scv,
                        c_loop,
                        c_loop-c_scv
                    )
                end

                flush(fc)
            end

            @printf(
                "  XX family finished in %.2f s\n",
                elapsed
            )
        end

        step == cfg.steps && break

        println("Evolving to step $(step+1)...")

        pair,evo = Z2.evolve_bp!(
            pair,circuit,cfg
        )

        for (name,e) in (
            ("scv",evo.scv),
            ("loop",evo.loop)
        )

            truncmax =
                isempty(e.errors) ? 0.0 : maximum(e.errors)

            truncsum = sum(e.errors)

            residual =
                isempty(e.bp) ? NaN :
                maximum(s.final_residual for s in e.bp)

            iterations =
                sum(s.iterations for s in e.bp)

            converged =
                all(s.converged for s in e.bp)

            @printf(
                fd,
                "%d,%.8f,%s,%.8f,%.12e,%.12e,%d,%d,%.12e,%s\n",
                step+1,
                (step+1)*cfg.dt,
                name,
                e.runtime,
                truncmax,
                truncsum,
                length(e.bp),
                iterations,
                residual,
                string(converged)
            )
        end

        flush(fd)
    end

    close(fc)
    close(fd)

    println()
    println("DONE: $outdir")
end

main()
