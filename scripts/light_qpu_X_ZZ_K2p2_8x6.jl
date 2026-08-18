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
    sites = Z2.dual_vertices(cfg)

    targets = Tuple{String,Int,Tuple{Int,Int}}[]

    for (dir,dx,dy) in [
        ("right",  1, 0),
        ("left",  -1, 0),
        ("up",     0, 1),
        ("down",   0,-1)
    ]
        for d in 1:3
            v = (defect[1] + d*dx, defect[2] + d*dy)
            if 1 <= v[1] <= cfg.nx && 1 <= v[2] <= cfg.ny
                push!(targets,(dir,d,v))
            end
        end
    end

    xobs = [("X",[v]) for v in sites]
    zobs = [("ZZ",[defect,v]) for (_,_,v) in targets]

    outdir = "results/light_qpu_X_ZZ_K2p2_8x6"
    rm(outdir; recursive=true, force=true)
    mkpath(outdir)

    fx = open(joinpath(outdir,"local_X.csv"),"w")
    fz = open(joinpath(outdir,"electric_strings_ZZ.csv"),"w")
    fd = open(joinpath(outdir,"diagnostics.csv"),"w")

    println(fx,"step,time,x,y,scv_x,loop_x,delta_x")
    println(fz,"step,time,id,direction,distance,target_x,target_y,scv_zz,loop_zz,delta_zz")
    println(fd,"step,time,state,evolution_seconds,trunc_max,trunc_sum,bp_updates,bp_iterations,residual_max,all_converged")

    println("Initializing BP states...")
    pair, _ = Z2.initialize_bp(cfg)
    circuit = build_trotter_layer(cfg,pair.graph)

    for step in 0:cfg.steps

        t = step * cfg.dt

        if iseven(step)

            println()
            @printf("Measurement step %d/%d, t=%.2f\n",step,cfg.steps,t)

            tx = @elapsed begin
                xscv  = real.(Z2.TNQS.expect(pair.scv,xobs))
                xloop = real.(Z2.TNQS.expect(pair.loop,xobs))

                for i in eachindex(sites)
                    v = sites[i]
                    @printf(
                        fx,"%d,%.8f,%d,%d,%.12f,%.12f,%.12f\n",
                        step,t,v[1],v[2],
                        xscv[i],xloop[i],xloop[i]-xscv[i]
                    )
                end
                flush(fx)
            end

            @printf("  local X finished in %.2f s\n",tx)

            tz = @elapsed begin
                zscv  = real.(Z2.TNQS.expect(pair.scv,zobs))
                zloop = real.(Z2.TNQS.expect(pair.loop,zobs))

                for i in eachindex(targets)
                    dir,d,v = targets[i]
                    @printf(
                        fz,"%d,%.8f,E_%s_d%d,%s,%d,%d,%d,%.12f,%.12f,%.12f\n",
                        step,t,dir,d,dir,d,v[1],v[2],
                        zscv[i],zloop[i],zloop[i]-zscv[i]
                    )
                end
                flush(fz)
            end

            @printf("  short ZZ strings finished in %.2f s\n",tz)
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
                truncmax,truncsum,length(e.bp),iterations,
                residual,string(converged)
            )
        end

        flush(fd)
    end

    close(fx)
    close(fz)
    close(fd)

    println()
    println("DONE: $outdir")
end

main()
