using Pkg
Pkg.activate(".")
using Z2GaugeBenchmarks
using Printf
using Statistics

const Z2 = Z2GaugeBenchmarks

function front_metrics(delta, radii)
    w = abs.(delta)
    total = sum(w)
    total < 1e-14 && return (NaN, NaN)

    rms = sqrt(sum(w .* radii.^2) / total)

    order = sortperm(radii)
    cw = cumsum(w[order])
    idx = searchsortedfirst(cw, 0.8 * total)
    idx = min(idx, length(order))
    q80 = radii[order[idx]]

    return rms, q80
end

function main()
    cfg = SimulationConfig(
        nx=8, ny=6,
        K=3.0, Gamma=1.0,
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

    links = dual_links(cfg)
    defect = resolved_defect(cfg)

    pair, _ = Z2.initialize_bp(cfg)
    circuit = build_trotter_layer(cfg, pair.graph)

    observables = [
        ("ZZ", [l.v1, something(l.v2)])
        for l in links
    ]

    function fluxes(cache)
        vals = Z2.TNQS.expect(cache, observables)
        return [0.5 * (1.0 - real(v)) for v in vals]
    end

    radii = [
        hypot(l.x - defect[1], l.y - defect[2])
        for l in links
    ]

    outdir = "results/fast_spatial_8x6"
    rm(outdir; recursive=true, force=true)
    mkpath(outdir)

    open(joinpath(outdir, "global.csv"), "w") do fg
        open(joinpath(outdir, "links.csv"), "w") do fl
            open(joinpath(outdir, "profiles.csv"), "w") do fp

                println(fg,
                    "step,time,scv_flux,loop_flux,delta_flux,front_rms,front_q80")

                println(fl,
                    "step,time,link_id,v1_x,v1_y,v2_x,v2_y,mid_x,mid_y,radius,scv_flux,loop_flux,delta_flux")

                println(fp,
                    "step,time,radius_bin,mean_delta_flux,mean_abs_delta_flux,nlinks")

                for step in 0:cfg.steps
                    t = step * cfg.dt

                    println("Measuring step $step / $(cfg.steps), t=$t ...")

                    fs = fluxes(pair.scv)
                    floop = fluxes(pair.loop)
                    delta = floop .- fs

                    Phi_s = sum(fs)
                    Phi_l = sum(floop)
                    DeltaPhi = sum(delta)

                    front_rms, front_q80 =
                        front_metrics(delta, radii)

                    @printf(
                        "step %2d  t=%.2f  DeltaPhi=%.8f  front_rms=%.5f  q80=%.5f\n",
                        step, t, DeltaPhi, front_rms, front_q80
                    )

                    @printf(
                        fg,
                        "%d,%.8f,%.12f,%.12f,%.12f,%.12f,%.12f\n",
                        step, t, Phi_s, Phi_l, DeltaPhi,
                        front_rms, front_q80
                    )

                    for (i, l) in enumerate(links)
                        v2 = something(l.v2)
                        @printf(
                            fl,
                            "%d,%.8f,%d,%d,%d,%d,%d,%.6f,%.6f,%.12f,%.12f,%.12f,%.12f\n",
                            step, t, l.id,
                            l.v1[1], l.v1[2],
                            v2[1], v2[2],
                            l.x, l.y, radii[i],
                            fs[i], floop[i], delta[i]
                        )
                    end

                    rbins = sort(unique(round.(radii ./ 0.5) .* 0.5))

                    for rb in rbins
                        inds = findall(
                            i -> abs(round(radii[i] / 0.5) * 0.5 - rb) < 1e-10,
                            eachindex(radii)
                        )

                        isempty(inds) && continue

                        md = mean(delta[inds])
                        mad = mean(abs.(delta[inds]))

                        @printf(
                            fp,
                            "%d,%.8f,%.6f,%.12f,%.12f,%d\n",
                            step, t, rb, md, mad, length(inds)
                        )
                    end

                    flush(fg)
                    flush(fl)
                    flush(fp)

                    step == cfg.steps && break

                    println("Evolving to step $(step+1)...")
                    pair, _ = Z2.evolve_bp!(pair, circuit, cfg)
                end
            end
        end
    end

    println("DONE: $outdir")
end

main()
