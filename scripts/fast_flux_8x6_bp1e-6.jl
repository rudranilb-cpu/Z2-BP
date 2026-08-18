using Pkg
Pkg.activate(".")
using Z2GaugeBenchmarks
using Printf

const Z2 = Z2GaugeBenchmarks

cfg = SimulationConfig(
    nx=8, ny=6,
    K=3.0, Gamma=1.0,
    dt=0.05, steps=20,
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

println("Initializing BP states...")
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

outdir = "results/fast_flux_8x6_bp1e-6"
mkpath(outdir)

open(joinpath(outdir, "global.csv"), "w") do fg
    println(fg, "step,time,scv_flux,loop_flux,delta_flux")

    for step in 0:cfg.steps
        println("Measuring step $step / $(cfg.steps)...")
        fs = fluxes(pair.scv)
        fl = fluxes(pair.loop)

        Φs = sum(fs)
        Φl = sum(fl)
        ΔΦ = Φl - Φs

        @printf("step %2d  t=%.2f  SCV=%.6f  LOOP=%.6f  Delta=%.6f\n",
                step, step*cfg.dt, Φs, Φl, ΔΦ)

        @printf(fg, "%d,%.8f,%.12f,%.12f,%.12f\n",
                step, step*cfg.dt, Φs, Φl, ΔΦ)
        flush(fg)

        step == cfg.steps && break

        println("Evolving to step $(step+1)...")
        Z2.evolve_bp!(pair, circuit, cfg)
    end
end

println("DONE: $outdir/global.csv")
