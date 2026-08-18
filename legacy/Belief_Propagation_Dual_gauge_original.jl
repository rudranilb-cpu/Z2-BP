# Z2 Gauge Theory - Total Electric Flux (Simple Update w/ BP)

using Pkg
Pkg.activate("C:\\Users\\dell\\TNQS_Project")

using TensorNetworkQuantumSimulator
using NamedGraphs
using NamedGraphs.NamedGraphGenerators: named_grid
using Graphs
using Printf
using LinearAlgebra

# Sim parameters
const Nx = 8
const Ny = 6
const K_coup = 3
const Gamma = 1
const dt = 0.05
const num_steps = 20
const CUTOFF = 1e-10

function build_trotter_layer(g, dt)
    layer = []
    
    # H = -K_coup*X - Gamma*ZZ
    append!(layer, ("Rx", [v], -2 * K_coup * dt) for v in vertices(g))

    # Two-site ZZ rotations (edge-colored for parallel BP efficiency)
    for colored_edges in edge_color(g, 4)
        append!(layer, ("Rzz", pair, -2 * Gamma * dt) for pair in colored_edges)
    end

    return layer
end

function measure_total_flux(ψ_bpc, g)
    # Total flux = Σ_edges (1 - <ZZ>)/2
    observables = [("ZZ", [src(e), dst(e)]) for e in edges(g)]
    zz_vals = expect(ψ_bpc, observables)
    
    return sum(0.5 * (1.0 - real(zz)) for zz in zz_vals)
end

function run_simulation(dbond)
    println("\n--- Starting run: $(Nx)x$(Ny), K=$K_coup, Γ=$Gamma, DBOND=$dbond ---")

    g = named_grid((Nx, Ny))

    # Init product states: SCV (all up) and MID (flipped spin at 4,3)
    psi_scv = tensornetworkstate(ComplexF64, v -> "Up", g, "S=1/2")
    psi_mid = tensornetworkstate(ComplexF64, v -> v == (4, 3) ? "Dn" : "Up", g, "S=1/2")

    layer = build_trotter_layer(g, dt)

    # Setup BP caches
    psi_scv_bpc = BeliefPropagationCache(psi_scv)
    psi_mid_bpc = BeliefPropagationCache(psi_mid)

    apply_kwargs = (; maxdim = dbond, cutoff = CUTOFF, normalize_tensors = true)

    # Storage
    steps, flux_scv_vals, flux_mid_vals, flux_diff_vals, time_vals = Int[], Float64[], Float64[], Float64[], Float64[]

    # t = 0 Measurement
    step_time = @elapsed begin
        f_scv = measure_total_flux(psi_scv_bpc, g)
        f_mid = measure_total_flux(psi_mid_bpc, g)
    end

    push!(steps, 0); push!(flux_scv_vals, f_scv); push!(flux_mid_vals, f_mid)
    push!(flux_diff_vals, f_mid - f_scv); push!(time_vals, step_time)

    @printf("Step %2d | SCV: %5.2f | MID: %5.2f | Diff: %5.2f | Time: %6.2fs\n", 
            0, f_scv, f_mid, f_mid - f_scv, step_time)

    # Time evolution
    for step in 1:num_steps
        step_time = @elapsed begin
            psi_scv_bpc, err_scv = apply_gates(layer, psi_scv_bpc; apply_kwargs)
            psi_mid_bpc, err_mid = apply_gates(layer, psi_mid_bpc; apply_kwargs)

            f_scv = measure_total_flux(psi_scv_bpc, g)
            f_mid = measure_total_flux(psi_mid_bpc, g)
        end

        f_diff = f_mid - f_scv

        push!(steps, step); push!(flux_scv_vals, f_scv); push!(flux_mid_vals, f_mid)
        push!(flux_diff_vals, f_diff); push!(time_vals, step_time)

        max_err = max(maximum(err_scv), maximum(err_mid))
        max_dim = max(maxvirtualdim(psi_scv_bpc), maxvirtualdim(psi_mid_bpc))

        @printf("Step %2d | SCV: %5.2f | MID: %5.2f | Diff: %5.2f | Time: %6.2fs | MaxErr: %.2e | MaxDim: %d\n",
                step, f_scv, f_mid, f_diff, step_time, max_err, max_dim)
    end

    return steps, flux_scv_vals, flux_mid_vals, flux_diff_vals, time_vals
end

# Run sweeps
dbonds_to_run = [16, 18, 20, 22, 24, 26] 
out_dir = "C:\\Users\\dell\\Desktop"

for dbond in dbonds_to_run
    steps, flux_scv, flux_mid, flux_diff, time_vals = run_simulation(dbond)

    save_path = joinpath(out_dir, "total_flux_results_dbond_$(dbond)_k_$(K_coup)_gamma_$(Gamma).csv")
    open(save_path, "w") do f
        println(f, "step,scv_flux,mid_flux,diff_flux,time_seconds")
        for i in eachindex(steps)
            @printf(f, "%d,%.6f,%.6f,%.6f,%.6f\n",
                    steps[i], flux_scv[i], flux_mid[i], flux_diff[i], time_vals[i])
        end
    end
    println("-> Saved DBOND $dbond to $save_path")
end

println("\nAll DBOND simulations complete.")
