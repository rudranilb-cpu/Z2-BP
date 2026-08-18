###############################################################################
#  Z2 Gauge Theory — Total Electric Flux (SCV, MID, and Difference)
###############################################################################

using ITensors
using ITensorMPS
using ITensorNetworks
using NamedGraphs
using NamedGraphs.NamedGraphGenerators: named_grid 
using Observers             
using LinearAlgebra: norm
using Graphs                
using Printf 

# Parameters 
Nx      = 8       
Ny      = 6          
K_coup  = 3.0      
Gamma   = 1.0        
DBOND   = 8    
dt      = 0.05       
num_steps = 20      
CUTOFF  = 1e-10      

# Lattice Setup 
function build_lattice(Nx, Ny)
    g     = named_grid((Nx, Ny))          
    sites = siteinds("S=1/2", g)          
    return g, sites
end

# Trotter Gates
function make_tfim_gates(g, sites, dt)
    gates = ITensor[]
    for e in edges(g)
        s1, s2 = sites[src(e)][1], sites[dst(e)][1]
        hj = -Gamma * op("Z", s1) * op("Z", s2)
        push!(gates, exp(-im * dt * hj))
    end
    for v in vertices(g)
        sv = sites[v][1]
        hj = -K_coup * op("X", sv)
        push!(gates, exp(-im * dt * hj))
    end
    return gates
end

# Measure TOTAL Flux 
function measure_total_flux(psi, g, sites)
    total_flux = 0.0
    
    norm_psi = inner(psi, psi; cache_update_kwargs=(; maxiter=100))
    
    # Sum the flux (1 - ZZ)/2 
    for e in edges(g)
        s1, s2 = sites[src(e)][1], sites[dst(e)][1]
        ZZ_op = op("Z", s1) * op("Z", s2)
        ZZ_ket = apply([ZZ_op], psi)
        
        zz_inner = inner(psi, ZZ_ket; cache_update_kwargs=(; maxiter=100))
        zz_val = real(zz_inner / norm_psi)
        
        total_flux += 0.5 * (1.0 - zz_val)
    end
    return total_flux
end

# Main Simulation 
function run_simulation()
    println("="^60)
    println("  Calculating Total Electric Flux (SCV vs MID)")
    println("  Lattice: $(Nx)x$(Ny) | K: $K_coup | Gamma: $Gamma | DBOND: $DBOND")
    println("="^60)

    g, sites = build_lattice(Nx, Ny)
    
    # Initialize States 
    psi_scv = ITensorNetwork(v -> "Up", sites)
    psi_mid = ITensorNetwork(v -> v == (4, 3) ? "Dn" : "Up", sites)
    
    gates = make_tfim_gates(g, sites, dt)

    # Storage
    steps = Int[]
    flux_scv_vals = Float64[]
    flux_mid_vals = Float64[]
    flux_diff_vals = Float64[]

    # t = 0 Measurement
    println("Measuring t=0... (This takes a moment)")
    f_scv = measure_total_flux(psi_scv, g, sites)
    f_mid = measure_total_flux(psi_mid, g, sites)
    
    push!(steps, 0)
    push!(flux_scv_vals, f_scv)
    push!(flux_mid_vals, f_mid)
    push!(flux_diff_vals, f_mid - f_scv)
    
    @printf("  Step %2d | SCV: %5.2f | MID: %5.2f | Diff: %5.2f\n", 0, f_scv, f_mid, f_mid - f_scv)

    println("\n--- Evolving ---")
    for step in 1:num_steps
        # Evolve
        psi_scv = apply(gates, psi_scv; maxdim=DBOND, cutoff=CUTOFF)
        psi_mid = apply(gates, psi_mid; maxdim=DBOND, cutoff=CUTOFF)

        # Measure
        f_scv = measure_total_flux(psi_scv, g, sites)
        f_mid = measure_total_flux(psi_mid, g, sites)
        f_diff = f_mid - f_scv
        
        push!(steps, step)
        push!(flux_scv_vals, f_scv)
        push!(flux_mid_vals, f_mid)
        push!(flux_diff_vals, f_diff)

        @printf("  Step %2d | SCV: %5.2f | MID: %5.2f | Diff: %5.2f\n", step, f_scv, f_mid, f_diff)
    end

    return steps, flux_scv_vals, flux_mid_vals, flux_diff_vals
end

# Execution & Export 
steps, flux_scv, flux_mid, flux_diff = run_simulation()


save_path = "total_flux_results_dbond_$(DBOND).csv"
open(save_path, "w") do f
    println(f, "step,scv_flux,mid_flux,diff_flux")
    for i in eachindex(steps)
        @printf(f, "%d,%.6f,%.6f,%.6f\n", steps[i], flux_scv[i], flux_mid[i], flux_diff[i])
    end
end
println("\n Data saved to: ", save_path)
