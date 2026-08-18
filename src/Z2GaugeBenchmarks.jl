module Z2GaugeBenchmarks

using Dates
using Graphs
using LinearAlgebra
using NamedGraphs
using Printf
using Statistics
using TOML
using TensorNetworkQuantumSimulator
using UUIDs

const TNQS = TensorNetworkQuantumSimulator
const TNQS_COMMIT = "b5d4089849de1cc23806aa8325e8db56a55f2e0b"
const SCHEMA_VERSION = "1.1.0"
const PACKAGE_ROOT = normpath(joinpath(@__DIR__, ".."))

include("config.jl")
include("lattice.jl")
include("circuits.jl")
include("measurement_plan.jl")
include("exact_backend.jl")
include("bp_backend.jl")
include("output.jl")
include("simulation.jl")

export
    SimulationConfig,
    load_config,
    validate_config,
    resolved_defect,
    DualLink,
    dual_links,
    internal_bond_count,
    full_open_link_count,
    build_trotter_layer,
    run_simulation,
    run_exact,
    run_bp,
    inspect_run,
    SCHEMA_VERSION

end
