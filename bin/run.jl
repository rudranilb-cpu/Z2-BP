#!/usr/bin/env julia
using Pkg

const ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(ROOT)
using Z2GaugeBenchmarks

length(ARGS) == 1 || error("usage: julia --project=. bin/run.jl CONFIG.toml")
config = load_config(abspath(ARGS[1]))
run_simulation(config)
