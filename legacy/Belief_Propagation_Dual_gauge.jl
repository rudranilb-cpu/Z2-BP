#!/usr/bin/env julia
# Compatibility entry point for the former ad-hoc BP script.
using Pkg
const ROOT = @__DIR__
Pkg.activate(ROOT)
using Z2GaugeBenchmarks

config_path = isempty(ARGS) ? joinpath(ROOT, "configs", "reproduce_8x6_bp.toml") : abspath(ARGS[1])
@warn "This compatibility entry point is deprecated; prefer `bin/run.jl CONFIG.toml`."
run_simulation(load_config(config_path))
