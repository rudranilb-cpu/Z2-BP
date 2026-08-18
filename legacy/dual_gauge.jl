#!/usr/bin/env julia
# Compatibility entry point for the former environment-free PEPS script.
using Pkg
const ROOT = @__DIR__
Pkg.activate(ROOT)
using Z2GaugeBenchmarks

config_path = isempty(ARGS) ? joinpath(ROOT, "configs", "smoke_exact.toml") : abspath(ARGS[1])
@warn "The old script was not an independent exact PEPS benchmark: it used environment-free simple update and BP contraction. This wrapper runs the controlled backend named in the supplied config."
run_simulation(load_config(config_path))
