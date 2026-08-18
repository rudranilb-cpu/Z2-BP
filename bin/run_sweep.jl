#!/usr/bin/env julia
using Pkg

const ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(ROOT)
using Z2GaugeBenchmarks

isempty(ARGS) && error("usage: julia --project=. bin/run_sweep.jl CONFIG1.toml [CONFIG2.toml ...]")
failures = Pair{String, Any}[]
for path in ARGS
    println("\n=== $path ===")
    try
        run_simulation(load_config(abspath(path)))
    catch exception
        push!(failures, path => (exception, catch_backtrace()))
        showerror(stderr, exception, catch_backtrace())
        println(stderr)
    end
end
isempty(failures) || error("$(length(failures)) sweep runs failed")
