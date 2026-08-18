#!/usr/bin/env julia
using Pkg

const ROOT = normpath(joinpath(@__DIR__, ".."))
const TNQS_COMMIT = "b5d4089849de1cc23806aa8325e8db56a55f2e0b"
const TNQS_URL = "https://github.com/JoeyT1994/TensorNetworkQuantumSimulator.jl.git"

Pkg.activate(ROOT)
Pkg.add(PackageSpec(url = TNQS_URL, rev = TNQS_COMMIT))
Pkg.instantiate()
Pkg.precompile()

println("Environment ready at $ROOT")
println("TensorNetworkQuantumSimulator pinned to $TNQS_COMMIT (package version 0.4.4).")
