Base.@kwdef struct SimulationConfig
    nx::Int = 8
    ny::Int = 6
    K::Float64 = 3.0
    Gamma::Float64 = 1.0
    dt::Float64 = 0.05
    steps::Int = 20
    trotter::Symbol = :lie_x_zz
    boundary::Symbol = :open_dual
    backend::Symbol = :bp
    defect_x::Int = 0
    defect_y::Int = 0

    maxdim::Int = 16
    cutoff::Float64 = 1.0e-10
    normalize_tensors::Bool = true
    bp_tolerance::Float64 = 1.0e-8
    bp_maxiter::Int = 50
    bp_fail_on_nonconvergence::Bool = true

    measure_every::Int = 1
    radial_bin_width::Float64 = 0.5
    angular_bins::Int = 16
    front_quantile::Float64 = 0.90
    correlation_max_distance::Int = 4
    wilson_max_area::Int = 4

    output_root::String = "results"
    run_label::String = "z2_loop"
    run_id::String = ""
    overwrite::Bool = false
end

_section(d::AbstractDict, name::String) = get(d, name, Dict{String, Any}())
_value(d::AbstractDict, key::String, default) = get(d, key, default)

function _symbol_value(d::AbstractDict, key::String, default::Symbol)
    return Symbol(lowercase(String(_value(d, key, String(default)))))
end

"""Load a simulation configuration from the documented TOML schema."""
function load_config(path::AbstractString)
    raw = TOML.parsefile(path)
    sim = _section(raw, "simulation")
    peps = _section(raw, "peps")
    bp = _section(raw, "bp")
    obs = _section(raw, "observables")
    out = _section(raw, "output")

    cfg = SimulationConfig(
        nx = Int(_value(sim, "nx", 8)),
        ny = Int(_value(sim, "ny", 6)),
        K = Float64(_value(sim, "K", 3.0)),
        Gamma = Float64(_value(sim, "Gamma", 1.0)),
        dt = Float64(_value(sim, "dt", 0.05)),
        steps = Int(_value(sim, "steps", 20)),
        trotter = _symbol_value(sim, "trotter", :lie_x_zz),
        boundary = _symbol_value(sim, "boundary", :open_dual),
        backend = _symbol_value(sim, "backend", :bp),
        defect_x = Int(_value(sim, "defect_x", 0)),
        defect_y = Int(_value(sim, "defect_y", 0)),
        maxdim = Int(_value(peps, "maxdim", 16)),
        cutoff = Float64(_value(peps, "cutoff", 1.0e-10)),
        normalize_tensors = Bool(_value(peps, "normalize_tensors", true)),
        bp_tolerance = Float64(_value(bp, "tolerance", 1.0e-8)),
        bp_maxiter = Int(_value(bp, "maxiter", 50)),
        bp_fail_on_nonconvergence = Bool(_value(bp, "fail_on_nonconvergence", true)),
        measure_every = Int(_value(obs, "measure_every", 1)),
        radial_bin_width = Float64(_value(obs, "radial_bin_width", 0.5)),
        angular_bins = Int(_value(obs, "angular_bins", 16)),
        front_quantile = Float64(_value(obs, "front_quantile", 0.90)),
        correlation_max_distance = Int(_value(obs, "correlation_max_distance", 4)),
        wilson_max_area = Int(_value(obs, "wilson_max_area", 4)),
        output_root = String(_value(out, "root", "results")),
        run_label = String(_value(out, "label", "z2_loop")),
        run_id = String(_value(out, "run_id", "")),
        overwrite = Bool(_value(out, "overwrite", false)),
    )
    validate_config(cfg)
    return cfg
end

function resolved_defect(cfg::SimulationConfig)
    x = iszero(cfg.defect_x) ? cld(cfg.nx, 2) : cfg.defect_x
    y = iszero(cfg.defect_y) ? cld(cfg.ny, 2) : cfg.defect_y
    return (x, y)
end

function validate_config(cfg::SimulationConfig)
    cfg.nx >= 2 || throw(ArgumentError("nx must be at least 2"))
    cfg.ny >= 2 || throw(ArgumentError("ny must be at least 2"))
    cfg.dt > 0 || throw(ArgumentError("dt must be positive"))
    cfg.steps >= 0 || throw(ArgumentError("steps must be nonnegative"))
    cfg.K >= 0 || throw(ArgumentError("K must be nonnegative in the frozen convention"))
    cfg.Gamma >= 0 || throw(ArgumentError("Gamma must be nonnegative in the frozen convention"))
    cfg.maxdim >= 1 || throw(ArgumentError("maxdim must be positive"))
    cfg.cutoff >= 0 || throw(ArgumentError("cutoff cannot be negative"))
    cfg.bp_tolerance > 0 || throw(ArgumentError("BP tolerance must be positive"))
    cfg.bp_maxiter >= 1 || throw(ArgumentError("BP maxiter must be positive"))
    cfg.measure_every >= 1 || throw(ArgumentError("measure_every must be positive"))
    cfg.radial_bin_width > 0 || throw(ArgumentError("radial_bin_width must be positive"))
    cfg.angular_bins >= 4 || throw(ArgumentError("angular_bins must be at least 4"))
    0 < cfg.front_quantile <= 1 || throw(ArgumentError("front_quantile must lie in (0,1]"))
    cfg.correlation_max_distance >= 0 || throw(ArgumentError("correlation_max_distance cannot be negative"))
    cfg.wilson_max_area >= 1 || throw(ArgumentError("wilson_max_area must be positive"))
    cfg.trotter in (:lie_x_zz, :lie_zz_x, :strang) ||
        throw(ArgumentError("trotter must be lie_x_zz, lie_zz_x, or strang"))
    cfg.boundary in (:open_dual, :fixed_exterior) ||
        throw(ArgumentError("boundary must be open_dual or fixed_exterior"))
    cfg.backend in (:bp, :exact) || throw(ArgumentError("backend must be bp or exact"))
    d = resolved_defect(cfg)
    1 <= d[1] <= cfg.nx || throw(ArgumentError("defect_x lies outside the lattice"))
    1 <= d[2] <= cfg.ny || throw(ArgumentError("defect_y lies outside the lattice"))
    if cfg.backend == :exact && cfg.nx * cfg.ny > 22
        throw(ArgumentError("exact backend is capped at 22 sites; use a smaller validation lattice"))
    end
    return cfg
end

function config_dictionary(cfg::SimulationConfig)
    defect = resolved_defect(cfg)
    return Dict(
        "simulation" => Dict(
            "nx" => cfg.nx,
            "ny" => cfg.ny,
            "K" => cfg.K,
            "Gamma" => cfg.Gamma,
            "dt" => cfg.dt,
            "steps" => cfg.steps,
            "trotter" => String(cfg.trotter),
            "boundary" => String(cfg.boundary),
            "backend" => String(cfg.backend),
            "defect_x" => defect[1],
            "defect_y" => defect[2],
        ),
        "peps" => Dict(
            "maxdim" => cfg.maxdim,
            "cutoff" => cfg.cutoff,
            "normalize_tensors" => cfg.normalize_tensors,
        ),
        "bp" => Dict(
            "tolerance" => cfg.bp_tolerance,
            "maxiter" => cfg.bp_maxiter,
            "fail_on_nonconvergence" => cfg.bp_fail_on_nonconvergence,
        ),
        "observables" => Dict(
            "measure_every" => cfg.measure_every,
            "radial_bin_width" => cfg.radial_bin_width,
            "angular_bins" => cfg.angular_bins,
            "front_quantile" => cfg.front_quantile,
            "correlation_max_distance" => cfg.correlation_max_distance,
            "wilson_max_area" => cfg.wilson_max_area,
        ),
        "output" => Dict(
            "root" => cfg.output_root,
            "label" => cfg.run_label,
            "run_id" => cfg.run_id,
            "overwrite" => cfg.overwrite,
        ),
    )
end
