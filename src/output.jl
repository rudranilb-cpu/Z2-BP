mutable struct RunWriter
    run_id::String
    directory::String
    global_io::IO
    links_io::IO
    sites_io::IO
    profiles_io::IO
    wilson_io::IO
    correlations_io::IO
    diagnostics_io::IO
    bp_updates_io::IO
    pauli_io::IO
end

function _git_revision()
    try
        return readchomp(`git -C $PACKAGE_ROOT rev-parse HEAD`)
    catch
        return "unknown"
    end
end

function _git_dirty()
    try
        return !isempty(readchomp(`git -C $PACKAGE_ROOT status --porcelain`))
    catch
        return "unknown"
    end
end

function _default_run_id(cfg::SimulationConfig)
    timestamp = Dates.format(now(UTC), dateformat"yyyymmddTHHMMSS")
    ratio = iszero(cfg.Gamma) ? "inf" : @sprintf("%.6g", cfg.K / cfg.Gamma)
    return "$(cfg.run_label)_$(cfg.backend)_$(cfg.nx)x$(cfg.ny)_r$(ratio)_D$(cfg.maxdim)_$(timestamp)_$(string(uuid4())[1:8])"
end

function _csv_escape(value)
    value === missing && return ""
    value isa AbstractFloat && isnan(value) && return "nan"
    text = string(value)
    if occursin(',', text) || occursin('"', text) || occursin('\n', text)
        return "\"" * replace(text, "\"" => "\"\"") * "\""
    end
    return text
end

function _write_row(io::IO, values...)
    println(io, join(_csv_escape.(values), ','))
end

function _prepare_run_directory(cfg::SimulationConfig)
    run_id = isempty(cfg.run_id) ? _default_run_id(cfg) : cfg.run_id
    directory = abspath(joinpath(cfg.output_root, run_id))
    if isdir(directory) && !cfg.overwrite
        throw(ArgumentError("run directory already exists: $directory"))
    end
    mkpath(directory)
    return run_id, directory
end

function _metadata(cfg::SimulationConfig, run_id::String, links::Vector{DualLink}, plan::MeasurementPlan)
    defect = resolved_defect(cfg)
    expected_initial_flux = length(star_links(links, defect))
    metadata = config_dictionary(cfg)
    metadata["schema"] = Dict(
        "name" => "z2-gauge-benchmark",
        "version" => SCHEMA_VERSION,
        "run_id" => run_id,
        "created_utc" => string(now(UTC)),
        "git_revision" => _git_revision(),
        "git_dirty" => _git_dirty(),
    )
    metadata["software"] = Dict(
        "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
        "cpu_threads" => Sys.CPU_THREADS,
        "architecture" => string(Sys.ARCH),
        "kernel" => string(Sys.KERNEL),
        "word_size" => Sys.WORD_SIZE,
        "tnqs_version" => string(Base.pkgversion(TNQS)),
        "tnqs_commit" => TNQS_COMMIT,
    )
    metadata["conventions"] = Dict(
        "hamiltonian" => "H=-K*sum(X)-Gamma*sum(ZZ) with optional fixed-exterior boundary Z fields",
        "duality" => "B_p=X_r; sigma_x(link)=Z_r*Z_rprime on internal links",
        "link_flux" => "n_b=(1-ZZ)/2; boundary n_b=(1-Z)/2 only for fixed_exterior",
        "differential_sign" => "loop_minus_scv",
        "lie_reference" => "sequential X then diagonal gates, U_step=U_ZZ*U_X",
        "qiskit_angles" => "Rx(-2*K*dt), Rzz(-2*Gamma*dt)",
        "initial_loop" => "one flipped dual Z spin",
        "expected_delta_phi_t0" => expected_initial_flux,
        "internal_bonds" => internal_bond_count(cfg.nx, cfg.ny),
        "measured_links" => length(links),
        "raw_pauli_observables" => length(plan.keys),
        "qpu_measurement_settings" => 2,
        "boundary_distance_sites" => distance_to_open_boundary(defect, cfg),
    )
    metadata["error_channels"] = Dict(
        "statistical" => "stderr columns; NaN for deterministic BP/exact runs",
        "physical_finite_size" => "compare relative-coordinate central observables across sizes in a pre-boundary window",
        "bond_dimension" => "compare separately across maxdim; not inferred from one run",
        "bp_environment" => "bp_updates.csv residual and iteration convergence",
        "trotter" => "compare dt and lie/strang runs separately",
        "local_truncation" => "diagnostics.csv max/sum discarded weight; not a global error bar",
    )
    return metadata
end

function open_run_writer(cfg::SimulationConfig, links::Vector{DualLink}, plan::MeasurementPlan)
    run_id, directory = _prepare_run_directory(cfg)
    open(joinpath(directory, "metadata.toml"), "w") do io
        TOML.print(io, _metadata(cfg, run_id, links, plan); sorted = true)
    end
    global_io = open(joinpath(directory, "global.csv"), "w")
    links_io = open(joinpath(directory, "links.csv"), "w")
    sites_io = open(joinpath(directory, "sites.csv"), "w")
    profiles_io = open(joinpath(directory, "profiles.csv"), "w")
    wilson_io = open(joinpath(directory, "wilson.csv"), "w")
    correlations_io = open(joinpath(directory, "correlations.csv"), "w")
    diagnostics_io = open(joinpath(directory, "diagnostics.csv"), "w")
    bp_updates_io = open(joinpath(directory, "bp_updates.csv"), "w")
    pauli_io = open(joinpath(directory, "pauli.csv"), "w")

    _write_row(global_io,
        "run_id", "backend", "step", "time", "phi_scv", "phi_loop", "delta_phi",
        "phi_scv_stderr", "phi_loop_stderr", "delta_phi_stderr",
        "flux_density_scv", "flux_density_loop", "delta_flux_density",
        "flux_density_scv_stderr", "flux_density_loop_stderr", "delta_flux_density_stderr",
        "loop_survival_scv", "loop_survival_loop", "delta_loop_survival",
        "loop_survival_scv_stderr", "loop_survival_loop_stderr", "delta_loop_survival_stderr",
        "central_star_delta_flux", "central_star_delta_flux_stderr",
        "front_mean", "front_rms", "front_quantile", "front_weight",
        "front_participation", "front_anisotropy",
        "front_mean_stderr", "front_rms_stderr", "front_quantile_stderr",
        "front_weight_stderr", "front_participation_stderr", "front_anisotropy_stderr")
    _write_row(links_io,
        "run_id", "backend", "step", "time", "link_id", "orientation",
        "v1_x", "v1_y", "v2_x", "v2_y", "mid_x", "mid_y", "radius", "theta",
        "n_scv", "n_loop", "delta_n", "n_scv_stderr", "n_loop_stderr", "delta_n_stderr")
    _write_row(sites_io,
        "run_id", "backend", "step", "time", "x", "y", "radius", "theta",
        "X_scv", "X_loop", "delta_X", "Z_scv", "Z_loop", "delta_Z",
        "local_flux_scv", "local_flux_loop", "delta_local_flux",
        "X_scv_stderr", "X_loop_stderr", "delta_X_stderr",
        "Z_scv_stderr", "Z_loop_stderr", "delta_Z_stderr",
        "local_flux_scv_stderr", "local_flux_loop_stderr", "delta_local_flux_stderr")
    _write_row(profiles_io,
        "run_id", "backend", "step", "time", "profile_kind", "bin_index", "bin_lo", "bin_hi",
        "bin_center", "count", "n_scv_mean", "n_loop_mean", "delta_n_mean", "abs_delta_n_mean",
        "n_scv_mean_stderr", "n_loop_mean_stderr", "delta_n_mean_stderr", "abs_delta_n_mean_stderr")
    _write_row(wilson_io,
        "run_id", "backend", "step", "time", "observable_id", "x0", "y0", "width", "height",
        "value_scv", "value_loop", "delta_value", "stderr_scv", "stderr_loop", "stderr_delta")
    _write_row(correlations_io,
        "run_id", "backend", "step", "time", "kind", "observable_id", "anchor_id", "target_id",
        "distance", "value_scv", "value_loop", "delta_value", "stderr_scv", "stderr_loop", "stderr_delta")
    _write_row(diagnostics_io,
        "run_id", "backend", "step", "time", "state", "evolution_seconds", "measurement_seconds",
        "max_bond_dimension", "truncation_max", "truncation_sum", "truncation_count",
        "bp_update_count", "bp_iterations", "bp_final_residual_max", "bp_all_converged",
        "norm_drift", "physical_link_violations", "live_memory_mb", "maxrss_mb")
    _write_row(bp_updates_io,
        "run_id", "backend", "step", "time", "state", "update_id", "iterations", "final_residual", "converged")
    _write_row(pauli_io,
        "run_id", "backend", "step", "time", "observable_id", "basis", "qpu_setting",
        "pauli", "support", "value_scv", "value_loop", "delta_value",
        "stderr_scv", "stderr_loop", "stderr_delta")
    return RunWriter(
        run_id,
        directory,
        global_io,
        links_io,
        sites_io,
        profiles_io,
        wilson_io,
        correlations_io,
        diagnostics_io,
        bp_updates_io,
        pauli_io,
    )
end

function Base.close(writer::RunWriter)
    for io in (
        writer.global_io,
        writer.links_io,
        writer.sites_io,
        writer.profiles_io,
        writer.wilson_io,
        writer.correlations_io,
        writer.diagnostics_io,
        writer.bp_updates_io,
        writer.pauli_io,
    )
        Base.close(io)
    end
    return nothing
end

function _weighted_quantile(radii::Vector{Float64}, weights::Vector{Float64}, q::Float64)
    total = sum(weights)
    total <= eps(Float64) && return NaN
    order = sortperm(radii)
    threshold = q * total
    cumulative = 0.0
    for index in order
        cumulative += weights[index]
        cumulative >= threshold && return radii[index]
    end
    return radii[order[end]]
end

function _front_metrics(cfg::SimulationConfig, links, scv::StateSnapshot, loop::StateSnapshot)
    defect = resolved_defect(cfg)
    radii = Float64[first(link_radius_angle(link, defect)) for link in links]
    dx = Float64[link.x - defect[1] for link in links]
    dy = Float64[link.y - defect[2] for link in links]
    weights = abs.(loop.link_flux .- scv.link_flux)
    total = sum(weights)
    if total <= eps(Float64)
        return (mean = NaN, rms = NaN, quantile = NaN, weight = total,
            participation = NaN, anisotropy = NaN)
    end
    mean_radius = sum(weights .* radii) / total
    rms = sqrt(sum(weights .* radii .^ 2) / total)
    quantile = _weighted_quantile(radii, weights, cfg.front_quantile)
    weight_square = sum(abs2, weights)
    participation = iszero(weight_square) ? NaN : total^2 / weight_square
    mxx = sum(weights .* dx .^ 2) / total
    myy = sum(weights .* dy .^ 2) / total
    mxy = sum(weights .* dx .* dy) / total
    trace = mxx + myy
    discriminant = hypot(mxx - myy, 2mxy)
    anisotropy = trace <= eps(Float64) ? NaN : discriminant / trace
    return (mean = mean_radius, rms = rms, quantile = quantile, weight = total,
        participation = participation, anisotropy = anisotropy)
end

function _local_flux(links, values::Vector{Float64})
    density = Dict{Vertex, Float64}()
    for (link, value) in zip(links, values)
        if isnothing(link.v2)
            density[link.v1] = get(density, link.v1, 0.0) + value
        else
            density[link.v1] = get(density, link.v1, 0.0) + value / 2
            density[link.v2] = get(density, link.v2, 0.0) + value / 2
        end
    end
    return density
end

function _write_profiles!(writer, cfg, step, time, links, scv, loop)
    defect = resolved_defect(cfg)
    radial = Dict{Int, Vector{Int}}()
    angular = Dict{Int, Vector{Int}}()
    angular_width = 2pi / cfg.angular_bins
    for (index, link) in enumerate(links)
        radius, theta = link_radius_angle(link, defect)
        rbin = floor(Int, radius / cfg.radial_bin_width)
        abin = clamp(floor(Int, (theta + pi) / angular_width), 0, cfg.angular_bins - 1)
        push!(get!(radial, rbin, Int[]), index)
        push!(get!(angular, abin, Int[]), index)
    end
    delta = loop.link_flux .- scv.link_flux
    for (kind, bins, width, origin) in (("radial", radial, cfg.radial_bin_width, 0.0), ("angular", angular, angular_width, -pi))
        for bin in sort(collect(keys(bins)))
            indices = bins[bin]
            lo, hi = origin + bin * width, origin + (bin + 1) * width
            _write_row(writer.profiles_io,
                writer.run_id, cfg.backend, step, time, kind, bin, lo, hi, (lo + hi) / 2,
                length(indices), mean(scv.link_flux[indices]), mean(loop.link_flux[indices]),
                mean(delta[indices]), mean(abs.(delta[indices])), NaN, NaN, NaN, NaN)
        end
    end
end

function write_observables!(writer::RunWriter, cfg::SimulationConfig, plan, links, step, scv, loop)
    time = step * cfg.dt
    phi_scv, phi_loop = sum(scv.link_flux), sum(loop.link_flux)
    delta_phi = phi_loop - phi_scv
    front = _front_metrics(cfg, links, scv, loop)
    central_ids = [link.id for link in star_links(links, resolved_defect(cfg))]
    central_star_delta = sum((loop.link_flux .- scv.link_flux)[central_ids])
    nlinks = length(links)
    _write_row(writer.global_io,
        writer.run_id, cfg.backend, step, time, phi_scv, phi_loop, delta_phi,
        NaN, NaN, NaN, phi_scv / nlinks, phi_loop / nlinks, delta_phi / nlinks,
        NaN, NaN, NaN,
        scv.loop_survival, loop.loop_survival, loop.loop_survival - scv.loop_survival,
        NaN, NaN, NaN, central_star_delta, NaN,
        front.mean, front.rms, front.quantile, front.weight, front.participation, front.anisotropy,
        NaN, NaN, NaN, NaN, NaN, NaN)

    defect = resolved_defect(cfg)
    for (index, link) in enumerate(links)
        radius, theta = link_radius_angle(link, defect)
        v2x = isnothing(link.v2) ? missing : link.v2[1]
        v2y = isnothing(link.v2) ? missing : link.v2[2]
        _write_row(writer.links_io,
            writer.run_id, cfg.backend, step, time, link.id, link.orientation,
            link.v1[1], link.v1[2], v2x, v2y, link.x, link.y, radius, theta,
            scv.link_flux[index], loop.link_flux[index], loop.link_flux[index] - scv.link_flux[index],
            NaN, NaN, NaN)
    end

    local_scv = _local_flux(links, scv.link_flux)
    local_loop = _local_flux(links, loop.link_flux)
    for v in dual_vertices(cfg)
        dx, dy = v[1] - defect[1], v[2] - defect[2]
        _write_row(writer.sites_io,
            writer.run_id, cfg.backend, step, time, v[1], v[2], hypot(dx, dy), atan(dy, dx),
            scv.site_x[v], loop.site_x[v], loop.site_x[v] - scv.site_x[v],
            scv.site_z[v], loop.site_z[v], loop.site_z[v] - scv.site_z[v],
            get(local_scv, v, 0.0), get(local_loop, v, 0.0), get(local_loop, v, 0.0) - get(local_scv, v, 0.0),
            NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
    end
    _write_profiles!(writer, cfg, step, time, links, scv, loop)

    for (index, key) in enumerate(plan.keys)
        observable_id = @sprintf("P%05d", index)
        pauli = repeat(string(key.axis), length(key.vertices))
        support = join(("$(v[1]):$(v[2])" for v in key.vertices), ';')
        a, b = scv.pauli[key], loop.pauli[key]
        _write_row(writer.pauli_io,
            writer.run_id, cfg.backend, step, time, observable_id, key.axis,
            key.axis == 'Z' ? "global_Z" : "global_X", pauli, support,
            a, b, b - a, NaN, NaN, NaN)
    end

    for spec in plan.wilson_specs
        a, b = scv.wilson[spec.id], loop.wilson[spec.id]
        _write_row(writer.wilson_io,
            writer.run_id, cfg.backend, step, time, spec.id, spec.x0, spec.y0, spec.width, spec.height,
            a, b, b - a, NaN, NaN, NaN)
    end
    for spec in plan.string_specs
        a, b = scv.strings[spec.id], loop.strings[spec.id]
        _write_row(writer.correlations_io,
            writer.run_id, cfg.backend, step, time, "electric_string", spec.id, "defect", "$(spec.target[1]):$(spec.target[2])",
            spec.distance, a, b, b - a, NaN, NaN, NaN)
    end
    for spec in plan.magnetic_specs
        a, b = scv.magnetic_connected[spec.id], loop.magnetic_connected[spec.id]
        _write_row(writer.correlations_io,
            writer.run_id, cfg.backend, step, time, "magnetic_connected", spec.id, "defect", "$(spec.target[1]):$(spec.target[2])",
            spec.distance, a, b, b - a, NaN, NaN, NaN)
    end
    for spec in plan.flux_specs
        key = (spec.anchor_link_id, spec.target_link_id)
        a, b = scv.flux_connected[key], loop.flux_connected[key]
        _write_row(writer.correlations_io,
            writer.run_id, cfg.backend, step, time, "flux_connected", "C_n_$(key[1])_$(key[2])", key[1], key[2], NaN,
            a, b, b - a, NaN, NaN, NaN)
    end
    return nothing
end

function _memory_mb()
    live = try
        Base.gc_live_bytes() / 2.0^20
    catch
        NaN
    end
    rss = try
        Sys.maxrss() / 2.0^20
    catch
        NaN
    end
    return live, rss
end

function write_diagnostics!(writer, cfg, step, state, evolution_seconds, measurement_seconds;
        maxdim = missing, errors = Float64[], bp_stats = BPUpdateStat[], norm_drift = NaN,
        physical_violations = 0)
    time = step * cfg.dt
    trunc_max = isempty(errors) ? 0.0 : maximum(errors)
    trunc_sum = sum(errors)
    bp_iterations = sum(stat.iterations for stat in bp_stats; init = 0)
    residual_max = isempty(bp_stats) ? NaN : maximum(stat.final_residual for stat in bp_stats)
    all_converged = isempty(bp_stats) ? true : all(stat.converged for stat in bp_stats)
    live, rss = _memory_mb()
    _write_row(writer.diagnostics_io,
        writer.run_id, cfg.backend, step, time, state, evolution_seconds, measurement_seconds,
        maxdim, trunc_max, trunc_sum, length(errors), length(bp_stats), bp_iterations,
        residual_max, all_converged, norm_drift, physical_violations, live, rss)
    for (update_id, stat) in enumerate(bp_stats)
        _write_row(writer.bp_updates_io,
            writer.run_id, cfg.backend, step, time, state, update_id, stat.iterations, stat.final_residual, stat.converged)
    end
    return nothing
end

function inspect_run(directory::AbstractString)
    metadata = TOML.parsefile(joinpath(directory, "metadata.toml"))
    return (
        run_id = metadata["schema"]["run_id"],
        schema_version = metadata["schema"]["version"],
        backend = metadata["simulation"]["backend"],
        nx = metadata["simulation"]["nx"],
        ny = metadata["simulation"]["ny"],
    )
end
