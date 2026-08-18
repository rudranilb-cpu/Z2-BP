function _physical_violations(snapshot::StateSnapshot; tolerance = 1.0e-8)
    return count(value -> value < -tolerance || value > 1.0 + tolerance, snapshot.link_flux)
end

function _check_t0(cfg, links, scv, loop)
    expected = length(star_links(links, resolved_defect(cfg)))
    phi_scv = sum(scv.link_flux)
    delta = sum(loop.link_flux) - phi_scv
    atol = cfg.backend == :exact ? 1.0e-12 : max(1.0e-8, 10 * cfg.bp_tolerance)
    isapprox(phi_scv, 0.0; atol) || error("t=0 SCV flux check failed: $phi_scv")
    isapprox(delta, expected; atol) || error("t=0 differential flux check failed: got $delta, expected $expected")
    return nothing
end

function _require_bp_convergence(cfg::SimulationConfig, stats, stage::AbstractString)
    cfg.bp_fail_on_nonconvergence || return nothing
    failed = [stat for stat in stats if !stat.converged]
    isempty(failed) && return nothing
    residual = maximum(stat.final_residual for stat in failed)
    error("BP did not converge during $stage: $(length(failed)) update(s) failed; " *
          "largest final residual=$residual, tolerance=$(cfg.bp_tolerance)")
end

function _measure_pair_exact(pair, cfg, plan, links)
    start = time_ns()
    scv = measure_exact_state(pair.scv, cfg, plan, links)
    scv_time = (time_ns() - start) / 1.0e9
    start = time_ns()
    loop = measure_exact_state(pair.loop, cfg, plan, links)
    loop_time = (time_ns() - start) / 1.0e9
    return scv, loop, (scv = scv_time, loop = loop_time)
end

function _measure_pair_bp(pair, plan, links)
    start = time_ns()
    scv = measure_bp_state(pair.scv, plan, links)
    scv_time = (time_ns() - start) / 1.0e9
    start = time_ns()
    loop = measure_bp_state(pair.loop, plan, links)
    loop_time = (time_ns() - start) / 1.0e9
    return scv, loop, (scv = scv_time, loop = loop_time)
end

function run_exact(cfg::SimulationConfig)
    cfg.backend == :exact || throw(ArgumentError("run_exact requires backend=exact"))
    return _run(cfg)
end

function run_bp(cfg::SimulationConfig)
    cfg.backend == :bp || throw(ArgumentError("run_bp requires backend=bp"))
    return _run(cfg)
end

run_simulation(cfg::SimulationConfig) = _run(validate_config(cfg))

function _run(cfg::SimulationConfig)
    links = dual_links(cfg)
    plan = build_measurement_plan(cfg, links)
    writer = open_run_writer(cfg, links, plan)
    println("Run $(writer.run_id)")
    println("Output: $(writer.directory)")
    println("Convention: H=-K sum(X)-Gamma sum(ZZ), deltaPhi=loop-SCV, trotter=$(cfg.trotter)")

    try
        if cfg.backend == :exact
            pair = initialize_exact(cfg)
            scv, loop, measurement_time = _measure_pair_exact(pair, cfg, plan, links)
            _check_t0(cfg, links, scv, loop)
            write_observables!(writer, cfg, plan, links, 0, scv, loop)
            write_diagnostics!(writer, cfg, 0, "scv", 0.0, measurement_time.scv;
                norm_drift = exact_norm_drift(pair.scv), physical_violations = _physical_violations(scv))
            write_diagnostics!(writer, cfg, 0, "loop", 0.0, measurement_time.loop;
                norm_drift = exact_norm_drift(pair.loop), physical_violations = _physical_violations(loop))

            for step in 1:cfg.steps
                pair, runtime = evolve_exact!(pair, cfg)
                should_measure = iszero(step % cfg.measure_every) || step == cfg.steps
                if should_measure
                    scv, loop, measurement_time = _measure_pair_exact(pair, cfg, plan, links)
                    write_observables!(writer, cfg, plan, links, step, scv, loop)
                else
                    measurement_time = (scv = 0.0, loop = 0.0)
                    scv = loop = nothing
                end
                write_diagnostics!(writer, cfg, step, "scv", runtime / 2, measurement_time.scv;
                    norm_drift = exact_norm_drift(pair.scv),
                    physical_violations = isnothing(scv) ? 0 : _physical_violations(scv))
                write_diagnostics!(writer, cfg, step, "loop", runtime / 2, measurement_time.loop;
                    norm_drift = exact_norm_drift(pair.loop),
                    physical_violations = isnothing(loop) ? 0 : _physical_violations(loop))
                @printf("step %d/%d t=%.6g exact %.3fs\n", step, cfg.steps, step * cfg.dt, runtime)
            end
        else
            pair, initial_bp = initialize_bp(cfg)
            circuit = build_trotter_layer(cfg, pair.graph)
            scv, loop, measurement_time = _measure_pair_bp(pair, plan, links)
            _check_t0(cfg, links, scv, loop)
            write_observables!(writer, cfg, plan, links, 0, scv, loop)
            write_diagnostics!(writer, cfg, 0, "scv", 0.0, measurement_time.scv;
                maxdim = bp_maxdim(pair.scv), bp_stats = initial_bp.scv, norm_drift = bp_norm_drift(pair.scv),
                physical_violations = _physical_violations(scv))
            write_diagnostics!(writer, cfg, 0, "loop", 0.0, measurement_time.loop;
                maxdim = bp_maxdim(pair.loop), bp_stats = initial_bp.loop, norm_drift = bp_norm_drift(pair.loop),
                physical_violations = _physical_violations(loop))
            _require_bp_convergence(cfg, vcat(initial_bp.scv, initial_bp.loop), "initialization")

            for step in 1:cfg.steps
                pair, evolution = evolve_bp!(pair, circuit, cfg)
                should_measure = iszero(step % cfg.measure_every) || step == cfg.steps
                if should_measure
                    scv, loop, measurement_time = _measure_pair_bp(pair, plan, links)
                    write_observables!(writer, cfg, plan, links, step, scv, loop)
                else
                    measurement_time = (scv = 0.0, loop = 0.0)
                    scv = loop = nothing
                end
                write_diagnostics!(writer, cfg, step, "scv", evolution.scv.runtime, measurement_time.scv;
                    maxdim = bp_maxdim(pair.scv), errors = evolution.scv.errors, bp_stats = evolution.scv.bp,
                    norm_drift = bp_norm_drift(pair.scv),
                    physical_violations = isnothing(scv) ? 0 : _physical_violations(scv))
                write_diagnostics!(writer, cfg, step, "loop", evolution.loop.runtime, measurement_time.loop;
                    maxdim = bp_maxdim(pair.loop), errors = evolution.loop.errors, bp_stats = evolution.loop.bp,
                    norm_drift = bp_norm_drift(pair.loop),
                    physical_violations = isnothing(loop) ? 0 : _physical_violations(loop))
                _require_bp_convergence(cfg, vcat(evolution.scv.bp, evolution.loop.bp), "step $step")
                @printf("step %d/%d t=%.6g BP %.3fs D=%d residual<=%.2e: %s\n",
                    step, cfg.steps, step * cfg.dt, evolution.scv.runtime + evolution.loop.runtime,
                    max(bp_maxdim(pair.scv), bp_maxdim(pair.loop)), cfg.bp_tolerance,
                    all(vcat(evolution.scv.bp, evolution.loop.bp)) do stat
                        stat.converged
                    end)
            end
        end
    finally
        close(writer)
    end
    println("Completed: $(writer.directory)")
    return writer.directory
end
