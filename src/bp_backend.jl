struct BPUpdateStat
    iterations::Int
    final_residual::Float64
    converged::Bool
end

mutable struct BPPair
    graph
    scv
    loop
end

function _message_distance(a, b)
    na, nb = norm(a), norm(b)
    (iszero(na) || iszero(nb)) && return Inf
    fidelity = abs2(dot(a, b) / (na * nb))
    return max(0.0, 1.0 - real(fidelity))
end

function _average_message_change(previous, current)
    current_keys = collect(keys(TNQS.messages(current)))
    isempty(current_keys) && return Inf
    residual = 0.0
    for edge in current_keys
        residual += _message_distance(TNQS.message(previous, edge), TNQS.message(current, edge))
    end
    return residual / length(current_keys)
end

"""
Warm-start BP and explicitly return the average projective message change. TNQS
0.4.4 does not expose its iteration count/residual, so this performs one public
`update` iteration at a time and computes the same message-distance definition
used internally. This is deliberately version-pinned and covered by tests.
"""
function converge_bp(cache, cfg::SimulationConfig)
    residual = Inf
    for iteration in 1:cfg.bp_maxiter
        next_cache = TNQS.update(cache; maxiter = 1, tolerance = nothing, verbose = false)
        residual = _average_message_change(cache, next_cache)
        cache = next_cache
        if residual <= cfg.bp_tolerance
            return cache, BPUpdateStat(iteration, residual, true)
        end
    end
    return cache, BPUpdateStat(cfg.bp_maxiter, residual, false)
end

function initialize_bp(cfg::SimulationConfig)
    g = named_grid((cfg.nx, cfg.ny))
    defect = resolved_defect(cfg)
    scv_tn = tensornetworkstate(ComplexF64, _ -> "Up", g, "S=1/2")
    loop_tn = tensornetworkstate(ComplexF64, v -> v == defect ? "Dn" : "Up", g, "S=1/2")
    scv, scv_stat = converge_bp(BeliefPropagationCache(scv_tn), cfg)
    loop, loop_stat = converge_bp(BeliefPropagationCache(loop_tn), cfg)
    return BPPair(g, scv, loop), (scv = [scv_stat], loop = [loop_stat])
end

function _apply_layer_controlled(cache, circuit::Vector, cfg::SimulationConfig)
    graph_cache = TNQS.graph(cache)
    converted = TNQS.toitensor(circuit, graph_cache, TNQS.siteinds(TNQS.network(cache)))
    gate_vertices = [gate[2] for gate in converted]
    gates = [gate[1] for gate in converted]
    cache = copy(cache)
    affected_vertices = Set{Any}()
    truncation_errors = zeros(Float64, length(gates))
    bp_stats = BPUpdateStat[]
    apply_kwargs = (
        maxdim = cfg.maxdim,
        cutoff = cfg.cutoff,
        normalize_tensors = cfg.normalize_tensors,
    )

    for (index, gate) in enumerate(gates)
        support = gate_vertices[index]
        update_required = length(support) >= 2 && any(v in affected_vertices for v in support)
        if update_required
            cache, stat = converge_bp(cache, cfg)
            push!(bp_stats, stat)
            empty!(affected_vertices)
        end

        adapted_gate = TNQS.adapt_gate(gate, cache)
        cache, truncation_errors[index] = TNQS.apply_gate!(
            adapted_gate,
            cache;
            v⃗ = support,
            apply_kwargs,
        )
        union!(affected_vertices, support)
    end

    cache, stat = converge_bp(cache, cfg)
    push!(bp_stats, stat)
    return cache, truncation_errors, bp_stats
end

function evolve_bp!(pair::BPPair, circuit::Vector, cfg::SimulationConfig)
    scv_start = time_ns()
    pair.scv, scv_errors, scv_bp = _apply_layer_controlled(pair.scv, circuit, cfg)
    scv_seconds = (time_ns() - scv_start) / 1.0e9

    loop_start = time_ns()
    pair.loop, loop_errors, loop_bp = _apply_layer_controlled(pair.loop, circuit, cfg)
    loop_seconds = (time_ns() - loop_start) / 1.0e9

    return pair, (
        scv = (runtime = scv_seconds, errors = scv_errors, bp = scv_bp),
        loop = (runtime = loop_seconds, errors = loop_errors, bp = loop_bp),
    )
end

function measure_bp_state(cache, plan::MeasurementPlan, links)
    observables = [pauli_observable(key) for key in plan.keys]
    raw = TNQS.expect(cache, observables)
    values = Dict{PauliKey, Float64}()
    for (key, value) in zip(plan.keys, raw)
        abs(imag(value)) <= 1.0e-8 || @warn "observable has a non-negligible imaginary part" key value
        values[key] = real(value)
    end
    return make_snapshot(plan, links, values)
end

function bp_norm_drift(cache)
    return abs(real(TNQS.norm_sqr(cache; alg = "bp")) - 1.0)
end

bp_maxdim(cache) = TNQS.maxvirtualdim(cache)
