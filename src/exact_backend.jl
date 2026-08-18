mutable struct ExactPair
    scv::Vector{ComplexF64}
    loop::Vector{ComplexF64}
end

function initialize_exact(cfg::SimulationConfig)
    n = cfg.nx * cfg.ny
    dim = Int(1) << n
    scv = zeros(ComplexF64, dim)
    loop = zeros(ComplexF64, dim)
    scv[1] = 1.0 + 0.0im
    defect_q = site_index(resolved_defect(cfg), cfg.nx) - 1
    loop[(Int(1) << defect_q) + 1] = 1.0 + 0.0im
    return ExactPair(scv, loop)
end

function _apply_x_site!(state::Vector{ComplexF64}, q::Int, alpha::Float64)
    stride = Int(1) << q
    block = stride << 1
    c = cos(alpha)
    s = im * sin(alpha)
    for base in 0:block:(length(state) - 1)
        for offset in 0:(stride - 1)
            i0 = base + offset + 1
            i1 = i0 + stride
            a0, a1 = state[i0], state[i1]
            state[i0] = c * a0 + s * a1
            state[i1] = s * a0 + c * a1
        end
    end
    return state
end

function _apply_x_layer!(state::Vector{ComplexF64}, cfg::SimulationConfig, scale::Float64)
    alpha = cfg.K * cfg.dt * scale
    for q in 0:(cfg.nx * cfg.ny - 1)
        _apply_x_site!(state, q, alpha)
    end
    return state
end

_zsign(bits::Int, q::Int) = iszero(bits & (Int(1) << q)) ? 1.0 : -1.0

function _apply_diagonal_layer!(state::Vector{ComplexF64}, cfg::SimulationConfig, scale::Float64)
    alpha = cfg.Gamma * cfg.dt * scale
    for bits in 0:(length(state) - 1)
        exponent = 0.0
        for y in 1:cfg.ny, x in 1:(cfg.nx - 1)
            q1 = site_index(x, y, cfg.nx) - 1
            q2 = site_index(x + 1, y, cfg.nx) - 1
            exponent += _zsign(bits, q1) * _zsign(bits, q2)
        end
        for y in 1:(cfg.ny - 1), x in 1:cfg.nx
            q1 = site_index(x, y, cfg.nx) - 1
            q2 = site_index(x, y + 1, cfg.nx) - 1
            exponent += _zsign(bits, q1) * _zsign(bits, q2)
        end
        if cfg.boundary == :fixed_exterior
            for v in dual_vertices(cfg)
                m = boundary_z_multiplicity(v, cfg)
                iszero(m) && continue
                q = site_index(v, cfg.nx) - 1
                exponent += m * _zsign(bits, q)
            end
        end
        state[bits + 1] *= cis(alpha * exponent)
    end
    return state
end

function exact_step!(state::Vector{ComplexF64}, cfg::SimulationConfig)
    if cfg.trotter == :lie_x_zz
        _apply_x_layer!(state, cfg, 1.0)
        _apply_diagonal_layer!(state, cfg, 1.0)
    elseif cfg.trotter == :lie_zz_x
        _apply_diagonal_layer!(state, cfg, 1.0)
        _apply_x_layer!(state, cfg, 1.0)
    elseif cfg.trotter == :strang
        _apply_x_layer!(state, cfg, 0.5)
        _apply_diagonal_layer!(state, cfg, 1.0)
        _apply_x_layer!(state, cfg, 0.5)
    end
    return state
end

function evolve_exact!(pair::ExactPair, cfg::SimulationConfig)
    t0 = time_ns()
    exact_step!(pair.scv, cfg)
    exact_step!(pair.loop, cfg)
    runtime = (time_ns() - t0) / 1.0e9
    return pair, runtime
end

function _expect_z(state::Vector{ComplexF64}, vertices::Tuple{Vararg{Vertex}}, cfg::SimulationConfig)
    total = 0.0
    qs = [site_index(v, cfg.nx) - 1 for v in vertices]
    for bits in 0:(length(state) - 1)
        sign = 1.0
        for q in qs
            sign *= _zsign(bits, q)
        end
        total += abs2(state[bits + 1]) * sign
    end
    return total
end

function _expect_x(state::Vector{ComplexF64}, vertices::Tuple{Vararg{Vertex}}, cfg::SimulationConfig)
    flipmask = 0
    for v in vertices
        flipmask ⊻= Int(1) << (site_index(v, cfg.nx) - 1)
    end
    total = 0.0 + 0.0im
    for bits in 0:(length(state) - 1)
        total += conj(state[bits + 1]) * state[(bits ⊻ flipmask) + 1]
    end
    return real(total)
end

function measure_exact_state(state::Vector{ComplexF64}, cfg::SimulationConfig, plan::MeasurementPlan, links)
    values = Dict{PauliKey, Float64}()
    for key in plan.keys
        value = key.axis == 'Z' ? _expect_z(state, key.vertices, cfg) : _expect_x(state, key.vertices, cfg)
        values[key] = value
    end
    return make_snapshot(plan, links, values)
end

function exact_norm_drift(state::Vector{ComplexF64})
    return abs(sum(abs2, state) - 1.0)
end
