struct PauliKey
    axis::Char
    vertices::Tuple{Vararg{Vertex}}
end

struct LinearTerm
    coefficient::Float64
    key::Union{Nothing, PauliKey}
end

struct WilsonSpec
    id::String
    x0::Int
    y0::Int
    width::Int
    height::Int
    key::PauliKey
end

struct StringSpec
    id::String
    direction::Symbol
    distance::Int
    target::Vertex
    key::PauliKey
end

struct MagneticCorrelationSpec
    id::String
    target::Vertex
    distance::Int
    x_pair_key::PauliKey
end

struct FluxCorrelationSpec
    anchor_link_id::Int
    target_link_id::Int
    terms::Vector{LinearTerm}
end

mutable struct MeasurementPlan
    nx::Int
    keys::Vector{PauliKey}
    key_index::Dict{PauliKey, Int}
    link_keys::Vector{PauliKey}
    x_keys::Dict{Vertex, PauliKey}
    z_keys::Dict{Vertex, PauliKey}
    loop_survival_terms::Vector{LinearTerm}
    wilson_specs::Vector{WilsonSpec}
    string_specs::Vector{StringSpec}
    magnetic_specs::Vector{MagneticCorrelationSpec}
    flux_specs::Vector{FluxCorrelationSpec}
    anchor_link_ids::Vector{Int}
end

function _ordered_vertices(vs, nx::Int)
    unique_vs = unique(Vertex[v for v in vs])
    sort!(unique_vs; by = v -> site_index(v, nx))
    return unique_vs
end

function register_pauli!(plan::MeasurementPlan, axis::Char, vertices)
    ordered = _ordered_vertices(vertices, plan.nx)
    isempty(ordered) && throw(ArgumentError("the identity is represented by a constant, not a PauliKey"))
    key = PauliKey(axis, Tuple(ordered))
    if !haskey(plan.key_index, key)
        push!(plan.keys, key)
        plan.key_index[key] = length(plan.keys)
    end
    return key
end

function _z_support(links::AbstractVector{DualLink})
    parity = Dict{Vertex, Bool}()
    for link in links
        parity[link.v1] = !get(parity, link.v1, false)
        if !isnothing(link.v2)
            parity[link.v2] = !get(parity, link.v2, false)
        end
    end
    return [v for (v, odd) in parity if odd]
end

"""Compile a product of link projectors n_b=(1-E_b)/2 into Pauli-Z terms."""
function compile_flux_projector_product!(plan::MeasurementPlan, links::Vector{DualLink})
    m = length(links)
    terms = LinearTerm[]
    for mask in 0:(Int(1) << m) - 1
        selected = DualLink[]
        parity = 0
        for j in 1:m
            if !iszero(mask & (Int(1) << (j - 1)))
                push!(selected, links[j])
                parity += 1
            end
        end
        coefficient = (isodd(parity) ? -1.0 : 1.0) / (2.0^m)
        support = _z_support(selected)
        key = isempty(support) ? nothing : register_pauli!(plan, 'Z', support)
        push!(terms, LinearTerm(coefficient, key))
    end
    return terms
end

function _empty_plan(nx::Int)
    return MeasurementPlan(
        nx,
        PauliKey[],
        Dict{PauliKey, Int}(),
        PauliKey[],
        Dict{Vertex, PauliKey}(),
        Dict{Vertex, PauliKey}(),
        LinearTerm[],
        WilsonSpec[],
        StringSpec[],
        MagneticCorrelationSpec[],
        FluxCorrelationSpec[],
        Int[],
    )
end

function _axial_targets(cfg::SimulationConfig, defect::Vertex)
    targets = Tuple{Symbol, Int, Vertex}[]
    directions = ((:right, 1, 0), (:left, -1, 0), (:up, 0, 1), (:down, 0, -1))
    for (name, dx, dy) in directions
        for d in 1:cfg.correlation_max_distance
            v = (defect[1] + d * dx, defect[2] + d * dy)
            1 <= v[1] <= cfg.nx || continue
            1 <= v[2] <= cfg.ny || continue
            push!(targets, (name, d, v))
        end
    end
    return targets
end

function build_measurement_plan(cfg::SimulationConfig, links::Vector{DualLink})
    plan = _empty_plan(cfg.nx)
    defect = resolved_defect(cfg)

    for link in links
        push!(plan.link_keys, register_pauli!(plan, 'Z', link_vertices(link)))
    end
    for v in dual_vertices(cfg)
        plan.x_keys[v] = register_pauli!(plan, 'X', [v])
        plan.z_keys[v] = register_pauli!(plan, 'Z', [v])
    end

    central_links = star_links(links, defect)
    plan.anchor_link_ids = [link.id for link in central_links]
    plan.loop_survival_terms = compile_flux_projector_product!(plan, central_links)

    # Small direct-lattice Wilson contours map to products of X over the
    # enclosed dual plaquettes. Keep areas small enough for a realistic QPU
    # principal set; all of them share one global X-basis measurement setting.
    for height in 1:cfg.wilson_max_area, width in 1:cfg.wilson_max_area
        width * height <= cfg.wilson_max_area || continue
        width <= cfg.nx || continue
        height <= cfg.ny || continue
        x0 = clamp(defect[1] - fld(width - 1, 2), 1, cfg.nx - width + 1)
        y0 = clamp(defect[2] - fld(height - 1, 2), 1, cfg.ny - height + 1)
        rect = [(x, y) for y in y0:(y0 + height - 1) for x in x0:(x0 + width - 1)]
        key = register_pauli!(plan, 'X', rect)
        id = "W_$(width)x$(height)_x$(x0)_y$(y0)"
        push!(plan.wilson_specs, WilsonSpec(id, x0, y0, width, height, key))
    end

    for (direction, distance, target) in _axial_targets(cfg, defect)
        zkey = register_pauli!(plan, 'Z', [defect, target])
        id = "E_$(direction)_d$(distance)"
        push!(plan.string_specs, StringSpec(id, direction, distance, target, zkey))

        xxkey = register_pauli!(plan, 'X', [defect, target])
        cid = "C_B_$(direction)_d$(distance)"
        push!(plan.magnetic_specs, MagneticCorrelationSpec(cid, target, distance, xxkey))
    end

    by_id = Dict(link.id => link for link in links)
    for anchor_id in plan.anchor_link_ids, target in links
        pair_terms = compile_flux_projector_product!(plan, [by_id[anchor_id], target])
        push!(plan.flux_specs, FluxCorrelationSpec(anchor_id, target.id, pair_terms))
    end
    return plan
end

function pauli_observable(key::PauliKey)
    pauli = repeat(string(key.axis), length(key.vertices))
    return (pauli, collect(key.vertices))
end

function evaluate_linear(terms::Vector{LinearTerm}, values::Dict{PauliKey, Float64})
    total = 0.0
    for term in terms
        total += term.coefficient * (isnothing(term.key) ? 1.0 : values[term.key])
    end
    return total
end

struct StateSnapshot
    pauli::Dict{PauliKey, Float64}
    link_flux::Vector{Float64}
    site_x::Dict{Vertex, Float64}
    site_z::Dict{Vertex, Float64}
    loop_survival::Float64
    wilson::Dict{String, Float64}
    strings::Dict{String, Float64}
    magnetic_connected::Dict{String, Float64}
    flux_connected::Dict{Tuple{Int, Int}, Float64}
end

function make_snapshot(plan::MeasurementPlan, links::Vector{DualLink}, values::Dict{PauliKey, Float64})
    link_flux = [(1.0 - values[key]) / 2.0 for key in plan.link_keys]
    site_x = Dict(v => values[key] for (v, key) in plan.x_keys)
    site_z = Dict(v => values[key] for (v, key) in plan.z_keys)
    loop_survival = evaluate_linear(plan.loop_survival_terms, values)
    wilson = Dict(spec.id => values[spec.key] for spec in plan.wilson_specs)
    strings = Dict(spec.id => values[spec.key] for spec in plan.string_specs)

    # Every magnetic specification is constructed from one common center. Recover
    # it as the non-target endpoint of the first XX key.
    magnetic_connected = Dict{String, Float64}()
    for spec in plan.magnetic_specs
        endpoints = collect(spec.x_pair_key.vertices)
        center = only(filter(v -> v != spec.target, endpoints))
        magnetic_connected[spec.id] = values[spec.x_pair_key] - site_x[center] * site_x[spec.target]
    end

    flux_connected = Dict{Tuple{Int, Int}, Float64}()
    for spec in plan.flux_specs
        joint = evaluate_linear(spec.terms, values)
        flux_connected[(spec.anchor_link_id, spec.target_link_id)] =
            joint - link_flux[spec.anchor_link_id] * link_flux[spec.target_link_id]
    end
    return StateSnapshot(
        values,
        link_flux,
        site_x,
        site_z,
        loop_survival,
        wilson,
        strings,
        magnetic_connected,
        flux_connected,
    )
end
