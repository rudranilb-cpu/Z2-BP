const Vertex = NTuple{2, Int}

struct DualLink
    id::Int
    v1::Vertex
    v2::Union{Nothing, Vertex}
    orientation::Symbol
    x::Float64
    y::Float64
end

site_index(v::Vertex, nx::Int) = v[1] + (v[2] - 1) * nx
site_index(x::Int, y::Int, nx::Int) = x + (y - 1) * nx

internal_bond_count(nx::Int, ny::Int) = (nx - 1) * ny + nx * (ny - 1)
full_open_link_count(nx::Int, ny::Int) = internal_bond_count(nx, ny) + 2 * nx + 2 * ny

function dual_vertices(cfg::SimulationConfig)
    return [(x, y) for y in 1:cfg.ny for x in 1:cfg.nx]
end

"""
Return a stable link ordering. Internal x bonds come first, then internal y bonds.
For `fixed_exterior`, the four physical boundary families follow. Boundary links
have `v2 = nothing` and map to a one-site Z operator with the exterior spin fixed +1.
"""
function dual_links(cfg::SimulationConfig)
    links = DualLink[]
    nextid = 1
    for y in 1:cfg.ny, x in 1:(cfg.nx - 1)
        push!(links, DualLink(nextid, (x, y), (x + 1, y), :x, x + 0.5, Float64(y)))
        nextid += 1
    end
    for y in 1:(cfg.ny - 1), x in 1:cfg.nx
        push!(links, DualLink(nextid, (x, y), (x, y + 1), :y, Float64(x), y + 0.5))
        nextid += 1
    end
    if cfg.boundary == :fixed_exterior
        for y in 1:cfg.ny
            push!(links, DualLink(nextid, (1, y), nothing, :left, 0.5, Float64(y)))
            nextid += 1
            push!(links, DualLink(nextid, (cfg.nx, y), nothing, :right, cfg.nx + 0.5, Float64(y)))
            nextid += 1
        end
        for x in 1:cfg.nx
            push!(links, DualLink(nextid, (x, 1), nothing, :bottom, Float64(x), 0.5))
            nextid += 1
            push!(links, DualLink(nextid, (x, cfg.ny), nothing, :top, Float64(x), cfg.ny + 0.5))
            nextid += 1
        end
    end
    return links
end

function link_vertices(link::DualLink)
    return isnothing(link.v2) ? Vertex[link.v1] : Vertex[link.v1, link.v2]
end

function internal_degree(v::Vertex, cfg::SimulationConfig)
    x, y = v
    return (x > 1) + (x < cfg.nx) + (y > 1) + (y < cfg.ny)
end

function star_links(links::Vector{DualLink}, v::Vertex)
    return [link for link in links if link.v1 == v || link.v2 == v]
end

function distance_to_open_boundary(v::Vertex, cfg::SimulationConfig)
    x, y = v
    return minimum((x - 1, cfg.nx - x, y - 1, cfg.ny - y))
end

function boundary_z_multiplicity(v::Vertex, cfg::SimulationConfig)
    cfg.boundary == :fixed_exterior || return 0
    x, y = v
    return (x == 1) + (x == cfg.nx) + (y == 1) + (y == cfg.ny)
end

function link_radius_angle(link::DualLink, defect::Vertex)
    dx = link.x - defect[1]
    dy = link.y - defect[2]
    return hypot(dx, dy), atan(dy, dx)
end
