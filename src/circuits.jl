function x_layer(cfg::SimulationConfig, g, scale::Float64 = 1.0)
    angle = -2.0 * cfg.K * cfg.dt * scale
    return [("Rx", [v], angle) for v in vertices(g)]
end

function diagonal_layer(cfg::SimulationConfig, g, scale::Float64 = 1.0)
    layer = Any[]
    angle = -2.0 * cfg.Gamma * cfg.dt * scale
    # Four colors are sufficient for an open square grid. Terms commute, so
    # coloring only schedules disjoint simple updates; it does not change U_ZZ.
    for colored_edges in edge_color(g, 4)
        append!(layer, [("Rzz", [src(pair), dst(pair)], angle) for pair in colored_edges])
    end
    if cfg.boundary == :fixed_exterior
        for v in vertices(g)
            m = boundary_z_multiplicity(v, cfg)
            iszero(m) || push!(layer, ("Rz", [v], angle * m))
        end
    end
    return layer
end

"""
Build one sequential gate list. In the reference `lie_x_zz` convention the list
applies X gates first and ZZ/boundary-Z gates second, hence the state is acted on
by U_ZZ(dt) U_X(dt). `strang` is X(dt/2), ZZ(dt), X(dt/2).
"""
function build_trotter_layer(cfg::SimulationConfig, g)
    if cfg.trotter == :lie_x_zz
        return Any[x_layer(cfg, g)..., diagonal_layer(cfg, g)...]
    elseif cfg.trotter == :lie_zz_x
        return Any[diagonal_layer(cfg, g)..., x_layer(cfg, g)...]
    elseif cfg.trotter == :strang
        return Any[x_layer(cfg, g, 0.5)..., diagonal_layer(cfg, g)..., x_layer(cfg, g, 0.5)...]
    end
    error("unreachable Trotter convention")
end
