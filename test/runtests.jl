using Test
using Z2GaugeBenchmarks

const Z2 = Z2GaugeBenchmarks

@testset "lattice and t=0 identities" begin
    cfg = SimulationConfig(nx = 3, ny = 3, backend = :exact, defect_x = 2, defect_y = 2)
    @test internal_bond_count(8, 6) == 82
    @test full_open_link_count(8, 6) == 110
    @test length(dual_links(cfg)) == 12

    pair = Z2.initialize_exact(cfg)
    links = dual_links(cfg)
    plan = Z2.build_measurement_plan(cfg, links)
    scv = Z2.measure_exact_state(pair.scv, cfg, plan, links)
    loop = Z2.measure_exact_state(pair.loop, cfg, plan, links)
    @test sum(scv.link_flux) ≈ 0 atol = 1.0e-14
    @test sum(loop.link_flux) - sum(scv.link_flux) ≈ 4 atol = 1.0e-14
    @test loop.loop_survival ≈ 1 atol = 1.0e-14
    front = Z2._front_metrics(cfg, links, scv, loop)
    @test front.mean ≈ 0.5 atol = 1.0e-14
    @test front.rms ≈ 0.5 atol = 1.0e-14
    @test front.quantile ≈ 0.5 atol = 1.0e-14
    @test front.participation ≈ 4 atol = 1.0e-14
    @test front.anisotropy ≈ 0 atol = 1.0e-14
end

@testset "gate sign and first Lie step" begin
    cfg = SimulationConfig(nx = 3, ny = 3, K = 3.0, Gamma = 1.0, dt = 0.05,
        backend = :exact, defect_x = 2, defect_y = 2, trotter = :lie_x_zz)
    pair = Z2.initialize_exact(cfg)
    Z2.exact_step!(pair.scv, cfg)
    Z2.exact_step!(pair.loop, cfg)
    links = dual_links(cfg)
    plan = Z2.build_measurement_plan(cfg, links)
    scv = Z2.measure_exact_state(pair.scv, cfg, plan, links)
    loop = Z2.measure_exact_state(pair.loop, cfg, plan, links)
    @test sum(loop.link_flux) - sum(scv.link_flux) ≈ 4cos(2cfg.K * cfg.dt)^2 atol = 1.0e-11
    @test Z2.exact_norm_drift(pair.scv) < 1.0e-12
    @test Z2.exact_norm_drift(pair.loop) < 1.0e-12

    one_qubit = ComplexF64[1, 0]
    Z2._apply_x_site!(one_qubit, 0, 0.17)
    @test one_qubit ≈ ComplexF64[cos(0.17), im * sin(0.17)] atol = 1.0e-14

    diagonal_cfg = SimulationConfig(nx = 2, ny = 2, K = 0.0, Gamma = 0.7, dt = 0.08,
        backend = :exact)
    diagonal_state = Z2.initialize_exact(diagonal_cfg).scv
    Z2._apply_diagonal_layer!(diagonal_state, diagonal_cfg, 1.0)
    @test diagonal_state[1] ≈ cis(4diagonal_cfg.Gamma * diagonal_cfg.dt) atol = 1.0e-14
end

@testset "fixed exterior counting" begin
    cfg = SimulationConfig(nx = 3, ny = 3, backend = :exact, boundary = :fixed_exterior,
        defect_x = 2, defect_y = 2)
    @test length(dual_links(cfg)) == full_open_link_count(3, 3) == 24
    pair = Z2.initialize_exact(cfg)
    links = dual_links(cfg)
    plan = Z2.build_measurement_plan(cfg, links)
    scv = Z2.measure_exact_state(pair.scv, cfg, plan, links)
    loop = Z2.measure_exact_state(pair.loop, cfg, plan, links)
    @test sum(scv.link_flux) ≈ 0 atol = 1.0e-14
    @test sum(loop.link_flux) - sum(scv.link_flux) ≈ 4 atol = 1.0e-14

    open_corner = SimulationConfig(nx = 3, ny = 3, backend = :exact, boundary = :open_dual,
        defect_x = 1, defect_y = 1)
    fixed_corner = SimulationConfig(nx = 3, ny = 3, backend = :exact, boundary = :fixed_exterior,
        defect_x = 1, defect_y = 1)
    @test length(Z2.star_links(dual_links(open_corner), (1, 1))) == 2
    @test length(Z2.star_links(dual_links(fixed_corner), (1, 1))) == 4
end

if get(ENV, "Z2_RUN_BP_TESTS", "0") == "1"
    @testset "BP minimal integration" begin
        cfg = SimulationConfig(nx = 3, ny = 3, backend = :bp, defect_x = 2, defect_y = 2,
            K = 3.0, Gamma = 1.0, dt = 0.05, trotter = :lie_x_zz,
            maxdim = 4, bp_tolerance = 1.0e-9, bp_maxiter = 20, steps = 1)
        pair, stats = Z2.initialize_bp(cfg)
        links = dual_links(cfg)
        observables = [(repeat("Z", length(Z2.link_vertices(link))), Z2.link_vertices(link)) for link in links]
        scv_flux = (1 .- real.(Z2.TNQS.expect(pair.scv, observables))) ./ 2
        loop_flux = (1 .- real.(Z2.TNQS.expect(pair.loop, observables))) ./ 2
        @test all(stat.converged for stat in vcat(stats.scv, stats.loop))
        @test sum(scv_flux) ≈ 0 atol = 1.0e-9
        @test sum(loop_flux) - sum(scv_flux) ≈ 4 atol = 1.0e-9

        circuit = Z2.build_trotter_layer(cfg, pair.graph)
        pair, evolution = Z2.evolve_bp!(pair, circuit, cfg)
        scv_flux = (1 .- real.(Z2.TNQS.expect(pair.scv, observables))) ./ 2
        loop_flux = (1 .- real.(Z2.TNQS.expect(pair.loop, observables))) ./ 2
        @test all(stat.converged for stat in vcat(evolution.scv.bp, evolution.loop.bp))
        @test sum(loop_flux) - sum(scv_flux) ≈ 4cos(2cfg.K * cfg.dt)^2 atol = 1.0e-7
        @test all(evolution.scv.errors .>= 0)
        @test all(evolution.loop.errors .>= 0)

        # Unlike the one-step Z-basis flux, local X after the full layer changes
        # sign if the ZZ rotation sign is reversed. This checks the actual TNQS
        # registry/gate path against the independent package state vector.
        exact_cfg = SimulationConfig(nx = cfg.nx, ny = cfg.ny, K = cfg.K, Gamma = cfg.Gamma,
            dt = cfg.dt, trotter = cfg.trotter, backend = :exact,
            defect_x = cfg.defect_x, defect_y = cfg.defect_y)
        exact_pair = Z2.initialize_exact(exact_cfg)
        Z2.exact_step!(exact_pair.scv, exact_cfg)
        Z2.exact_step!(exact_pair.loop, exact_cfg)
        center = resolved_defect(cfg)
        bp_x_scv = real(Z2.TNQS.expect(pair.scv, ("X", [center])))
        bp_x_loop = real(Z2.TNQS.expect(pair.loop, ("X", [center])))
        @test bp_x_scv ≈ Z2._expect_x(exact_pair.scv, (center,), exact_cfg) atol = 1.0e-10
        @test bp_x_loop ≈ Z2._expect_x(exact_pair.loop, (center,), exact_cfg) atol = 1.0e-10
    end
end
