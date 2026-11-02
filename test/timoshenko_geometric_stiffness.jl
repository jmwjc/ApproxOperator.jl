using ApproxOperator
using Test
import Gmsh: gmsh

function build_timoshenko_geo_model!(n_elem::Int, L::Float64)
    gmsh.model.add("timoshenko-geometric-stiffness")

    p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0)
    p2 = gmsh.model.geo.addPoint(L, 0.0, 0.0)
    Ω = gmsh.model.geo.addLine(p1, p2)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.setTransfiniteCurve(Ω, n_elem + 1)
    gmsh.model.addPhysicalGroup(1, [Ω], 1)
    gmsh.model.setPhysicalName(1, 1, "Ω")
    gmsh.model.mesh.generate(1)

    entities = ApproxOperator.GmshImport.getPhysicalGroups()
    nodes = ApproxOperator.GmshImport.get𝑿ᵢ()
    return entities, nodes
end

@testset "ApproxOperator.jl - Timoshenko Geometric Stiffness Operators" begin
    gmsh.initialize()
    try
        entities, nodes = build_timoshenko_geo_model!(2, 1.0)
        elements = ApproxOperator.GmshImport.getElements(nodes, entities["Ω"], 2)
        set∇𝝭!(elements)

        n = length(nodes)
        KwwG = zeros(n, n)
        KwφGH = zeros(n, n)
        KφwGH = zeros(n, n)
        KφφGH = zeros(n, n)

        (ApproxOperator.Timoshenko.∫wwGdΩ => elements)(KwwG)
        (ApproxOperator.Timoshenko.∫wφGHdΩ => elements)(KwφGH)
        (ApproxOperator.Timoshenko.∫φwGHdΩ => elements)(KφwGH)
        (ApproxOperator.Timoshenko.∫φφGHdΩ => elements)(KφφGH)

        @test all(isfinite, KwwG)
        @test all(isfinite, KwφGH)
        @test all(isfinite, KφwGH)
        @test all(isfinite, KφφGH)

        @test KwwG ≈ KwwG' atol=1.0e-12
        @test KφφGH ≈ KφφGH' atol=1.0e-12
        @test KwφGH ≈ KφwGH' atol=1.0e-12

        @test norm(KwwG) > 0.0
        @test norm(KwφGH) > 0.0
        @test norm(KφwGH) > 0.0
        @test norm(KφφGH) > 0.0
    finally
        gmsh.finalize()
    end
end

@testset "ApproxOperator.jl - Timoshenko Geometric Stiffness Spot Check (Seg2, 1 element)" begin
    gmsh.initialize()
    try
        entities, nodes = build_timoshenko_geo_model!(1, 1.0)
        elements = ApproxOperator.GmshImport.getElements(nodes, entities["Ω"], 1)
        set∇𝝭!(elements)

        n = length(nodes)
        KwwG = zeros(n, n)
        KwφGH = zeros(n, n)
        KφwGH = zeros(n, n)
        KφφGH = zeros(n, n)

        (ApproxOperator.Timoshenko.∫wwGdΩ => elements)(KwwG)
        (ApproxOperator.Timoshenko.∫wφGHdΩ => elements)(KwφGH)
        (ApproxOperator.Timoshenko.∫φwGHdΩ => elements)(KφwGH)
        (ApproxOperator.Timoshenko.∫φφGHdΩ => elements)(KφφGH)

        # Seg2, L=1, 1-point rule at center:
        # B = [-1, 1], N = [0.5, 0.5], dΩ = 1
        KwwG_expected = [1.0 -1.0; -1.0 1.0]
        Kwφ_expected = [-0.5 -0.5; 0.5 0.5]
        Kφw_expected = [-0.5 0.5; -0.5 0.5]
        Kφφ_expected = [0.25 0.25; 0.25 0.25]

        @test KwwG ≈ KwwG_expected atol=1.0e-12
        @test KwφGH ≈ Kwφ_expected atol=1.0e-12
        @test KφwGH ≈ Kφw_expected atol=1.0e-12
        @test KφφGH ≈ Kφφ_expected atol=1.0e-12
    finally
        gmsh.finalize()
    end
end
