using ApproxOperator
using Test
import Gmsh: gmsh

function build_timoshenko_ss_model!(n_elem::Int, L::Float64)
    gmsh.model.add("timoshenko-static-ss")

    p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0)
    p2 = gmsh.model.geo.addPoint(L, 0.0, 0.0)
    Ω = gmsh.model.geo.addLine(p1, p2)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.setTransfiniteCurve(Ω, n_elem + 1)
    gmsh.model.addPhysicalGroup(1, [Ω], 1)
    gmsh.model.setPhysicalName(1, 1, "Ω")
    gmsh.model.mesh.generate(1)

    return ApproxOperator.GmshImport.getPhysicalGroups(), ApproxOperator.GmshImport.get𝑿ᵢ()
end

@testset "ApproxOperator.jl - Timoshenko Static SS" begin
    E = 1.0e8
    ν = 0.3
    κ = 5.0 / 6.0
    q = 1.0
    L = 1.0
    h = 0.001
    n_elem = 100

    G = E / (2.0 * (1.0 + ν))
    A = h
    I = h^3 / 12.0
    EI = E * I
    kGA = κ * G * A

    gmsh.initialize()
    try
        entities, nodes = build_timoshenko_ss_model!(n_elem, L)

        elements = ApproxOperator.GmshImport.getElements(nodes, entities["Ω"], 1)
        prescribe!(elements, :EI=>EI, :kGA=>kGA, :q=>q)
        set∇𝝭!(elements)

        ndofs = 2 * length(nodes)
        K = zeros(ndofs, ndofs)
        F = zeros(ndofs)

        𝑎 = [
            ApproxOperator.Timoshenko.∫κEIκds=>elements,
            ApproxOperator.Timoshenko.∫γkGAγds=>elements
        ]
        𝑎(K)

        𝑓 = ApproxOperator.Timoshenko.∫vqds=>elements
        𝑓(F)

        fixed = Int[]
        for xᵢ in nodes
            if isapprox(xᵢ.x, 0.0; atol=1.0e-12) || isapprox(xᵢ.x, L; atol=1.0e-12)
                push!(fixed, 2 * xᵢ.𝐼 - 1)
            end
        end
        sort!(fixed)
        free = setdiff(collect(eachindex(F)), fixed)

        U = zeros(length(F))
        U[free] = K[free, free] \ F[free]

        d = U[1:2:end]
        φ = U[2:2:end]
        push!(nodes, :d => d, :φ => φ)

        ξ₅ = [
            -0.906179845938664,
            -0.5384693101056831,
            0.0,
            0.5384693101056831,
            0.906179845938664
        ]
        w₅ = [
            0.2369268850561891,
            0.4786286704993665,
            0.5688888888888889,
            0.4786286704993665,
            0.2369268850561891
        ]
        localcoord = Float64[]
        for ξ in ξ₅
            append!(localcoord, (ξ, 0.0, 0.0))
        end

        elements_l2 = ApproxOperator.GmshImport.getElements(nodes, entities["Ω"], (localcoord, w₅))
        prescribe!(
            elements_l2,
            :u => (x, y, z) -> ApproxOperator.Timoshenko.w_exact_ss(x, E, I, κ, G, A, L, q)
        )
        set𝝭!(elements_l2)
        L2_ratio = ApproxOperator.Timoshenko.L₂(elements_l2)
        L2_pct = 100.0 * L2_ratio

        xcoords = [xᵢ.x for xᵢ in nodes]
        midpoint = argmin(abs.(xcoords .- L / 2.0))
        left_node = (first(fixed) + 1) ÷ 2
        right_node = (last(fixed) + 1) ÷ 2
        @test abs(xcoords[midpoint] - L / 2.0) ≤ 1.0e-8
        @test d[left_node] ≈ 0.0 atol=1.0e-12
        @test d[right_node] ≈ 0.0 atol=1.0e-12
        @test d[midpoint] ≈ 1.5622539164827454 atol=1.0e-7
        @test L2_pct ≈ 0.0249617938201175 atol=5.0e-6
    finally
        gmsh.finalize()
    end
end
