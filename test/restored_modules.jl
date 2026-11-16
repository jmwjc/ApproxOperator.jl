# Smoke tests for the modules that S5 brought back into the build.
#
# `phasefield.jl` and `plasticity.jl` were not `include`d before S5: they still
# dispatched on `Operator{:tag}`, a struct that had been deleted, so any call was
# a `UndefVarError`. `grkgsi.jl` depended on `AbstractReproducingKernel`, a type
# that was never defined anywhere in the package, and has been removed.
# `rkwave.jl` only extends `ReproducingKernel{:Wave2D}`, which does exist, so it
# is now included.
#
# These tests do not check the physics — they check that every operator is
# reachable and callable on a real element, which is what was broken.

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
using Test
import Gmsh: gmsh
using LinearAlgebra: norm

@testset "S5: modules brought back into the build" begin

    @testset "operators are defined" begin
        for name in (:∫v²uₓuₓdx, :∫vₓvₓvvdx, :∫vₓvₓvvdx_hard_device,
                     :∫∫∇v∇vvvdxdy, :UPDATE_PFM_1D, :UPDATE_PFM_2D, :CRACK_NORMAL,
                     :∫vₓσdx, :∫vᵢσdΩ_tresca, :∫vᵢσdΩ_mohr_coulomb,
                     :∫vᵢσdΩ_frictional_contact, :∫vᵢσdΩ_frictional_contact_2,
                     :∫vᵢσdΩ_mc_phasefield)
            @test isdefined(ApproxOperator, name)
        end
    end

    @testset "rkwave adds the Wave2D basis" begin
        for tag in (:Wave2D, :Wave2DY)
            T = ApproxOperator.ReproducingKernel{tag,:Linear1D,:Linear1D}
            rk = T(ApproxOperator.𝑿ᵢ[], ApproxOperator.𝑿ₛ[])
            @test ApproxOperator.get𝑛𝒑(rk) > 0
            @test length(ApproxOperator.get𝒑(rk, (1.0, 2.0, 0.0))) ==
                  ApproxOperator.get𝑛𝒑(rk)
        end
        # The trigonometric enrichment really is trigonometric: at x = 0 the
        # sin-term vanishes and the cos-term is 1.
        T = ApproxOperator.ReproducingKernel{:Wave2D,:Linear1D,:Linear1D}
        rk = T(ApproxOperator.𝑿ᵢ[], ApproxOperator.𝑿ₛ[])
        @test ApproxOperator.get𝒑(rk, (0.0, 0.0, 0.0))[2] ≈ 0.0 atol = 1.0e-14
        @test ApproxOperator.get𝒑(rk, (0.0, 0.0, 0.0))[3] ≈ 1.0 atol = 1.0e-14
    end

    @testset "plasticity operators assemble on a real element" begin
        gmsh.initialize()
        try
            gmsh.model.add("plasticity")
            p = [gmsh.model.geo.addPoint(x, y, 0.0)
                 for (x,y) in ((0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0,1.0))]
            l = [gmsh.model.geo.addLine(p[1],p[2]), gmsh.model.geo.addLine(p[2],p[3]),
                 gmsh.model.geo.addLine(p[3],p[4]), gmsh.model.geo.addLine(p[4],p[1])]
            cl = gmsh.model.geo.addCurveLoop(l)
            s  = gmsh.model.geo.addPlaneSurface([cl])
            gmsh.model.geo.synchronize()
            for li in l; gmsh.model.mesh.setTransfiniteCurve(li, 3); end
            gmsh.model.mesh.setTransfiniteSurface(s)
            gmsh.model.mesh.setRecombine(2, s)
            gmsh.model.addPhysicalGroup(2, [s], 1)
            gmsh.model.setPhysicalName(2, 1, "Ω")
            gmsh.model.mesh.generate(2)

            entities = getPhysicalGroups()
            nodes    = get𝑿ᵢ()
            elms     = getElements(nodes, entities["Ω"], 2)
            set∇𝝭!(elms)

            # State the plasticity operators read and update. These fields are
            # normally created by `prescribe!` on a mesh built for plasticity;
            # here they are seeded directly so the operators can be exercised.
            data = getfield(elms[1].𝓖[1], :data)
            nq   = length(data[:x][2])
            nn   = maximum(xᵢ.𝐼 for ap in elms for xᵢ in ap.𝓒)
            # Fields read as `ξ.name` (one value per quadrature point).
            for name in (:σ₁₁, :σ₂₂, :σ₁₂, :εᵖ₁₁, :εᵖ₂₂, :εᵖ₁₂, :μ̄, :𝜙)
                data[name] = (2, zeros(nq))
            end
            # Fields read as `xᵢ.name` live in the *support nodes'* dictionary,
            # which is a different object from the quadrature dictionary.
            d𝒙 = getfield(elms[1].𝓒[1], :data)
            for name in (:d₁, :d₂, :Δd₁, :Δd₂, :v, :d)
                d𝒙[name] = (1, zeros(nn))
            end
            data[:μ̄][2] .= 0.3
            data[:𝜙][2] .= pi/3
            d𝒙[:Δd₁][2] .= 1.0e-6

            n = 2 * maximum(xᵢ.𝐼 for ap in elms for xᵢ in ap.𝓒)

            for (name, kwargs) in (
                    (:∫vᵢσdΩ_tresca, (; λ = 1.0, μ = 0.5, c = 0.1)),
                    (:∫vᵢσdΩ_mohr_coulomb, (; λ = 1.0, μ = 0.5, c = 0.1, tol = 1.0e-13)),
                    (:∫vᵢσdΩ_mc_phasefield, (; E = 2.0, ν = 0.3, tol = 1.0e-13)),
                )
                f = getfield(ApproxOperator, name)
                k = zeros(n, n); fint = zeros(n)
                f(elms[1], k, fint; kwargs...)
                @test all(isfinite, k)
                @test all(isfinite, fint)
            end
        finally
            gmsh.finalize()
        end
    end

    @testset "phase-field update operators run on a real element" begin
        gmsh.initialize()
        try
            gmsh.model.add("phasefield")
            p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0)
            p2 = gmsh.model.geo.addPoint(1.0, 0.0, 0.0)
            l  = gmsh.model.geo.addLine(p1, p2)
            gmsh.model.geo.synchronize()
            gmsh.model.mesh.setTransfiniteCurve(l, 3)
            gmsh.model.addPhysicalGroup(1, [l], 1)
            gmsh.model.setPhysicalName(1, 1, "Ω")
            gmsh.model.mesh.generate(1)

            entities = getPhysicalGroups()
            nodes    = get𝑿ᵢ()
            elms     = getElements(nodes, entities["Ω"], 2)
            set∇𝝭!(elms)
            # Seed the fields `UPDATE_PFM_1D` reads: nodal displacement `u`,
            # damage `v`, and the history variable `ℋ`.
            data = getfield(elms[1].𝓖[1], :data)
            d𝒙   = getfield(elms[1].𝓒[1], :data)
            nn   = maximum(xᵢ.𝐼 for ap in elms for xᵢ in ap.𝓒)
            nq   = length(data[:x][2])
            for name in (:d, :Δd, :v, :u)
                d𝒙[name] = (1, zeros(nn))
            end
            data[:ℋ] = (2, zeros(nq))
            data[:k] = (2, ones(nq))   # fracture toughness, read as `ξ.k`
            d𝒙[:d][2] .= 0.5

            # `UPDATE_PFM_1D` writes the history variable back into `ξ.ℋ`.
            for ap in elms
                ApproxOperator.UPDATE_PFM_1D(ap)
            end
            @test all(x.ℋ ≥ 0 for ap in elms for x in ap.𝓖)

            # The 1-D damage operator assembles a finite system.
            n = maximum(xᵢ.𝐼 for ap in elms for xᵢ in ap.𝓒)
            k = zeros(n, n); f = zeros(n)
            ApproxOperator.∫vₓvₓvvdx(elms, k, f; l = 0.01, kc = 1.0)
            @test all(isfinite, k)
            @test norm(k) > 0.0
        finally
            gmsh.finalize()
        end
    end
end
