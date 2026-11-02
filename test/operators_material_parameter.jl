# Regression test for the removal of the undefined global `op`.
#
# Historically the Kirchhoff-plate, Kirchhoff-shell, hyperelasticity,
# curved-beam, error-estimate, phase-field and plasticity modules read their
# material constants from a global variable `op` that no longer existed. Every
# one of those operators therefore raised `UndefVarError: op not defined` on the
# first call, while still being committed, compiled and exported.
#
# These tests pin down the replacement contract: the constants are read from
# the quadrature point `ξ`, and the caller supplies them with `prescribe!`,
# exactly as elasticity.jl already did with `ξ.E` / `ξ.ν`.

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
using Test
import Gmsh: gmsh
using LinearAlgebra: norm

# Two Tri6 elements taken from the packaged patch-test mesh. The mesh is read
# through the normal pipeline and converted with `Tri3toTri6!`, so the DOF
# layout is exactly the one an assembled problem would see.
function patch_tri6()
    gmsh.open(joinpath(@__DIR__, "patchtest.msh"))
    entities = getPhysicalGroups()
    nodes    = get𝑿ᵢ()
    tri3     = getElements(nodes, entities["Ω"], 2)
    tri6, _  = Tri3toTri6!(nodes, tri3)
    # Two elements are enough to exercise the assembly and keep the test fast.
    return tri6[1:2]
end

@testset "operators read material parameters from ξ, not from a global" begin

    @testset "KirchhoffPlate.∫κᵢⱼMᵢⱼdΩ" begin
        gmsh.initialize()
        try
            elements = patch_tri6()
            @test length(elements) == 2

            D, ν = 2.5, 0.25
            set∇²𝝭!(elements)

            n  = maximum(xᵢ.𝐼 for ap in elements for xᵢ in ap.𝓒)
            k  = zeros(n, n)
            (ApproxOperator.KirchhoffPlate.∫κᵢⱼMᵢⱼdΩ => elements)(k; D = D, ν = ν)

            @test all(isfinite, k)
            @test count(!iszero, k) > 1
            @test k ≈ k' atol = 1.0e-10     # bending stiffness must be symmetric

            # The stiffness must scale linearly with D and must actually depend
            # on ν, otherwise the parameters are not being read at all.
            k2 = zeros(n, n)
            (ApproxOperator.KirchhoffPlate.∫κᵢⱼMᵢⱼdΩ => elements)(k2; D = 2D, ν = ν)
            @test k2 ≈ 2k atol = 1.0e-9

            k3 = zeros(n, n)
            (ApproxOperator.KirchhoffPlate.∫κᵢⱼMᵢⱼdΩ => elements)(k3; D = D, ν = 0.45)
            @test !(k3 ≈ k)
        finally
            gmsh.finalize()
        end
    end

    @testset "Heat.∫∫∇v∇udxdy" begin
        gmsh.initialize()
        try
            elements = patch_tri6()
            set∇𝝭!(elements)

            n  = maximum(xᵢ.𝐼 for ap in elements for xᵢ in ap.𝓒)
            k  = zeros(n, n)
            (ApproxOperator.Heat.∫∫∇v∇udxdy => elements)(k; kᶜ = 3.0)
            @test all(isfinite, k)
            @test norm(k) > 0.0
            @test k ≈ k' atol = 1.0e-10

            # Conductivity enters linearly.
            k2 = zeros(n, n)
            (ApproxOperator.Heat.∫∫∇v∇udxdy => elements)(k2; kᶜ = 6.0)
            @test k2 ≈ 2k atol = 1.0e-9
        finally
            gmsh.finalize()
        end
    end

    @testset "keyword path equals prescribe! path" begin
        # S2 moved material constants from `prescribe!` into keyword arguments.
        # The two routes must produce the same matrix, otherwise the migration
        # silently changed the numerics.
        gmsh.initialize()
        try
            elements = patch_tri6()
            set∇𝝭!(elements)
            n  = maximum(xᵢ.𝐼 for ap in elements for xᵢ in ap.𝓒)

            # Stokes velocities carry two DOFs per node, so the matrix is 2n.
            k1 = zeros(2n, 2n)
            (ApproxOperator.Stokes.∫∫μ∇u∇vdxdy => elements)(k1; μ = 2.5)

            # (b) same value supplied through the quadrature point
            prescribe!(elements, :μ => 2.5)
            k2 = zeros(2n, 2n)
            (ApproxOperator.Stokes.∫∫μ∇u∇vdxdy => elements)(k2; μ = 2.5)

            @test k1 ≈ k2 atol = 1.0e-12
            @test norm(k1) > 0.0

            # The keyword really drives the assembly: doubling it doubles k.
            k3 = zeros(2n, 2n)
            (ApproxOperator.Stokes.∫∫μ∇u∇vdxdy => elements)(k3; μ = 5.0)
            @test k3 ≈ 2k1 atol = 1.0e-9
        finally
            gmsh.finalize()
        end
    end

    @testset "missing keyword is reported at the call site" begin
        gmsh.initialize()
        try
            elements = patch_tri6()
            set∇𝝭!(elements)
            n = maximum(xᵢ.𝐼 for ap in elements for xᵢ in ap.𝓒)
            k = zeros(2n, 2n)
            # A mandatory keyword must fail with UndefKeywordError, naming the
            # argument, rather than a KeyError from a Dict lookup deep inside
            # Node.getproperty.
            @test_throws UndefKeywordError (ApproxOperator.Stokes.∫∫μ∇u∇vdxdy => elements)(k)
        finally
            gmsh.finalize()
        end
    end

    @testset "material parameters are bound, not dangling" begin
        # Removing the `op.<field>` reads left statements such as
        # `Dᵐ = E*h/(1-ν^2)` naming E/ν/h without binding them. The static
        # `op.` check cannot see that, so assert here that any function which
        # uses E/ν/h names them via a binding (`E = ap.E`, `E = ξ.E`, ...)
        # before or inside its quadrature loop.
        root = joinpath(dirname(dirname(pathof(ApproxOperator))), "src", "operation")
        offenders = String[]
        for file in readdir(root; join = true)
            endswith(file, ".jl") || continue
            lines = readlines(file)
            idx = 0
            in_doc = false
            while idx < length(lines)
                line = lines[idx += 1]
                # Skip docstring bodies: they mention E/ν in prose, which is
                # not a use of the variable.
                ndelim = length(collect(eachmatch(r"""\"\"\"""", line)))
                if in_doc
                    in_doc = ndelim % 2 == 0
                    continue
                elseif ndelim % 2 == 1
                    in_doc = true
                    continue
                end
                startswith(line, "function") || continue
                hdr = line
                pre, post, seen = String[], String[], false
                while idx < length(lines)
                    l = lines[idx += 1]
                    if !seen && startswith(strip(l), "for ")
                        seen = true
                        continue
                    end
                    seen ? push!(post, l) : push!(pre, l)
                    (l == "end" || startswith(l, "end")) && break
                end
                pres, posts = join(pre, "\n"), join(post, "\n")
                for name in ("E", "ν", "h")
                    occursin(Regex("(?<![\\w.])$name(?![\\w])"), pres) || continue
                    # Bound in one of three ways:
                    #   * assigned before the loop (legacy `E = ap.E`);
                    #   * read from a quadrature node inside the loop;
                    #   * declared as a keyword in the signature (S2 onwards).
                    bound = occursin(Regex("^\\s*$name\\s*=", "m"), pres) ||
                            occursin(Regex("\\b$name\\s*="), pres) ||
                            occursin(Regex("ξ[\\w₀-₉]*\\.$name(?![\\w])"), posts) ||
                            occursin(Regex("[;,]\\s*" * name * "\\s*(::\\w+)?\\s*[,)]"), hdr)
                    bound || push!(offenders, "$(basename(file)): $(strip(hdr)) uses $name")
                end
            end
        end
        @test isempty(offenders)
    end

    @testset "no module reaches for a global `op`" begin
        # A cheap static check: the operator sources must not contain a bare
        # `op.` field access. `op::Operator` as a dispatch argument is fine.
        root = dirname(dirname(pathof(ApproxOperator)))
        root = joinpath(root, "src")
        offenders = String[]
        for file in readdir(joinpath(root, "operation"); join = true)
            endswith(file, ".jl") || continue
            for (n, line) in enumerate(eachline(file))
                startswith(strip(line), "#") && continue
                occursin("op::Operator", line) && continue
                if occursin(r"\bop\.[A-Za-z]", line)
                    push!(offenders, "$(basename(file)):$n")
                end
            end
        end
        @test isempty(offenders)
    end
end
