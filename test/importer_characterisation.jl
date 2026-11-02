# Characterisation tests for the mesh importer.
#
# These tests pin the *current* output of `GmshImport.getElements` so that the
# S3 rewrite (which replaces the `@eval`/`quote` machinery with a plain
# function) can be proved behaviour-preserving. They deliberately record
# observables rather than assert "correctness":
#
#   * how many elements come back,
#   * which keys the shared `data` dictionary holds after the call,
#   * the sum of the integration weights (`Σ ξ.𝑤`), which changes if the
#     quadrature rule or the Jacobian handling changes,
#   * the quadrature-point coordinates `ξ.x`, `ξ.y`, `ξ.z`.
#
# A change in any of these is a signal that the rewrite altered behaviour —
# either a bug, or a deliberate change that must be called out.

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements, getCurvedElements,
                                  addCentroidNodes!, corner_count
using Test
import Gmsh: gmsh

# ---------------------------------------------------------------------------
# Mesh under test: a unit square of 2x2 Quad4 elements.
# ---------------------------------------------------------------------------
function characterisation_mesh()
    gmsh.model.add("characterisation")
    p = [gmsh.model.geo.addPoint(x, y, 0.0)
         for (x,y) in ((0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0,1.0))]
    l = [gmsh.model.geo.addLine(p[1], p[2]),
         gmsh.model.geo.addLine(p[2], p[3]),
         gmsh.model.geo.addLine(p[3], p[4]),
         gmsh.model.geo.addLine(p[4], p[1])]
    cl = gmsh.model.geo.addCurveLoop(l)
    s  = gmsh.model.geo.addPlaneSurface([cl])
    gmsh.model.geo.synchronize()
    for li in l
        gmsh.model.mesh.setTransfiniteCurve(li, 3)
    end
    gmsh.model.mesh.setTransfiniteSurface(s)
    gmsh.model.mesh.setRecombine(2, s)
    gmsh.model.addPhysicalGroup(2, [s], 1)
    gmsh.model.setPhysicalName(2, 1, "Ω")
    gmsh.model.addPhysicalGroup(1, l, 2)
    gmsh.model.setPhysicalName(1, 2, "Γ")
    gmsh.model.mesh.generate(2)
end

# The observables that the rewrite must preserve.
#
# `𝑔` and `:ξ` are included deliberately: during the S3 port the `𝑔` counter
# was briefly advanced across elements (it must stay 1 for every element, since
# it indexes the quadrature point *within* one element) and `:ξ`/`:η` stopped
# being written at all. Both broke `set∇𝝭!` in ways the element count alone
# could not see.
function signature(elements)
    isempty(elements) && return (n=0, keys=Symbol[], w=0.0, centroid=(0.0,0.0))
    data = getfield(elements[1].𝓖[1], :data)
    w = sum(sum(x.𝑤) for ap in elements for x in ap.𝓖)
    xs = [x.x for ap in elements for x in ap.𝓖]
    ys = [x.y for ap in elements for x in ap.𝓖]
    return (
        n        = length(elements),
        keys     = sort(collect(keys(data))),
        w        = w,
        centroid = (sum(xs)/length(xs), sum(ys)/length(ys)),
        g_first  = elements[1].𝓖[1].𝑔,
        g_last   = elements[end].𝓖[end].𝑔,
        ξ_len    = length(data[:ξ][2]),
        x_len    = length(data[:x][2]),
        η_len    = haskey(data, :η) ? length(data[:η][2]) : 0,
    )
end

@testset "importer characterisation" begin

    @testset "Quad4, Gmsh quadrature" begin
        gmsh.initialize()
        try
            characterisation_mesh()
            entities = getPhysicalGroups()
            nodes    = get𝑿ᵢ()
            elms     = getElements(nodes, entities["Ω"], 2)
            s        = signature(elms)

            @test s.n == 4
            @test s.keys == [
                :w, :x, :y, :z, :η, :ξ,
                :∂η∂x, :∂η∂y, :∂ξ∂x, :∂ξ∂y,
                :𝐽, :𝑤,
            ]
            # Gmsh's quadrature weights sum to the element area; the four
            # 0.25-area quads give 1.0 in total. Any change here means the
            # Jacobian or the weight convention moved.
            @test s.w ≈ 1.0 atol = 1.0e-12
            @test s.centroid[1] ≈ 0.5 atol = 1.0e-12
            @test s.centroid[2] ≈ 0.5 atol = 1.0e-12

            # `:𝑔` indexes the quadrature point within one element, so it runs
            # 1:ng afresh for every element rather than being a global counter.
            @test s.g_first == 1
            @test s.g_last == 4          # Quad4 with a 4-point Gmsh rule
            # Local coordinates are recorded once per quadrature point of the
            # entity (not once per element), matching the historical layout:
            # `:ξ` has ng entries while `:x` has ng*ne.
            @test s.ξ_len == 4
            @test s.η_len == 4
            @test s.x_len == 16
        finally
            gmsh.finalize()
        end
    end

    @testset "Quad4, explicit quadrature" begin
        gmsh.initialize()
        try
            characterisation_mesh()
            entities = getPhysicalGroups()
            nodes    = get𝑿ᵢ()
            # One-point rule at the element centre.
            # A one-point rule at the element centre. The weights are the raw
            # quadrature weights, so sum(w) = 4 (one per element) while the
            # Jacobian-weighted `𝑤` sums to the area, 1.0.
            elms = getElements(nodes, entities["Ω"], ([0.0, 0.0, 0.0], [4.0]))
            s    = signature(elms)

            @test s.n == 4
            @test s.w ≈ 1.0 atol = 1.0e-12
        finally
            gmsh.finalize()
        end
    end

    @testset "explicit element type" begin
        gmsh.initialize()
        try
            characterisation_mesh()
            entities = getPhysicalGroups()
            nodes    = get𝑿ᵢ()
            elms     = getElements(nodes, entities["Ω"],
                                   ApproxOperator.Element{:Quad4}, 2)
            @test all(typeof(a) == ApproxOperator.Element{:Quad4} for a in elms)
            @test length(elms) == 4
        finally
            gmsh.finalize()
        end
    end

    @testset "1-D boundary entities" begin
        gmsh.initialize()
        try
            characterisation_mesh()
            entities = getPhysicalGroups()
            nodes    = get𝑿ᵢ()
            elms     = getElements(nodes, entities["Γ"], 2)
            @test all(typeof(a) == ApproxOperator.Element{:Seg2} for a in elms)
            # Eight boundary segments of the 2x2 quad mesh, total length 4.
            @test length(elms) == 8
            s = signature(elms)
            @test s.w ≈ 4.0 atol = 1.0e-12
        finally
            gmsh.finalize()
        end
    end

    @testset "normal flag adds the normal fields" begin
        gmsh.initialize()
        try
            characterisation_mesh()
            entities = getPhysicalGroups()
            nodes    = get𝑿ᵢ()
            plain    = getElements(nodes, entities["Γ"], 2)
            with_n   = getElements(nodes, entities["Γ"], 2; normal = true)
            kplain   = signature(plain).keys
            kwithn   = signature(with_n).keys
            @test :n₁ ∉ kplain
            @test :n₁ ∈ kwithn
            @test :s₁ ∈ kwithn
        finally
            gmsh.finalize()
        end
    end

    @testset "getCurvedElements uses cs.𝐽, not the flat determinant" begin
        # `cs` carries the surface as *fields*, so a NamedTuple of functions is
        # the documented form. A flat map (cs.𝐽 = 1) must reproduce the flat
        # area; scaling the map must scale every weight.
        flat_map(scale) = (
            𝐽   = x -> scale,
            a₁₁ = x -> 1.0, a₁₂ = x -> 0.0, a₂₂ = x -> 1.0,
            a¹¹ = x -> 1.0, a¹² = x -> 0.0, a²² = x -> 1.0,
            𝒂₁  = x -> (1.0, 0.0, 0.0),
            𝒂₂  = x -> (0.0, 1.0, 0.0),
            𝒂₃  = x -> (0.0, 0.0, 1.0),
        )

        gmsh.initialize()
        try
            characterisation_mesh()                 # 2x2 Quad4 on the unit square
            entities = getPhysicalGroups()
            nodes    = get𝑿ᵢ()

            elms = getCurvedElements(nodes, entities["Ω"], flat_map(1.0), 2)
            @test length(elms) == 4
            @test all(typeof(a) == ApproxOperator.Element{:Quad4} for a in elms)

            # The map is flat, so the total area is the unit square.
            total = sum(sum(x.𝑤) for ap in elms for x in ap.𝓖)
            @test total ≈ 1.0 atol = 1.0e-12

            # `cs.𝐽` really drives the weight: doubling the area scale doubles
            # every quadrature weight. If `fill_curved_coordinates!` were still
            # using Gmsh's flat determinant this would stay at 1.0.
            doubled = getCurvedElements(nodes, entities["Ω"], flat_map(2.0), 2)
            total2  = sum(sum(x.𝑤) for ap in doubled for x in ap.𝓖)
            @test total2 ≈ 2.0 atol = 1.0e-12
            @test total2 ≈ 2total atol = 1.0e-12

            # The per-point determinant recorded in `:𝐽` is Gmsh's flat one; the
            # curvature correction lives in `:𝑤` alone, which is the documented
            # contract.
            d = getfield(elms[1].𝓖[1], :data)
            @test :𝐽 ∈ keys(d)
            @test :𝑤 ∈ keys(d)
        finally
            gmsh.finalize()
        end
    end

    @testset "addCentroidNodes! appends one node per element" begin
        gmsh.initialize()
        try
            characterisation_mesh()                 # 2x2 Quad4 on the unit square
            entities = getPhysicalGroups()
            nodes    = get𝑿ᵢ()
            n0       = length(nodes)
            𝐼_before = [x.𝐼 for x in nodes]

            addCentroidNodes!(nodes, entities["Ω"])

            @test length(nodes) == n0 + 4
            # The original nodes are untouched and keep their 𝐼.
            @test [nodes[i].𝐼 for i in 1:n0] == 𝐼_before
            # New nodes get consecutive 𝐼 values after the existing ones.
            @test [nodes[i].𝐼 for i in n0+1:length(nodes)] == n0+1:n0+4
            # ... and share the same `data` dictionary, so `prescribe!` and
            # `getElements` keep working on the extended list.
            @test getfield(nodes[end], :data) === getfield(nodes[1], :data)

            # The four centroids of a 2x2 quad mesh.
            got = sort([(nodes[i].x, nodes[i].y) for i in n0+1:length(nodes)])
            expected = [(0.25,0.25), (0.25,0.75), (0.75,0.25), (0.75,0.75)]
            for (g, e) in zip(got, expected)
                @test g[1] ≈ e[1] atol = 1.0e-10
                @test g[2] ≈ e[2] atol = 1.0e-10
            end

            # Elements can still be built from the extended node list.
            elms = getElements(nodes, entities["Ω"], 2)
            @test length(elms) == 4
        finally
            gmsh.finalize()
        end
    end

    @testset "addCentroidNodes! uses the corner mean for quadratic elements" begin
        gmsh.initialize()
        try
            # A single order-2 quad: the mean of *all* its nodes would be pulled
            # towards the mid-side points, so the centroid must average only the
            # four corner nodes.
            gmsh.model.add("quadratic")
            p = [gmsh.model.geo.addPoint(x, y, 0.0)
                 for (x,y) in ((0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0,1.0))]
            l = [gmsh.model.geo.addLine(p[1],p[2]), gmsh.model.geo.addLine(p[2],p[3]),
                 gmsh.model.geo.addLine(p[3],p[4]), gmsh.model.geo.addLine(p[4],p[1])]
            cl = gmsh.model.geo.addCurveLoop(l)
            s  = gmsh.model.geo.addPlaneSurface([cl])
            gmsh.model.geo.synchronize()
            for li in l; gmsh.model.mesh.setTransfiniteCurve(li, 2); end
            gmsh.model.mesh.setTransfiniteSurface(s)
            gmsh.model.mesh.setRecombine(2, s)
            gmsh.model.addPhysicalGroup(2, [s], 1)
            gmsh.model.setPhysicalName(2, 1, "Ω")
            gmsh.model.mesh.generate(2)
            gmsh.model.mesh.setOrder(2)

            entities = getPhysicalGroups()
            nodes    = get𝑿ᵢ()
            n0       = length(nodes)
            addCentroidNodes!(nodes, entities["Ω"])

            @test length(nodes) == n0 + 1
            @test nodes[end].x ≈ 0.5 atol = 1.0e-10
            @test nodes[end].y ≈ 0.5 atol = 1.0e-10
        finally
            gmsh.finalize()
        end
    end

    @testset "get𝑿ᵢ(centroid = …) matches get𝑿ᵢ() + addCentroidNodes!" begin
        gmsh.initialize()
        try
            characterisation_mesh()                 # 2x2 Quad4 on the unit square
            entities = getPhysicalGroups()

            # Default: mesh nodes only. This must stay the behaviour every
            # existing script relies on.
            plain = get𝑿ᵢ()
            @test length(plain) == 9

            # Keyword form.
            withc = get𝑿ᵢ(centroid = entities["Ω"])
            @test length(withc) == 9 + 4

            # Explicit two-step form must agree exactly.
            stepwise = get𝑿ᵢ()
            addCentroidNodes!(stepwise, entities["Ω"])
            @test length(stepwise) == length(withc)
            for (a, b) in zip(withc, stepwise)
                @test a.𝐼 == b.𝐼
                @test a.x ≈ b.x atol = 1.0e-12
                @test a.y ≈ b.y atol = 1.0e-12
            end

            # Mesh nodes keep their 𝐼, so code that indexes `nodes` is unaffected.
            @test all(withc[i].𝐼 == i for i in 1:9)
            # One shared dictionary, as `prescribe!` and `getElements` require.
            @test getfield(withc[end], :data) === getfield(withc[1], :data)

            # Restricting to an element type that the entity does not have
            # yields no extra nodes.
            tri = get𝑿ᵢ(centroid = entities["Ω"], elementType = 2)   # Tri3
            @test length(tri) == 9
        finally
            gmsh.finalize()
        end
    end

    @testset "corner_count" begin
        @test corner_count(1)  == 2    # Seg2
        @test corner_count(2)  == 3    # Tri3
        @test corner_count(9)  == 3    # Tri6  — corner nodes only
        @test corner_count(11) == 4    # Tet10 — corner nodes only
        @test corner_count(5)  == 8    # Hex8
        @test corner_count(16) == 4    # Quad8
        @test corner_count(99) == -1   # unknown: caller averages all nodes
    end

    @testset "data dictionary is shared across all elements" begin
        gmsh.initialize()
        try
            characterisation_mesh()
            entities = getPhysicalGroups()
            nodes    = get𝑿ᵢ()
            elms     = getElements(nodes, entities["Ω"], 2)
            # One `data` dictionary serves every quadrature point; the S3 rewrite
            # must keep this invariant, since `push!(elements, …)` and
            # `prescribe!` both rely on it.
            first_data = getfield(elms[1].𝓖[1], :data)
            @test all(getfield(x, :data) === first_data for ap in elms for x in ap.𝓖)
        finally
            gmsh.finalize()
        end
    end
end
