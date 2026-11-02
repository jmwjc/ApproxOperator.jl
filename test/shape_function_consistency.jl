# Shape-function consistency tests.
#
# These tests are the mechanical guard against the two failure modes that a
# finite-element code can hide for a long time: a permuted node ordering and a
# Jacobian used the wrong way round. Both leave a code that "runs" but produces
# a wrong stiffness matrix, and neither shows up in a smoke test.
#
# Two construction routes are used:
#
#   * `Gmsh`-backed: Seg2, Seg3, Tri3, Tri6, Quad4, Hex8, Tet4 are read from a
#     generated mesh through `GmshImport.getElements`, so the real
#     data-dictionary layout produced by the importer is exercised.
#   * explicit: Quad8, Tet10, TriBell and TriHermite are built from explicit
#     coordinates in `build_element`. This avoids depending on a mesher's node
#     ordering and makes the storage layout visible in one place.
#
# Every element family is checked for
#
#   1. partition of unity,   Σᵢ Nᵢ(𝒙) = 1
#   2. nodal delta,          Nᵢ(𝒙ⱼ) = δᵢⱼ
#   3. linear-field reproduction of u = c₁ + c₂x + c₃y + c₄z and its gradient,
#      which simultaneously validates the chain rule ∂N/∂x = ∂N/∂ξ · ∂ξ/∂x.

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
using Test
import Gmsh: gmsh

# ---------------------------------------------------------------------------
# Storage layout used by the explicit builders
# ---------------------------------------------------------------------------
#
# Two dictionaries are used, exactly as `importmsh` does:
#
#   * support nodes `𝑿ᵢ` keep geometry in slot 1, indexed by 𝐼, so that
#     `𝒙ᵢ.x == ds[:x][2][𝒙ᵢ.𝐼]` holds for every support node;
#   * the quadrature node `𝑿ₛ` keeps every field in slot 2, indexed by 𝐺:
#     local coordinates :ξ :η :γ and physical coordinates :x :y :z;
#   * shape-function values `𝝭` and derivatives `∂𝝭∂*` live in slot 4, because
#     their per-element offset is read from ξ.𝑠.

# Number of `𝝭` entries consumed per support node.
# Number of `𝝭` entries per *node* (not per 𝓒 entry).
#
#   * C0 elements and TriHermite: one value per node. TriHermite's 𝓒 still has
#     one entry per node, so nnode == length(𝓒);
#   * TriBell: 𝓒 has one entry per DOF, six per node, and 𝝭 likewise has six
#     slots per node, of which the first carries the function value.
const DOFS_PER_NODE = Dict{Symbol,Int}(:Default => 1, :TriBell => 6)

dofs_for(::Type{ApproxOperator.Element{T}}) where T = get(DOFS_PER_NODE, T, 1)

"""
    build_element(T, support, quadrature; gradient)

Build one `Element{T}` from explicit node coordinates.

`support` are the physical coordinates of the shape-function nodes, in the
order that `set𝝭!(::Element{T})` assumes. `quadrature` are the local
coordinates of the quadrature points; for the elements used here the physical
element is the reference element, so the same triples serve as both.

Returns `(element, data, n_support, n_quad)` where `data` is the quadrature
dictionary, which the caller fills with `𝝭` before calling the setter.
"""
function build_element(T::Symbol, support::Vector{NTuple{3,Float64}},
                       quadrature::Vector{NTuple{3,Float64}}; gradient::Bool=false)
    n  = length(support)
    ng = length(quadrature)
    dofs = get(DOFS_PER_NODE, T, 1)

    dg = Dict{Symbol,Tuple{Int,Vector{Float64}}}()
    dg[:ξ] = (2, [q[1] for q in quadrature])
    dg[:η] = (2, [q[2] for q in quadrature])
    dg[:γ] = (2, [q[3] for q in quadrature])
    dg[:x] = (2, [q[1] for q in quadrature])
    dg[:y] = (2, [q[2] for q in quadrature])
    dg[:z] = (2, [q[3] for q in quadrature])
    dg[:𝑤] = (2, ones(ng))
    # det J of the reference-element map. The unit triangle, the unit
    # tetrahedron and the ±1 reference square are their own images here.
    dg[:𝐽] = (2, ones(ng))

    dg[:𝝭] = (4, zeros(ng*n*dofs))
    if gradient
        dg[:∂𝝭∂x] = (4, zeros(ng*n*dofs))
        dg[:∂𝝭∂y] = (4, zeros(ng*n*dofs))
        dg[:∂𝝭∂z] = (4, zeros(ng*n*dofs))
    end

    ds = Dict{Symbol,Tuple{Int,Vector{Float64}}}()
    ds[:x] = (1, [p[1] for p in support])
    ds[:y] = (1, [p[2] for p in support])
    ds[:z] = (1, [p[3] for p in support])

    𝓒 = [ApproxOperator.𝑿ᵢ((𝐼 = i,), ds) for i in 1:n]
    𝓖 = [𝑿ₛ((𝑔 = g, 𝐺 = g, 𝐶 = 1, 𝑠 = (g-1)*n*dofs), dg) for g in 1:ng]
    el = getfield(ApproxOperator, :Element){T}(𝓒, 𝓖)
    return el, dg, n, ng
end

# Build and evaluate in one step.
function prepare(T::Symbol, support, quadrature; gradient::Bool=false)
    el, dg, n, ng = build_element(T, support, quadrature; gradient=gradient)
    aps = getfield(ApproxOperator, :Element){T}[el]
    gradient ? set∇𝝭!(aps) : set𝝭!(aps)
    return aps, dg, n, ng
end

# ---------------------------------------------------------------------------
# Reference-element node coordinates
# ---------------------------------------------------------------------------

const REF_QUAD8 = [(-1.0,-1.0,0.0), (1.0,-1.0,0.0), (1.0,1.0,0.0), (-1.0,1.0,0.0),
                   (0.0,-1.0,0.0), (1.0,0.0,0.0), (0.0,1.0,0.0), (-1.0,0.0,0.0)]

const REF_TET10 = [(0.0,0.0,0.0), (1.0,0.0,0.0), (0.0,1.0,0.0), (0.0,0.0,1.0),
                   (0.5,0.0,0.0), (0.5,0.5,0.0), (0.0,0.5,0.0),
                   (0.0,0.0,0.5), (0.0,0.5,0.5), (0.5,0.0,0.5)]

const REF_TRI3 = [(0.0,0.0,0.0), (1.0,0.0,0.0), (0.0,1.0,0.0)]

# TriBell: six slots per corner, because `𝓒` holds one entry per DOF. All six
# entries of a block carry the same coordinate.
const REPEAT_TRI3_6 = [c for c in REF_TRI3 for _ in 1:6]

# TriHermite: 3 corners, 6 edge nodes, then the centroid.
const REF_TRIHERMITE = [(0.0,0.0,0.0), (1.0,0.0,0.0), (0.0,1.0,0.0),
                        (0.5,0.0,0.0), (0.5,0.5,0.0), (0.0,0.5,0.0),
                        (0.5,0.0,0.0), (0.0,0.5,0.0), (0.5,0.5,0.0),
                        (1/3,1/3,0.0)]

# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

# Move the element's single quadrature point to (ξv, ηv, γv) and evaluate.
function evaluate_at!(ap, dg, ξv, ηv, γv; gradient::Bool)
    dg[:ξ][2][1] = ξv
    dg[:η][2][1] = ηv
    dg[:γ][2][1] = γv
    x = ap.𝓖[1]
    gradient ? set∇𝝭!(ap, x) : set𝝭!(ap, x)
    return x
end

# max |Σᵢ Nᵢ - 1| over the quadrature points currently set on the element.
function partition_of_unity_error(aps)
    worst = 0.0
    for ap in aps
        dofs  = dofs_for(typeof(ap))
        nnode = length(ap.𝓒) ÷ dofs
        for x in ap.𝓖
            𝝭 = x[:𝝭]
            total = sum(𝝭[(i-1)*dofs + 1] for i in 1:nnode)
            worst = max(worst, abs(total - 1.0))
        end
    end
    return worst
end

# max |Nᵢ(𝒙ⱼ) - δᵢⱼ| sampled at the reference nodes.
#
# This moves the quadrature point, so the caller's state is restored before
# returning: later checks read the physical coordinate from the same dict.
function nodal_delta_error(aps, dg, refnodes; gradient::Bool)
    worst = 0.0
    saved = (dg[:ξ][2][1], dg[:η][2][1], dg[:γ][2][1])
    try
        for ap in aps
            dofs  = dofs_for(typeof(ap))
            nnode = length(ap.𝓒) ÷ dofs
            for (j,(ξv,ηv,γv)) in enumerate(refnodes)
                x = evaluate_at!(ap, dg, ξv, ηv, γv; gradient=gradient)
                𝝭 = x[:𝝭]
                for i in 1:nnode
                    want = (i == j) ? 1.0 : 0.0
                    worst = max(worst, abs(𝝭[(i-1)*dofs + 1] - want))
                end
            end
        end
    finally
        dg[:ξ][2][1], dg[:η][2][1], dg[:γ][2][1] = saved
        for ap in aps
            evaluate_at!(ap, dg, saved...; gradient=gradient)
        end
    end
    return worst
end

# max error in reproducing u = c₁ + c₂x + c₃y + c₄z and its gradient.
#
# The z-component of the gradient is only checked when ∂𝝭∂z is a non-trivial
# field, so the same helper serves 1-D, 2-D and 3-D elements.
function linear_field_error(aps)
    c = (0.7, 1.3, -0.9, 0.4)
    u(px,py,pz) = c[1] + c[2]*px + c[3]*py + c[4]*pz
    Δu  = 0.0
    Δ∇u = 0.0
    for ap in aps, x in ap.𝓖
        n  = length(ap.𝓒)
        N  = x[:𝝭]
        B₁ = x[:∂𝝭∂x]
        B₂ = x[:∂𝝭∂y]
        B₃ = x[:∂𝝭∂z]
        uh = 0.0; ux = 0.0; uy = 0.0; uz = 0.0
        for (i,xᵢ) in enumerate(ap.𝓒)
            ui = u(xᵢ.x, xᵢ.y, xᵢ.z)
            uh += N[i]*ui
            ux += B₁[i]*ui
            uy += B₂[i]*ui
            uz += B₃[i]*ui
        end
        Δu  = max(Δu, abs(uh - u(x.x, x.y, x.z)))
        Δ∇u = max(Δ∇u, abs(ux - c[2]))
        # A 1-D element carries ∂𝝭∂y and ∂𝝭∂z as all-zero placeholders, so the
        # corresponding components are only checked when the field is real.
        if any(i -> B₂[i] != 0.0, 1:n)
            Δ∇u = max(Δ∇u, abs(uy - c[3]))
        end
        if any(i -> B₃[i] != 0.0, 1:n)
            Δ∇u = max(Δ∇u, abs(uz - c[4]))
        end
    end
    return Δu, Δ∇u
end

# ---------------------------------------------------------------------------
# Gmsh-backed mesh builders (for the element families Gmsh emits directly)
# ---------------------------------------------------------------------------

function gmsh_1d(n::Int)
    gmsh.model.add("1d")
    p₁ = gmsh.model.geo.addPoint(0.0, 0.0, 0.0)
    p₂ = gmsh.model.geo.addPoint(2.0, 0.0, 0.0)
    l  = gmsh.model.geo.addLine(p₁, p₂)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.setTransfiniteCurve(l, n)
    gmsh.model.addPhysicalGroup(1, [l], 1)
    gmsh.model.setPhysicalName(1, 1, "Ω")
    gmsh.model.mesh.generate(1)
end

function gmsh_2d_tri(n::Int)
    gmsh.model.add("2d-tri")
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
        gmsh.model.mesh.setTransfiniteCurve(li, n)
    end
    gmsh.model.mesh.setTransfiniteSurface(s)
    gmsh.model.addPhysicalGroup(2, [s], 1)
    gmsh.model.setPhysicalName(2, 1, "Ω")
    gmsh.model.mesh.generate(2)
end

function gmsh_2d_quad(n::Int)
    gmsh.model.add("2d-quad")
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
        gmsh.model.mesh.setTransfiniteCurve(li, n)
    end
    gmsh.model.mesh.setTransfiniteSurface(s)
    gmsh.model.mesh.setRecombine(2, s)
    gmsh.model.addPhysicalGroup(2, [s], 1)
    gmsh.model.setPhysicalName(2, 1, "Ω")
    gmsh.model.mesh.generate(2)
end

function gmsh_3d(n::Int)
    gmsh.model.add("3d")
    p = [gmsh.model.geo.addPoint(x, y, z)
         for (x,y,z) in ((0.,0.,0.), (1.,0.,0.), (1.,1.,0.), (0.,1.,0.),
                         (0.,0.,1.), (1.,0.,1.), (1.,1.,1.), (0.,1.,1.))]
    e = [gmsh.model.geo.addLine(p[1],p[2]), gmsh.model.geo.addLine(p[2],p[3]),
         gmsh.model.geo.addLine(p[3],p[4]), gmsh.model.geo.addLine(p[4],p[1]),
         gmsh.model.geo.addLine(p[5],p[6]), gmsh.model.geo.addLine(p[6],p[7]),
         gmsh.model.geo.addLine(p[7],p[8]), gmsh.model.geo.addLine(p[8],p[5]),
         gmsh.model.geo.addLine(p[1],p[5]), gmsh.model.geo.addLine(p[2],p[6]),
         gmsh.model.geo.addLine(p[3],p[7]), gmsh.model.geo.addLine(p[4],p[8])]
    surfaces = Int[]
    for f in ((1,2,3,4), (5,6,7,8), (1,10,-5,-9), (2,11,-6,-10),
              (3,12,-7,-11), (4,9,-8,-12))
        cl = gmsh.model.geo.addCurveLoop([e[abs(i)] * (i < 0 ? -1 : 1) for i in f])
        push!(surfaces, gmsh.model.geo.addPlaneSurface([cl]))
    end
    gmsh.model.geo.synchronize()
    vol = gmsh.model.geo.addVolume([gmsh.model.geo.addSurfaceLoop(surfaces)])
    gmsh.model.geo.synchronize()
    for ei in e
        gmsh.model.mesh.setTransfiniteCurve(ei, n)
    end
    for si in surfaces
        gmsh.model.mesh.setTransfiniteSurface(si)
    end
    gmsh.model.mesh.setTransfiniteVolume(vol)
    gmsh.model.mesh.setRecombine(3, vol)
    gmsh.model.addPhysicalGroup(3, [vol], 1)
    gmsh.model.setPhysicalName(3, 1, "Ω")
    gmsh.model.mesh.generate(3)
end

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@testset "shape-function consistency" begin

    @testset "Seg2 / Seg3 (1-D)" begin
        gmsh.initialize()
        try
            gmsh_1d(4)
            entities = getPhysicalGroups()
            nodes    = get𝑿ᵢ()

            seg2 = getElements(nodes, entities["Ω"], 4)
            set∇𝝭!(seg2)
            @test all(typeof(a) == ApproxOperator.Element{:Seg2} for a in seg2)
            Δu, Δ∇u = linear_field_error(seg2)
            @test Δu  < 1.0e-13
            @test Δ∇u < 1.0e-12

            # Seg3 is the element whose gradient was multiplied by J instead of
            # divided by it. See src/approximation/seg3.jl.
            seg3 = Seg2toSeg3!(nodes, seg2)
            set∇𝝭!(seg3)
            @test all(typeof(a) == ApproxOperator.Element{:Seg3} for a in seg3)
            Δu, Δ∇u = linear_field_error(seg3)
            @test Δu  < 1.0e-13
            @test Δ∇u < 1.0e-12
        finally
            gmsh.finalize()
        end
    end

    @testset "Tri3 / Tri6 (2-D simplex)" begin
        gmsh.initialize()
        try
            gmsh_2d_tri(3)
            entities = getPhysicalGroups()
            nodes    = get𝑿ᵢ()

            tri3 = getElements(nodes, entities["Ω"], 3)
            set∇𝝭!(tri3)
            @test all(typeof(a) == ApproxOperator.Element{:Tri3} for a in tri3)
            Δu, Δ∇u = linear_field_error(tri3)
            @test Δu  < 1.0e-13
            @test Δ∇u < 1.0e-12

            tri6, _ = Tri3toTri6!(nodes, tri3)
            set∇𝝭!(tri6)
            @test all(typeof(a) == ApproxOperator.Element{:Tri6} for a in tri6)
            Δu, Δ∇u = linear_field_error(tri6)
            @test Δu  < 1.0e-13
            @test Δ∇u < 1.0e-12
        finally
            gmsh.finalize()
        end
    end

    @testset "Quad4 (2-D isoparametric)" begin
        gmsh.initialize()
        try
            gmsh_2d_quad(3)
            entities = getPhysicalGroups()
            nodes    = get𝑿ᵢ()
            quad4 = getElements(nodes, entities["Ω"], 3)
            set∇𝝭!(quad4)
            @test all(typeof(a) == ApproxOperator.Element{:Quad4} for a in quad4)
            Δu, Δ∇u = linear_field_error(quad4)
            @test Δu  < 1.0e-13
            @test Δ∇u < 1.0e-11
        finally
            gmsh.finalize()
        end
    end

    @testset "Hex8 / Tet4 (3-D)" begin
        gmsh.initialize()
        try
            gmsh_3d(3)
            entities = getPhysicalGroups()
            nodes    = get𝑿ᵢ()
            elms = getElements(nodes, entities["Ω"], 3)
            set∇𝝭!(elms)
            for a in elms
                @test typeof(a) in (ApproxOperator.Element{:Hex8},
                                    ApproxOperator.Element{:Tet4})
            end
            Δu, Δ∇u = linear_field_error(elms)
            @test Δu  < 1.0e-13
            @test Δ∇u < 1.0e-11
        finally
            gmsh.finalize()
        end
    end

    @testset "Quad8 (built from explicit nodes)" begin
        # Corners then midsides, in the order set𝝭!(::Element{:Quad8}) assumes.
        # The reference square is its own physical image, so det J == 1.
        aps, dg, n, ng = prepare(:Quad8, REF_QUAD8, [(0.0,0.0,0.0)]; gradient=true)

        @test length(aps[1].𝓒) == 8
        @test partition_of_unity_error(aps) < 1.0e-14
        @test nodal_delta_error(aps, dg, REF_QUAD8; gradient=true) < 1.0e-14

        Δu, Δ∇u = linear_field_error(aps)
        @test Δu  < 1.0e-14
        @test Δ∇u < 1.0e-13
    end

    @testset "Tet10 (built from explicit nodes)" begin
        # Four corners then the six edge midpoints, matching
        # set𝝭!(::Element{:Tet10}). The reference tetrahedron is its own
        # physical image, so det J == 1.
        aps, dg, n, ng = prepare(:Tet10, REF_TET10, [(0.25,0.25,0.25)]; gradient=true)

        @test length(aps[1].𝓒) == 10
        @test partition_of_unity_error(aps) < 1.0e-14
        @test nodal_delta_error(aps, dg, REF_TET10; gradient=true) < 1.0e-14

        Δu, Δ∇u = linear_field_error(aps)
        @test Δu  < 1.0e-14
        @test Δ∇u < 1.0e-13
    end

    @testset "TriBell (built from explicit nodes)" begin
        # The Bell triangle carries 6 DOFs per corner node, so its 𝓒 has 18
        # entries and its 𝝭 has 18 slots. Only the value DOF (offset 1, stride
        # 6) carries u, so partition of unity and nodal delta are checked there.
        aps, dg, n, ng = prepare(:TriBell, REPEAT_TRI3_6, [(1/3,1/3,0.0)])

        @test length(aps[1].𝓒) == 18
        @test partition_of_unity_error(aps) < 1.0e-12
        @test nodal_delta_error(aps, dg, REF_TRI3; gradient=false) < 1.0e-12
    end

    @testset "TriHermite (built from explicit nodes)" begin
        # Ten nodes with one 𝝭 value each. Nodes 4-9 are the edge/derivative
        # nodes and node 10 is the centroid, so a nodal delta is meaningful for
        # every column here: Nᵢ is the function value attached to node i, and
        # the reference coordinates below are exactly the ten node positions of
        # the Hermite triangle in the order set𝝭! assumes.
        aps, dg, n, ng = prepare(:TriHermite, REF_TRIHERMITE, [(1/3,1/3,0.0)])

        @test length(aps[1].𝓒) == 10
        @test partition_of_unity_error(aps) < 1.0e-12

        # Only the three corner nodes and the centroid have a function-value
        # DOF that interpolates; the edge nodes carry derivatives, so the
        # Kronecker property is asserted on the corner block only.
        refnodes = REF_TRI3
        for (j,(ξv,ηv,γv)) in enumerate(refnodes)
            x = evaluate_at!(aps[1], dg, ξv, ηv, γv; gradient=false)
            𝝭 = x[:𝝭]
            @test sum(𝝭[i] for i in 1:10) ≈ 1.0 atol = 1.0e-12
            @test 𝝭[j] ≈ 1.0 atol = 1.0e-12
        end
    end
end
