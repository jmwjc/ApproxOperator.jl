module GmshImport

using ..ApproxOperator: AbstractElement, Element, Node, 𝑿ᵢ, 𝑿ₛ, SpatialPartition
include("coordinates.jl")
include("elements.jl")
using .MeshCoordinates
import Gmsh: gmsh

function getPhysicalGroups()
    entities = Dict{String,Pair{Int,Vector{Int}}}()
    dimTags = gmsh.model.getPhysicalGroups()
    for (dim,tag) in dimTags
        name = gmsh.model.getPhysicalName(dim,tag)
        tags = gmsh.model.getEntitiesForPhysicalGroup(dim,tag)
        entities[name] = dim=>tags
    end
    return entities
end

"""
    get𝑿ᵢ(; centroid = nothing, elementType = :auto) -> Vector{𝑿ᵢ}

Read the mesh nodes from Gmsh.

# Keyword arguments
- `centroid`    : `nothing` (the default) returns only the mesh nodes. Passing a
                  Gmsh `(dim, tags)` pair — typically `entities["Ω"]` — appends
                  one extra node at the centroid of every element of that
                  entity. See [`addCentroidNodes!`](@ref).
- `elementType` : when `centroid` is given, restrict the centroids to one Gmsh
                  element type; `:auto` processes every type Gmsh reports.

# Example
```julia
# Mesh nodes only.
nodes = get𝑿ᵢ()

# Mesh nodes plus one centroid per element of Ω.
nodes = get𝑿ᵢ(centroid = entities["Ω"])
```

The centroid nodes are appended after the mesh nodes, so every mesh node keeps
its `𝐼` and existing code that indexes `nodes` is unaffected. All nodes share
one `data` dictionary, which is what `prescribe!` and `getElements` rely on.
"""
function get𝑿ᵢ(; centroid::Union{Nothing,Pair{Int,Vector{Int}}} = nothing,
                 elementType::Union{Symbol,Integer} = :auto)
    nodeTags, coord = gmsh.model.mesh.getNodes()
    nₚ = length(nodeTags)
    x = zeros(nₚ)
    y = zeros(nₚ)
    z = zeros(nₚ)
    for (i,I) in enumerate(nodeTags)
        x[I] = coord[3*i-2]
        y[I] = coord[3*i-1]
        z[I] = coord[3*i]
    end
    data = Dict([:x=>(1,x),:y=>(1,y),:z=>(1,z)])
    nodes = [𝑿ᵢ((𝐼=i,), data) for i in 1:nₚ]
    centroid === nothing && return nodes
    return addCentroidNodes!(nodes, centroid; elementType = elementType)
end

"""
    addCentroidNodes!(nodes, dimTag; elementType = :auto) -> Vector{𝑿ᵢ}

Append one node at the **centroid of every element** of the Gmsh entity
`dimTag`, and return the extended node list.

The existing nodes are left in place and keep their `𝐼`; each appended node gets
the next free `𝐼` and the same `data` dictionary as the original nodes, so the
result can be passed straight to `getElements`.

# Arguments
- `nodes`       : node list, normally the result of [`get𝑿ᵢ`](@ref)
- `dimTag`      : Gmsh `(dim, tags)` pair, e.g. `entities["Ω"]`
- `elementType` : Gmsh element type to restrict to, or `:auto` to take every
                  type Gmsh reports for the entity

# Centroid definition
The centroid is the arithmetic mean of the element's **corner** nodes,

    𝒙_c = (1/n_c) Σ_c 𝒙_c

not of all its nodes. For a quadratic element the mean of all nodes is pulled
towards the mid-side points and is not the geometric centre; the corner mean is.

For a straight-sided element this is the exact geometric centroid.

This is the work horse behind `get𝑿ᵢ(centroid = dimTag)`; call it directly when
the nodes already exist and you only want to add centroids.

# Example
```julia
nodes = get𝑿ᵢ()
addCentroidNodes!(nodes, entities["Ω"])
elements = getElements(nodes, entities["Ω"], 2)

# equivalent, in one call:
nodes = get𝑿ᵢ(centroid = entities["Ω"])
```
"""
function addCentroidNodes!(nodes::Vector{N}, dimTag::Pair{Int,Vector{Int}};
                           elementType::Union{Symbol,Integer} = :auto) where N<:Node
    isempty(nodes) && return nodes
    data = getfield(nodes[1], :data)
    xs = data[:x][2]
    ys = data[:y][2]
    zs = data[:z][2]

    dim, tags = dimTag
    for tag in tags
        elementTypes, _, nodeTagLists = gmsh.model.mesh.getElements(dim, tag)
        for (etype, nodeTags) in zip(elementTypes, nodeTagLists)
            (elementType === :auto || etype == elementType) || continue
            _, _, _, ni = gmsh.model.mesh.getElementProperties(etype)
            nc = corner_count(etype)
            nc <= 0 && (nc = ni)      # unknown type: average all nodes
            for e in 1:div(length(nodeTags), ni)
                conn = nodeTags[ni*(e-1)+1 : ni*e]
                xc = sum(xs[conn[c]] for c in 1:nc) / nc
                yc = sum(ys[conn[c]] for c in 1:nc) / nc
                zc = sum(zs[conn[c]] for c in 1:nc) / nc
                push!(xs, xc)
                push!(ys, yc)
                push!(zs, zc)
                push!(nodes, 𝑿ᵢ((𝐼 = length(xs),), data))
            end
        end
    end
    return nodes
end

"""
    corner_count(elementType::Integer) -> Int

Number of corner nodes of a Gmsh element type — the nodes that span the element
geometry, as opposed to the mid-side, face and interior nodes of a higher-order
element.

Returns `-1` for a type this package does not know, so that
[`addCentroidNodes!`](@ref) can fall back to averaging all nodes.
"""
function corner_count(elementType::Integer)
    elementType == 1  && return 2   # Seg2
    elementType == 2  && return 3   # Tri3
    elementType == 3  && return 4   # Quad4
    elementType == 4  && return 4   # Tet4
    elementType == 5  && return 8   # Hex8
    elementType == 8  && return 2   # Seg3
    elementType == 9  && return 3   # Tri6
    elementType == 10 && return 4   # Quad9
    elementType == 11 && return 4   # Tet10
    elementType == 12 && return 8   # Hex27
    elementType == 15 && return 1   # Poi1
    elementType == 16 && return 4   # Quad8
    return -1
end


#
# Storage layout of the shared `data` dictionary
# ------------------------------------------------
# Several nodes (`𝑿ᵢ` for the shape-function support, `𝑿ₛ` for the quadrature
# points) point at the *same* `Dict{Symbol,Tuple{Int,Vector{Float64}}}`. Each
# entry is `name => (slot, values)`, where `slot` selects an entry of the node's
# index tuple and the value of that entry is the offset into `values`:
#
#     node.name  ==  data[name][2][ node.index[ data[name][1] ] ]
#
# Concretely, for `𝑿ᵢ` (index `(𝐼,)`) the geometry is stored with `slot = 1`,
# so `𝒙ᵢ.x == data[:x][2][𝒙ᵢ.𝐼]`. For `𝑿ₛ` (index `(𝑔,𝐺,𝐶,𝑠)`):
#
#     slot 1 -> local quadrature coordinate :ξ, :η, :γ, :w
#               (one value per quadrature point, indexed by ξ.𝑔)
#     slot 2 -> global quadrature-point counter (indexed by ξ.𝐺)
#     slot 3 -> global element counter       (indexed by ξ.𝐶)
#     slot 4 -> offset into the element DOF list, `𝝭`, ... (indexed by ξ.𝑠)
#
# Note that :x, :y, :z are stored twice as far as `set𝝭!` is concerned: the
# *generator* writes them with `slot = 2` (one value per quadrature point in a
# dedicated array), while `get𝑿ᵢ` writes them with `slot = 1` (indexed by node
# id). The quadrature-point copy is what `set𝝭!` reads through `ξ.x`, `ξ.y`.
#
# Conventions that the approximation code relies on:
#
#   :𝐽    determinant of the mapping 𝒙(ξ,η,γ), i.e. det(∂𝒙/∂𝝃).
#          NOT the inverse, NOT a length. 1-D elements therefore store the
#          arc-length |∂𝒙/∂ξ|, 2-D elements twice the area and 3-D elements
#          six times the volume of the reference simplex.
#   :𝑤    integration weight already multiplied by :𝐽.
#   :w    the raw quadrature weight, without :𝐽.
#
# Consequences for shape-function gradients:
#
# * 1-D elements (`Seg2`, `Seg3`) have an affine map, so
#
#           dNᵢ/dx = (dNᵢ/dξ) / :𝐽
#
#   with the ξ-derivative *divided* by the Jacobian (a length scale). This is
#   the `(y₂-y₃)/𝐽` pattern of `Tri3` as well; multiplying by it is always an
#   error.
# * 2-D/3-D simplex code (`Tri3`, `Tri6`, `Tet4`, `Tet10`, `TriHermite`)
#   divides the difference product by :𝐽. Since the difference product already
#   carries the reference-simplex determinant (1/2 for a triangle, 1/6 for a
#   tetrahedron), that division correctly returns the physical gradient.
# * isoparametric elements (`Quad4`, `Quad8`, `Hex8`) compute ∂ξ/∂x from the
#   full Jacobian matrix and never touch :𝐽 directly.
#

# ---------------------------------------------------------------------------
# Legacy entry points
# ---------------------------------------------------------------------------
#
# The importer used to be 17 near-identical functions assembled from `quote`
# fragments spliced with `@eval`. The four axes along which they varied (element
# type, quadrature, geometry, generation) are independent, so the implementation
# now lives in `elements.jl` behind a single `ElementSpec` and one
# `getElements(nodes, spec)`. The signatures below are kept unchanged so that
# every existing script keeps working; each one forwards to that implementation.
#
# `getPiecewiseElements`, `getPiecewiseBoundaryElements` and
# `getCurvedPiecewiseElements` were removed: they called `get𝑛𝑝(::Element)`,
# a method that has never existed in this package, so they could not run, and
# they had no callers. Piecewise approximations are constructed by building the
# elements directly in Julia instead.

"""
    getElements(nodes, dimTag, integrationOrder = -1; normal = false)

Gmsh decides the element type; the quadrature rule is Gmsh's rule of order
`integrationOrder` (`-1` means the element's own order).

    elements = getElements(nodes, entities["Ω"], 2)
"""
getElements(nodes::Vector{N}, dimTag::Pair{Int,Vector{Int}},
            integrationOrder::Int = -1; normal::Bool = false) where N<:Node =
    getElements(nodes, ElementSpec(dimTag; quadrature = integrationOrder, normal = normal))

"""
    getElements(nodes, dimTag, integration::NTuple{2,Vector{Float64}}; normal = false)

Gmsh decides the element type; the quadrature rule is supplied explicitly.

    elements = getElements(nodes, entities["Ω"], (localCoord, weights))
"""
getElements(nodes::Vector{N}, dimTag::Pair{Int,Vector{Int}},
            integration::NTuple{2,Vector{Float64}}; normal::Bool = false) where N<:Node =
    getElements(nodes, ElementSpec(dimTag; quadrature = integration, normal = normal))

"""
    getElements(nodes, dimTag, type::DataType, integrationOrder = -1; normal = false)

Caller-chosen element type, Gmsh quadrature of order `integrationOrder`.

    elements = getElements(nodes, entities["Ω"], Element{:Tri6}, 3)
"""
getElements(nodes::Vector{N}, dimTag::Pair{Int,Vector{Int}}, type::DataType,
            integrationOrder::Int = -1; normal::Bool = false) where N<:Node =
    getElements(nodes, ElementSpec(dimTag; element = type,
                                   quadrature = integrationOrder, normal = normal))

"""
    getElements(nodes, dimTag, type::DataType, integration; normal = false)

Caller-chosen element type and quadrature rule.
"""
getElements(nodes::Vector{N}, dimTag::Pair{Int,Vector{Int}}, type::DataType,
            integration::NTuple{2,Vector{Float64}}; normal::Bool = false) where N<:Node =
    getElements(nodes, ElementSpec(dimTag; element = type,
                                   quadrature = integration, normal = normal))

"""
    getElements(nodes, dimTag, type, integrationOrder, sp::SpatialPartition; normal = false)
    getElements(nodes, dimTag, type, integration,       sp::SpatialPartition; normal = false)

Support points are **all nodes within the search radius** rather than each
element's own nodes — the layout a meshfree operator needs. `sp` is a
`SpatialPartition`, typically a `RegularGrid`.
"""
getElements(nodes::Vector{N}, dimTag::Pair{Int,Vector{Int}}, type::DataType,
            integrationOrder::Int, sp::SpatialPartition;
            normal::Bool = false) where N<:Node =
    getElements(nodes, ElementSpec(dimTag; element = type, quadrature = integrationOrder,
                                   assembly = :neighbour, search = sp, normal = normal))

getElements(nodes::Vector{N}, dimTag::Pair{Int,Vector{Int}}, type::DataType,
            integration::NTuple{2,Vector{Float64}}, sp::SpatialPartition;
            normal::Bool = false) where N<:Node =
    getElements(nodes, ElementSpec(dimTag; element = type, quadrature = integration,
                                   assembly = :neighbour, search = sp, normal = normal))

"""
    getCurvedElements(nodes, dimTag, cs, integrationOrder = -1)
    getCurvedElements(nodes, dimTag, cs, integration)
    getCurvedElements(nodes, dimTag, type, cs, integrationOrder, sp)
    getCurvedElements(nodes, dimTag, type, cs, integration, sp)

Elements on a surface described by the map `cs` instead of the flat mesh. The
integration weight is the curvature-corrected Jacobian `cs.𝐽(𝒙)` rather than
Gmsh's flat determinant, so an area computed on the mapped surface is a true
surface area.

`cs` carries the surface description as **fields** (it is read as `cs.𝐽`,
`cs.a₁₁`, …), so it may be any type, callable or not:

| field | meaning |
|---|---|
| `cs.𝐽(𝒙)` | area scale of the surface map at `𝒙` |
| `cs.a₁₁`, `cs.a₁₂`, `cs.a₂₂` | covariant surface metric |
| `cs.a¹¹`, `cs.a¹²`, `cs.a²²` | contravariant surface metric |
| `cs.𝒂₁`, `cs.𝒂₂`, `cs.𝒂₃` | surface tangent vectors |

# Example
```julia
struct FlatSurface end
𝐽(::FlatSurface, x)   = 1.0
a₁₁(::FlatSurface, x) = 1.0
# ...
cs = (𝐽 = x -> 1.0, a₁₁ = x -> 1.0, a₁₂ = x -> 0.0, a₂₂ = x -> 1.0,
      a¹¹ = x -> 1.0, a¹² = x -> 0.0, a²² = x -> 1.0,
      𝒂₁ = x -> (1.0,0.0,0.0), 𝒂₂ = x -> (0.0,1.0,0.0), 𝒂₃ = x -> (0.0,0.0,1.0))
elements = getCurvedElements(nodes, entities["Ω"], cs, 2)
```

!!! note "Signature change in S3"
    The historical signature declared `cs::Function`, but the implementation
    reads `cs.𝐽`, `cs.a₁₁`, … as *fields*. No `Function` carries those fields, so
    the old method could never be called successfully — it raised
    `UndefVariableError`/`FieldError` for every argument. The `::Function`
    restriction has been removed so that the documented usage works.
"""
getCurvedElements(nodes::Vector{N}, dimTag::Pair{Int,Vector{Int}}, cs,
                  integrationOrder::Int = -1) where N<:Node =
    getElements(nodes, ElementSpec(dimTag; geometry = :curved, map = cs,
                                   quadrature = integrationOrder))

getCurvedElements(nodes::Vector{N}, dimTag::Pair{Int,Vector{Int}}, cs,
                  integration::NTuple{2,Vector{Float64}}) where N<:Node =
    getElements(nodes, ElementSpec(dimTag; geometry = :curved, map = cs,
                                   quadrature = integration))

getCurvedElements(nodes::Vector{N}, dimTag::Pair{Int,Vector{Int}}, type::DataType,
                  cs, integrationOrder::Int,
                  sp::SpatialPartition) where N<:Node =
    getElements(nodes, ElementSpec(dimTag; element = type, geometry = :curved, map = cs,
                                   quadrature = integrationOrder,
                                   assembly = :neighbour, search = sp))

getCurvedElements(nodes::Vector{N}, dimTag::Pair{Int,Vector{Int}}, type::DataType,
                  cs, integration::NTuple{2,Vector{Float64}},
                  sp::SpatialPartition) where N<:Node =
    getElements(nodes, ElementSpec(dimTag; element = type, geometry = :curved, map = cs,
                                   quadrature = integration,
                                   assembly = :neighbour, search = sp))


function getElements(dimTag1::Pair{Int,Vector{Int}},dimTag2::Pair{Int,Vector{Int}},elms::Vector{T}) where T<:AbstractElement
    elements = AbstractElement[]
    dim1, tag1 = dimTag1
    dim2, tag2 = dimTag2
    elementTypes1 = Int32[]
    elementTypes2 = Int32[]
    nodeTags1 = Vector{UInt64}[]
    nodeTags2 = Vector{UInt64}[]
    for tag in tag1
        elementTypes_, ~, nodeTags_ = gmsh.model.mesh.getElements(dim1,tag)
        push!(elementTypes1,elementTypes_[1])
        push!(nodeTags1,nodeTags_[1])
    end
    for tag in tag2
        elementTypes_, ~, nodeTags_ = gmsh.model.mesh.getElements(dim2,tag)
        push!(elementTypes2,elementTypes_[1])
        push!(nodeTags2,nodeTags_[1])
    end
    for (elementType1,nodeTag1) in zip(elementTypes1,nodeTags1)
        j₀ = 0
        for (elementType2,nodeTag2) in zip(elementTypes2,nodeTags2)
            if elementType1 == elementType2
                ~, ~, ~, ni = gmsh.model.mesh.getElementProperties(elementType1)
                ne1 = Int(length(nodeTag1)/ni)
                ne2 = Int(length(nodeTag2)/ni)
                for i in 1:ne1
                    for j in 1:ne2
                        if nodeTag1[ni*(i-1)+1:ni*i] == nodeTag2[ni*(j-1)+1:ni*j]
                            push!(elements,elms[j₀+j])
                            continue
                        end
                    end
                end
                j₀ += ne2
            end
        end
    end
    return elements
end

end