"""
    ElementSpec

Everything that distinguishes one `getElements` call from another.

Before this type existed the importer was 17 near-identical functions
assembled from `quote` fragments spliced with `@eval`; the four axes along
which they varied were element type, quadrature rule, geometry and element
generation. Those axes are independent, so one function with four dispatch
points replaces all 17.

# Fields
- `dimTag`     : Gmsh `(dim, tags)` pair identifying the entities to read
- `element`    : `:auto` to use the element type Gmsh reports, or a concrete
                 `DataType` such as `Element{:Tri6}` or
                 `PiecewisePolynomial{:Linear2D}`
- `quadrature` : `:auto` for the Gmsh rule of the element's own order, an `Int`
                 for a Gmsh rule of that order, or an explicit
                 `(localCoord, weights)` tuple
- `geometry`   : `:straight`, `:curved` (requires `map`) or `:boundary`
                 (requires `dimTagΩ`)
- `assembly`   : `:fem` or `:neighbour` (requires `search`)
- `normal`     : also fill the outward normal `ξ.nᵢ` and tangent `ξ.sᵢ`
- `map`        : surface map for `geometry = :curved`
- `dimTagΩ`    : ambient entity pairs for `geometry = :boundary`
- `search`     : spatial partition for `assembly = :neighbour`

# Example
```julia
# The common case: Gmsh decides everything, order-2 rule.
getElements(nodes, entities["Ω"], 2)

# Explicit element type.
getElements(nodes, entities["Ω"], Element{:Tri6}, 3)

# Explicit quadrature rule.
getElements(nodes, entities["Ω"], (localCoord, weights))

# Boundary normals.
getElements(nodes, entities["Γᵍ"], 2; normal = true)
```
"""
struct ElementSpec
    dimTag    :: Pair{Int,Vector{Int}}
    element   :: Any
    quadrature:: Any
    geometry  :: Symbol
    assembly  :: Symbol
    normal    :: Bool
    map       :: Any
    dimTagΩ   :: Any
    search    :: Any
end

# Convenience constructor: everything optional except the entities.
ElementSpec(dimTag::Pair{Int,Vector{Int}};
            element = :auto, quadrature = :auto,
            geometry::Symbol = :straight, assembly::Symbol = :fem,
            normal::Bool = false, map = nothing,
            dimTagΩ = nothing, search = nothing) =
    ElementSpec(dimTag, element, quadrature, geometry, assembly, normal,
                map, dimTagΩ, search)

# Gmsh element-type tag -> element family name. The tags are the Gmsh/VTK
# identifiers, so this table is fixed by the file format.
const GMSH_ELEMENT_TYPES = Dict(
    1  => :Seg2,  2  => :Tri3,  3  => :Quad4,  4  => :Tet4,
    5  => :Hex8,  8  => :Seg3,  9  => :Tri6,  10  => :Quad9,
    11 => :Tet10, 12 => :Hex27, 15 => :Poi1,  16 => :Quad8,
)

"""
    resolve_element(spec, elementType) -> DataType

Return the element type to construct. `spec.element === :auto` maps the Gmsh
element tag through [`GMSH_ELEMENT_TYPES`](@ref); otherwise the caller's type is
used unchanged.
"""
resolve_element(spec::ElementSpec, elementType) =
    spec.element === :auto ? Element{GMSH_ELEMENT_TYPES[elementType]} : spec.element

"""
    resolve_quadrature(spec, elementType) -> (localCoord, weights, ni)

Return the quadrature rule to use and `ni`, the number of nodes per Gmsh
element.

  * `spec.quadrature === :auto` — Gmsh's rule of the element's own order;
  * `spec.quadrature::Int`      — Gmsh's rule of that order;
  * otherwise                   — the supplied `(localCoord, weights)` tuple.

Quadrilaterals and hexahedra (Gmsh tags 3, 5, 10, 12, 16) use `CompositeGauss`,
which is what makes a non-tensor-product rule usable on them.
"""
function resolve_quadrature(spec::ElementSpec, elementType)
    _, _, order, ni = gmsh.model.mesh.getElementProperties(elementType)
    if spec.quadrature isa Int
        integrationOrder = spec.quadrature < 0 ? order : spec.quadrature
        kind = elementType ∈ (3, 5, 10, 12, 16) ? "CompositeGauss" : "Gauss"
        localCoord, weights = gmsh.model.mesh.getIntegrationPoints(
            elementType, kind * string(integrationOrder))
    elseif spec.quadrature === :auto
        kind = elementType ∈ (3, 5, 10, 12, 16) ? "CompositeGauss" : "Gauss"
        localCoord, weights = gmsh.model.mesh.getIntegrationPoints(
            elementType, kind * string(order))
    else
        localCoord, weights = spec.quadrature
    end
    return localCoord, weights, ni
end

"""
    entity_list(dimTag::Pair{Int,Vector{Int}}) -> Vector{(elementType, nodeTags, tag)}

Flatten a Gmsh physical group into the `(element type, node tags, entity tag)`
triples that [`getElements`](@ref) iterates over.

A physical group may span several entities and, in principle, several element
types. Only the elements of each entity are returned together with that
entity's tag, because the quadrature rule and the Jacobian are queried per
entity. A physical group containing mixed element types is therefore read one
type at a time, which is the historical behaviour.
"""
function entity_list(dimTag::Pair{Int,Vector{Int}})
    dim, tags = dimTag
    out = Tuple{Int32,Vector{UInt64},Int}[]
    for tag in tags
        elementTypes, _, nodeTags = gmsh.model.mesh.getElements(dim, tag)
        for (elementType, nodeTags_) in zip(elementTypes, nodeTags)
            push!(out, (elementType, nodeTags_, tag))
        end
    end
    return out
end

"""
    generate_fem!(elements, nodes, type, nodeTag, ni, ng, data, 𝐺, 𝐶, 𝑠)
        -> (𝐺, 𝐶, 𝑠)

One element per Gmsh element, with the support points taken from `nodes`.

The four counters index the shared `data` dictionary and are threaded through
every element of a call, which is why they are passed in and returned rather
than reset per element:

| counter | index component | meaning |
|---|---|---|
| `𝐺` | `:𝐺` | quadrature point, global |
| `𝐶` | `:𝐶` | element |
| `𝑠` | `:𝑠` | offset into the per-element DOF list |
"""
function generate_fem!(elements, nodes, type, nodeTag, ni, ng, data, 𝐺, 𝐶, 𝑠)
    ne = div(length(nodeTag), ni)
    for C in 1:ne
        𝐶 += 1
        𝓒 = nodes[nodeTag[ni*(C-1)+1 : ni*C]]
        # `:𝑔` is the quadrature-point index *within* one element, so it always
        # runs 1:ng; `:𝐺` is the global counter.
        𝓖 = [𝑿ₛ((𝑔 = g, 𝐺 = 𝐺 + g, 𝐶 = 𝐶, 𝑠 = 𝑠 + (g-1)*ni), data) for g in 1:ng]
        𝐺 += ng
        𝑠 += ng*ni
        push!(elements, type(𝓒, 𝓖))
    end
    return 𝐺, 𝐶, 𝑠
end

"""
    generate_neighbour!(elements, nodes, type, ni, ng, data, sp, x, y, z, 𝐺, 𝐶, 𝑠)
        -> (𝐺, 𝐶, 𝑠)

One element per Gmsh element, whose support points are **all nodes within the
search radius** rather than the element's own nodes.

This is what a meshfree / reproducing-kernel operator needs: each quadrature
point gets the cloud of nodes inside its support, taken from the spatial
partition `sp`. The same counters as in [`generate_fem!`](@ref) are threaded
through.
"""
function generate_neighbour!(elements, nodes, type, ni, ng, data, sp, x, y, z, 𝐺, 𝐶, 𝑠)
    ne = div(length(x), ng)
    for C in 1:ne
        𝐶 += 1
        indices = Set{Int}()
        for g in 1:ng
            union!(indices, sp(x[ng*(C-1)+g], y[ng*(C-1)+g], z[ng*(C-1)+g]))
        end
        ni_ = length(indices)
        𝓒 = [nodes[i] for i in indices]
        𝓖 = [𝑿ₛ((𝑔 = g, 𝐺 = 𝐺 + g, 𝐶 = 𝐶, 𝑠 = 𝑠 + (g-1)*ni_), data) for g in 1:ng]
        𝐺 += ng
        𝑠 += ng*ni_
        push!(elements, type(𝓒, 𝓖))
    end
    return 𝐺, 𝐶, 𝑠
end

"""
    log_summary(elements_so_far, type, ng)

Print the one-line summary the importer has always emitted. Kept because
several scripts grep for it.
"""
function log_summary(elements, type, ng)
    println("Info    : Generate $(length(elements)) elements of $type with $ng integration points.")
end

"""
    getElements(nodes, spec::ElementSpec) -> Vector{<:AbstractElement}

Build integration elements from a Gmsh physical group.

This is the single implementation behind every `getElements` /
`getCurvedElements` overload. It reads the entities of `spec.dimTag`, and for
each one fills the shared `data` dictionary and generates the elements, in four
steps whose behaviour is selected by `spec`:

1. **element type** — [`resolve_element`](@ref)
2. **quadrature** — [`resolve_quadrature`](@ref)
3. **geometry** — straight (`fill_coordinates!`), boundary
   (`fill_edge_coordinates!`) or curved (`fill_curved_coordinates!`), followed
   by the Jacobian and, if requested, the normal
4. **generation** — [`generate_fem!`](@ref) or [`generate_neighbour!`](@ref)

# Storage
Every element of the result points at one `Dict{Symbol,Tuple{Int,Vector{Float64}}}`.
Each entry is `name => (slot, values)` and a node reads

    node.name == data[name][2][ node.index[ data[name][1] ] ]

so the slot number decides which component of `(𝑔, 𝐺, 𝐶, 𝑠)` indexes the value.
Support nodes `𝑿ᵢ` use slot 1 (indexed by `𝐼`); quadrature points `𝑿ₛ` use
slot 2 for physical/derived data, slot 4 for the per-element DOF offset.

# Example
```julia
elements = getElements(nodes, entities["Ω"], 2)
prescribe!(elements, :E => 210e9, :ν => 0.3)
set∇𝝭!(elements)
```
"""
function getElements(nodes::Vector{N}, spec::ElementSpec) where N<:Node
    dim = spec.dimTag.first
    elements = AbstractElement[]
    𝐺 = 0; 𝐶 = 0; 𝑠 = 0
    # One `data` dictionary serves every entity of the call. The original code
    # created it once, before the entity loop, so all quadrature points share
    # the same vectors; recreating it per entity would give each entity its own
    # dictionary and break `prescribe!`, which writes into `elements[1]`'s dict.
    data = initial_data(dim; normal = spec.normal)
    for (elementType, nodeTag, tag) in entity_list(spec.dimTag)
        type = resolve_element(spec, elementType)
        localCoord, weights, ni = resolve_quadrature(spec, elementType)
        ng = length(weights)

        if spec.geometry === :curved
            fill_curved_coordinates!(data, dim, elementType, tag, localCoord, weights, spec.map)
            log_summary(elements, type, ng)   # the curved path has always summarised early
        else
            if spec.geometry === :boundary
                fill_edge_coordinates!(data, dim, elementType, spec.dimTagΩ.first,
                                       tag, spec.dimTagΩ.second, 𝐶, nodeTag,
                                       localCoord, weights)
            else
                fill_coordinates!(data, elementType, tag, localCoord, weights)
            end
            jacobians = gmsh.model.mesh.getJacobians(elementType, localCoord, tag)[1]
            determinants = gmsh.model.mesh.getJacobians(elementType, localCoord, tag)[2]
            fill_jacobian!(data, dim, elementType, tag, localCoord, jacobians)
            spec.normal && fill_normal!(data, dim, elementType, tag, determinants, ng,
                                        div(length(nodeTag), ni))
        end

        if spec.assembly === :neighbour
            x, y, z = data[:x][2], data[:y][2], data[:z][2]
            𝐺, 𝐶, 𝑠 = generate_neighbour!(elements, nodes, type, ni, ng, data,
                                          spec.search, x, y, z, 𝐺, 𝐶, 𝑠)
        else
            𝐺, 𝐶, 𝑠 = generate_fem!(elements, nodes, type, nodeTag, ni, ng, data,
                                    𝐺, 𝐶, 𝑠)
        end
        spec.geometry === :curved || log_summary(elements, type, ng)
    end
    return elements
end
