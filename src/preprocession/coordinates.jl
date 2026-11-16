"""
Helpers that fill the shared quadrature-point dictionary from a Gmsh entity.

Every element produced by `GmshImport.getElements` points at one
`Dict{Symbol,Tuple{Int,Vector{Float64}}}`; the storage layout is described in
`elements.jl`. This module provides the five filling steps:

  * [`initial_data`](@ref)             — create the dictionary with the keys a dimension needs
  * [`fill_coordinates!`](@ref)        — physical coordinates and Jacobian-weighted weights
  * [`fill_edge_coordinates!`](@ref)   — the same plus `:Δ` and ambient local coordinates
  * [`fill_curved_coordinates!`](@ref) — the same on a mapped (curvilinear) surface
  * [`fill_jacobian!`](@ref)           — `:𝐽` and the inverse-Jacobian entries
  * [`fill_normal!`](@ref)             — outward normal `:nᵢ` and tangent `:sᵢ`

`Gmsh` returns three local coordinates per quadrature point regardless of the
element dimension, so `localCoord` is always ordered `[ξ₁, η₁, γ₁, ξ₂, …]`.
"""
module MeshCoordinates

import Gmsh: gmsh

export initial_data, fill_coordinates!, fill_edge_coordinates!,
       fill_curved_coordinates!, fill_jacobian!, fill_normal!

"""
    initial_data(dim; normal=false)

Create the shared `data` dictionary for elements of dimension `dim`.

Each entry is `name => (slot, values)`, where `slot` selects a component of the
node's index tuple and its value is the offset into `values`; see the layout
description in `elements.jl`. All vectors start empty.

With `normal = true` the outward-normal fields `:n₁ … :s₃` are created as well.
"""
function initial_data(dim::Int; normal::Bool=false)
    data = Dict{Symbol,Tuple{Int,Vector{Float64}}}()
    data[:w]    = (1, Float64[])
    data[:ξ]    = (1, Float64[])
    data[:x]    = (2, Float64[])
    data[:y]    = (2, Float64[])
    data[:z]    = (2, Float64[])
    data[:𝑤]    = (2, Float64[])
    data[:𝐽]    = (2, Float64[])
    data[:∂ξ∂x] = (2, Float64[])
    if normal
        data[:n₁] = (3, Float64[]); data[:n₂] = (3, Float64[])
        data[:s₁] = (3, Float64[]); data[:s₂] = (3, Float64[])
    end
    if dim >= 2
        data[:η]    = (1, Float64[])
        data[:∂ξ∂y] = (2, Float64[])
        data[:∂η∂x] = (2, Float64[])
        data[:∂η∂y] = (2, Float64[])
        if normal
            data[:n₃] = (3, Float64[]); data[:s₃] = (3, Float64[])
        end
    end
    if dim >= 3
        data[:γ]    = (1, Float64[])
        data[:∂ξ∂z] = (2, Float64[])
        data[:∂η∂z] = (2, Float64[])
        data[:∂γ∂x] = (2, Float64[])
        data[:∂γ∂y] = (2, Float64[])
        data[:∂γ∂z] = (2, Float64[])
    end
    return data
end

"""
    fill_coordinates!(data, elementType, tag, localCoord, weights)

Append the physical coordinates of every quadrature point of the Gmsh entity
`tag`, together with the raw weight `:w` and the Jacobian-weighted weight `:𝑤`.

`gmsh.model.mesh.getJacobians` returns, for every quadrature point, the image of
`localCoord` (that is, the physical `x`, `y`, `z`) and the determinant of the
mapping. `:𝑤 = det · w` is what an integral is actually summed with.
"""
function fill_coordinates!(data, elementType, tag, localCoord, weights)
    ng = length(weights)
    _, determinants, coord = gmsh.model.mesh.getJacobians(elementType, localCoord, tag)
    append!(data[:x][2], coord[1:3:end])
    append!(data[:y][2], coord[2:3:end])
    append!(data[:z][2], coord[3:3:end])
    append!(data[:w][2], weights)
    # Local coordinates of the quadrature points. Only the components that
    # `initial_data` created for this dimension are written.
    append!(data[:ξ][2], localCoord[1:3:end])
    haskey(data, :η) && append!(data[:η][2], localCoord[2:3:end])
    haskey(data, :γ) && append!(data[:γ][2], localCoord[3:3:end])
    for i in 1:div(length(determinants), ng), (j, w) in enumerate(weights)
        push!(data[:𝑤][2], determinants[ng*(i-1)+j]*w)
    end
    return data
end

"""
    fill_edge_coordinates!(data, dim, elementType, elementTypeΩ, tag, tagsΩ, CΩ,
                           nodeTag, localCoord, weights)

Coordinates and quadrature points for a boundary entity of a 2-D or 3-D domain.

On top of what [`fill_coordinates!`](@ref) does, two boundary-specific pieces of
data are recorded:

  * `:Δ` — the sign of the end point hit by each quadrature point (+1, −1, or 0
    in the interior). The corner-force terms of the Kirchhoff plate read it.
  * `:ξ`, `:η`, `:γ` — the local coordinates of every boundary quadrature point
    **inside the ambient element**, which the piecewise-boundary operators need.

`ne` (the number of boundary elements of this entity) is derived from the node
tags, so the caller does not have to supply it.
"""
function fill_edge_coordinates!(data, dim, elementType, elementTypeΩ, tag, tagsΩ, CΩ,
                                nodeTag, localCoord, weights)
    ng  = length(weights)
    nb  = elementTypeΩ ∈ (2, 9) ? 3 : 4          # edges per ambient element
    _, _, _, ni = gmsh.model.mesh.getElementProperties(elementType)
    ne  = Int(length(nodeTag)/ni)                # number of boundary elements
    _, determinants, coord = gmsh.model.mesh.getJacobians(elementType, localCoord, tag)
    x = coord[1:3:end]; y = coord[2:3:end]; z = coord[3:3:end]
    append!(data[:x][2], x)
    append!(data[:y][2], y)
    append!(data[:z][2], z)
    append!(data[:w][2], weights)
    for i in 1:div(length(determinants), ng), (j, w) in enumerate(weights)
        push!(data[:𝑤][2], determinants[ng*(i-1)+j]*w)
    end
    for g in 1:ng
        ξg = localCoord[3*g-2]
        push!(data[:Δ][2], ξg ≈ 1.0 ? 1.0 : (ξg ≈ -1.0 ? -1.0 : 0.0))
    end
    # Map each boundary quadrature point back into its ambient element.
    for C in 1:ne
        tagΩ = tagsΩ[CΩ + div(ne, nb)]
        for g in 1:ng
            G = ng*(C-1) + g
            ξ, η, γ = gmsh.model.mesh.getLocalCoordinatesInElement(tagΩ, x[G], y[G], z[G])
            push!(data[:ξ][2], ξ)
            push!(data[:η][2], η)
            haskey(data, :γ) && push!(data[:γ][2], γ)
        end
    end
    return data
end

"""
    fill_curved_coordinates!(data, dim, elementType, tag, localCoord, weights, cs)

Coordinates and weights for elements on a surface described by the map `cs`
instead of by the flat mesh.

Instead of Gmsh's determinant the weight uses the curvature-corrected Jacobian
`cs.𝐽(x)`, so areas computed on the mapped surface are true surface areas.

For `dim == 1` the surface metric at every quadrature point is recorded too:
`:𝒂₁ :𝒂₂ :𝒂₃` (tangents), `:a¹¹ :a¹² :a²²` (contravariant metric),
`:nᵢ :sᵢ :nⁱ :sⁱ` (normal, tangent, and their contravariant forms) and the
finite-difference derivatives `:∂₁nᵢ :∂₂nᵢ :∂₁sᵢ :∂₂sᵢ`.
"""
function fill_curved_coordinates!(data, dim, elementType, tag, localCoord, weights, cs)
    ng = length(weights)
    _, determinants, coord = gmsh.model.mesh.getJacobians(elementType, localCoord, tag)
    jacobians = gmsh.model.mesh.getJacobians(elementType, localCoord, tag)[1]
    x = coord[1:3:end]; y = coord[2:3:end]; z = coord[3:3:end]
    ne = div(length(determinants), ng)
    𝑤 = zeros(length(determinants))

    if dim == 2
        for i in 1:ne, (j, w) in enumerate(weights)
            G  = ng*(i-1) + j
            x_ = (x[G], y[G], z[G])
            𝑤[G] = determinants[G]*cs.𝐽(x_)*w
        end
        data[:w] = (1, weights)
        data[:x] = (2, x); data[:y] = (2, y); data[:z] = (2, z)
        data[:𝑤] = (2, 𝑤)
    elseif dim == 1
        ∂x∂ξ = jacobians[1:9:end]
        ∂y∂ξ = jacobians[2:9:end]
        ∂z∂ξ = jacobians[3:9:end]
        nodeTags = gmsh.model.mesh.getElementEdgeNodes(elementType, tag, true)
        n₁ = zeros(ne*ng); n₂ = zeros(ne*ng)
        n¹ = zeros(ne*ng); n² = zeros(ne*ng)
        s₁ = zeros(ne*ng); s₂ = zeros(ne*ng)
        s¹ = zeros(ne*ng); s² = zeros(ne*ng)
        ∂₁n₁ = zeros(ne*ng); ∂₁n₂ = zeros(ne*ng)
        ∂₂n₁ = zeros(ne*ng); ∂₂n₂ = zeros(ne*ng)
        ∂₁s₁ = zeros(ne*ng); ∂₁s₂ = zeros(ne*ng)
        ∂₂s₁ = zeros(ne*ng); ∂₂s₂ = zeros(ne*ng)
        for C in 1:ne
            𝐿 = 2*determinants[C*ng]
            c₁, = gmsh.model.mesh.getNode(nodeTags[2*C-1])
            c₂, = gmsh.model.mesh.getNode(nodeTags[2*C])
            x₁, y₁ = c₁[1], c₁[2]
            x₂, y₂ = c₂[1], c₂[2]
            t¹ = (x₂-x₁)/𝐿
            t² = (y₂-y₁)/𝐿
            t₁(x_) = cs.a₁₁(x_)*t¹ + cs.a₁₂(x_)*t²
            t₂(x_) = cs.a₁₂(x_)*t¹ + cs.a₂₂(x_)*t²
            t(x_)  = (t₁(x_)*t¹ + t₂(x_)*t²)^0.5
            s¹_(x_) = t¹/t(x_)
            s²_(x_) = t²/t(x_)
            s₁_(x_) = t₁(x_)/t(x_)
            s₂_(x_) = t₂(x_)/t(x_)
            deta(x_) = (cs.a₁₁(x_)*cs.a₂₂(x_) - cs.a₁₂(x_)^2)^0.5
            n₁_(x_) =  s²_(x_)*deta(x_)
            n₂_(x_) = -s¹_(x_)*deta(x_)
            n¹_(x_) = cs.a¹¹(x_)*n₁_(x_) + cs.a¹²(x_)*n₂_(x_)
            n²_(x_) = cs.a¹²(x_)*n₁_(x_) + cs.a²²(x_)*n₂_(x_)
            for (j, w) in enumerate(weights)
                G  = ng*(C-1) + j
                x_ = (x[G], y[G], z[G])
                𝒂₁_, 𝒂₂_, 𝒂₃_ = cs.𝒂₁(x_), cs.𝒂₂(x_), cs.𝒂₃(x_)
                J = ((𝒂₁_[1]*∂x∂ξ[G] + 𝒂₂_[1]*∂y∂ξ[G] + 𝒂₃_[1]*∂z∂ξ[G])^2
                   + (𝒂₁_[2]*∂x∂ξ[G] + 𝒂₂_[2]*∂y∂ξ[G] + 𝒂₃_[2]*∂z∂ξ[G])^2
                   + (𝒂₁_[3]*∂x∂ξ[G] + 𝒂₂_[3]*∂y∂ξ[G] + 𝒂₃_[3]*∂z∂ξ[G])^2)^0.5
                s₁[G] = s₁_(x_); s₂[G] = s₂_(x_)
                s¹[G] = s¹_(x_); s²[G] = s²_(x_)
                n₁[G] = n₁_(x_); n₂[G] = n₂_(x_)
                n¹[G] = n¹_(x_); n²[G] = n²_(x_)
                ∂₁n₁[G] = gradient(n₁_, x_)[1]; ∂₂n₁[G] = gradient(n₁_, x_)[2]
                ∂₁n₂[G] = gradient(n₂_, x_)[1]; ∂₂n₂[G] = gradient(n₂_, x_)[2]
                ∂₁s₁[G] = gradient(s₁_, x_)[1]; ∂₂s₁[G] = gradient(s₁_, x_)[2]
                ∂₁s₂[G] = gradient(s₂_, x_)[1]; ∂₂s₂[G] = gradient(s₂_, x_)[2]
                𝑤[G] = J*w
            end
        end
        Δ = [localCoord[3g-2] ≈ 1.0 ? 1.0 : (localCoord[3g-2] ≈ -1.0 ? -1.0 : 0.0)
             for g in 1:ng]
        data[:w] = (1, weights)
        data[:x] = (2, x); data[:y] = (2, y); data[:z] = (2, z)
        data[:𝑤] = (2, 𝑤)
        data[:n₁] = (2, n₁); data[:n₂] = (2, n₂)
        data[:n¹] = (2, n¹); data[:n²] = (2, n²)
        data[:s₁] = (2, s₁); data[:s₂] = (2, s₂)
        data[:s¹] = (2, s¹); data[:s²] = (2, s²)
        data[:∂₁n₁] = (2, ∂₁n₁); data[:∂₁n₂] = (2, ∂₁n₂)
        data[:∂₂n₁] = (2, ∂₂n₁); data[:∂₂n₂] = (2, ∂₂n₂)
        data[:∂₁s₁] = (2, ∂₁s₁); data[:∂₁s₂] = (2, ∂₁s₂)
        data[:∂₂s₁] = (2, ∂₂s₁); data[:∂₂s₂] = (2, ∂₂s₂)
        data[:Δ]  = (1, Δ)
    end
    data[:ξ] = (1, localCoord[1:3:end])
    dim >= 2 && (data[:η] = (1, localCoord[2:3:end]))
    return data
end

"""
    fill_jacobian!(data, dim, elementType, tag, localCoord, jacobians)

Record `data[:𝐽]` (the determinant) and the inverse-Jacobian entries
`:∂ξ∂x`, `:∂ξ∂y`, … read by the isoparametric shape functions.

`jacobians` is the flat 3×3-per-point block returned by
`gmsh.model.mesh.getJacobians`; its inverse is exactly the `∂ξ/∂x` block.

!!! note "Preserved quirk"
    The original implementation appended the `z`-direction inverse-Jacobian
    entries in an `elseif dim >= 3` branch that can never be reached, because
    `dim >= 2` is tested first. The arrays `:∂ξ∂z`, `:∂η∂z`, `:∂γ∂x/y/z` are
    therefore **never filled** for 3-D elements. This rewrite keeps that
    behaviour so that results do not change; see S3_DESIGN.md §2.2.
"""
function fill_jacobian!(data, dim, elementType, tag, localCoord, jacobians)
    ng = div(length(localCoord), 3)
    ne = div(length(jacobians), 9*ng)
    determinants = gmsh.model.mesh.getJacobians(elementType, localCoord, tag)[2]
    append!(data[:𝐽][2], determinants)
    ∂ξ∂x = zeros(ne*ng); ∂ξ∂y = zeros(ne*ng); ∂ξ∂z = zeros(ne*ng)
    ∂η∂x = zeros(ne*ng); ∂η∂y = zeros(ne*ng); ∂η∂z = zeros(ne*ng)
    ∂γ∂x = zeros(ne*ng); ∂γ∂y = zeros(ne*ng); ∂γ∂z = zeros(ne*ng)
    J = zeros(3, 3)
    for C in 1:ne, g in 1:ng
        b = 9*(ng*(C-1)+g) - 9
        for a in 1:9
            J[a] = jacobians[b + a]
        end
        J⁻¹ = inv(J)
        k = ng*(C-1) + g
        ∂ξ∂x[k] = J⁻¹[1,1]; ∂ξ∂y[k] = J⁻¹[1,2]; ∂ξ∂z[k] = J⁻¹[1,3]
        ∂η∂x[k] = J⁻¹[2,1]; ∂η∂y[k] = J⁻¹[2,2]; ∂η∂z[k] = J⁻¹[2,3]
        ∂γ∂x[k] = J⁻¹[3,1]; ∂γ∂y[k] = J⁻¹[3,2]; ∂γ∂z[k] = J⁻¹[3,3]
    end
    append!(data[:∂ξ∂x][2], ∂ξ∂x)
    if dim >= 2
        append!(data[:∂ξ∂y][2], ∂ξ∂y)
        append!(data[:∂η∂x][2], ∂η∂x)
        append!(data[:∂η∂y][2], ∂η∂y)
    elseif dim >= 3
        # Unreachable, kept to preserve the historical output exactly.
        append!(data[:∂ξ∂z][2], ∂ξ∂z)
        append!(data[:∂η∂z][2], ∂η∂z)
        append!(data[:∂γ∂x][2], ∂γ∂x)
        append!(data[:∂γ∂y][2], ∂γ∂y)
        append!(data[:∂γ∂z][2], ∂γ∂z)
    end
    return data
end

"""
    fill_normal!(data, dim, elementType, tag, determinants, ng, ne)

Outward unit normal `:nᵢ` and surface tangent `:sᵢ` of a boundary entity.

For a 1-D boundary the normal follows from the difference of the two end points.
For a 2-D boundary it is accumulated from the cross products of successive
edge-point pairs, the standard way to obtain an outward normal of a
quadrilateral or triangle face in 3-D.
"""
function fill_normal!(data, dim, elementType, tag, determinants, ng, ne)
    nodeTags = gmsh.model.mesh.getElementEdgeNodes(elementType, tag, true)
    if dim == 1
        for C in 1:ne
            𝐿 = 2*determinants[C*ng]
            c₁, = gmsh.model.mesh.getNode(nodeTags[2*C-1])
            c₂, = gmsh.model.mesh.getNode(nodeTags[2*C])
            x₁, y₁ = c₁[1], c₁[2]
            x₂, y₂ = c₂[1], c₂[2]
            push!(data[:n₁][2], (y₂-y₁)/𝐿)
            push!(data[:n₂][2], (x₁-x₂)/𝐿)
            push!(data[:s₁][2], (x₂-x₁)/𝐿)
            push!(data[:s₂][2], (y₂-y₁)/𝐿)
        end
    end
    if dim == 2
        nₙ = Int(length(nodeTags)/ne)
        for C in 1:ne
            𝐽 = determinants[C*ng]
            n₁ = 0.0; n₂ = 0.0; n₃ = 0.0
            for i in 1:2:nₙ
                c₁, = gmsh.model.mesh.getNode(nodeTags[nₙ*(C-1)+i])
                c₂, = gmsh.model.mesh.getNode(nodeTags[nₙ*(C-1)+i+1])
                x₁, y₁, z₁ = c₁[1], c₁[2], c₁[3]
                x₂, y₂, z₂ = c₂[1], c₂[2], c₂[3]
                n₁ += y₁*z₂ - y₂*z₁
                n₂ += z₁*x₂ - z₂*x₁
                n₃ += x₁*y₂ - x₂*y₁
            end
            elementType == 3 && (𝐽 *= 8)   # Quad4 repeats each edge node pair
            push!(data[:n₁][2], n₁/𝐽)
            push!(data[:n₂][2], n₂/𝐽)
            push!(data[:n₃][2], n₃/𝐽)
        end
    end
    return data
end

end # module MeshCoordinates
