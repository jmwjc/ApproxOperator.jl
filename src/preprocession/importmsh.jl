module GmshImport
    
using ..ApproxOperator: AbstractElement, Element, Node, 𝑿ᵢ, 𝑿ₛ, SpatialPartition, get𝑛𝑝, NodeData, PerNodeValue
import Gmsh: gmsh

function _vec(dat::NodeData, s::Symbol)
    v = dat.value[s]
    return v.values
end

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

function get𝑿ᵢ()
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
    dat = NodeData(x=(1,PerNodeValue(x)), y=(1,PerNodeValue(y)), z=(1,PerNodeValue(z)))
    return [𝑿ᵢ((𝐼=i,), dat) for i in 1:nₚ ]
end

prequote = quote
    types = Dict([1=>:Seg2, 2=>:Tri3, 3=>:Quad, 4=>:Tet4, 5=>:Hex8, 8=>:Seg3, 9=>:Tri6, 10=>:Quad9, 11=>:Tet10, 12=>:Hex27, 15=>:Poi1, 16=>:Quad8])
    dim, tags = dimTag
    elementTypes = Int32[]
    nodeTags = Vector{UInt64}[]
    for tag in tags
        elementTypes_, ~, nodeTags_ = gmsh.model.mesh.getElements(dim,tag)
        push!(elementTypes,elementTypes_[1])
        push!(nodeTags,nodeTags_[1])
    end
    elements = AbstractElement[]

    𝑔 = 0; 𝐺 = 0; 𝐶 = 0;𝑠 = 0;
    data = NodeData()
    data.value[:w] = PerNodeValue(Float64[]); data.index[:w] = 1
    data.value[:ξ] = PerNodeValue(Float64[]); data.index[:ξ] = 1
    data.value[:x] = PerNodeValue(Float64[]); data.index[:x] = 2
    data.value[:y] = PerNodeValue(Float64[]); data.index[:y] = 2
    data.value[:z] = PerNodeValue(Float64[]); data.index[:z] = 2
    data.value[:𝑤] = PerNodeValue(Float64[]); data.index[:𝑤] = 2
    data.value[:𝐽] = PerNodeValue(Float64[]); data.index[:𝐽] = 2
    data.value[:∂ξ∂x] = PerNodeValue(Float64[]); data.index[:∂ξ∂x] = 2
    if normal
        data.value[:n₁] = PerNodeValue(Float64[]); data.index[:n₁] = 3
        data.value[:n₂] = PerNodeValue(Float64[]); data.index[:n₂] = 3
        data.value[:s₁] = PerNodeValue(Float64[]); data.index[:s₁] = 3
        data.value[:s₂] = PerNodeValue(Float64[]); data.index[:s₂] = 3
    end
    if dim >= 2
        data.value[:η] = PerNodeValue(Float64[]); data.index[:η] = 1

        data.value[:∂ξ∂y] = PerNodeValue(Float64[]); data.index[:∂ξ∂y] = 2
        data.value[:∂η∂x] = PerNodeValue(Float64[]); data.index[:∂η∂x] = 2
        data.value[:∂η∂y] = PerNodeValue(Float64[]); data.index[:∂η∂y] = 2
        if normal
            data.value[:n₃] = PerNodeValue(Float64[]); data.index[:n₃] = 3
            data.value[:s₃] = PerNodeValue(Float64[]); data.index[:s₃] = 3
        end
    end
    if dim >= 3
        data.value[:γ] = PerNodeValue(Float64[]); data.index[:γ] = 1

        data.value[:∂ξ∂z] = PerNodeValue(Float64[]); data.index[:∂ξ∂z] = 2
        data.value[:∂η∂z] = PerNodeValue(Float64[]); data.index[:∂η∂z] = 2
        data.value[:∂γ∂x] = PerNodeValue(Float64[]); data.index[:∂γ∂x] = 2
        data.value[:∂γ∂y] = PerNodeValue(Float64[]); data.index[:∂γ∂y] = 2
        data.value[:∂γ∂z] = PerNodeValue(Float64[]); data.index[:∂γ∂z] = 2
    end
end

preForEdge = quote
    dimΩ,tagΩ = dimTagΩ
    tagsΩ = UInt64[]
    elementTypesΩ = Int32[]
    CΩ = 0
    for tagΩ_ in tagΩ
        elementTypesΩ_, tagsΩ_ = gmsh.model.mesh.getElements(dimΩ,tagΩ_)
        push!(elementTypesΩ,elementTypesΩ_[1])
        push!(tagsΩ,tagsΩ_[1]...)
    end

    data.value[:w] = PerNodeValue(Float64[]); data.index[:w] = 1
    data.value[:Δ] = PerNodeValue(Float64[]); data.index[:Δ] = 1
    data.value[:ξ] = PerNodeValue(Float64[]); data.index[:ξ] = 2
    data.value[:η] = PerNodeValue(Float64[]); data.index[:η] = 2
    data.value[:n₁] = PerNodeValue(Float64[]); data.index[:n₁] = 3
    data.value[:n₂] = PerNodeValue(Float64[]); data.index[:n₂] = 3
    data.value[:s₁] = PerNodeValue(Float64[]); data.index[:s₁] = 3
    data.value[:s₂] = PerNodeValue(Float64[]); data.index[:s₂] = 3
    if dim > 1
        data.value[:γ] = PerNodeValue(Float64[]); data.index[:γ] = 2
        data.value[:n₃] = PerNodeValue(Float64[]); data.index[:n₃] = 3
        data.value[:s₃] = PerNodeValue(Float64[]); data.index[:s₃] = 3
    end   
end

coordinates = quote
    ng = length(weights)
    ne = Int(length(nodeTag)/ni)

    append!(_vec(data,:w),weights)
    haskey(data,:ξ) ? append!(_vec(data,:ξ),localCoord[1:3:end]) : nothing
    haskey(data,:η) ? append!(_vec(data,:η),localCoord[2:3:end]) : nothing
    haskey(data,:γ) ? append!(_vec(data,:γ),localCoord[3:3:end]) : nothing
    jacobians, determinants, coord = gmsh.model.mesh.getJacobians(elementType, localCoord, tag)
    x = coord[1:3:end]
    y = coord[2:3:end]
    z = coord[3:3:end]
    append!(_vec(data,:x),x)
    append!(_vec(data,:y),y)
    append!(_vec(data,:z),z)
    for i in 1:Int(length(determinants)/ng)
        for (j,w) in enumerate(weights)
            G = ng*(i-1)+j
            push!(_vec(data,:𝑤), determinants[G]*w)
        end
    end
end

coordinatesForEdges = quote
    ng = length(weights)
    ne = Int(length(nodeTag)/ni)
    if elementTypeΩ ∈ (2,9)
        nb = 3
    elseif elementTypeΩ ∈ (3,10,16)
        nb = 4
    end

    nodeTag = gmsh.model.mesh.getElementEdgeNodes(elementType,tag,true)

    append!(_vec(data,:w),weights)
    jacobians, determinants, coord = gmsh.model.mesh.getJacobians(elementType, localCoord, tag)
    x = coord[1:3:end]
    y = coord[2:3:end]
    z = coord[3:3:end]
    append!(_vec(data,:x),x)
    append!(_vec(data,:y),y)
    append!(_vec(data,:z),z)
    for i in 1:Int(length(determinants)/ng)
        for (j,w) in enumerate(weights)
            G = ng*(i-1)+j
            push!(_vec(data,:𝑤), determinants[G]*w)
        end
    end

    for g in 1:ng
        ξg = localCoord[3*g-2]
        if ξg ≈ 1.0
            push!(_vec(data,:Δ), 1.0)
        elseif ξg ≈ -1.0
            push!(_vec(data,:Δ), -1.0)
        else
            push!(_vec(data,:Δ), 0.0)
        end
    end
    for CΩ_ in 1:Int(ne/nb)
        tagΩ = tagsΩ[CΩ+CΩ_]
        for C in nb*(CΩ_-1)+1:nb*CΩ_
            𝐿 = 2*determinants[C*ng]
            coord, = gmsh.model.mesh.getNode(nodeTag[2*C-1])
            x₁ = coord[1]
            y₁ = coord[2]
            coord, = gmsh.model.mesh.getNode(nodeTag[2*C])
            x₂ = coord[1]
            y₂ = coord[2]
            push!(_vec(data,:n₁), (y₂-y₁)/𝐿)
            push!(_vec(data,:n₂), (x₁-x₂)/𝐿)
            push!(_vec(data,:s₁), (x₂-x₁)/𝐿)
            push!(_vec(data,:s₂), (y₂-y₁)/𝐿)
            for g in 1:ng
                G = ng*(C-1)+g
                ξ, η, γ = gmsh.model.mesh.getLocalCoordinatesInElement(tagΩ, x[G], y[G], z[G])
                push!(_vec(data,:ξ), ξ)
                push!(_vec(data,:η), η)
                haskey(data,:γ) ? push!(_vec(data,:γ), γ) : nothing
            end
        end
    end
end

curvilinearCoordinates = quote
    ng = length(weights)
    ne = Int(length(nodeTag)/ni)

    ξ = localCoord[1:3:end]
    η = localCoord[2:3:end]
    γ = localCoord[3:3:end]
    jacobians, determinants, coord = gmsh.model.mesh.getJacobians(elementType, localCoord, tag)
    x = coord[1:3:end]
    y = coord[2:3:end]
    z = coord[3:3:end]
    𝑤 = zeros(length(determinants))
    if dim == 2
        for i in 1:Int(length(determinants)/ng)
            for (j,w) in enumerate(weights)
                G = ng*(i-1)+j
                x_ = Vec{3}((x[G],y[G],z[G]))
                𝑤[G] = determinants[G]*cs.𝐽(x_)*w
            end
        end
        data = NodeData(w=(1,PerNodeValue(weights)),x=(2,PerNodeValue(x)),y=(2,PerNodeValue(y)),z=(2,PerNodeValue(z)),𝑤=(2,PerNodeValue(𝑤)))
    elseif dim == 1
        Δ = zeros(ng)
        ∂x∂ξ = jacobians[1:9:end]
        ∂y∂ξ = jacobians[2:9:end]
        ∂z∂ξ = jacobians[3:9:end]
        n₁ = zeros(ne*ng)
        n₂ = zeros(ne*ng)
        n¹ = zeros(ne*ng)
        n² = zeros(ne*ng)
        s₁ = zeros(ne*ng)
        s₂ = zeros(ne*ng)
        s¹ = zeros(ne*ng)
        s² = zeros(ne*ng)
        ∂₁n₁ = zeros(ne*ng)
        ∂₁n₂ = zeros(ne*ng)
        ∂₂n₁ = zeros(ne*ng)
        ∂₂n₂ = zeros(ne*ng)
        ∂₁s₁ = zeros(ne*ng)
        ∂₁s₂ = zeros(ne*ng)
        ∂₂s₁ = zeros(ne*ng)
        ∂₂s₂ = zeros(ne*ng)
        nodeTags = gmsh.model.mesh.getElementEdgeNodes(elementType, tag, true)
        for C in 1:ne
            𝐿 = 2*determinants[C*ng]
            coord, = gmsh.model.mesh.getNode(nodeTags[2*C-1])
            x₁ = coord[1]
            y₁ = coord[2]
            coord, = gmsh.model.mesh.getNode(nodeTags[2*C])
            x₂ = coord[1]
            y₂ = coord[2]
            t¹ = (x₂-x₁)/𝐿
            t² = (y₂-y₁)/𝐿
            t₁(x) = cs.a₁₁(x)*t¹ + cs.a₁₂(x)*t²
            t₂(x) = cs.a₁₂(x)*t¹ + cs.a₂₂(x)*t²
            t(x) = (t₁(x)*t¹ + t₂(x)*t²)^0.5
            s¹_(x) = t¹/t(x)
            s²_(x) = t²/t(x)
            s₁_(x) = t₁(x)/t(x)
            s₂_(x) = t₂(x)/t(x)
            deta(x) = (cs.a₁₁(x)*cs.a₂₂(x) - cs.a₁₂(x)^2)^0.5
            n₁_(x) = s²_(x)*deta(x)
            n₂_(x) =-s¹_(x)*deta(x)
            n¹_(x) = cs.a¹¹(x)*n₁_(x) + cs.a¹²(x)*n₂_(x)
            n²_(x) = cs.a¹²(x)*n₁_(x) + cs.a²²(x)*n₂_(x)
            ∂₁n₁_(x) = gradient(n₁_,x)[1]
            ∂₂n₁_(x) = gradient(n₁_,x)[2]
            ∂₁n₂_(x) = gradient(n₂_,x)[1]
            ∂₂n₂_(x) = gradient(n₂_,x)[2]
            ∂₁s₁_(x) = gradient(s₁_,x)[1]
            ∂₂s₁_(x) = gradient(s₁_,x)[2]
            ∂₁s₂_(x) = gradient(s₂_,x)[1]
            ∂₂s₂_(x) = gradient(s₂_,x)[2]
            for (j,w) in enumerate(weights)
                G = ng*(C-1)+j
                x_ = Vec{3}((x[G],y[G],z[G]))
                𝒂₁_ = cs.𝒂₁(x_)
                𝒂₂_ = cs.𝒂₂(x_)
                𝒂₃_ = cs.𝒂₃(x_)
                J = ((𝒂₁_[1]*∂x∂ξ[G] + 𝒂₂_[1]*∂y∂ξ[G] + 𝒂₃_[1]*∂z∂ξ[G])^2
                  +  (𝒂₁_[2]*∂x∂ξ[G] + 𝒂₂_[2]*∂y∂ξ[G] + 𝒂₃_[2]*∂z∂ξ[G])^2
                  +  (𝒂₁_[3]*∂x∂ξ[G] + 𝒂₂_[3]*∂y∂ξ[G] + 𝒂₃_[3]*∂z∂ξ[G])^2)^0.5
                s₁[G] = s₁_(x_)
                s₂[G] = s₂_(x_)
                s¹[G] = s¹_(x_)
                s²[G] = s²_(x_)
                n₁[G] = n₁_(x_)
                n₂[G] = n₂_(x_)
                n¹[G] = n¹_(x_)
                n²[G] = n²_(x_)
                ∂₁n₁[G] = ∂₁n₁_(x_)
                ∂₁n₂[G] = ∂₁n₂_(x_)
                ∂₂n₁[G] = ∂₂n₁_(x_)
                ∂₂n₂[G] = ∂₂n₂_(x_)
                ∂₁s₁[G] = ∂₁s₁_(x_)
                ∂₁s₂[G] = ∂₁s₂_(x_)
                ∂₂s₁[G] = ∂₂s₁_(x_)
                ∂₂s₂[G] = ∂₂s₂_(x_)
                𝑤[G] = J*w
            end
        end
        for g in 1:ng
            ξg = localCoord[3*g-2]
            if ξg ≈ 1.0
                Δ[g] = 1.0
            elseif ξg ≈ -1.0
                Δ[g] = -1.0
            else
                Δ[g] = 0.0
            end
        end
        data = NodeData(
            w=(1,PerNodeValue(weights)),
            x=(2,PerNodeValue(x)),
            y=(2,PerNodeValue(y)),
            z=(2,PerNodeValue(z)),
            𝑤=(2,PerNodeValue(𝑤)),
            n₁=(2,PerNodeValue(n₁)),
            n₂=(2,PerNodeValue(n₂)),
            n¹=(2,PerNodeValue(n¹)),
            n²=(2,PerNodeValue(n²)),
            s₁=(2,PerNodeValue(s₁)),
            s₂=(2,PerNodeValue(s₂)),
            s¹=(2,PerNodeValue(s¹)),
            s²=(2,PerNodeValue(s²)),
            ∂₁n₁=(2,PerNodeValue(∂₁n₁)),
            ∂₁n₂=(2,PerNodeValue(∂₁n₂)),
            ∂₂n₁=(2,PerNodeValue(∂₂n₁)),
            ∂₂n₂=(2,PerNodeValue(∂₂n₂)),
            ∂₁s₁=(2,PerNodeValue(∂₁s₁)),
            ∂₁s₂=(2,PerNodeValue(∂₁s₂)),
            ∂₂s₁=(2,PerNodeValue(∂₂s₁)),
            ∂₂s₂=(2,PerNodeValue(∂₂s₂)),
            Δ=(1,PerNodeValue(Δ)),
        )
    end
    if dim == 2
        data.value[:ξ] = PerNodeValue(ξ); data.index[:ξ] = 1
        data.value[:η] = PerNodeValue(η); data.index[:η] = 1
    else
        data.value[:ξ] = PerNodeValue(ξ); data.index[:ξ] = 1
    end
end

cal_jacobe = quote
    append!(_vec(data,:𝐽),determinants)
    J = zeros(3,3)
    ∂ξ∂x = zeros(ne*ng)
    ∂ξ∂y = zeros(ne*ng)
    ∂ξ∂z = zeros(ne*ng)
    ∂η∂x = zeros(ne*ng)
    ∂η∂y = zeros(ne*ng)
    ∂η∂z = zeros(ne*ng)
    ∂γ∂x = zeros(ne*ng)
    ∂γ∂y = zeros(ne*ng)
    ∂γ∂z = zeros(ne*ng)
    for C in 1:ne
        for g in 1:ng
            J[1,1] = jacobians[9*(ng*(C-1)+g)-8]
            J[1,2] = jacobians[9*(ng*(C-1)+g)-7]
            J[1,3] = jacobians[9*(ng*(C-1)+g)-6]
            J[2,1] = jacobians[9*(ng*(C-1)+g)-5]
            J[2,2] = jacobians[9*(ng*(C-1)+g)-4]
            J[2,3] = jacobians[9*(ng*(C-1)+g)-3]
            J[3,1] = jacobians[9*(ng*(C-1)+g)-2]
            J[3,2] = jacobians[9*(ng*(C-1)+g)-1]
            J[3,3] = jacobians[9*(ng*(C-1)+g)]
            J⁻¹ = inv(J)
            ∂ξ∂x[ng*(C-1)+g] = J⁻¹[1,1]
            ∂ξ∂y[ng*(C-1)+g] = J⁻¹[1,2]
            ∂ξ∂z[ng*(C-1)+g] = J⁻¹[1,3]
            ∂η∂x[ng*(C-1)+g] = J⁻¹[2,1]
            ∂η∂y[ng*(C-1)+g] = J⁻¹[2,2]
            ∂η∂z[ng*(C-1)+g] = J⁻¹[2,3]
            ∂γ∂x[ng*(C-1)+g] = J⁻¹[3,1]
            ∂γ∂y[ng*(C-1)+g] = J⁻¹[3,2]
            ∂γ∂z[ng*(C-1)+g] = J⁻¹[3,3]
        end
    end
    append!(_vec(data,:∂ξ∂x),∂ξ∂x)
    if dim >= 2
        append!(_vec(data,:∂ξ∂y),∂ξ∂y)
        append!(_vec(data,:∂η∂x),∂η∂x)
        append!(_vec(data,:∂η∂y),∂η∂y)
    end
    if dim >= 3
        append!(_vec(data,:∂ξ∂z),∂ξ∂z)
        append!(_vec(data,:∂η∂z),∂η∂z)
        append!(_vec(data,:∂γ∂x),∂γ∂x)
        append!(_vec(data,:∂γ∂y),∂γ∂y)
        append!(_vec(data,:∂γ∂z),∂γ∂z)
    end
end

typeForFEM = quote
    type = Element{types[elementType]}
end

cal_normal = quote
    if normal
        nodeTags = gmsh.model.mesh.getElementEdgeNodes(elementType,tag,true)
        if dim == 1
            for C in 1:ne
                𝐿 = 2*determinants[C*ng]
                coord, = gmsh.model.mesh.getNode(nodeTags[2*C-1])
                x₁ = coord[1]
                y₁ = coord[2]
                coord, = gmsh.model.mesh.getNode(nodeTags[2*C])
                x₂ = coord[1]
                y₂ = coord[2]
                push!(data.value[:n₁], (y₂-y₁)/𝐿)
                push!(data.value[:n₂], (x₁-x₂)/𝐿)
                push!(data.value[:s₁], (x₂-x₁)/𝐿)
                push!(data.value[:s₂], (y₂-y₁)/𝐿)
            end
        end
        if dim == 2
            nₙ = Int(length(nodeTags)/ne)
            for C in 1:ne
                𝐽 = determinants[C*ng]
                n₁ = 0.0
                n₂ = 0.0
                n₃ = 0.0
                for i in 1:2:nₙ
                    coord, = gmsh.model.mesh.getNode(nodeTags[nₙ*(C-1)+i])
                    x₁ = coord[1]
                    y₁ = coord[2]
                    z₁ = coord[3]
                    coord, = gmsh.model.mesh.getNode(nodeTags[nₙ*(C-1)+i+1])
                    x₂ = coord[1]
                    y₂ = coord[2]
                    z₂ = coord[3]

                    n₁ += y₁*z₂-y₂*z₁
                    n₂ += z₁*x₂-z₂*x₁
                    n₃ += x₁*y₂-x₂*y₁
                end
                if elementType == 3
                    𝐽 *= 8
                end
                push!(data.value[:n₁], n₁/𝐽)
                push!(data.value[:n₂], n₂/𝐽)
                push!(data.value[:n₃], n₃/𝐽)
            end
        end
    end
end

integrationByGmsh = quote
    ~, ~, order, ni = gmsh.model.mesh.getElementProperties(elementType)
    if integrationOrder < 0 integrationOrder = order end
    integrationType = "Gauss"*string(integrationOrder)
    localCoord, weights = gmsh.model.mesh.getIntegrationPoints(elementType,integrationType)
end

integrationByManual = quote
    ~, ~, ~, ni = gmsh.model.mesh.getElementProperties(elementType)
    localCoord, weights = integration
end

generateForFEM = quote
    for C in 1:ne
        𝐶 += 1
        𝓒 = nodes[nodeTag[ni*(C-1)+1:ni*C]]
        𝓖 = [𝑿ₛ((𝑔 = 𝑔+g, 𝐺 = 𝐺+g, 𝐶 = 𝐶, 𝑠 = 𝑠+(g-1)*ni), data) for g in 1:ng]
        𝐺 += ng
        𝑠 += ng*ni
        push!(elements,type(𝓒,𝓖))
    end
    𝑔 += ng
end

generateForNeighbor = quote
    for C in 1:ne
        𝐶 += 1
        indices = Set{Int}()
        for g in 1:ng
            xᵢ = x[ng*(C-1)+g]
            yᵢ = y[ng*(C-1)+g]
            zᵢ = z[ng*(C-1)+g]
            union!(indices,sp(xᵢ,yᵢ,zᵢ))
        end
        ni = length(indices)
        𝓒 = [nodes[i] for i in indices]
        𝓖 = [𝑿ₛ((𝑔 = 𝑔+g, 𝐺 = 𝐺+g, 𝐶 = 𝐶, 𝑠 = 𝑠+(g-1)*ni), data) for g in 1:ng]
        𝐺 += ng
        𝑠 += ng*ni
        push!(elements,type(𝓒,𝓖))
    end
    𝑔 += ng
end

generateForMarco = quote
    elements = Vector{type}(undef,ne)
    data𝓒 = Dict{Symbol,Tuple{Int,Vector{Float64}}}()
    ni = get𝑛𝑝(type(𝑿ᵢ[],𝑿ₛ[]))
    n₁ = Int(round(n/nₕ))
    n₂ = Int(round(ne/nₐ/n₁/nₕ^2))
    for j in 1:n₂
        for i in 1:n₁
            𝓒 = [𝑿ᵢ((𝐼=n₁*ni*(j-1)+ni*(i-1)+k,),data𝓒) for k in 1:ni]
            for k in 1:nₕ
                for l in 1:nₐ*nₕ
                    C = nₐ*nₕ*n₁*(nₕ*(j-1)+k-1)+nₐ*nₕ*(i-1)+l
                    G = ng*(C-1)
                    s = G*ni
                    𝓖 = [𝑿ₛ((𝑔 = g, 𝐺 = G+g, 𝐶 = C, 𝑠 = s+(g-1)*ni), data) for g in 1:ng]
                    elements[C] = type(𝓒,𝓖)
                end
            end
        end
    end
end

generateForPiecewise = quote
    data𝓒 = Dict{Symbol,Tuple{Int,Vector{Float64}}}()
    ni = get𝑛𝑝(type(𝑿ᵢ[],𝑿ₛ[]))
    for C in 1:ne
        𝐶 += 1
        𝓒 = [𝑿ᵢ((𝐼=ni*(𝐶-1)+j,),data𝓒) for j in 1:ni]
        𝓖 = [𝑿ₛ((𝑔 = 𝑔+g, 𝐺 = 𝐺+g, 𝐶 = 𝐶, 𝑠 = 𝑠+(g-1)*ni), data) for g in 1:ng]
        𝐺 += ng
        𝑠 += ng*ni
        push!(elements,type(𝓒,𝓖))
    end
    𝑔 += ng
end

generateForPiecewiseBoundary = quote
    data𝓒 = Dict{Symbol,Tuple{Int,Vector{Float64}}}()
    ni = get𝑛𝑝(type(𝑿ᵢ[],𝑿ₛ[]))
    for CΩ_ in 1:Int(ne/nb)
        tagΩ = tagsΩ[CΩ+CΩ_]
        for C in nb*(CΩ_-1)+1:nb*CΩ_
            𝐶 += 1
            𝓒 = [𝑿ᵢ((𝐼=ni*(CΩ+CΩ_-1)+j,),data𝓒) for j in 1:ni]
            𝓖 = [𝑿ₛ((𝑔 = 𝑔+g, 𝐺 = 𝐺+g, 𝐶 = 𝐶, 𝑠 = 𝑠+(g-1)*ni), data) for g in 1:ng]
            𝐺 += ng
            𝑠 += ng*ni
            push!(elements,type(𝓒,𝓖))
        end
    end
    𝑔 += ng
    CΩ += Int(ne/nb)
end

generateSummary = quote
    println("Info    : Generate $ne elements of $type with $ng integration points.")
end

@eval begin

# function getElements(nodes::Vector{N},dimTag::Pair{Int,Vector{Int}};
#                         type::Union{Int,DataType} = -1,
#                         integration::Union{Int,NTuple{2,Vector{Float64}}} = -1,
#                         searching::Union{Int,SpatialPartition} = -1,
#                         coordinate::Union{Int,Function} = -1,
#                         normal::Bool=false
#                     ) where N<:Node
#     $prequote
#     for (elementType,nodeTag,tag) in zip(elementTypes,nodeTags,tags)
#         if isa(type,Int)
#             type = Element{types[elementType]}
#         end
#         ~, ~, order, ni = gmsh.model.mesh.getElementProperties(elementType)
#         if isa(integration,Int)
#             integrationOrder = integration < 0 : order : integration
#             integrationType = "Gauss"*string(integrationOrder)
#             integration = gmsh.model.mesh.getIntegrationPoints(elementType,integrationType)
#         end
#         localCoord, weights = integration
#         if isa(coordinate,Int)
#             $coordiantes
#         else
#             $curvilinearCoordinates
#         end
#         if isa(searching,Int)
#             if searching == 0
#                 $
#             else
#                 $generateForFEM
#             end
#         else
#             $generateForNeighbor
#         end
#         println("Info: Generate $ne elements of $type with $ng integration points.")
#     end
#     return elements
# end

function getElements(nodes::Vector{N},dimTag::Pair{Int,Vector{Int}},integrationOrder::Int = -1;normal::Bool=false) where N<:Node
    $prequote
    for (elementType,nodeTag,tag) in zip(elementTypes,nodeTags,tags)
        ## element type
        $typeForFEM
        ## integration rule
        $integrationByGmsh
        ## coordinates
        $coordinates
        ## special variables
        $cal_jacobe
        $cal_normal # unit outernal normal
        ## generate element
        $generateForFEM
        ## summary
        $generateSummary
    end
    return elements
end

function getElements(nodes::Vector{N},dimTag::Pair{Int,Vector{Int}},integration::NTuple{2,Vector{Float64}};normal::Bool=false) where N<:Node
    $prequote
    for (elementType,nodeTag,tag) in zip(elementTypes,nodeTags,tags)
        ## element type
        $typeForFEM
        ## integration rule
        $integrationByManual
        ## coordiantes
        $coordinates
        ## special variables
        $cal_jacobe
        $cal_normal # unit outernal normal
        ## generate element
        $generateForFEM
        ## summary
        $generateSummary
    end
    return elements
end

function getElements(nodes::Vector{N},dimTag::Pair{Int,Vector{Int}},type::DataType,integrationOrder::Int = -1;normal::Bool=false) where N<:Node
    $prequote
    for (elementType,nodeTag,tag) in zip(elementTypes,nodeTags,tags)
        ## integration rule
        $integrationByGmsh
        ## coordinates
        $coordinates
        ## special variables
        $cal_jacobe
        $cal_normal # unit outernal normal
        ## generate element
        $generateForFEM
        ## summary
        $generateSummary
    end
    return elements
end

function getElements(nodes::Vector{N},dimTag::Pair{Int,Vector{Int}},type::DataType,integration::NTuple{2,Vector{Float64}};normal::Bool=false) where N<:Node
    $prequote
    for (elementType,nodeTag,tag) in zip(elementTypes,nodeTags,tags)
        ## integration rule
        $integrationByManual
        ## coordiantes
        $coordinates
        ## special variables
        $cal_jacobe
        $cal_normal # unit outernal normal
        ## generate element
        $generateForFEM
        ## summary
        $generateSummary
    end
    return elements
end

function getElements(nodes::Vector{N},dimTag::Pair{Int,Vector{Int}},type::DataType,integrationOrder::Int,sp::SpatialPartition;normal::Bool=false) where N<:Node
    $prequote
    for (elementType,nodeTag,tag) in zip(elementTypes,nodeTags,tags)
        ## integration rule
        $integrationByGmsh
        ## coordinates
        $coordinates
        ## special variables
        $cal_jacobe
        $cal_normal # unit outernal normal
        ## generate element
        $generateForNeighbor
        ## summary
        $generateSummary
    end
    return elements
end

function getElements(nodes::Vector{N},dimTag::Pair{Int,Vector{Int}},type::DataType,integration::NTuple{2,Vector{Float64}},sp::SpatialPartition;normal::Bool=false) where N<:Node
    $prequote
    for (elementType,nodeTag,tag) in zip(elementTypes,nodeTags,tags)
        ## integration rule
        $integrationByManual
        ## coordinates
        $coordinates
        ## special variables
        $cal_jacobe
        $cal_normal # unit outernal normal
        ## generate element
        $generateForNeighbor
        ## summary
        $generateSummary
    end
    return elements
end

function getPiecewiseElements(dimTag::Pair{Int,Vector{Int}},type::DataType,integrationOrder::Int;normal::Bool=false)
    $prequote
    for (elementType,nodeTag,tag) in zip(elementTypes,nodeTags,tags)
        ## integration rule
        $integrationByGmsh
        ## coordinates
        $coordinates
        ## special variables
        $cal_jacobe
        ## generate element
        $generateForPiecewise
        ## summary
        $generateSummary
    end
    return elements
end

function getPiecewiseElements(dimTag::Pair{Int,Vector{Int}},type::DataType,integration::NTuple{2,Vector{Float64}};normal::Bool=false)
    $prequote
    for (elementType,nodeTag,tag) in zip(elementTypes,nodeTags,tags)
        ## integration rule
        $integrationByManual
        ## coordinates
        $coordinates
        ## special variables
        $cal_jacobe
        ## generate element
        $generateForPiecewise
        ## summary
        $generateSummary
    end
    return elements
end

function getPiecewiseBoundaryElements(dimTag::Pair{Int,Vector{Int}},dimTagΩ::Pair{Int,Vector{Int}},type::DataType,integrationOrder::Int)
    normal = false
    $prequote
    $preForEdge
    for (elementType,elementTypeΩ,nodeTag,tag) in zip(elementTypes,elementTypesΩ,nodeTags,tags)
        ## integration rule
        $integrationByGmsh
        ## coordinates
        $coordinatesForEdges
        ## special variables
        $cal_jacobe 
        ## generate element
        $generateForPiecewiseBoundary
        ## summary
        $generateSummary
    end
    return elements
end

function getPiecewiseBoundaryElements(dimTag::Pair{Int,Vector{Int}},dimTagΩ::Pair{Int,Vector{Int}},type::DataType,integration::NTuple{2,Vector{Float64}})
    normal = false
    $prequote
    $preForEdge
    for (elementType,elementTypeΩ,nodeTag,tag) in zip(elementTypes,elementTypesΩ,nodeTags,tags)
        ## integration rule
        $integrationByManual
        ## coordinates
        $coordinatesForEdges
        ## special variables
        $cal_jacobe 
        ## generate element
        $generateForPiecewiseBoundary
        ## summary
        $generateSummary
    end
    return elements
end

function getCurvedElements(nodes::Vector{N},dimTag::Pair{Int,Vector{Int}},cs::Function,integrationOrder::Int = -1) where N<:Node
    $prequote
    for (elementType,nodeTag,tag) in zip(elementTypes,nodeTags,tags)
        ## element type
        $typeForFEM
        ## integration rule
        $integrationByGmsh
        ## coordinates
        $curvilinearCoordinates
        ## special variables
        ## generate element
        $generateForFEM
        ## summary
        $generateSummary
    end
    return elements
end

function getCurvedElements(nodes::Vector{N},dimTag::Pair{Int,Vector{Int}},cs::Function,integration::NTuple{2,Vector{Float64}}) where N<:Node
    $prequote
    for (elementType,nodeTag,tag) in zip(elementTypes,nodeTags,tags)
        ## element type
        $typeForFEM
        ## integration rule
        $integrationByManual
        ## coordinates
        $curvilinearCoordinates
        ## special variables
        ## generate element
        $generateForFEM
        ## summary
        $generateSummary
    end
    return elements
end

function getCurvedElements(nodes::Vector{N},dimTag::Pair{Int,Vector{Int}},type::DataType,cs::Function,integrationOrder::Int,sp::SpatialPartition) where N<:Node
    $prequote
    for (elementType,nodeTag,tag) in zip(elementTypes,nodeTags,tags)
        ## integration rule
        $integrationByGmsh
        ## coordinates
        $curvilinearCoordinates
        ## special variables
        ## generate element
        $generateForNeighbor
        ## summary
        $generateSummary
    end
    return elements
end

function getCurvedElements(nodes::Vector{N},dimTag::Pair{Int,Vector{Int}},type::DataType,cs::Function,integration::NTuple{2,Vector{Float64}},sp::SpatialPartition) where N<:Node
    $prequote
    for (elementType,nodeTag,tag) in zip(elementTypes,nodeTags,tags)
        ## integration rule
        $integrationByManual
        ## coordinates
        $curvilinearCoordinates
        ## special variables
        ## generate element
        $generateForNeighbor
        ## summary
        $generateSummary
    end
    return elements
end

function getCurvedPiecewiseElements(dimTag::Pair{Int,Vector{Int}},type::DataType,cs::Function,integrationOrder::Int,nb::Int=1)
    $prequote
    for (elementType,nodeTag,tag) in zip(elementTypes,nodeTags,tags)
        ## integration rule
        $integrationByGmsh
        ## coordinates
        $curvilinearCoordinates
        ## special variables
        ## generate element
        $generateForPiecewise
        ## summary
        $generateSummary
    end
    return elements
end

function getCurvedPiecewiseElements(dimTag::Pair{Int,Vector{Int}},type::DataType,cs::Function,integration::NTuple{2,Vector{Float64}},nb::Int=1)
    $prequote
    for (elementType,nodeTag,tag) in zip(elementTypes,nodeTags,tags)
        ## integration rule
        $integrationByManual
        ## coordinates
        $curvilinearCoordinates
        ## special variables
        ## generate element
        $generateForPiecewise
        ## summary
        $generateSummary
    end
    return elements
end


# function getMacroCurvedElements(dimTag::Tuple{Int,Int},type::DataType,integrationOrder::Int,n::Int;nₕ::Int=1,nₐ::Int=2)
#     $prequote
#     for (elementType,nodeTag) in zip(elementTypes,nodeTags)
#         ## integration rule
#         $integrationByGmsh
#         ## coordinates
#         $coordinates
#         ## special variables
#         $cal_length_area_volume # length area and volume
#         ## generate element
#         $generateForPiecewise
#         ## summary
#         $generateSummary
#     end
#     return elements
# end

# function getMacroCurvedElements(dimTag::Tuple{Int,Int},type::DataType,integration::NTuple{2,Vector{Float64}},n::Int;nₕ::Int=1,nₐ::Int=2)
#     $prequote
#     for (elementType,nodeTag) in zip(elementTypes,nodeTags)
#         ## integration rule
#         $integrationByManual
#         ## coordinates
#         $coordinates
#         ## special variables
#         $cal_length_area_volume # length area and volume
#         ## generate element
#         $generateForPiecewise
#         ## summary
#         $generateSummary
#     end
#     return elements
# end

# function getMacroBoundaryElements(dimTag::Tuple{Int,Int},dimTagΩ::Tuple{Int,Int},type::DataType,integrationOrder::Int,n::Int;nₕ::Int=1,nₐ::Int=6)
#     $prequote
#     for (elementType,nodeTag) in zip(elementTypes,nodeTags)
#         ## integration rule
#         $integrationByGmsh
#         ## coordinates
#         $coordinatesForEdges
#         ## special variables
#         $cal_length_area_volume # length area and volume
#         ## generate element
#         $generateForPiecewise
#         ## summary
#         $generateSummary
#     end
#     return elements
# end

# function getMacroBoundaryElements(dimTag::Tuple{Int,Int},dimTagΩ::Tuple{Int,Int},type::DataType,integration::NTuple{2,Vector{Float64}},n::Int;nₕ::Int=1,nₐ::Int=6)
#     $prequote
#     for (elementType,nodeTag) in zip(elementTypes,nodeTags)
#         ## integration rule
#         $integrationByManual
#         ## coordinates
#         $coordinatesForEdges
#         ## special variables
#         $cal_length_area_volume # length area and volume
#         ## generate element
#         $generateForPiecewise
#         ## summary
#         $generateSummary
#     end
#     return elements
# end

end 

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