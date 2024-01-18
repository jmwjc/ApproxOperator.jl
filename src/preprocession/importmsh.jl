
function getPhysicalGroups()
    entities = Dict{String,Tuple{Int,Int}}()
    dimTags = gmsh.model.getPhysicalGroups()
    for (dim,tag) in dimTags
        name = gmsh.model.getPhysicalName(dim,tag)
        tags = gmsh.model.getEntitiesForPhysicalGroup(dim,tag)
        entities[name] = (dim,tags[1])
    end
    return entities
end

function get𝑿ᵢ()
    nodeTags, coord = gmsh.model.mesh.getNodes()
    nₚ = length(nodeTags)
    x = coord[1:3:3*nₚ]
    y = coord[2:3:3*nₚ]
    z = coord[3:3:3*nₚ]
    data = Dict([:x=>(1,x),:y=>(1,y),:z=>(1,z)])
    return [𝑿ᵢ((𝐼=i,),data) for i in 1:nₚ]
end

prequote = quote
    types = Dict([1=>:Seg2, 2=>:Tri3, 3=>:Quad, 4=>:Tet4, 8=>:Seg3, 9=>:Tri6, 10=>:Quad9, 11=>:Tet10, 15=>:Poi1, 16=>Quad8])
    dim, tag = dimTag
    elementTypes, ~, nodeTags = gmsh.model.mesh.getElements(dim,tag)
    elements = AbstractElement[]
end

coordinates = quote
    ng = length(weights)
    ne = Int(length(nodeTag)/ni)

    ξ = localCoord[1:3:end]
    η = localCoord[2:3:end]
    γ = localCoord[3:3:end]
    jacobians, determinants, coord = gmsh.model.mesh.getJacobians(elementType, localCoord, tag)
    x = coord[1:3:end]
    y = coord[2:3:end]
    z = coord[3:3:end]
    𝑤 = [weight*determinant for determinant in determinants for weight in weights]
    data = Dict([
        :w=>(1,weights),
        :x=>(2,x),
        :y=>(2,y),
        :z=>(2,z),
        :𝑤=>(2,𝑤),
    ])
    if dim == 3
        push!(data, :ξ=>(1,ξ), :η=>(1,η), :γ=>(1,γ))
    elseif dim == 2
        push!(data, :ξ=>(1,ξ), :η=>(1,η))
    else
        push!(data, :ξ=>(1,ξ))
    end
end

coordinatesForEdges = quote
    ng = length(weights)
    ne = Int(length(nodeTag)/ni)

    ξ = zeros(ne*ng)
    η = zeros(ne*ng)
    γ = zeros(ne*ng)
    n₁ = zeros(ne)
    n₂ = zeros(ne)
    s₁ = zeros(ne)
    s₂ = zeros(ne)
    Δ = zeros(ng)
    jacobians, determinants, coord = gmsh.model.mesh.getJacobians(elementType, localCoord, tag)
    x = coord[1:3:end]
    y = coord[2:3:end]
    z = coord[3:3:end]
    𝑤 = [weight*determinant for determinant in determinants for weight in weights]

    for g in 1:ng
        ξg = localCoord[3*g-2]
        if abs(ξg) ≈ 1.0 Δ[g] = 1.0 end
    end

    nodeTags = gmsh.model.mesh.getElementEdgeNodes(elementType,tag,true)
    dimΩ,tagΩ = dimTagΩ
    ~, tagsΩ = gmsh.model.mesh.getElements(dimΩ,tagΩ)
    for (CΩ,tagΩ) in enumerate(tagsΩ[1])
        for C in 3*CΩ-2:3*CΩ
            𝐿 = 2*determinants[C*ng]
            coord, = gmsh.model.mesh.getNode(nodeTags[2*C-1])
            x₁ = coord[1]
            y₁ = coord[2]
            coord, = gmsh.model.mesh.getNode(nodeTags[2*C])
            x₂ = coord[1]
            y₂ = coord[2]
            n₁[C] = (y₂-y₁)/𝐿
            n₂[C] = (x₁-x₂)/𝐿
            s₁[C] = -n₂[C]
            s₂[C] =  n₁[C]
            for g in 1:ng
                G = ng*(C-1)+g
                ξ[G], η[G], γ[G] = gmsh.model.mesh.getLocalCoordinatesInElement(tagΩ, x[G], y[G], z[G])
            end
        end
    end
    data = Dict([
        :w=>(1,weights),
        :x=>(2,x),
        :y=>(2,y),
        :z=>(2,z),
        :𝑤=>(2,𝑤),
        :n₁=>(3,n₁),
        :n₂=>(3,n₂),
        :s₁=>(3,s₁),
        :s₂=>(3,s₂),
        :Δ=>(1,Δ),
    ])
    if dim == 2
        push!(data, :ξ=>(1,ξ), :η=>(1,η), :γ=>(1,γ))
    else
        push!(data, :ξ=>(1,ξ), :η=>(1,η))
    end
end

cal_length_area_volume = quote
    if elementType == 1
        𝐿 = [2*determinants[C*ng] for C in 1:ne]
        push!(data, :𝐿=>(3,𝐿))
    elseif elementType == 2
        𝐴 = [determinants[C*ng]/2 for C in 1:ne]
        push!(data, :𝐴=>(3,𝐴))
    end
end

typeForFEM = quote
    type = Element{types[elementType]}
end

cal_normal = quote
    if normal
        nodeTags = gmsh.model.mesh.getElementEdgeNodes(elementType,tag,true)
        if dim == 1
            n₁ = zeros(ne)
            n₂ = zeros(ne)
            for C in 1:ne
                𝐿 = 2*determinants[C*ng]
                coord, = gmsh.model.mesh.getNode(nodeTags[2*C-1])
                x₁ = coord[1]
                y₁ = coord[2]
                coord, = gmsh.model.mesh.getNode(nodeTags[2*C])
                x₂ = coord[1]
                y₂ = coord[2]
                n₁[C] = (y₂-y₁)/𝐿
                n₂[C] = (x₁-x₂)/𝐿
            end
            push!(data,:n₁=>(3,n₁),:n₂=>(3,n₂))
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
    G = 0
    s = 0
    for C in 1:ne
        𝓒 = nodes[nodeTag[ni*(C-1)+1:ni*C]]
        𝓖 = [𝑿ₛ((𝑔 = g, 𝐺 = G+g, 𝐶 = C, 𝑠 = s+(g-1)*ni), data) for g in 1:ng]
        G += ng
        s += ng*ni
        push!(elements,type(𝓒,𝓖))
    end
end

generateForNeighbor = quote
    G = 0
    s = 0
    for C in 1:ne
        indices = Set{Int}()
        for g in 1:ng
            xᵢ = x[G+g]
            yᵢ = y[G+g]
            zᵢ = z[G+g]
            union!(indices,sp(xᵢ,yᵢ,zᵢ))
        end
        ni = length(indices)
        𝓒 = [nodes[i] for i in indices]
        𝓖 = [𝑿ₛ((𝑔 = g, 𝐺 = G+g, 𝐶 = C, 𝑠 = s+(g-1)*ni), data) for g in 1:ng]
        G += ng
        s += ng*ni
        push!(elements,type(𝓒,𝓖))
    end
end

generateForPiecewise = quote
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

generateSummary = quote
    println("Info: Generate $ne elements of $type with $ng integration points.")
end

@eval begin

function getElements(nodes::Vector{N},dimTag::Pair{Int,Int},integrationOrder::Int = -1;normal::Bool=false) where N<:Node
    $prequote
    for (elementType,nodeTag) in zip(elementTypes,nodeTags)
        ## element type
        $typeForFEM
        ## integration rule
        $integrationByGmsh
        ## coordinates
        $coordinates
        ## special variables
        $cal_length_area_volume # length area and volume
        $cal_normal # unit outernal normal
        ## generate element
        $generateForFEM
        ## summary
        $generateSummary
    end
    return elements
end

function getElements(nodes::Vector{N},dimTag::Tuple{Int,Int},integration::NTuple{2,Vector{Float64}};normal::Bool=false) where N<:Node
    $prequote
    for (elementType,nodeTag) in zip(elementTypes,nodeTags)
        ## element type
        $typeForFEM
        ## integration rule
        $integrationByManual
        ## coordiantes
        $coordinates
        ## special variables
        $cal_length_area_volume # length area and volume
        $cal_normal # unit outernal normal
        ## generate element
        $generateForFEM
        ## summary
        $generateSummary
    end
    return elements
end

function getElements(nodes::Vector{N},dimTag::Tuple{Int,Int},type::DataType,integrationOrder::Int = -1;normal::Bool=false) where N<:Node
    $prequote
    for (elementType,nodeTag) in zip(elementTypes,nodeTags)
        ## integration rule
        $integrationByGmsh
        ## coordinates
        $coordinates
        ## special variables
        $cal_length_area_volume # length area and volume
        $cal_normal # unit outernal normal
        ## generate element
        $generateForFEM
        ## summary
        $generateSummary
    end
    return elements
end

function getElements(nodes::Vector{N},dimTag::Tuple{Int,Int},type::DataType,integration::NTuple{2,Vector{Float64}};normal::Bool=false) where N<:Node
    $prequote
    for (elementType,nodeTag) in zip(elementTypes,nodeTags)
        ## integration rule
        $integrationByManual
        ## coordiantes
        $coordinates
        ## special variables
        $cal_length_area_volume # length area and volume
        $cal_normal # unit outernal normal
        ## generate element
        $generateForFEM
        ## summary
        $generateSummary
    end
    return elements
end

function getElements(nodes::Vector{N},dimTag::Tuple{Int,Int},type::DataType,integrationOrder::Int,sp::SpatialPartition;normal::Bool=false) where N<:Node
    $prequote
    for (elementType,nodeTag) in zip(elementTypes,nodeTags)
        ## integration rule
        $integrationByGmsh
        ## coordinates
        $coordinates
        ## special variables
        $cal_length_area_volume # length area and volume
        $cal_normal # unit outernal normal
        ## generate element
        $generateForNeighbor
        ## summary
        $generateSummary
    end
    return elements
end

function getElements(nodes::Vector{N},dimTag::Tuple{Int,Int},type::DataType,integration::NTuple{2,Vector{Float64}},sp::SpatialPartition;normal::Bool=false) where N<:Node
    $prequote
    for (elementType,nodeTag) in zip(elementTypes,nodeTags)
        ## integration rule
        $integrationByManual
        ## coordinates
        $coordinates
        ## special variables
        $cal_length_area_volume # length area and volume
        $cal_normal # unit outernal normal
        ## generate element
        $generateForNeighbor
        ## summary
        $generateSummary
    end
    return elements
end

function getMacroElementsForTriangles(dimTag::Tuple{Int,Int},type::DataType,integration::NTuple{2,Vector{Float64}},n::Int)
    $prequote
    nₕ = 2; nₐ = 2
    for (elementType,nodeTag) in zip(elementTypes,nodeTags)
        ## integration rule
        $integrationByManual
        ## coordinates
        $coordinates
        ## special variables
        $cal_length_area_volume # length area and volume
        ## generate element
        $generateForPiecewise
        ## summary
        $generateSummary
    end
    return elements
end

function getMacroBoundaryElementsForTriangles(dimTag::Tuple{Int,Int},dimTagΩ::Tuple{Int,Int},type::DataType,integration::NTuple{2,Vector{Float64}},n::Int)
    $prequote
    nₕ = 2; nₐ = 6
    for (elementType,nodeTag) in zip(elementTypes,nodeTags)
        ## integration rule
        $integrationByManual
        ## coordinates
        $coordinatesForEdges
        ## special variables
        $cal_length_area_volume # length area and volume
        ## generate element
        $generateForPiecewise
        ## summary
        $generateSummary
    end
    return elements
end

end 
