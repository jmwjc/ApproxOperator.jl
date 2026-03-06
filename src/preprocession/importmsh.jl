

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
    data = Dict([:x=>(1,x),:y=>(1,y),:z=>(1,z)])
    return [𝑿ᵢ((𝐼=i,), data) for i in 1:nₚ ]
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
    data = Dict{Symbol,Tuple{Int,Vector{Float64}}}()
    data[:w] = (1,Float64[])
    data[:ξ] = (1,Float64[])
    data[:x] = (2,Float64[])
    data[:y] = (2,Float64[])
    data[:z] = (2,Float64[])
    data[:𝑤] = (2,Float64[])
    data[:𝐽] = (2,Float64[])
    data[:∂ξ∂x] = (2,Float64[])
    if normal
        data[:n₁] = (3,Float64[])
        data[:n₂] = (3,Float64[])
        data[:s₁] = (3,Float64[])
        data[:s₂] = (3,Float64[])
    end
    if dim == 2
        data[:η] = (1,Float64[])

        data[:∂ξ∂y] = (2,Float64[])
        data[:∂η∂x] = (2,Float64[])
        data[:∂η∂y] = (2,Float64[])
        if normal
            data[:n₃] = (3,Float64[])
            data[:s₃] = (3,Float64[])
        end
    end
    if dim == 3
        data[:γ] = (1,Float64[])
        data[:η] = (1,Float64[])

        
        data[:∂ξ∂y] = (2,Float64[])
        data[:∂η∂x] = (2,Float64[])
        data[:∂η∂y] = (2,Float64[])

        data[:∂ξ∂z] = (2,Float64[])
        data[:∂η∂z] = (2,Float64[])
        data[:∂γ∂x] = (2,Float64[])
        data[:∂γ∂y] = (2,Float64[])
        data[:∂γ∂z] = (2,Float64[])
    end

    if normal
        
        data[:n₁] = (3,Float64[])
        data[:n₂] = (3,Float64[])
        data[:n₃] = (3,Float64[])
        data[:s₁] = (3,Float64[])
        data[:s₂] = (3,Float64[])
        data[:s₃] = (3,Float64[])
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

    data[:w] = (1,Float64[])
    data[:Δ] = (1,Float64[])
    data[:ξ] = (2,Float64[])
    data[:η] = (2,Float64[])
    data[:n₁] = (3,Float64[])
    data[:n₂] = (3,Float64[])
    data[:s₁] = (3,Float64[])
    data[:s₂] = (3,Float64[])
    if dim > 1
        data[:γ] = (2,Float64[])
        data[:n₃] = (3,Float64[])
        data[:s₃] = (3,Float64[])
    end   
end

coordinates = quote
    ng = length(weights)
    ne = Int(length(nodeTag)/ni)

    append!(data[:w][2],weights)
    haskey(data,:ξ) ? append!(data[:ξ][2],localCoord[1:3:end]) : nothing
    haskey(data,:η) ? append!(data[:η][2],localCoord[2:3:end]) : nothing
    haskey(data,:γ) ? append!(data[:γ][2],localCoord[3:3:end]) : nothing
    jacobians, determinants, coord = gmsh.model.mesh.getJacobians(elementType, localCoord, tag)
    x = coord[1:3:end]
    y = coord[2:3:end]
    z = coord[3:3:end]
    append!(data[:x][2],x)
    append!(data[:y][2],y)
    append!(data[:z][2],z)
    for i in 1:Int(length(determinants)/ng)
        for (j,w) in enumerate(weights)
            G = ng*(i-1)+j
            push!(data[:𝑤][2], determinants[G]*w)
        end
    end
end

coordinatesForFaces2 = quote
    ng = length(weights)
    ne = Int(length(nodeTag)/ni)
    if elementTypeΩ ∈ (4)
        face_type = 2
        nb = 4  # Number of faces per element
        nf = 3  # Nodes per face
    elseif elementTypeΩ ∈ (5)
        face_type = 3
        nb = 6  # Number of faces per element
        nf = 4  # Nodes per face
    end

    append!(data[:w][2],weights)
    haskey(data,:ξ) ? append!(data[:ξ][2],localCoord[1:3:end]) : nothing
    haskey(data,:η) ? append!(data[:η][2],localCoord[2:3:end]) : nothing
    haskey(data,:γ) ? append!(data[:γ][2],localCoord[3:3:end]) : nothing
    jacobians, determinants, coord = gmsh.model.mesh.getJacobians(elementType, localCoord, tag)
    x = coord[1:3:end]
    y = coord[2:3:end]
    z = coord[3:3:end]
    append!(data[:x][2],x)
    append!(data[:y][2],y)
    append!(data[:z][2],z)
    for i in 1:Int(length(determinants)/ng)
        for (j,w) in enumerate(weights)
            G = ng*(i-1)+j
            push!(data[:𝑤][2], determinants[G]*w)
        end
    end
end


# coordinatesForFaces = quote
#     if elementTypeΩ ∈ (4)
#         face_type = 2
#         nb = 4  # Number of faces per element
#         nf = 3  # Nodes per face
#     elseif elementTypeΩ ∈ (5)
#         face_type = 3
#         nb = 6  # Number of faces per element
#         nf = 4  # Nodes per face
#     end

#     ng = length(weights)
#     ne = Int(length(nodeTag)/ni)
#     println(ne)
#     println(ni)
#     println(nb)
#     # nodeTag = gmsh.model.mesh.getElementFaceNodes(elementTypeΩ,nf,tagΩ_,true)
   
#     append!(data[:w][2], weights)

#     # Get face Jacobians using face element type
#     jacobians, determinants, coord = gmsh.model.mesh.getJacobians(elementType, localCoord, tag)
#     x = coord[1:3:end]
#     y = coord[2:3:end]
#     z = coord[3:3:end]
#     append!(data[:x][2],x)
#     append!(data[:y][2],y)
#     append!(data[:z][2],z)

#     # Calculate weighted determinants
#     for i in 1:Int(length(determinants)/ng)
#         for (j,w) in enumerate(weights)
#             G = ng*(i-1)+j
#             push!(data[:𝑤][2], determinants[G]*w)
#         end
#     end

#     # for g in 1:ng
#     #     u = localCoord[2g-1]  # Surface parametric coordinates
#     #     v = localCoord[2g]
#     #     # Example boundary detection (modify according to your needs)
#     #     Δ = ifelse(u ≈ 1.0 || v ≈ 1.0 || u+v ≈ 1.0, 1.0, 0.0)
#     #     push!(data[:Δ][2], Δ)
#     # end
    
#     tolerance = 1e-6
#     for g in 1:ng
#         ξg = localCoord[3*g-2]
#         ηg = localCoord[3*g-1]
#         if abs(ξg - 1.0) < tolerance && abs(ηg - 0.0) < tolerance
#             push!(data[:Δ][2], 1.0)
#         elseif abs(ξg + 1.0) < tolerance && abs(ηg - 0.0) < tolerance
#             push!(data[:Δ][2], -1.0)
#         else
#             push!(data[:Δ][2], 0.0)
#         end
#     end

#     # Face normal calculation
#     for CΩ_ in 1:Int(ne/nb)
#         tagΩ = tagsΩ[CΩ+CΩ_]
#         for C in nb*(CΩ_-1)+1:nb*CΩ_
#         @assert length(nodeTag) >= C*nf "nodeTag 长度不足"
#             face_nodes = nodeTag[(C-1)*nf+1:C*nf]

#             coord1, = gmsh.model.mesh.getNode(face_nodes[1])
#             coord2, = gmsh.model.mesh.getNode(face_nodes[2])
#             coord3, = gmsh.model.mesh.getNode(face_nodes[3])
#             v1 = [coord2[1]-coord1[1], coord2[2]-coord1[2], coord2[3]-coord1[3]]  # 边向量1
#             v2 = [coord3[1]-coord1[1], coord3[2]-coord1[2], coord3[3]-coord1[3]]  # 边向量2
            
        
#             # 手动计算叉积（法向量）
#             n = [
#                 v1[2] * v2[3] - v1[3] * v2[2],  # x 分量
#                 v1[3] * v2[1] - v1[1] * v2[3],  # y 分量
#                 v1[1] * v2[2] - v1[2] * v2[1]   # z 分量
#             ]
        
#             # 归一化（单位法向量）
#             n_length = sqrt(n[1]^2 + n[2]^2 + n[3]^2)
#             normal = [n[1] / n_length, n[2] / n_length, n[3] / n_length]
            
#             # normals = gmsh.model.mesh.getNormals(face_type, face_nodes)
           
#             # append!(data[:n₁][2], normals[1:3:end])
#             # append!(data[:n₂][2], normals[2:3:end])
#             # append!(data[:n₃][2], normals[3:3:end])

#             # push!(data[:n₁][2], normals[1:3:end])
#             # push!(data[:n₂][2], normals[2:3:end])
#             # push!(data[:n₃][2], normals[3:3:end])
#             push!(data[:n₁][2], normal[1])
#             push!(data[:n₂][2], normal[2])
#             push!(data[:n₃][2], normal[3])

#               # Get local coordinates in parent volume element
#             for g in 1:ng
#                 G = ng*(C-1)+g
#                 ξ, η, γ = gmsh.model.mesh.getLocalCoordinatesInElement(tagΩ, x[G], y[G], z[G])  
#                 push!(data[:ξ][2], ξ)
#                 push!(data[:η][2], η)
#                 push!(data[:γ][2], γ)
#                 # haskey(data,:γ) ? push!(data[:γ][2], γ) : nothing
#             end
#         end
#     end
# end


# coordinatesForFaces = quote
#     # -------------------------
#     # 基本数量
#     # -------------------------
#     ng = length(weights)
#     ne = Int(length(nodeTag) / ni)   # 2D face elements count (from discrete surface entity)

#     # -------------------------
#     # 体单元类型 -> 每个体单元的面数 nb
#     # （Tet4=4 faces, Hex8=6 faces）
#     # -------------------------
#     if elementTypeΩ == 4
#         nb = 4
#     elseif elementTypeΩ == 5
#         nb = 6
#     else
#         error("coordinatesForFaces2: unsupported elementTypeΩ = $elementTypeΩ")
#     end

#     # -------------------------
#     # 写入参考权重、参考坐标
#     # （对 2D 单元：通常 localCoord 是 (u,v,0) 三元组，取前两列即可）
#     # -------------------------
#     append!(data[:w][2], weights)
#     haskey(data, :ξ) ? append!(data[:ξ][2], localCoord[1:3:end]) : nothing
#     haskey(data, :η) ? append!(data[:η][2], localCoord[2:3:end]) : nothing
#     haskey(data, :γ) ? append!(data[:γ][2], localCoord[3:3:end]) : nothing  # 可能全 0，留着不影响

#     # -------------------------
#     # 从 gmsh 一次性拿到：所有面单元、所有高斯点的 Jacobian & 物理坐标
#     # determinants 对 2D 面单元就是 |xu × xv| （面积尺度）
#     # -------------------------
#     jacobians, determinants, coord = gmsh.model.mesh.getJacobians(elementType, localCoord, tag)
#     x = coord[1:3:end]
#     y = coord[2:3:end]
#     z = coord[3:3:end]
#     append!(data[:x][2], x)
#     append!(data[:y][2], y)
#     append!(data[:z][2], z)

#     # -------------------------
#     # 物理权重：det * w
#     # -------------------------
#     for i in 1:Int(length(determinants) / ng)
#         for (j, w) in enumerate(weights)
#             G = ng*(i-1) + j
#             push!(data[:𝑤][2], determinants[G] * w)
#         end
#     end

#     # -------------------------
#     # 法向：优先用 gmsh.getNormals（最稳），失败则用 jacobians 做叉积
#     # 注意：这里 push 的法向是“每个高斯点一个”
#     # -------------------------
#     normals_ok = true
#     normals = Float64[]
#     try
#         normals = gmsh.model.mesh.getNormals(elementType, localCoord, tag)
#     catch
#         normals_ok = false
#     end

#     if normals_ok
#         # normals 是 [nx1, ny1, nz1, nx2, ny2, nz2, ...]
#         for G in 1:(Int(length(coord)/3))  # total gauss points across all face elems
#             push!(data[:n₁][2], normals[3*G-2])
#             push!(data[:n₂][2], normals[3*G-1])
#             push!(data[:n₃][2], normals[3*G])
#         end
#     else
#         # fallback: use surface tangents from jacobians (3x3 per GP)
#         # take first two columns as xu and xv: (∂x/∂u, ∂x/∂v)
#         nGP = Int(length(coord)/3)
#         for G in 1:nGP
#             # jacobians layout: 9 per GP
#             base = 9*(G-1)
#             xu = (jacobians[base+1], jacobians[base+2], jacobians[base+3])
#             xv = (jacobians[base+4], jacobians[base+5], jacobians[base+6])
#             cx = xu[2]*xv[3] - xu[3]*xv[2]
#             cy = xu[3]*xv[1] - xu[1]*xv[3]
#             cz = xu[1]*xv[2] - xu[2]*xv[1]
#             nlen = sqrt(cx*cx + cy*cy + cz*cz)
#             if nlen < 1e-14
#                 push!(data[:n₁][2], 0.0); push!(data[:n₂][2], 0.0); push!(data[:n₃][2], 0.0)
#             else
#                 push!(data[:n₁][2], cx/nlen); push!(data[:n₂][2], cy/nlen); push!(data[:n₃][2], cz/nlen)
#             end
#         end
#     end

#     # -------------------------
#     # （可选）把面上每个高斯点映射回所属体单元的局部坐标 (ξ,η,γ)
#     # 这一步需要 face->parent volume element 的映射。
#     #
#     # 如果你的面单元集合来自 addFaceElements 且未去重：
#     #   ne(face) == nb * ne(volume)
#     #   第 iFace 个面单元属于第 ceil(iFace/nb) 个体单元
#     # -------------------------
#     # 体单元 tags（已由 preForEdge 构造 tagsΩ; 也可在这里现取）
#     # tagsΩ 是所有体单元 element tags 的平铺序列
#     neΩ = length(tagsΩ)
#     @assert ne == nb * neΩ "coordinatesForFaces2: mapping requires ne(face)=nb*ne(volume), got ne=$ne, nb=$nb, neΩ=$neΩ"

#     # 逐面单元逐高斯点映射
#     for iFace in 1:ne
#         iVol = Int(ceil(iFace/nb))
#         tagΩ = tagsΩ[iVol]
#         for g in 1:ng
#             G = ng*(iFace-1) + g
#             ξv, ηv, γv = gmsh.model.mesh.getLocalCoordinatesInElement(tagΩ, x[G], y[G], z[G])
#             push!(data[:ξ][2], ξv)
#             push!(data[:η][2], ηv)
#             push!(data[:γ][2], γv)
#         end
#     end
# end


coordinatesForFaces = quote
    ng = length(weights)
    ne = Int(length(nodeTag)/ni)   # number of 2D face elements in THIS surface entity

    # -------------------------
    # Determine faces-per-volume (nb) and face node count (nf) from volume type
    # -------------------------
    if elementTypeΩ == 4
        nb = 4
        nf = 3
        faceTypeId = 2   # Tri3
    elseif elementTypeΩ == 5
        nb = 6
        nf = 4
        faceTypeId = 3   # Quad4
    else
        error("coordinatesForFaces: unsupported elementTypeΩ = $elementTypeΩ")
    end

    # -------------------------
    # (1) surface integration meta (same style as edges2)
    # -------------------------
    append!(data[:w][2], weights)
    haskey(data,:ξ) ? append!(data[:ξ][2], localCoord[1:3:end]) : nothing
    haskey(data,:η) ? append!(data[:η][2], localCoord[2:3:end]) : nothing
    haskey(data,:γ) ? append!(data[:γ][2], localCoord[3:3:end]) : nothing  # often 0 for surface, keep for consistency

    jacobians, determinants, coord = gmsh.model.mesh.getJacobians(elementType, localCoord, tag)
    x = coord[1:3:end]
    y = coord[2:3:end]
    z = coord[3:3:end]
    append!(data[:x][2], x)
    append!(data[:y][2], y)
    append!(data[:z][2], z)

    for i in 1:Int(length(determinants)/ng)
        for (j,w) in enumerate(weights)
            G = ng*(i-1)+j
            push!(data[:𝑤][2], determinants[G]*w)   # surface weight
        end
    end

    # -------------------------
    # (2) normals (per Gauss point): prefer gmsh normals, else cross product from jacobians
    # -------------------------
    normals_ok = true
    normals = Float64[]
    try
        normals = gmsh.model.mesh.getNormals(elementType, localCoord, tag)
    catch
        normals_ok = false
    end

    nGP = Int(length(coord)/3)  # total Gauss points = ne*ng

    if normals_ok
        for G in 1:nGP
            push!(data[:n₁][2], normals[3*G-2])
            push!(data[:n₂][2], normals[3*G-1])
            push!(data[:n₃][2], normals[3*G])
        end
    else
        # jacobians are 3x3 per GP; take first two columns as tangents xu, xv
        for G in 1:nGP
            base = 9*(G-1)
            xu = (jacobians[base+1], jacobians[base+2], jacobians[base+3])
            xv = (jacobians[base+4], jacobians[base+5], jacobians[base+6])
            cx = xu[2]*xv[3] - xu[3]*xv[2]
            cy = xu[3]*xv[1] - xu[1]*xv[3]
            cz = xu[1]*xv[2] - xu[2]*xv[1]
            nlen = sqrt(cx*cx + cy*cy + cz*cz)
            if nlen < 1e-14
                push!(data[:n₁][2], 0.0); push!(data[:n₂][2], 0.0); push!(data[:n₃][2], 0.0)
            else
                push!(data[:n₁][2], cx/nlen); push!(data[:n₂][2], cy/nlen); push!(data[:n₃][2], cz/nlen)
            end
        end
    end

    # -------------------------
    # (3) map each face Gauss point back to its parent volume element local coords
    # same pattern as your edges code: group faces into blocks of nb
    # -------------------------
    for CΩ_ in 1:Int(ne/nb)
        tagΩ = tagsΩ[CΩ + CΩ_]

        for C in nb*(CΩ_-1)+1 : nb*CΩ_
            for g in 1:ng
                G = ng*(C-1) + g
                ξv, ηv, γv = gmsh.model.mesh.getLocalCoordinatesInElement(tagΩ, x[G], y[G], z[G])
                push!(data[:ξ][2], ξv)
                push!(data[:η][2], ηv)
                push!(data[:γ][2], γv)
            end
        end
    end
end



coordinatesForEdges2 = quote
    ng = length(weights)
    ne = Int(length(nodeTag)/ni)
    if elementTypeΩ ∈ (2,9)
        nb = 3
    elseif elementTypeΩ ∈ (3,4,10,16)
        nb = 4
    end

    append!(data[:w][2],weights)
    haskey(data,:ξ) ? append!(data[:ξ][2],localCoord[1:3:end]) : nothing
    haskey(data,:η) ? append!(data[:η][2],localCoord[2:3:end]) : nothing
    haskey(data,:γ) ? append!(data[:γ][2],localCoord[3:3:end]) : nothing
    jacobians, determinants, coord = gmsh.model.mesh.getJacobians(elementType, localCoord, tag)
    x = coord[1:3:end]
    y = coord[2:3:end]
    z = coord[3:3:end]
    append!(data[:x][2],x)
    append!(data[:y][2],y)
    append!(data[:z][2],z)
    for i in 1:Int(length(determinants)/ng)
        for (j,w) in enumerate(weights)
            G = ng*(i-1)+j
            push!(data[:𝑤][2], determinants[G]*w)
        end
    end
end

coordinatesForEdges = quote
    ng = length(weights)
    ne = Int(length(nodeTag)/ni)
    if elementTypeΩ ∈ (2,9)
        nb = 3
    elseif elementTypeΩ ∈ (3,4,10,16)
        
        nb = 4
    end
    nodeTag = gmsh.model.mesh.getElementEdgeNodes(elementTypeΩ,tagΩ_,true)
    
    println(ne)
    println(ni)
    println(nb)
    append!(data[:w][2],weights)
    jacobians, determinants, coord = gmsh.model.mesh.getJacobians(elementType, localCoord, tag)
    x = coord[1:3:end]
    y = coord[2:3:end]
    z = coord[3:3:end]
    append!(data[:x][2],x)
    append!(data[:y][2],y)
    append!(data[:z][2],z)

    for i in 1:Int(length(determinants)/ng)
        for (j,w) in enumerate(weights)
            G = ng*(i-1)+j
            push!(data[:𝑤][2], determinants[G]*w)
        end
    end

    for g in 1:ng
        ξg = localCoord[3*g-2]
        if ξg ≈ 1.0
            push!(data[:Δ][2], 1.0)
        elseif ξg ≈ -1.0
            push!(data[:Δ][2], -1.0)
        else
            push!(data[:Δ][2], 0.0)
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
            
            push!(data[:n₁][2], (y₂-y₁)/𝐿)
            push!(data[:n₂][2], (x₁-x₂)/𝐿)
            push!(data[:s₁][2], (x₂-x₁)/𝐿)
            push!(data[:s₂][2], (y₂-y₁)/𝐿)
           
            for g in 1:ng
                G = ng*(C-1)+g
                ξ, η, γ = gmsh.model.mesh.getLocalCoordinatesInElement(tagΩ, x[G], y[G], z[G])
                push!(data[:ξ][2], ξ)
                push!(data[:η][2], η)
                haskey(data,:γ) ? push!(data[:γ][2], γ) : nothing
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
        data = Dict([
            :w=>(1,weights),
            :x=>(2,x),
            :y=>(2,y),
            :z=>(2,z),
            :𝑤=>(2,𝑤),
        ])
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
        data = Dict([
            :w=>(1,weights),
            :x=>(2,x),
            :y=>(2,y),
            :z=>(2,z),
            :𝑤=>(2,𝑤),
            :n₁=>(2,n₁),
            :n₂=>(2,n₂),
            :n¹=>(2,n¹),
            :n²=>(2,n²),
            :s₁=>(2,s₁),
            :s₂=>(2,s₂),
            :s¹=>(2,s¹),
            :s²=>(2,s²),
            :∂₁n₁=>(2,∂₁n₁),
            :∂₁n₂=>(2,∂₁n₂),
            :∂₂n₁=>(2,∂₂n₁),
            :∂₂n₂=>(2,∂₂n₂),
            :∂₁s₁=>(2,∂₁s₁),
            :∂₁s₂=>(2,∂₁s₂),
            :∂₂s₁=>(2,∂₂s₁),
            :∂₂s₂=>(2,∂₂s₂),
            :Δ=>(1,Δ),
        ])
    end
    if dim == 2
        push!(data, :ξ=>(1,ξ), :η=>(1,η))
    else
        push!(data, :ξ=>(1,ξ))
    end
end

cal_jacobe = quote
    append!(data[:𝐽][2],determinants)
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
    append!(data[:∂ξ∂x][2],∂ξ∂x)
    if dim == 2
        append!(data[:∂ξ∂y][2],∂ξ∂y)
        append!(data[:∂η∂x][2],∂η∂x)
        append!(data[:∂η∂y][2],∂η∂y)
    elseif dim == 3
        append!(data[:∂ξ∂y][2],∂ξ∂y)
        append!(data[:∂η∂x][2],∂η∂x)
        append!(data[:∂η∂y][2],∂η∂y)

        append!(data[:∂ξ∂z][2],∂ξ∂z)
        append!(data[:∂η∂z][2],∂η∂z)
        append!(data[:∂γ∂x][2],∂γ∂x)
        append!(data[:∂γ∂y][2],∂γ∂y)
        append!(data[:∂γ∂z][2],∂γ∂z)
    end
end

typeForFEM = quote
    type = Element{types[elementType]}
end

# cal_normal = quote
#     if normal
#         nodeTags = gmsh.model.mesh.getElementEdgeNodes(elementType,tag,true)
       
#         if dim == 1 
            
#             for C in 1:ne
#                 𝐿 = 2*determinants[C*ng]
#                 coord, = gmsh.model.mesh.getNode(nodeTags[2*C-1])
#                 x₁ = coord[1]
#                 y₁ = coord[2]
#                 coord, = gmsh.model.mesh.getNode(nodeTags[2*C])
#                 x₂ = coord[1]
#                 y₂ = coord[2]
#                 push!(data[:n₁][2], (y₂-y₁)/𝐿)
#                 push!(data[:n₂][2], (x₁-x₂)/𝐿)
#                 push!(data[:s₁][2], (x₂-x₁)/𝐿)
#                 push!(data[:s₂][2], (y₂-y₁)/𝐿)
#             end
#         end
#         if dim == 2
       
#             nₙ = Int(length(nodeTags)/ne)
           
#             for C in 1:ne
#                 𝐽 = determinants[C*ng]
#                 n₁ = 0.0
#                 n₂ = 0.0
#                 n₃ = 0.0

               

#                 for i in 1:2:nₙ


#                     coord, = gmsh.model.mesh.getNode(nodeTags[nₙ*(C-1)+i])
#                     x₁ = coord[1]
#                     y₁ = coord[2]
#                     z₁ = coord[3]
#                     coord, = gmsh.model.mesh.getNode(nodeTags[nₙ*(C-1)+i+1])
#                     x₂ = coord[1]
#                     y₂ = coord[2]
#                     z₂ = coord[3]

#                     n₁ += y₁*z₂-y₂*z₁
#                     n₂ += z₁*x₂-z₂*x₁
#                     n₃ += x₁*y₂-x₂*y₁
#                 end
                
#                 if elementType == 3
#                     𝐽 *= 8
#                 end
               
#                 # push!(data[:n₁][2], n[1]/𝐽)
#                 # push!(data[:n₂][2], n[2]/𝐽)
#                 # push!(data[:n₃][2], n[3]/𝐽)
            
#                 push!(data[:n₁][2], n₁/𝐽)
#                 push!(data[:n₂][2], n₂/𝐽)
#                 push!(data[:n₃][2], n₃/𝐽)

                
#                 # push!(data[:n₁][2], n₁)
#                 # push!(data[:n₂][2], n₂)
#                 # push!(data[:n₃][2], n₃)
#             end
#         end

        
#     end
# end

cal_normal = quote
    if normal
        nodeTags = gmsh.model.mesh.getElementEdgeNodes(elementType, tag, true)

        # -------------------------
        # dim == 1 : 线单元外法向/切向
        # -------------------------
        if dim == 1
            for C in 1:ne
                𝐿 = 2 * determinants[C*ng]

                coord, = gmsh.model.mesh.getNode(nodeTags[2*C-1])
                x₁ = coord[1]
                y₁ = coord[2]

                coord, = gmsh.model.mesh.getNode(nodeTags[2*C])
                x₂ = coord[1]
                y₂ = coord[2]

                # 单位外法向（平面内旋转90度）
                push!(data[:n₁][2], (y₂ - y₁) / 𝐿)
                push!(data[:n₂][2], (x₁ - x₂) / 𝐿)

                # 单位切向
                push!(data[:s₁][2], (x₂ - x₁) / 𝐿)
                push!(data[:s₂][2], (y₂ - y₁) / 𝐿)
            end
        end

        # -------------------------
        # dim == 2 : 面单元法向（修正版）
        # -------------------------
        if dim == 2
            nₙ = Int(length(nodeTags) / ne)  # 每个面单元的边节点数

            for C in 1:ne
                n₁ = 0.0
                n₂ = 0.0
                n₃ = 0.0

                base = nₙ * (C - 1)

                # 遍历所有边，并闭合：i -> i2（最后一个连回第一个）
                for i in 1:nₙ
                    i2 = (i == nₙ) ? 1 : (i + 1)

                    coord1, = gmsh.model.mesh.getNode(nodeTags[base + i])
                    x₁ = coord1[1]; y₁ = coord1[2]; z₁ = coord1[3]

                    coord2, = gmsh.model.mesh.getNode(nodeTags[base + i2])
                    x₂ = coord2[1]; y₂ = coord2[2]; z₂ = coord2[3]

                    # ∑ (x_i × x_{i+1}) 形成面积法向（方向与节点顺序一致）
                    n₁ += y₁*z₂ - y₂*z₁
                    n₂ += z₁*x₂ - z₂*x₁
                    n₃ += x₁*y₂ - x₂*y₁
                end

                # 归一化为单位法向
                nrm = sqrt(n₁*n₁ + n₂*n₂ + n₃*n₃)
                if nrm == 0.0
                    # 退化面：给一个零向量，避免 NaN
                    push!(data[:n₁][2], 0.0)
                    push!(data[:n₂][2], 0.0)
                    push!(data[:n₃][2], 0.0)
                else
                    push!(data[:n₁][2], n₁ / nrm)
                    push!(data[:n₂][2], n₂ / nrm)
                    push!(data[:n₃][2], n₃ / nrm)
                end
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
    println("Info: Generate $ne elements of $type with $ng integration points.")
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

function getRKGSIElements(nodes::Vector{N},dimTag::Pair{Int,Vector{Int}},type::DataType,integrationOrder::Int,sp::SpatialPartition;normal::Bool=false) where N<:Node
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
    for (elementType,elementTypeΩ,nodeTag,tag,tagΩ_) in zip(elementTypes,elementTypesΩ,nodeTags,tags,tagΩ)
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


function getPiecewiseBoundaryElements(dimTag::Pair{Int,Vector{Int}},dimTagΩ::Pair{Int,Vector{Int}},type::DataType,integrationOrder::Int,Dim::Int)
    normal = false
    $prequote
    $preForEdge
    for (elementType,elementTypeΩ,nodeTag,tag,tagΩ_) in zip(elementTypes,elementTypesΩ,nodeTags,tags,tagΩ)
        ## integration rule
        $integrationByGmsh
        ## coordinates
        $coordinatesForFaces
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

function getElements(
    nodes::Vector{N},
    dimTag::Pair{Int,Vector{Int}},
    ::Type{KGS},                      # KGradientSmoothing{...,:Tri3}
    scheme_smooth::Dict{Symbol,<:Any},# e.g. quadraturerule(:TriGI3)
    scheme_main::Dict{Symbol,<:Any},  # e.g. quadraturerule(:TriRK6) 或 Gauss高阶
    sp::SpatialPartition;
    normal::Bool=false
) where {N<:Node, KGS<:AbstractElement}

    # 1) gmsh 读 element connectivity（跟你原 prequote 同源）
    dim, tags = dimTag
    elementTypes = Int32[]
    nodeTagsByEntity = Vector{UInt64}[]
    for tag in tags
        elementTypes_, ~, nodeTags_ = gmsh.model.mesh.getElements(dim, tag)
        push!(elementTypes, elementTypes_[1])
        push!(nodeTagsByEntity, nodeTags_[1])
    end

    elements = AbstractElement[]

    # 2) 两套积分规则
    localCoord̃ = _localCoord3(scheme_smooth)
    weights̃    = scheme_smooth[:w]
    localCoord  = _localCoord3(scheme_main)
    weights     = scheme_main[:w]

    ng̃ = length(weights̃)
    ng  = length(weights)

    # 3) 全局计数（与旧版 import_rkgsi 对齐的思想）
    c  = 0
    g̃ = 0; s̃ = 0
    g  = 0; s  = 0
    Cglobal = 0

    for (elementType, nodeTag, tag) in zip(elementTypes, nodeTagsByEntity, tags)
        ~, ~, ~, ni = gmsh.model.mesh.getElementProperties(elementType)
        ne = Int(length(nodeTag) / ni)

        # 4) 计算两套积分点在物理空间的坐标（gmsh给出）
        jac_m, det_m, coord_m = gmsh.model.mesh.getJacobians(elementType, localCoord,  tag)
        x_m = coord_m[1:3:end]; y_m = coord_m[2:3:end]; z_m = coord_m[3:3:end]

        jac_s, det_s, coord_s = gmsh.model.mesh.getJacobians(elementType, localCoord̃, tag)
        x_s = coord_s[1:3:end]; y_s = coord_s[2:3:end]; z_s = coord_s[3:3:end]

        # 5) 两套 data（全元素共享：节点 Node 里引用同一份数组）
        data_main = Dict{Symbol,Tuple{Int,Vector{Float64}}}(
            :w => (1, collect(weights)),
            :ξ => (1, collect(scheme_main[:ξ])),
            :η => (1, haskey(scheme_main, :η) ? collect(scheme_main[:η]) : Float64[]),
            :x => (2, Float64[]), :y => (2, Float64[]), :z => (2, Float64[]),
            :𝑤 => (2, Float64[]),
        )

        data_smooth = Dict{Symbol,Tuple{Int,Vector{Float64}}}(
            :w => (1, collect(weights̃)),
            :ξ => (1, collect(scheme_smooth[:ξ])),
            :η => (1, haskey(scheme_smooth, :η) ? collect(scheme_smooth[:η]) : Float64[]),
            :x => (2, Float64[]), :y => (2, Float64[]), :z => (2, Float64[]),
            :𝑤 => (2, Float64[]),
        )

        # 预填坐标与加权雅可比（按 element 展平 append）
        append!(data_main[:x][2], x_m); append!(data_main[:y][2], y_m); append!(data_main[:z][2], z_m)
        for i in 1:length(det_m)
            wloc = weights[(i - 1) % ng + 1]
            push!(data_main[:𝑤][2], det_m[i] * wloc)
        end

        append!(data_smooth[:x][2], x_s); append!(data_smooth[:y][2], y_s); append!(data_smooth[:z][2], z_s)
        for i in 1:length(det_s)
            wloc = weights̃[(i - 1) % ng̃ + 1]
            push!(data_smooth[:𝑤][2], det_s[i] * wloc)
        end

        # 6) 逐单元生成 KGradientSmoothing：((c,nc,𝓒),(g̃,ng̃,𝓖̃),(g,ng,𝓖))
        for Ce in 1:ne
            Cglobal += 1

            # patch 支撑域：用主积分点坐标 union sp
            idxset = Set{Int}()
            for ig in 1:ng
                Gidx = ng * (Ce - 1) + ig
                union!(idxset, sp(x_m[Gidx], y_m[Gidx], z_m[Gidx]))
            end
            nc = length(idxset)
            𝓒 = [nodes[i] for i in idxset]

            # 记录该单元的起始全局积分点编号（更稳健）
            g̃0 = g̃
            g0  = g

            # smooth 积分点集合 𝓖̃
            𝓖̃ = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
            for ig in 1:ng̃
                g̃ += 1
                push!(𝓖̃,
                    Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}(
                        (𝑔=ig, 𝐺=g̃, 𝐶=Cglobal, 𝑠=s̃),   # ← 关键修复：NamedTuple
                        data_smooth
                    )
                )
                s̃ += nc
            end

            # main 积分点集合 𝓖
            𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
            for ig in 1:ng
                g += 1
                push!(𝓖,
                    Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}(
                        (𝑔=ig, 𝐺=g, 𝐶=Cglobal, 𝑠=s),    # ← 关键修复：NamedTuple
                        data_main
                    )
                )
                s += nc
            end

            push!(elements, KGS((c, nc, 𝓒), (g̃0, ng̃, 𝓖̃), (g0, ng, 𝓖)))
            c += nc
        end

        println("Info: Generate $ne elements of $KGS with $ng(main) and $ng̃(smooth) integration points.")
    end

    return elements
end




function _localCoord3(scheme::Dict{Symbol,<:Any})
    ξ = scheme[:ξ]
    η = haskey(scheme, :η) ? scheme[:η] : zeros(length(ξ))
    n = length(ξ)
    lc = Vector{Float64}(undef, 3n)
    @inbounds for i in 1:n
        lc[3i-2] = ξ[i]
        lc[3i-1] = η[i]
        lc[3i]   = 0.0
    end
    return lc
end