
"""
getnβ
"""
getnβ(ap::T) where T<:AbstractElement = length(getfield(ap.π[1],:data)[:x][2])
@inline getnβ(aps::Vector{T}) where T<:AbstractElement = getnβ(aps[1])

"""
# Element
"""
struct Element{T}<:AbstractElement{T}
    π::Vector{Node}
    π::Vector{SNode}
end
Element{T}(π::Vector{Node}) where T = Element{T}(π,SNode[])
Element{T}(a::S) where {T,S<:AbstractElement} = Element{T}(a.π)

# function Element{T}(as::Vector{S};renumbering::Bool=false) where {T,S<:AbstractElement}
#     aps = Element{T}[]
#     if renumbering
#         index, data = renumber(as)
#         for a in as
#             π = [Node(index[x.id],data) for x in a.π]
#             π = Node[]
#             push!(aps,Element{T}(π,π))
#         end
#     else
#         for a in as
#             push!(aps,Element{T}(a))
#         end
#     end
#     return aps
# end

# function Element{T}(a::AbstractElement,b::AbstractElement) where T
#     π = a.π
#     π = getπ(a,b)
#     π β  nothing ? Element{T}(π,π) : nothing
# end

# function Element{T}(as::Vector{A},bs::Vector{B}) where {T,A<:AbstractElement,B<:AbstractElement}
#     aps = Element{T}[]
#     for a in as
#         for b in bs
#             ap = Element{T}(a,b)
#             ap β  nothing ? push!(aps,ap) : nothing
#         end
#     end
#     return aps
# end

# function renumber(aps::Vector{T}) where T<:AbstractElement
#     index = Dict{Int,Int}()
#     n = 0
#     for ap in aps
#         for x in ap.π
#             I = x.id
#             if ~haskey(index,I)
#                 n += 1
#                 index[I] = n
#             end
#         end
#     end
#     data_ = aps[1].π[1].data
#     data = Dict(:x=>zeros(n),:y=>zeros(n),:z=>zeros(n))
#     for (j,i) in index
#         data[:x][i] = data_[:x][j]
#         data[:y][i] = data_[:y][j]
#         data[:z][i] = data_[:z][j]
#     end
#     return index, data
# end

"""
getπ(ap::T,x::SNode) where T<:AbstractElement
getπ(ap::T,ΞΎ::Float64...) where T<:AbstractElement
"""
@inline getπ(ap::T,::Any) where T<:AbstractElement{:Poi1} = (ap.π[1].x,ap.π[1].y,ap.π[1].z)
@inline getπ(ap::T,ΞΎ::SNode) where T<:AbstractElement{:Seg2} = getπ(ap,ΞΎ.ΞΎ)
@inline getπ(ap::T,ΞΎ::SNode) where T<:AbstractElement{:Tri3} = getπ(ap,ΞΎ.ΞΎ,ΞΎ.Ξ·)
@inline getπ(ap::T,ΞΎ::SNode) where T<:AbstractElement{:Quad} = getπ(ap,ΞΎ.ΞΎ,ΞΎ.Ξ·)

function getπ(ap::T,ΞΎ::Float64) where T<:AbstractElement{:Seg2}
    xβ = ap.π[1].x
    yβ = ap.π[1].y
    zβ = ap.π[1].z
    xβ = ap.π[2].x
    yβ = ap.π[2].y
    zβ = ap.π[2].z
    Nβ = 0.5*(1-ΞΎ)
    Nβ = 0.5*(1+ΞΎ)
    return (xβ*Nβ+xβ*Nβ,yβ*Nβ+yβ*Nβ,zβ*Nβ+zβ*Nβ)
end
function getπ(ap::T,ΞΎ::Float64,Ξ·::Float64) where T<:AbstractElement{:Tri3}
    xβ = ap.π[1].x
    yβ = ap.π[1].y
    zβ = ap.π[1].z
    xβ = ap.π[2].x
    yβ = ap.π[2].y
    zβ = ap.π[2].z
    xβ = ap.π[3].x
    yβ = ap.π[3].y
    zβ = ap.π[3].z
    Nβ = ΞΎ
    Nβ = Ξ·
    Nβ = 1.0-ΞΎ-Ξ·
    return (xβ*Nβ+xβ*Nβ+xβ*Nβ,yβ*Nβ+yβ*Nβ+yβ*Nβ,zβ*Nβ+zβ*Nβ+zβ*Nβ)
end

function getπ(ap::T,ΞΎ::Float64,Ξ·::Float64) where T<:AbstractElement{:Quad}
    xβ = ap.π[1].x
    yβ = ap.π[1].y
    zβ = ap.π[1].z
    xβ = ap.π[2].x
    yβ = ap.π[2].y
    zβ = ap.π[2].z
    xβ = ap.π[3].x
    yβ = ap.π[3].y
    zβ = ap.π[3].z
    xβ = ap.π[4].x
    yβ = ap.π[4].y
    zβ = ap.π[4].z
    Nβ,Nβ,Nβ,Nβ = getπ­(ap,ΞΎ,Ξ·)
    return (xβ*Nβ+xβ*Nβ+xβ*Nβ+xβ*Nβ,yβ*Nβ+yβ*Nβ+yβ*Nβ+yβ*Nβ,zβ*Nβ+zβ*Nβ+zβ*Nβ+zβ*Nβ)
end

function getπ±(ap::T,ΞΎ::SNode) where T<:AbstractElement{:Quad}
    xβ = ap.π[1].x
    xβ = ap.π[2].x
    xβ = ap.π[3].x
    xβ = ap.π[4].x
    yβ = ap.π[1].y
    yβ = ap.π[2].y
    yβ = ap.π[3].y
    yβ = ap.π[4].y
    βNββΞΎ,βNββΞΎ,βNββΞΎ,βNββΞΎ = getβπ­βΞΎ(ap,ΞΎ)
    βNββΞ·,βNββΞ·,βNββΞ·,βNββΞ· = getβπ­βΞ·(ap,ΞΎ)
    Jββ = βNββΞΎ*xβ + βNββΞΎ*xβ + βNββΞΎ*xβ + βNββΞΎ*xβ
    Jββ = βNββΞ·*xβ + βNββΞ·*xβ + βNββΞ·*xβ + βNββΞ·*xβ
    Jββ = βNββΞΎ*yβ + βNββΞΎ*yβ + βNββΞΎ*yβ + βNββΞΎ*yβ
    Jββ = βNββΞ·*yβ + βNββΞ·*yβ + βNββΞ·*yβ + βNββΞ·*yβ
    return Jββ,Jββ,Jββ,Jββ
end

"""
getπ½(ap::T,x::SNode) where T<:AbstractElement
"""
@inline getπ½(  ::T,::Any) where T<:AbstractElement{:Poi1} = 1.0
@inline getπ½(ap::T,::Any) where T<:AbstractElement{:Seg2} = 0.5*getπΏ(ap)
@inline getπ½(ap::T,::Any) where T<:AbstractElement{:Tri3} = 2.0*getπ΄(ap)
@inline function getπ½(ap::T,ΞΎ::SNode) where T<:AbstractElement{:Quad}
    Jββ,Jββ,Jββ,Jββ = getπ±(ap,ΞΎ)
    return Jββ*Jββ-Jββ*Jββ
end

"""
getπ€(ap::T,x::SNode) where T<:AbstractElement
"""
@inline getπ€(  ::T,::Any) where T<:AbstractElement{:Poi1} = 1.0
@inline getπ€(ap::T,ΞΎ::SNode) where T<:AbstractElement{:Seg2} = 0.5*getπΏ(ap)*ΞΎ.w
@inline getπ€(ap::T,ΞΎ::SNode) where T<:AbstractElement{:Tri3} = getπ΄(ap)*ΞΎ.w
@inline getπ€(ap::T,ΞΎ::SNode) where T<:AbstractElement{:Quad} = getπ½(ap,ΞΎ)*ΞΎ.w

"""
getπΏ,getπ΄,getπ
"""
@inline function getπΏ(ap::T) where T<:AbstractElement{:Seg2}
    xβ = ap.π[1].x
    yβ = ap.π[1].y
    zβ = ap.π[1].z
    xβ = ap.π[2].x
    yβ = ap.π[2].y
    zβ = ap.π[2].z
    return ((xβ-xβ)^2+(yβ-yβ)^2+(zβ-zβ)^2)^0.5
end
function getπ΄(ap::T) where T<:AbstractElement{:Tri3}
    xβ = ap.π[1].x
    yβ = ap.π[1].y
    zβ = ap.π[1].z
    xβ = ap.π[2].x
    yβ = ap.π[2].y
    zβ = ap.π[2].z
    xβ = ap.π[3].x
    yβ = ap.π[3].y
    zβ = ap.π[3].z
    π΄β = 0.5*(yβ*zβ+yβ*zβ+yβ*zβ-yβ*zβ-yβ*zβ-yβ*zβ)
    π΄β = 0.5*(zβ*xβ+zβ*xβ+zβ*xβ-zβ*xβ-zβ*xβ-zβ*xβ)
    π΄β = 0.5*(xβ*yβ+xβ*yβ+xβ*yβ-xβ*yβ-xβ*yβ-xβ*yβ)
    return (π΄β^2 + π΄β^2 + π΄β^2)^0.5
end

"""
setπ!
"""
function setπ!(aps::Vector{T}) where T<:AbstractElement{:Seg2}
    data = getfield(aps[1].π[1],:data)
    n = length(data[:x][2])
    push!(data,:nβ=>(2,zeros(n)))
    push!(data,:nβ=>(2,zeros(n)))
    for ap in aps
        setπ!(ap)
    end
end

function setπ!(ap::T) where T<:AbstractElement{:Seg2}
    xβ = ap.π[1].x
    yβ = ap.π[1].y
    xβ = ap.π[2].x
    yβ = ap.π[2].y
    πΏ = getπΏ(ap)
    nβ = (yβ-yβ)/πΏ
    nβ = (xβ-xβ)/πΏ
    for ΞΎ in ap.π
        ΞΎ.nβ = nβ
        ΞΎ.nβ = nβ
    end
end

# function getπ(ap::T,ΞΎ::SNode) where T<:AbstractElement


# ## setπ!
# function setπ!(aps::Vector{T}) where T<:AbstractElement
#     for ap in aps
#         setπ!(ap)
#     end
# end

# function setπ!(ap::T) where T<:AbstractElement{:Seg2}
#     π = ap.π
#     for ΞΎ in π
#         ΞΎ.nβ = getπ(ap,ΞΎ)
#     end
# end

# function setπ!(ap::T) where T<:AbstractElement{:Tri3}
#     π = ap.π
#     for ΞΎ in π
#         ΞΎ.nβ, ΞΎ.nβ = getπ(ap,ΞΎ)
#     end
# end

## shape functions
"""
setπ­!
"""
function setπ­!(aps::Vector{T}) where T<:AbstractElement
    for ap in aps
        setπ­!(ap)
    end
end

function setπ­!(ap::Element{S}) where S
    π = ap.π
    for ΞΎ in π
        N = getπ­(ap,ΞΎ)
        for i in 1:length(ap.π)
            π­ = ΞΎ[:π­]
            π­[i] = N[i]
        end
    end
end

"""
getπ­(ap::Element,ΞΎ::SNode)
"""
# ------------- Poi1 ---------------
@inline getπ­(::Element{:Poi1},::Any) = 1.0

# ------------- Seg2 ---------------
@inline getπ­(ap::Element{:Seg2},ΞΎ::SNode) = (0.5*(1.0-ΞΎ.ΞΎ),0.5*(1.0+ΞΎ.ΞΎ))
@inline function getβπ­βx(ap::Element{:Seg2},::Any)
    πΏ = getπΏ(ap)
    return (-1.0/πΏ,1.0/πΏ)
end

# ------------- Tri3 ---------------
@inline getπ­(ap::Element{:Tri3},ΞΎ::SNode) = (ΞΎ.ΞΎ,ΞΎ.Ξ·,1.0-ΞΎ.ΞΎ-ΞΎ.Ξ·)
@inline function getβπ­βx(ap::Element{:Tri3},ΞΎ::SNode)
    yβ = ap.π[1].y
    yβ = ap.π[2].y
    yβ = ap.π[3].y
    π΄ = getπ΄(ap)
    return (yβ-yβ)/2.0/π΄,(yβ-yβ)/2.0/π΄,(yβ-yβ)/2.0/π΄
end
@inline function getβπ­βy(ap::Element{:Tri3},ΞΎ::SNode)
    xβ = ap.π[1].x
    xβ = ap.π[2].x
    xβ = ap.π[3].x
    π΄ = getπ΄(ap)
    return (xβ-xβ)/2.0/π΄,(xβ-xβ)/2.0/π΄,(xβ-xβ)/2.0/π΄
end

# ------------- Quad ---------------
@inline getπ­(ap::Element{:Quad},ΞΎ::SNode) = getπ­(ap,ΞΎ.ΞΎ,ΞΎ.Ξ·)
@inline getβπ­βΞΎ(ap::Element{:Quad},ΞΎ::SNode) = getβπ­βΞΎ(ap,ΞΎ.Ξ·)
@inline getβπ­βΞ·(ap::Element{:Quad},ΞΎ::SNode) = getβπ­βΞ·(ap,ΞΎ.ΞΎ)

function getπ­(ap::Element{:Quad},ΞΎ::Float64,Ξ·::Float64)
    Nβ = 0.25*(1.0-ΞΎ)*(1.0-Ξ·)
    Nβ = 0.25*(1.0+ΞΎ)*(1.0-Ξ·)
    Nβ = 0.25*(1.0+ΞΎ)*(1.0+Ξ·)
    Nβ = 0.25*(1.0-ΞΎ)*(1.0+Ξ·)
    return (Nβ,Nβ,Nβ,Nβ)
end
function getβπ­βΞΎ(ap::Element{:Quad},Ξ·::Float64)
    βNββΞΎ = - 0.25*(1-Ξ·)
    βNββΞΎ =   0.25*(1-Ξ·)
    βNββΞΎ =   0.25*(1+Ξ·)
    βNββΞΎ = - 0.25*(1+Ξ·)
    return (βNββΞΎ,βNββΞΎ,βNββΞΎ,βNββΞΎ)
end
function getβπ­βΞ·(ap::Element{:Quad},ΞΎ::Float64)
    βNββΞ· = - 0.25*(1-ΞΎ)
    βNββΞ· = - 0.25*(1+ΞΎ)
    βNββΞ· =   0.25*(1+ΞΎ)
    βNββΞ· =   0.25*(1-ΞΎ)
    return (βNββΞ·,βNββΞ·,βNββΞ·,βNββΞ·)
end
function getβπ­βxβπ­βy(ap::Element{:Quad},ΞΎ::SNode)
    xβ = ap.π[1].x
    xβ = ap.π[2].x
    xβ = ap.π[3].x
    xβ = ap.π[4].x
    yβ = ap.π[1].y
    yβ = ap.π[2].y
    yβ = ap.π[3].y
    yβ = ap.π[4].y
    βNββΞΎ,βNββΞΎ,βNββΞΎ,βNββΞΎ = getβπ­βΞΎ(ap,ΞΎ)
    βNββΞ·,βNββΞ·,βNββΞ·,βNββΞ· = getβπ­βΞ·(ap,ΞΎ)
    βxβΞΎ = βNββΞΎ*xβ + βNββΞΎ*xβ + βNββΞΎ*xβ + βNββΞΎ*xβ
    βxβΞ· = βNββΞ·*xβ + βNββΞ·*xβ + βNββΞ·*xβ + βNββΞ·*xβ
    βyβΞΎ = βNββΞΎ*yβ + βNββΞΎ*yβ + βNββΞΎ*yβ + βNββΞΎ*yβ
    βyβΞ· = βNββΞ·*yβ + βNββΞ·*yβ + βNββΞ·*yβ + βNββΞ·*yβ
    detJ = βxβΞΎ*βyβΞ· - βxβΞ·*βyβΞΎ
    βΞΎβx =   βyβΞ·/detJ
    βΞ·βx = - βyβΞΎ/detJ
    βΞΎβy = - βxβΞ·/detJ
    βΞ·βy =   βxβΞΎ/detJ
    βNββx = βNββΞΎ*βΞΎβx + βNββΞ·*βΞ·βx
    βNββx = βNββΞΎ*βΞΎβx + βNββΞ·*βΞ·βx
    βNββx = βNββΞΎ*βΞΎβx + βNββΞ·*βΞ·βx
    βNββx = βNββΞΎ*βΞΎβx + βNββΞ·*βΞ·βx
    βNββy = βNββΞΎ*βΞΎβy + βNββΞ·*βΞ·βy
    βNββy = βNββΞΎ*βΞΎβy + βNββΞ·*βΞ·βy
    βNββy = βNββΞΎ*βΞΎβy + βNββΞ·*βΞ·βy
    βNββy = βNββΞΎ*βΞΎβy + βNββΞ·*βΞ·βy
    return (βNββx,βNββx,βNββx,βNββx),(βNββy,βNββy,βNββy,βNββy)
end
# @inline getβπ­(ap::Element{:Quad},ΞΎ::π) where π<:AbstractNode = getπ­(ap,ΞΎ),getβπ­βxβπ­βy(ap,ΞΎ)...,(0.0,0.0,0.0,0.0)

"""
β,β©
"""
function issubset(a::T,b::S) where {T<:AbstractElement{:Poi1},S<:AbstractElement{:Seg2}}
    i = findfirst(x->x==a.π[1],b.π)
    return i β  nothing && i β€ 2
end

@inline intersect(a::T,b::T) where T<:AbstractElement = a.π == b.π ? a : nothing
@inline function intersect(a::T,b::S) where {T<:AbstractElement{:Seg2},S<:AbstractElement{:Poi1}}
    i = findfirst(x->x==b.π[1],a.π)
    return i β  nothing && i β€ 2 ? a : nothing
end
@inline function intersect(a::T,b::S) where {T<:AbstractElement{:Tri3},S<:AbstractElement{:Poi1}}
    i = findfirst(x->x==b.π[1],a.π)
    return i β  nothing && i β€ 3 ? a : nothing
end
@inline function intersect(a::T,b::S) where {T<:AbstractElement{:Tri3},S<:AbstractElement{:Seg2}}
    i = findfirst(x->x==b.π[1],a.π)
    j = findfirst(x->x==b.π[2],a.π)
    return i β  nothing && j β  nothing && i β€ 3 && j β€ 3 ? a : nothing
end
function intersect(as::Vector{T},bs::Vector{S}) where {T<:AbstractElement,S<:AbstractElement}
    aps = T[]
    for b in bs
        for a in as
            ap = aβ©b
            ap β  nothing ? push!(aps,ap) : nothing
        end
    end
    return aps
end