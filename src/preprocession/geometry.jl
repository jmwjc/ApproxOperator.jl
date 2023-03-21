
abstract type AbstractGeometry end
abstract type AbstractSegment<:AbstractGeometry end
abstract type AbstractSurface<:AbstractGeometry end
abstract type AbstractVolume<:AbstractGeometry end
abstract type AbstractTriangle<:AbstractSurface end
abstract type AbstracQuadrilateralt<:AbstractSurface end
abstract type AbstractTetrahedron<:AbstractVolume end

struct Point
    i::Int
    x::Float64
    y::Float64
    z::Float64
end

struct Seg2<:AbstractSegment
    vertices::NTuple{2,Point}
end

struct Seg3<:AbstractSegment
    vertices::NTuple{3,Point}
end

struct Tri3<:AbstractTriangle
    vertices::NTuple{3,Point}
    edges::NTuple{3,Seg2}
end

struct Tri6<:AbstractTriangle
    vertices::NTuple{6,Point}
    edges::NTuple{3,Seg3}
end

struct Quad4<:AbstracQuadrilateralt
    vertices::NTuple{4,Point}
    edges::NTuple{4,Seg2}
end

struct Quad9<:AbstracQuadrilateralt
    vertices::NTuple{9,Point}
    edges::NTuple{4,Seg3}
end

struct Tet4<:AbstractTetrahedron
    vertices::NTuple{4,Point}
    edges::NTuple{6,Seg2}
    surfaces::NTuple{4,Tri3}
end

function (a::Seg2)(ξ::Float64)
    v₁ = a.vertices[1]
    v₂ = a.vertices[2]
    x₁ = v₁.x
    y₁ = v₁.y
    z₁ = v₁.z
    x₂ = v₂.x
    y₂ = v₂.y
    z₂ = v₂.z
    N₁ = 0.5*(1.0-ξ)
    N₂ = 0.5*(1.0+ξ)
    return N₁*x₁+N₂*x₂,
           N₁*y₁+N₂*y₂,
           N₁*z₁+N₂*z₂
end

function (a::Seg2)(v::Point)
    i = findfirst(x->x==v,a.vertices)
    if i == 1
        return -1.0
    elseif i == 2
        return 1.0
    else
        return nothing
    end
end

function (a::Seg3)(ξ::Float64)
    v₁ = a.vertices[1]
    v₂ = a.vertices[2]
    v₃ = a.vertices[3]
    x₁ = v₁.x
    y₁ = v₁.y
    z₁ = v₁.z
    x₂ = v₂.x
    y₂ = v₂.y
    z₂ = v₂.z
    x₃ = v₃.x
    y₃ = v₃.y
    z₃ = v₃.z
    N₁ = 0.5*ξ*(ξ-1.0)
    N₂ = 1.0-ξ^2
    N₃ = 0.5*ξ*(ξ+1.0)
    return N₁*x₁+N₂*x₂+N₃*x₃,
           N₁*y₁+N₂*y₂+N₃*y₃,
           N₁*z₁+N₂*z₂+N₃*z₃
end

function (a::Seg3)(v::Point)
    i = findfirst(x->x==v,a.vertices)
    if i == 1
        return -1.0
    elseif i == 3
        return 1.0
    else
        return nothing
    end
end

function (a::Tri3)(ξ::Float64,η::Float64)
    v₁ = a.vertices[1]
    v₂ = a.vertices[2]
    v₃ = a.vertices[3]
    x₁ = v₁.x
    y₁ = v₁.y
    z₁ = v₁.z
    x₂ = v₂.x
    y₂ = v₂.y
    z₂ = v₂.z
    x₃ = v₃.x
    y₃ = v₃.y
    z₃ = v₃.z
    N₁ = ξ
    N₂ = η
    N₃ = 1.0-ξ-η
    return N₁*x₁+N₂*x₂+N₃*x₃,
           N₁*y₁+N₂*y₂+N₃*y₃,
           N₁*z₁+N₂*z₂+N₃*z₃
end

function (a::Tri3)(v::Point)
    i = findfirst(x->x==v,a.vertices)
    if i == 1
        return 1.0,0.0
    elseif i == 2
        return 0.0,1.0
    elseif i == 3
        return 0.0,0.0
    else
        return nothing
    end
end

function (a::Tri3)(e::Seg2,ξ::Float64)
    i = findfirst(x->x==e,a.edges)
    if i == 1
        return 0.0,0.5*(1.0-ξ)
    elseif i == 2
        return 0.5*(1.0+ξ),0.0
    elseif i == 3
        return 0.5*(1.0-ξ),0.5*(1.0+ξ)
    else
        return nothing
    end
end

function (a::Tri6)(ξ::Float64,η::Float64)
    γ = 1.0-ξ-η
    x₁ = a.vertices[1].x;y₁ = a.vertices[1].y;z₁ = a.vertices[1].z
    x₂ = a.vertices[2].x;y₂ = a.vertices[2].y;z₂ = a.vertices[2].z
    x₃ = a.vertices[3].x;y₃ = a.vertices[3].y;z₃ = a.vertices[3].z
    x₄ = a.vertices[4].x;y₄ = a.vertices[4].y;z₄ = a.vertices[4].z
    x₅ = a.vertices[5].x;y₅ = a.vertices[5].y;z₅ = a.vertices[5].z
    x₆ = a.vertices[6].x;y₆ = a.vertices[6].y;z₆ = a.vertices[6].z
    N₁ = ξ*(2*ξ-1)
    N₂ = η*(2*η-1)
    N₃ = γ*(2*γ-1)
    N₄ = 4*ξ*η
    N₅ = 4*η*γ
    N₆ = 4*γ*ξ
    return x₁*N₁+x₂*N₂+x₃*N₃+x₄*N₄+x₅*N₅+x₆*N₆,
           y₁*N₁+y₂*N₂+y₃*N₃+y₄*N₄+y₅*N₅+y₆*N₆,
           z₁*N₁+z₂*N₂+z₃*N₃+z₄*N₄+z₅*N₅+z₆*N₆
end

function get𝐴(a::Tri3)
    x₁ = a.vertices[1].x
    x₂ = a.vertices[2].x
    x₃ = a.vertices[3].x
    y₁ = a.vertices[1].y
    y₂ = a.vertices[2].y
    y₃ = a.vertices[3].y
    z₁ = a.vertices[1].z
    z₂ = a.vertices[2].z
    z₃ = a.vertices[3].z

    return 0.5*(x₁*y₂+x₂*y₃+x₃*y₁-x₂*y₁-x₃*y₂-x₁*y₃)
end