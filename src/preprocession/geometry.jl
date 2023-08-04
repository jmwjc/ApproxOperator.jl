abstract type AbstractGeometry end
abstract type AbstractPoint<:AbstractGeometry end
abstract type AbstractSegment<:AbstractGeometry end
abstract type AbstractSurface<:AbstractGeometry end
abstract type AbstractVolume<:AbstractGeometry end
abstract type AbstractTriangle<:AbstractSurface end
abstract type AbstracQuadrilateralt<:AbstractSurface end
abstract type AbstractTetrahedron<:AbstractVolume end

struct Point
    i::Int
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
end

struct Poi1<:AbstractPoint
    i::Int
    vertex::Point
end

struct Seg2<:AbstractSegment
    i::Int
    vertex::NTuple{2,Point}
end

struct Seg3<:AbstractSegment
    i::Int
    vertex::NTuple{3,Point}
end

struct Tri3<:AbstractTriangle
    i::Int
    vertex::NTuple{3,Point}
end

struct Tri6<:AbstractTriangle
    i::Int
    vertex::NTuple{6,Point}
end

struct Quad4<:AbstracQuadrilateralt
    i::Int
    vertex::NTuple{4,Point}
end

struct Quad9<:AbstracQuadrilateralt
    i::Int
    vertex::NTuple{9,Point}
end

struct Tet4<:AbstractTetrahedron
    i::Int
    vertex::NTuple{4,Point}
end

function Tri3(vertex::NTuple{3,Point})
end

function Quad4(vertex::NTuple{4,Point})
end


function (a::Poi1)(::Any)
    v₁ = a.vertex[1]
    x = v₁.x
    y = v₁.y
    z = v₁.z
    return x,y,z
end

function (a::Seg2)(ξ::Float64)
    v₁ = a.vertex[1]
    v₂ = a.vertex[2]
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
    i = findfirst(x->x==v,a.vertex)
    if i == 1
        return -1.0
    elseif i == 2
        return 1.0
    else
        return nothing
    end
end

function (a::Seg3)(ξ::Float64)
    v₁ = a.vertex[1]
    v₂ = a.vertex[2]
    v₃ = a.vertex[3]
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
    i = findfirst(x->x==v,a.vertex)
    if i == 1
        return -1.0
    elseif i == 3
        return 1.0
    else
        return nothing
    end
end

function (a::Tri3)(ξ::Float64,η::Float64)
    v₁ = a.vertex[1]
    v₂ = a.vertex[2]
    v₃ = a.vertex[3]
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
    i = findfirst(x->x==v,a.vertex)
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
    x₁ = a.vertex[1].x;y₁ = a.vertex[1].y;z₁ = a.vertex[1].z
    x₂ = a.vertex[2].x;y₂ = a.vertex[2].y;z₂ = a.vertex[2].z
    x₃ = a.vertex[3].x;y₃ = a.vertex[3].y;z₃ = a.vertex[3].z
    x₄ = a.vertex[4].x;y₄ = a.vertex[4].y;z₄ = a.vertex[4].z
    x₅ = a.vertex[5].x;y₅ = a.vertex[5].y;z₅ = a.vertex[5].z
    x₆ = a.vertex[6].x;y₆ = a.vertex[6].y;z₆ = a.vertex[6].z
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

function (a::Quad4)(ξ::Float64,η::Float64)
    x₁ = a.vertex[1].x
    y₁ = a.vertex[1].y
    z₁ = a.vertex[1].z
    x₂ = a.vertex[2].x
    y₂ = a.vertex[2].y
    z₂ = a.vertex[2].z
    x₃ = a.vertex[3].x
    y₃ = a.vertex[3].y
    z₃ = a.vertex[3].z
    x₄ = a.vertex[4].x
    y₄ = a.vertex[4].y
    z₄ = a.vertex[4].z
    N₁,N₂,N₃,N₄ = get𝝭(a,ξ,η)
    return (x₁*N₁+x₂*N₂+x₃*N₃+x₄*N₄,y₁*N₁+y₂*N₂+y₃*N₃+y₄*N₄,z₁*N₁+z₂*N₂+z₃*N₃+z₄*N₄)
end
function get𝐴(a::Tri3)
    x₁ = a.vertex[1].x
    x₂ = a.vertex[2].x
    x₃ = a.vertex[3].x
    y₁ = a.vertex[1].y
    y₂ = a.vertex[2].y
    y₃ = a.vertex[3].y
    z₁ = a.vertex[1].z
    z₂ = a.vertex[2].z
    z₃ = a.vertex[3].z

    return 0.5*(x₁*y₂+x₂*y₃+x₃*y₁-x₂*y₁-x₃*y₂-x₁*y₃)
end

function get𝐿(a::Seg2)
    x₁ = a.vertex[1].x
    x₂ = a.vertex[2].x
    y₁ = a.vertex[1].y
    y₂ = a.vertex[2].y

    return ((x₂-x₁)^2+(y₂-y₁)^2)^0.5
end


function get𝑱(a::Quad4,ξ::Float64,η::Float64)
    x₁ = a.vertex[1].x
    y₁ = a.vertex[1].y
    x₂ = a.vertex[2].x
    y₂ = a.vertex[2].y
    x₃ = a.vertex[3].x
    y₃ = a.vertex[3].y
    x₄ = a.vertex[4].x
    y₄ = a.vertex[4].y
    ∂N₁∂ξ,∂N₂∂ξ,∂N₃∂ξ,∂N₄∂ξ = get∂𝝭∂ξ(a,ξ)
    ∂N₁∂η,∂N₂∂η,∂N₃∂η,∂N₄∂η = get∂𝝭∂η(a,η)
    J₁₁ = ∂N₁∂ξ*x₁ + ∂N₂∂ξ*x₂ + ∂N₃∂ξ*x₃ + ∂N₄∂ξ*x₄
    J₁₂ = ∂N₁∂η*x₁ + ∂N₂∂η*x₂ + ∂N₃∂η*x₃ + ∂N₄∂η*x₄
    J₂₁ = ∂N₁∂ξ*y₁ + ∂N₂∂ξ*y₂ + ∂N₃∂ξ*y₃ + ∂N₄∂ξ*y₄
    J₂₂ = ∂N₁∂η*y₁ + ∂N₂∂η*y₂ + ∂N₃∂η*y₃ + ∂N₄∂η*y₄
    return J₁₁,J₂₁,J₁₂,J₂₂
end

function get𝐽(a::Quad4,ξ::Float64,η::Float64)
    J₁₁,J₂₁,J₁₂,J₂₂ = get𝑱(a,ξ,η)
    return J₁₁*J₂₂-J₂₁*J₁₂
end

function get∂𝝭∂ξ(::Quad4,η::Float64)
    ∂N₁∂ξ = - 0.25*(1-η)
    ∂N₂∂ξ =   0.25*(1-η)
    ∂N₃∂ξ =   0.25*(1+η)
    ∂N₄∂ξ = - 0.25*(1+η)
    return (∂N₁∂ξ,∂N₂∂ξ,∂N₃∂ξ,∂N₄∂ξ)
end
function get∂𝝭∂η(::Quad4,ξ::Float64)
    ∂N₁∂η = - 0.25*(1-ξ)
    ∂N₂∂η = - 0.25*(1+ξ)
    ∂N₃∂η =   0.25*(1+ξ)
    ∂N₄∂η =   0.25*(1-ξ)
    return (∂N₁∂η,∂N₂∂η,∂N₃∂η,∂N₄∂η)
end

function get𝝭(::Quad4,ξ::Float64,η::Float64)
    N₁ = 0.25*(1.0-ξ)*(1.0-η)
    N₂ = 0.25*(1.0+ξ)*(1.0-η)
    N₃ = 0.25*(1.0+ξ)*(1.0+η)
    N₄ = 0.25*(1.0-ξ)*(1.0+η)
    return N₁,N₂,N₃,N₄
end

Base.issubset(a<:AbstractGeometry,b<:AbstractGeometry) = a.vertex ⊆ b.vertex