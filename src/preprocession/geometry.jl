abstract type AbstractGeometry end

struct Poi1<:AbstractGeometry
    i::Tuple{Int}
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
end

struct Seg2<:AbstractGeometry
    i::NTuple{2,Int}
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
end

struct Seg3<:AbstractGeometry
    i::NTuple{3,Int}
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
end

struct Tri3<:AbstractGeometry
    i::NTuple{3,Int}
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
end

struct Tri6<:AbstractGeometry
    i::NTuple{6,Int}
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
end

struct Quad4<:AbstractGeometry
    i::NTuple{4,Int}
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
end

struct Quad8<:AbstractGeometry
    i::NTuple{8,Int}
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
end

struct Quad9<:AbstractGeometry
    i::NTuple{9,Int}
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
end

struct Tet4<:AbstractGeometry
    i::NTuple{4,Int}
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
end
struct Hex8<:AbstractGeometry
    i::NTuple{8,Int}
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
end

struct UnStructGeo<:AbstractGeometry
    i::Set{Int}
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
end

function (a::Poi1)(::Any)
    x = a.x[a.i[1]]
    y = a.y[a.i[1]]
    z = a.z[a.i[1]]
    return x, y, z
end

function (a::Seg2)(ξ::Float64)
    x₁ = a.x[a.i[1]]
    y₁ = a.y[a.i[1]]
    z₁ = a.z[a.i[1]]
    x₂ = a.x[a.i[2]]
    y₂ = a.y[a.i[2]]
    z₂ = a.z[a.i[2]]
    N₁ = 0.5*(1.0-ξ)
    N₂ = 0.5*(1.0+ξ)
    return N₁*x₁+N₂*x₂,
           N₁*y₁+N₂*y₂,
           N₁*z₁+N₂*z₂
end

function (a::Seg2)(b::Poi1,::Any)
    i = findfirst(x->x==b.i[1],a.i)
    if i == 1
        ξ = -1.0
    elseif i == 2
        ξ = 1.0
    else
        return (nothing,nothing)
    end
    return (a(ξ),ξ)
end

function (a::Seg3)(ξ::Float64)
    x₁ = a.x[a.i[1]]
    y₁ = a.y[a.i[1]]
    z₁ = a.z[a.i[1]]
    x₂ = a.x[a.i[3]]
    y₂ = a.y[a.i[3]]
    z₂ = a.z[a.i[3]]
    x₃ = a.x[a.i[2]]
    y₃ = a.y[a.i[2]]
    z₃ = a.z[a.i[2]]
    N₁ = 0.5*ξ*(ξ-1.0)
    N₂ = 1.0-ξ^2
    N₃ = 0.5*ξ*(ξ+1.0)
    return N₁*x₁+N₂*x₂+N₃*x₃,
           N₁*y₁+N₂*y₂+N₃*y₃,
           N₁*z₁+N₂*z₂+N₃*z₃
end

function (a::Seg3)(b::Poi1,::Any)
    i = findfirst(x->x==b.i[1],a.i)
    if i == 1
        return -1.0
    elseif i == 2
        return 1.0
    else
        return nothing
    end
end

function (a::Tri3)(ξ::Float64,η::Float64)
    x₁ = a.x[a.i[1]]
    y₁ = a.y[a.i[1]]
    z₁ = a.z[a.i[1]]
    x₂ = a.x[a.i[2]]
    y₂ = a.y[a.i[2]]
    z₂ = a.z[a.i[2]]
    x₃ = a.x[a.i[3]]
    y₃ = a.y[a.i[3]]
    z₃ = a.z[a.i[3]]
    N₁ = ξ
    N₂ = η
    N₃ = 1.0-ξ-η
    return N₁*x₁+N₂*x₂+N₃*x₃,
           N₁*y₁+N₂*y₂+N₃*y₃,
           N₁*z₁+N₂*z₂+N₃*z₃
end

function (a::Tri3)(b::Poi1,::Any)
    i = findfirst(x->x==b.i[1],a.i)
    if i == 1
        ξ = 1.0; η = 0.0
    elseif i == 2
        ξ = 0.0; η = 1.0
    elseif i == 3
        ξ = 0.0; η = 0.0
    else
        return (nothing,nothing)
    end
    return (a(ξ,η),(ξ,η))
end

function (a::Tri3)(b::Seg2,ξ₀::Float64)
    i = findfirst(x->x==b.i[1],a.i)
    if i == 1
        ξ = 0.5*(1.0-ξ₀); η = 0.5*(1.0+ξ₀)
    elseif i == 2
        ξ = 0.0; η = 0.5*(1.0-ξ₀)
    elseif i == 3
        ξ = 0.5*(1.0+ξ₀); η = 0.0
    else
        return (nothing,nothing)
    end
    return (a(ξ,η),(ξ,η))
end

function (a::Tri6)(ξ::Float64,η::Float64)
    γ = 1.0-ξ-η
    x₁ = a.x[a.i[1]];y₁ = a.y[a.i[1]];z₁ = a.z[a.i[1]]
    x₂ = a.x[a.i[2]];y₂ = a.y[a.i[2]];z₂ = a.z[a.i[2]]
    x₃ = a.x[a.i[3]];y₃ = a.y[a.i[3]];z₃ = a.z[a.i[3]]
    x₄ = a.x[a.i[4]];y₄ = a.y[a.i[4]];z₄ = a.z[a.i[4]]
    x₅ = a.x[a.i[5]];y₅ = a.y[a.i[5]];z₅ = a.z[a.i[5]]
    x₆ = a.x[a.i[6]];y₆ = a.y[a.i[6]];z₆ = a.z[a.i[6]]
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
function (a::Tet4)(ξ::Float64,η::Float64,ζ::Float64)
    
    x₁ = a.x[a.i[1]];y₁ = a.y[a.i[1]];z₁ = a.z[a.i[1]]
    x₂ = a.x[a.i[2]];y₂ = a.y[a.i[2]];z₂ = a.z[a.i[2]]
    x₃ = a.x[a.i[3]];y₃ = a.y[a.i[3]];z₃ = a.z[a.i[3]]
    x₄ = a.x[a.i[4]];y₄ = a.y[a.i[4]];z₄ = a.z[a.i[4]]
    N₁ = 1.0-x.ξ-x.η-x.ζ
    N₂ = x.ξ
    N₃ = x.η
    N₄ = x.ζ
    return x₁*N₁+x₂*N₂+x₃*N₃+x₄*N₄,
           y₁*N₁+y₂*N₂+y₃*N₃+y₄*N₄,
           z₁*N₁+z₂*N₂+z₃*N₃+z₄*N₄
end

function (a::Quad4)(ξ::Float64,η::Float64)
    x₁ = a.x[a.i[1]]
    y₁ = a.y[a.i[1]]
    z₁ = a.z[a.i[1]]
    x₂ = a.x[a.i[2]]
    y₂ = a.y[a.i[2]]
    z₂ = a.z[a.i[2]]
    x₃ = a.x[a.i[3]]
    y₃ = a.y[a.i[3]]
    z₃ = a.z[a.i[3]]
    x₄ = a.x[a.i[4]]
    y₄ = a.y[a.i[4]]
    z₄ = a.z[a.i[4]]
    N₁,N₂,N₃,N₄ = get𝝭(a,ξ,η)
    return (x₁*N₁+x₂*N₂+x₃*N₃+x₄*N₄,y₁*N₁+y₂*N₂+y₃*N₃+y₄*N₄,z₁*N₁+z₂*N₂+z₃*N₃+z₄*N₄)
end

function (a::Quad8)(ξ::Float64,η::Float64)
    x₁ = a.x[a.i[1]]
    y₁ = a.y[a.i[1]]
    z₁ = a.z[a.i[1]]
    x₂ = a.x[a.i[2]]
    y₂ = a.y[a.i[2]]
    z₂ = a.z[a.i[2]]
    x₃ = a.x[a.i[3]]
    y₃ = a.y[a.i[3]]
    z₃ = a.z[a.i[3]]
    x₄ = a.x[a.i[4]]
    y₄ = a.y[a.i[4]]
    z₄ = a.z[a.i[4]]
    x₅ = a.x[a.i[5]]
    y₅ = a.y[a.i[5]]
    z₅ = a.z[a.i[5]]
    x₆ = a.x[a.i[6]]
    y₆ = a.y[a.i[6]]
    z₆ = a.z[a.i[6]]
    x₇ = a.x[a.i[7]]
    y₇ = a.y[a.i[7]]
    z₇ = a.z[a.i[7]]
    x₈ = a.x[a.i[8]]
    y₈ = a.y[a.i[8]]
    z₈ = a.z[a.i[8]]
    N₁,N₂,N₃,N₄,N₅,N₆,N₇,N₈ = get𝝭(a,ξ,η)
    return (x₁*N₁+x₂*N₂+x₃*N₃+x₄*N₄+x₅*N₅+x₆*N₆+x₇*N₇+x₈*N₈,
            y₁*N₁+y₂*N₂+y₃*N₃+y₄*N₄+y₅*N₅+y₆*N₆+y₇*N₇+y₈*N₈,
            z₁*N₁+z₂*N₂+z₃*N₃+z₄*N₄+z₅*N₅+z₆*N₆+z₇*N₇+z₈*N₈)
end
function (a::Hex8)(ξ::Float64,η::Float64,ζ::Float64)
    x₁ = a.x[a.i[1]]
    y₁ = a.y[a.i[1]]
    z₁ = a.z[a.i[1]]
    x₂ = a.x[a.i[2]]
    y₂ = a.y[a.i[2]]
    z₂ = a.z[a.i[2]]
    x₃ = a.x[a.i[3]]
    y₃ = a.y[a.i[3]]
    z₃ = a.z[a.i[3]]
    x₄ = a.x[a.i[4]]
    y₄ = a.y[a.i[4]]
    z₄ = a.z[a.i[4]]
    x₅ = a.x[a.i[5]]
    y₅ = a.y[a.i[5]]
    z₅ = a.z[a.i[5]]
    x₆ = a.x[a.i[6]]
    y₆ = a.y[a.i[6]]
    z₆ = a.z[a.i[6]]
    x₇ = a.x[a.i[7]]
    y₇ = a.y[a.i[7]]
    z₇ = a.z[a.i[7]]
    x₈ = a.x[a.i[8]]
    y₈ = a.y[a.i[8]]
    z₈ = a.z[a.i[8]]
    N₁,N₂,N₃,N₄,N₅,N₆,N₇,N₈ = get𝝭(a,ξ,η,ζ)
    return (x₁*N₁+x₂*N₂+x₃*N₃+x₄*N₄+x₅*N₅+x₆*N₆+x₇*N₇+x₈*N₈,
            y₁*N₁+y₂*N₂+y₃*N₃+y₄*N₄+y₅*N₅+y₆*N₆+y₇*N₇+y₈*N₈,
            z₁*N₁+z₂*N₂+z₃*N₃+z₄*N₄+z₅*N₅+z₆*N₆+z₇*N₇+z₈*N₈)
end
function get𝐴(a::Tri3)
    x₁ = a.x[a.i[1]]
    x₂ = a.x[a.i[2]]
    x₃ = a.x[a.i[3]]
    y₁ = a.y[a.i[1]]
    y₂ = a.y[a.i[2]]
    y₃ = a.y[a.i[3]]
    z₁ = a.z[a.i[1]]
    z₂ = a.z[a.i[2]]
    z₃ = a.z[a.i[3]]

    return 0.5*(x₁*y₂+x₂*y₃+x₃*y₁-x₂*y₁-x₃*y₂-x₁*y₃)
end
get𝐽(a::Tri3,::Float64,::Float64) = get𝐴(a)
function get𝐴(a::Tri6)
    x₁ = a.x[a.i[1]]
    x₂ = a.x[a.i[2]]
    x₃ = a.x[a.i[3]]
    y₁ = a.y[a.i[1]]
    y₂ = a.y[a.i[2]]
    y₃ = a.y[a.i[3]]
    z₁ = a.z[a.i[1]]
    z₂ = a.z[a.i[2]]
    z₃ = a.z[a.i[3]]

    return 0.5*(x₁*y₂+x₂*y₃+x₃*y₁-x₂*y₁-x₃*y₂-x₁*y₃)
end
get𝐽(a::Tri6,::Float64,::Float64) = get𝐴(a)

function get𝐿(a::Seg2)
    x₁ = a.x[a.i[1]]
    x₂ = a.x[a.i[2]]
    y₁ = a.y[a.i[1]]
    y₂ = a.y[a.i[2]]

    return ((x₂-x₁)^2+(y₂-y₁)^2)^0.5
end
get𝐽(a::Seg2,::Float64) = 0.5*get𝐿(a)

get𝐽(a::Poi1,::Float64) = 1

function get𝐿(a::Seg3)
    x₁ = a.x[a.i[1]]
    y₁ = a.y[a.i[1]]
    z₁ = a.z[a.i[1]]
    x₂ = a.x[a.i[2]]
    y₂ = a.y[a.i[2]]
    z₂ = a.z[a.i[2]]
    return ((x₂-x₁)^2+(y₂-y₁)^2+(z₂-z₁)^2)^0.5
end
get𝐽(a::Seg3,::Float64) = 0.5*get𝐿(a)

function get𝑱(a::Quad4,ξ::Float64,η::Float64)
    x₁ = a.x[a.i[1]]
    y₁ = a.y[a.i[1]]
    x₂ = a.x[a.i[2]]
    y₂ = a.y[a.i[2]]
    x₃ = a.x[a.i[3]]
    y₃ = a.y[a.i[3]]
    x₄ = a.x[a.i[4]]
    y₄ = a.y[a.i[4]]
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

function get𝑱(a::Quad8,ξ::Float64,η::Float64)
    x₁ = a.x[a.i[1]]
    y₁ = a.y[a.i[1]]
    x₂ = a.x[a.i[2]]
    y₂ = a.y[a.i[2]]
    x₃ = a.x[a.i[3]]
    y₃ = a.y[a.i[3]]
    x₄ = a.x[a.i[4]]
    y₄ = a.y[a.i[4]]
    x₅ = a.x[a.i[5]]
    y₅ = a.y[a.i[5]]
    x₆ = a.x[a.i[6]]
    y₆ = a.y[a.i[6]]
    x₇ = a.x[a.i[7]]
    y₇ = a.y[a.i[7]]
    x₈ = a.x[a.i[8]]
    y₈ = a.y[a.i[8]]
    ∂N₁∂ξ,∂N₂∂ξ,∂N₃∂ξ,∂N₄∂ξ,∂N₅∂ξ,∂N₆∂ξ,∂N₇∂ξ,∂N₈∂ξ = get∂𝝭∂ξ(a,ξ,η)
    ∂N₁∂η,∂N₂∂η,∂N₃∂η,∂N₄∂η,∂N₅∂η,∂N₆∂η,∂N₇∂η,∂N₈∂η = get∂𝝭∂η(a,ξ,η)
    J₁₁ = ∂N₁∂ξ*x₁ + ∂N₂∂ξ*x₂ + ∂N₃∂ξ*x₃ + ∂N₄∂ξ*x₄ + ∂N₅∂ξ*x₅ + ∂N₆∂ξ*x₆ + ∂N₇∂ξ*x₇ + ∂N₈∂ξ*x₈
    J₁₂ = ∂N₁∂η*x₁ + ∂N₂∂η*x₂ + ∂N₃∂η*x₃ + ∂N₄∂η*x₄ + ∂N₅∂η*x₅ + ∂N₆∂η*x₆ + ∂N₇∂η*x₇ + ∂N₈∂η*x₈
    J₂₁ = ∂N₁∂ξ*y₁ + ∂N₂∂ξ*y₂ + ∂N₃∂ξ*y₃ + ∂N₄∂ξ*y₄ + ∂N₅∂ξ*y₅ + ∂N₆∂ξ*y₆ + ∂N₇∂ξ*y₇ + ∂N₈∂ξ*y₈
    J₂₂ = ∂N₁∂η*y₁ + ∂N₂∂η*y₂ + ∂N₃∂η*y₃ + ∂N₄∂η*y₄ + ∂N₅∂η*y₅ + ∂N₆∂η*y₆ + ∂N₇∂η*y₇ + ∂N₈∂η*y₈
    return J₁₁,J₂₁,J₁₂,J₂₂
end

function get𝐽(a::Quad8,ξ::Float64,η::Float64)
    J₁₁,J₂₁,J₁₂,J₂₂ = get𝑱(a,ξ,η)
    return J₁₁*J₂₂-J₂₁*J₁₂
end

function get∂𝝭∂ξ(::Quad8,ξ::Float64,η::Float64)
    ∂N₁∂ξ =   0.25*(2*ξ+η)*(1-η)
    ∂N₂∂ξ =   0.25*(2*ξ-η)*(1-η)
    ∂N₃∂ξ =   0.25*(2*ξ+η)*(1+η)
    ∂N₄∂ξ =   0.25*(2*ξ-η)*(1+η)
    ∂N₅∂ξ = - ξ*(1-η)
    ∂N₆∂ξ =   0.5*(1-η^2)
    ∂N₇∂ξ = - ξ*(1+η)
    ∂N₈∂ξ = - 0.5*(1-η^2)
    return (∂N₁∂ξ,∂N₂∂ξ,∂N₃∂ξ,∂N₄∂ξ,∂N₅∂ξ,∂N₆∂ξ,∂N₇∂ξ,∂N₈∂ξ)
end

function get∂𝝭∂η(::Quad8,ξ::Float64,η::Float64)
    ∂N₁∂η =   0.25*(1-ξ)*(2*η+ξ)
    ∂N₂∂η =   0.25*(1+ξ)*(2*η-ξ)
    ∂N₃∂η =   0.25*(1+ξ)*(2*η+ξ)
    ∂N₄∂η =   0.25*(1-ξ)*(2*η-ξ)
    ∂N₅∂η = - 0.5*(1-ξ^2)
    ∂N₆∂η = - η*(1+ξ)
    ∂N₇∂η =   0.5*(1-ξ^2)
    ∂N₈∂η = - η*(1-ξ)
    return (∂N₁∂η,∂N₂∂η,∂N₃∂η,∂N₄∂η,∂N₅∂η,∂N₆∂η,∂N₇∂η,∂N₈∂η)
end

function get𝝭(::Quad8,ξ::Float64,η::Float64)
    N₁ = 0.25*(1.0-ξ)*(1.0-η)*(-ξ-η-1)
    N₂ = 0.25*(1.0+ξ)*(1.0-η)*(ξ-η-1)
    N₃ = 0.25*(1.0+ξ)*(1.0+η)*(ξ+η-1)
    N₄ = 0.25*(1.0-ξ)*(1.0+η)*(-ξ+η-1)
    N₅ = 0.5*(1.0-ξ^2)*(1.0-η)
    N₆ = 0.5*(1.0+ξ)*(1.0-η^2)
    N₇ = 0.5*(1.0-ξ^2)*(1.0+η)
    N₈ = 0.5*(1.0-ξ)*(1.0-η^2)
    return N₁,N₂,N₃,N₄,N₅,N₆,N₇,N₈
end

function get𝒏(a::Seg2)
    𝐿 = get𝐿(a)
    x₁ = a.x[a.i[1]]
    y₁ = a.y[a.i[1]]
    x₂ = a.x[a.i[2]]
    y₂ = a.y[a.i[2]]
    return (y₂-y₁)/𝐿, (x₁-x₂)/𝐿
end

function get𝑫(a::Seg2)
    return (-1.0,1.0)
end

Base.issubset(a::AbstractGeometry,b::AbstractGeometry) = a.i ⊆ b.i

function Base.intersect(as::Vector{Seg2},bs::Vector{Tri3})
    𝓢 = Tuple{Seg2,Tri3}[]
    for a in as
        for (i,b) in enumerate(bs)
            if a ⊆ b
                index = collect(findfirst(x->x==i,b.i) for i in a.i)
                if index == [1,3]
                    index = [2,1]
                elseif index == [3,1]
                    index = [1,2]
                else
                    index = sortperm(index)
                end
                c = Seg2(Tuple(a.i[index]),a.x,a.y,a.z)
                push!(𝓢,(c,b))
            end
        end
    end
    return 𝓢
end

function Base.intersect(as::Vector{T₁},bs::Vector{T₂}) where {T₁<:AbstractGeometry,T₂<:AbstractGeometry}
    𝓢 = Tuple{T₁,T₂}[]
    for a in as
        for (i,b) in enumerate(bs)
            if a ⊆ b
                index = sortperm(collect(findfirst(x->x==i,b.i) for i in a.i))
                c = T₁(Tuple(a.i[index]),a.x,a.y,a.z)
                push!(𝓢,(c,b))
            end
        end
    end
    return 𝓢
end

function getboundaries(as::Vector{Tri3})
    𝓑 = Set{Int}[]
    for a in as
        for i in a.i
            b = Set(setdiff(a.i,i))
            if b ∉ 𝓑
                push!(𝓑,b)
            end
        end
    end
    return 𝓑
end
