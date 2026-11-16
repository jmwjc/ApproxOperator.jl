# Quadratic segment (Seg3), nodes ordered as
#   1: ξ = -1,  2: ξ = +1,  3: ξ = 0
#
# Shape functions and their ξ-derivatives:
#   N₁ = ½ξ(ξ-1)   dN₁/dξ = ξ - ½
#   N₂ = ½ξ(ξ+1)   dN₂/dξ = ξ + ½
#   N₃ = 1 - ξ²    dN₃/dξ = -2ξ
#
# `𝐽` stores det(∂x/∂ξ) = |∂x/∂ξ| (arc-length scale), see importmsh.jl.
# For a 1-D element the chain rule therefore reads
#       ∂Nᵢ/∂x = (dNᵢ/dξ) / 𝐽,
# i.e. the ξ-derivative is *divided* by the Jacobian, not multiplied by it.

function set𝝭!(::Element{:Seg3},x::Node)
    ξ = x.ξ
    𝝭 = x[:𝝭]
    𝝭[1] = 0.5*ξ*(ξ-1.0)
    𝝭[2] = 0.5*ξ*(ξ+1.0)
    𝝭[3] = 1.0-ξ^2
end

function set∇𝝭!(ap::Element{:Seg3},x::Node)
    𝐽 = ap.𝐽
    ξ = x.ξ
    𝝭 = x[:𝝭]
    𝝭[1] = 0.5*ξ*(ξ-1.0)
    𝝭[2] = 0.5*ξ*(ξ+1.0)
    𝝭[3] = 1.0-ξ^2
    ∂𝝭∂x = x[:∂𝝭∂x]
    ∂𝝭∂x[1] = (ξ-0.5)/𝐽
    ∂𝝭∂x[2] = (ξ+0.5)/𝐽
    ∂𝝭∂x[3] = -2.0*ξ/𝐽
end
