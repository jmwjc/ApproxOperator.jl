
struct TRElement{T} <: AbstractElement
    𝓒::Tuple{Int,Int,Vector{Node{(:𝐼,:𝐽),2}}}
    𝓖::Tuple{Int,Int,Vector{Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}}}
end

function set𝝭!(::TRElement{:Tri3},x::Node)
    ξ₁ = x.ξ
    ξ₂ = x.η
    ξ₃ = 1.0-x.ξ-x.η
    N₁ = ξ₂+ξ₃-ξ₁
    N₂ = ξ₃+ξ₁-ξ₂
    N₃ = ξ₁+ξ₂-ξ₃
    𝝭 = x[:𝝭]
    𝝭[1] = N₁
    𝝭[2] = N₂
    𝝭[3] = N₃
end

function set∇𝝭!(ap::TRElement{:Tri3},x::Node)
    ξ₁ = x.ξ
    ξ₂ = x.η
    ξ₃ = 1.0-x.ξ-x.η
    N₁ = ξ₂+ξ₃-ξ₁
    N₂ = ξ₃+ξ₁-ξ₂
    N₃ = ξ₁+ξ₂-ξ₃
    𝝭 = x[:𝝭]
    𝝭[1] = N₁
    𝝭[2] = N₂
    𝝭[3] = N₃
    v₁,v₂,v₃ = ap.𝓒
    x₁ = v₁.x
    x₂ = v₂.x
    x₃ = v₃.x
    y₁ = v₁.y
    y₂ = v₂.y
    y₃ = v₃.y
    𝐴 = ap.𝐴
    ∂𝝭∂x = x[:∂𝝭∂x]
    ∂𝝭∂y = x[:∂𝝭∂y]
    ∂𝝭∂x[1] = (y₃-y₂)/𝐴
    ∂𝝭∂x[2] = (y₁-y₃)/𝐴
    ∂𝝭∂x[3] = (y₂-y₁)/𝐴
    ∂𝝭∂y[1] = (x₂-x₃)/𝐴
    ∂𝝭∂y[2] = (x₃-x₁)/𝐴
    ∂𝝭∂y[3] = (x₁-x₂)/𝐴
end