module Stokes

using ..ApproxOperator: AbstractElement

#===== 粘性项算子：μ∫2∇u:∇v dΩ → 对应矩阵 A =====#
"""
    ∫∫μ∇u∇vdxdy(aᵤ, k; μ)

Viscous term of the Stokes problem, `2μ ∫ ∇u:∇v dΩ`.

The factor 2 comes from `∇u:∇v` over the symmetric part rather than the
dot product of gradients alone:

    k[2I-1, 2J-1] += μ(2 ∂ₓNᵢ ∂ₓNⱼ + ∂ᵧNᵢ ∂ᵧNⱼ) w
    k[2I,   2J  ] += μ(∂ₓNᵢ ∂ₓNⱼ + 2 ∂ᵧNᵢ ∂ᵧNⱼ) w

Two velocity DOFs per node.

# Keyword arguments
- `μ` : dynamic viscosity

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫∫μ∇u∇vdxdy(aᵤ::T, k::AbstractMatrix{Float64}; μ::Real) where T<:AbstractElement
    𝓒 = aᵤ.𝓒; 𝓖 = aᵤ.𝓖
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]  # 速度形函数 x 导数
        B₂ = ξ[:∂𝝭∂y]  # 速度形函数 y 导数
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼   # 速度自由度全局索引（每个节点2自由度）
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                # 粘性项贡献：μ ∫ (∇u_x ⋅ ∇v_x + ∇u_y ⋅ ∇v_y) dΩ
                k[2I-1,2J-1] += μ * (2*B₁[i]*B₁[j] + B₂[i]*B₂[j]) * 𝑤
                k[2I,2J-1]   += μ * (B₁[i]*B₂[j]) * 𝑤
                k[2I-1,2J]   += μ * (B₂[i]*B₁[j]) * 𝑤
                k[2I,2J]     += μ * (B₁[i]*B₁[j] + 2*B₂[i]*B₂[j]) * 𝑤
            end
        end
    end
end

end