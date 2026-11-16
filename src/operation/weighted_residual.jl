module WeightedResidual

using ..ApproxOperator: AbstractElement

# ===== 二维（时空）加权残余法算子 =====
# 二维坐标 (x,y) 中 y 为时间方向
# 1. Kuu: k ∫∫ Nᵢ ∂Nⱼ/∂t dxdt
"""
    ∫kNṄdxdt(ap, k; kᶜ)

Space–time weighted-residual operator `k ∫∫ Nᵢ ∂Nⱼ/∂t dx dt`.

In the 2-D space–time mesh the local coordinate `y` is the **time** direction,
so the time derivative is `ξ[:∂𝝭∂y]`.

# Keyword arguments
- `kᶜ` : the coefficient of the first-order time term (kept as `kᶜ` because
  `k` is the stiffness matrix)

# Element requirements
`set∇𝝭!` — reads `ξ[:𝝭]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫kNṄdxdt(ap::T, k::AbstractMatrix{Float64}; kᶜ::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        Ṅ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += kᶜ * N[i] * Ṅ[j] * 𝑤
            end
        end
    end
end

# 2. Kuv: k ∫∫ Nᵢ Nⱼ dxdt 
"""
    ∫kNNdxdt(ap, k; kᶜ)

Space–time weighted-residual mass-like term `-k ∫∫ Nᵢ Nⱼ dx dt`.

Note the **minus** sign: this is the reaction term, assembled into the same
matrix as [`∫kNṄdxdt`](@ref) so that both form one space–time operator.

# Keyword arguments
- `kᶜ` : the coefficient (kept as `kᶜ` because `k` is the stiffness matrix)

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.𝑤`.
"""
function ∫kNNdxdt(ap::T, k::AbstractMatrix{Float64}; kᶜ::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] -= kᶜ * N[i] * N[j] * 𝑤
            end
        end
    end
end

# 3. Kvv: ∫∫ Nᵢ ∂Nⱼ/∂t dxdt
"""
    ∫NṄdxdt(ap, k)

Space–time weighted-residual operator `∫∫ Nᵢ ∂Nⱼ/∂t dx dt` with unit
coefficient.

# Element requirements
`set∇𝝭!` — reads `ξ[:𝝭]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫NṄdxdt(ap::T, k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        Ṅ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += N[i] * Ṅ[j] * 𝑤
            end
        end
    end
end

# 4. Kvu: c² ∫∫ (∂Nᵢ/∂x ∂Nⱼ/∂x) dxdt
"""
    ∫c²B₁B₁dxdt(ap, k; c)

Space–time weighted-residual stiffness, `c² ∫∫ (∂ₓNᵢ)(∂ₓNⱼ) dx dt`.

# Keyword arguments
- `c` : wave speed, squared internally

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ.𝑤`.
"""
function ∫c²B₁B₁dxdt(ap::T, k::AbstractMatrix{Float64}; c::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += c^2 * B₁[i] * B₁[j] * 𝑤
            end
        end
    end
end
end