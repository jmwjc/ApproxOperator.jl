"""
    ∫vᵢθᵢds(ap, k, f)

Curved-beam penalty/constraint term over the reference surface,

    ∫ α (𝒏 · v)(𝒏 · θ) ds

with `𝒏` the surface covector from `ξ.n₁₁`, `ξ.n₂₂`, `ξ.n₁₂` and the prescribed
vector `gᵢ` from `ξ.g₁`, `ξ.g₂`, `ξ.g₃`. Three DOFs per node.

# Keyword arguments
- `α` : penalty coefficient

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.n₁₁`, `ξ.n₂₂`, `ξ.n₁₂`, `ξ.g₁`, `ξ.g₂`, `ξ.g₃`, `ξ.𝑤`.
"""
function ∫vᵢθᵢds(ap::T, k::AbstractMatrix{Float64}, f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        n₁₁ = ξ.n₁₁
        n₂₂ = ξ.n₂₂
        n₁₂ = ξ.n₁₂
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        g₃ = ξ.g₃
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += α*N[i]*n₁₁*N[j]*𝑤
                k[3*I-2,3*J-1] += α*N[i]*n₁₂*N[j]*𝑤
                k[3*I-1,3*J-2] += α*N[i]*n₁₂*N[j]*𝑤
                k[3*I-1,3*J-1] += α*N[i]*n₂₂*N[j]*𝑤
                k[3*I,3*J]     += α*N[i]*N[j]*𝑤
            end
            f[3*I-2] += α*N[i]*(n₁₁*g₁+n₁₂*g₂)*𝑤
            f[3*I-1] += α*N[i]*(n₁₂*g₁+n₂₂*g₂)*𝑤
            f[3*I]   += α*N[i]*g₃*𝑤
        end
    end
end

"""
    ∫κεγds(ap, k; EA, EI)
    ∫κεds(ap, k; EA, EI)

Curved-beam stiffness from membrane strain `ε`, curvature `κ` and shear `γ`,

    a = ∫ (EA ε(v)ε(u) + EI κ(v)κ(u) + kGA γ(v)γ(u)) ds

`∫κεds` is the same operator with the shear term omitted (thin/Euler–Bernoulli
limit). Three DOFs per node: two translations and one rotation.

# Keyword arguments
- `EA` : axial stiffness
- `EI` : bending stiffness

# Element requirements
`set∇𝝭!` — reads `ξ[:𝝭]`, `ξ[:∂𝝭∂x]`, `ξ.𝑤`.
"""
function ∫κεγds(ap::T, k::AbstractMatrix{Float64}; EA::Real, EI::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        B = ξ[:∂𝝭∂x]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += (N[i]*kGA/R^2*N[j]+B[i]*EA*B[j])*𝑤
                k[3*I-2,3*J-1] += (N[i]*kGA/R*B[j]-B[i]*EA/R*N[j])*𝑤
                k[3*I-2,3*J]   += (-N[i]*kGA/R*N[j])*𝑤
                k[3*I-1,3*J-2] += (B[i]*kGA/R*N[j]-N[i]*EA/R*B[j])*𝑤
                k[3*I-1,3*J-1] += (B[i]*kGA*B[j]+N[i]*EA/R^2*N[j])*𝑤
                k[3*I-1,3*J]   += (-B[i]*kGA*N[j])*𝑤
                k[3*I,3*J-2]   += (-N[i]*kGA/R*N[j])*𝑤
                k[3*I,3*J-1]   += (-N[i]*kGA*B[j])*𝑤
                k[3*I,3*J]     += (B[i]*EI*B[j]+N[i]*kGA*N[j])*𝑤
            end
        end
    end
end

"""
    ∫δNEA⁻¹Nds(a::T,b::S;k::AbstractMatrix{Float64},g::AbstractMatrix{Float64})

Element-set overload of `∫δNEA⁻¹Nds`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function ∫δNEA⁻¹Nds(a::T,b::S, k::AbstractMatrix{Float64}, g::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒₁ = a.𝓒; 𝓖₁ = a.𝓖
    𝓒₂ = a.𝓒; 𝓖₂ = a.𝓖
    for (ξ₁,ξ₂) in zip(𝓖₁,𝓖₂)
        N = ξ₁[:𝝭]
        B = ξ₁[:∂𝝭∂x]
        Ñ = ξ₂[:𝝭]
        𝑤 = ξ₁.𝑤
        𝑤ᵇ = ξ₁.𝑤ᵇ
        for (i,xᵢ) in enumerate(𝓒₁)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒₁)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += N[i]*N[j]/EA*𝑤
                k[2*I,2*J]     += N[i]*N[j]/EA*𝑤
            end
            for (j,xⱼ) in enumerate(𝓒₂)
                J = xⱼ.𝐼
                g[2*I-1,2*J-1] += B[i]*Ñ[j]*𝑤 - N[i]*Ñ[j]*𝑤ᵇ
                g[2*I,2*J]     += B[i]*Ñ[j]*𝑤 - N[i]*Ñ[j]*𝑤ᵇ
            end
        end
    end
end

"""
    ∫δMEI⁻¹Mds(a::T,b::S;k::AbstractMatrix{Float64},g::AbstractMatrix{Float64})

Element-set overload of `∫δMEI⁻¹Mds`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function ∫δMEI⁻¹Mds(a::T,b::S, k::AbstractMatrix{Float64}, g::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒₁ = a.𝓒; 𝓖₁ = a.𝓖
    𝓒₂ = a.𝓒; 𝓖₂ = a.𝓖
    for (ξ₁,ξ₂) in zip(𝓖₁,𝓖₂)
        N = ξ₁[:𝝭]
        B = ξ₁[:∂𝝭∂x]
        C = ξ₁[:∂²𝝭∂x²]
        Ñ = ξ₂[:𝝭]
        B̃ = ξ₂[:∂𝝭∂x]
        𝑤 = ξ₁.𝑤
        𝑤ᵇ = ξ₁.𝑤ᵇ
        for (i,xᵢ) in enumerate(𝓒₁)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒₁)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += N[i]*N[j]/EI*𝑤
                k[2*I,2*J]     += N[i]*N[j]/EI*𝑤
            end
            for (j,xⱼ) in enumerate(𝓒₂)
                J = xⱼ.𝐼
                g[2*I-1,2*J-1] += C[i]*Ñ[j]*𝑤 - (B[i]*Ñ[j] - N[i]*B̃[j])*𝑤ᵇ
                g[2*I,2*J]     += C[i]*Ñ[j]*𝑤 - (B[i]*Ñ[j] - N[i]*B̃[j])*𝑤ᵇ
            end
        end
    end
end

"""
    Bᵢvᵢ(a::T,b::S;k::AbstractMatrix{Float64},f::AbstractVector{Float64})

Curved-beam strain matrix applied to the displacement, assembling
`∫ Bᵀ v ds` into the right-hand side.

# Element requirements
`set∇𝝭!` — reads `ξ[:𝝭]`, first derivatives and the surface metric.
"""
function Bᵢvᵢ(a::T,b::S, k::AbstractMatrix{Float64}, f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒₁ = a.𝓒; 𝓖₁ = a.𝓖
    𝓒₂ = a.𝓒; 𝓖₂ = a.𝓖
    for (ξ₁,ξ₂) in zip(𝓖₁,𝓖₂)
        N = ξ₁[:𝝭]
        B = ξ₁[:∂𝝭∂x]
        C = ξ₁[:∂²𝝭∂x²]
        Ñ = ξ₂[:𝝭]
        𝑤 = ξ₁.𝑤
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        n₁ = ξ.n₁    
        for (i,xᵢ) in enumerate(𝓒₁)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒₂)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += n₁*B[i]*Ñ[j]*𝑤
                k[2*I,2*J]     += n₁*B[i]*Ñ[j]*𝑤
            end
            f[2*I-1] += n₁*B[i]*g₁*𝑤
            f[2*I]   += n₁*B[i]*g₂*𝑤
        end
    end
end

"""
    ∫κεds(ap::T;k::AbstractMatrix{Float64}, EA, EI)

Element-set overload of `∫κεds`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function ∫κεds(ap::T, k::AbstractMatrix{Float64}; EA::Real, EI::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B = ξ[:∂𝝭∂x]
        C = ξ[:∂²𝝭∂x²]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += B[i]*(ξ.a₁[1]+ξ.a₁[2])*EA*(ξ.a₁[1]+ξ.a₁[2])*B[j]*𝑤
                k[2*I,2*J]     -= C[i]*(ξ.a₃[1]+ξ.a₃[2])*EI*(ξ.a₃[1]+ξ.a₃[2])*C[j]*𝑤
            end
        end
    end
end