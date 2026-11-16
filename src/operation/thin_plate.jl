module KirchhoffPlate
    
using ..ApproxOperator: AbstractElement

"""
    ∫κᵢⱼMᵢⱼdΩ(ap, k; D, ν)

Kirchhoff plate bending stiffness, `∫ κᵢⱼ(w) Mᵢⱼ(v) dΩ`.

    a(w,v) = ∫_Ω D[∇²w ∇²v + ν(∇²w v,ₓᵧ + w,ₓᵧ ∇²v) + 2(1-ν) w,ₓᵧ v,ₓᵧ] dΩ

The `w,ₓᵧ` terms mean the mixed second derivative, read from
`ξ[:∂²𝝭∂x∂y]`. Requires a C¹ element: `set∇²𝝭!`.

# Keyword arguments
- `D` : bending stiffness `E h³ / (12(1-ν²))`
- `ν` : Poisson ratio

# Element requirements
`set∇²𝝭!` — reads `ξ[:∂²𝝭∂x²]`, `ξ[:∂²𝝭∂x∂y]`, `ξ[:∂²𝝭∂y²]`, `ξ.𝑤`.
"""
function ∫κᵢⱼMᵢⱼdΩ(ap::T,k::AbstractMatrix{Float64}; ν::Real, D::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B₁₁ = ξ[:∂²𝝭∂x²]
        B₁₂ = ξ[:∂²𝝭∂x∂y]
        B₂₂ = ξ[:∂²𝝭∂y²]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += D*(B₁₁[i]*B₁₁[j] + ν*(B₁₁[i]*B₂₂[j] + B₂₂[i]*B₁₁[j]) + B₂₂[i]*B₂₂[j] + 2*(1-ν)*B₁₂[i]*B₁₂[j])*𝑤
            end
        end
    end
end

"""
    ∫κ̃ᵢⱼM̃ᵢⱼdΩ(ap, k; D, ν)

Kirchhoff plate bending stiffness with the **symmetric** mixed-derivative
treatment.

Identical to [`∫κᵢⱼMᵢⱼdΩ`](@ref) except that both `w,ₓᵧ` and `w,ᵧₓ` are kept
(`ξ[:∂²𝝭∂y∂x]` in addition to `ξ[:∂²𝝭∂x∂y]`), which is the form needed when
the approximation is not C¹ in the mixed derivative.

# Keyword arguments
- `D` : bending stiffness `E h³ / (12(1-ν²))`
- `ν` : Poisson ratio

# Element requirements
`set∇²𝝭!` — reads `ξ[:∂²𝝭∂x²]`, `ξ[:∂²𝝭∂x∂y]`, `ξ[:∂²𝝭∂y²]`, `ξ.𝑤`.
"""
function ∫κ̃ᵢⱼM̃ᵢⱼdΩ(ap::T,k::AbstractMatrix{Float64}; ν::Real, D::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B₁₁ = ξ[:∂²𝝭∂x²]
        B₁₂ = ξ[:∂²𝝭∂x∂y]
        B₂₁ = ξ[:∂²𝝭∂y∂x]
        B₂₂ = ξ[:∂²𝝭∂y²]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += D*(B₁₁[i]*B₁₁[j] + ν*(B₁₁[i]*B₂₂[j] + B₂₂[i]*B₁₁[j]) + B₂₂[i]*B₂₂[j] + (1-ν)*(B₁₂[i]*B₁₂[j]+B₂₁[i]*B₂₁[j]))*𝑤
            end
        end
    end
end

"""
    ∫wqdΩ(ap::T,f::AbstractVector{Float64})

Transverse load vector, `f[I] += ∫ q Nᵢ dΩ`.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.q`, `ξ.𝑤`.
"""
function ∫wqdΩ(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        q = ξ.q
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] += N[i]*q*𝑤
        end
    end
end

"""
    ∫wVdΓ(ap::T,f::AbstractVector{Float64})

Kirchhoff shear boundary load `∫ V Nᵢ dΓ`.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.V`, `ξ.𝑤`.
"""
function ∫wVdΓ(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        V = ξ.V
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] += N[i]*V*𝑤
        end
    end
end

"""
    ∫θₙMₙₙdΓ(ap::T,f::AbstractVector{Float64})

Rotation–moment boundary coupling `∫ θₙ Mₙₙ dΓ`.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.M`, `ξ.𝑤`.
"""
function ∫θₙMₙₙdΓ(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        M = ξ.M
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] -= (B₁[i]*n₁+B₂[i]*n₂)*M*𝑤
        end
    end
end

"""
    ∫∇wMdΓ(ap::T,f::AbstractVector{Float64})

Deflection-gradient boundary term `∫ (∂ₙw) M dΓ`.

# Element requirements
`set∇𝝭!` — reads gradients, `ξ.M`, `ξ.n₁`, `ξ.n₂`, `ξ.𝑤`.
"""
function ∫∇wMdΓ(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        M₁ = ξ.M₁
        M₂ = ξ.M₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] -= (B₁[i]*M₁+B₂[i]*M₂)*𝑤
        end
    end
end

"""
    ∫∂ₙwθdΓ(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}; α)

Boundary penalty on the normal derivative, `α ∫ (∂ₙw) θ dΓ`.

# Keyword arguments
- `α` : Nitsche/penalty parameter

# Element requirements
`set∇𝝭!` — reads gradients, `ξ.n₁`, `ξ.n₂`, `ξ.θ`, `ξ.𝑤`.
"""
function ∫∂ₙwθdΓ(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}; α::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        θ = ξ.θ
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            θᵢ = B₁[i]*n₁+B₂[i]*n₂
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                θⱼ = B₁[j]*n₁+B₂[j]*n₂
                k[I,J] += α*θᵢ*θⱼ*𝑤
            end
            f[I] += α*θᵢ*θ*𝑤
        end
    end
end

"""
    ∫MₙₙθdΓ(ap, k, f; D, ν, α)

Kirchhoff plate boundary term `∫ (Mₙₙ(w) + α θ) θ̄ dΓ` — normal bending moment
plus penalty, as used in the weak (Nitsche/penalty) free-edge condition.

    Mₙₙ = -D[(n₁² + νn₂²) w,ₓₓ + 2(1-ν)n₁n₂ w,ₓᵧ + (νn₁² + n₂²) w,ᵧᵧ]

# Keyword arguments
- `D` : bending stiffness
- `ν` : Poisson ratio
- `α` : Nitsche/penalty parameter

# Element requirements
`set∇²𝝭!` — reads the second derivatives, `ξ.n₁`, `ξ.n₂`, `ξ.θ`, `ξ.𝑤`.
"""
function ∫MₙₙθdΓ(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}; ν::Real, D::Real, α::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        D₁₁ = -D*(n₁^2+ν*n₂^2)
        D₁₂ = -2*D*n₁*n₂*(1-ν)
        D₂₂ = -D*(ν*n₁^2+n₂^2)
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₁₁ = ξ[:∂²𝝭∂x²]
        B₁₂ = ξ[:∂²𝝭∂x∂y]
        B₂₂ = ξ[:∂²𝝭∂y²]
        θ = ξ.θ
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            θᵢ = B₁[i]*n₁ + B₂[i]*n₂
            Mᵢ = D₁₁*B₁₁[i] + D₁₂*B₁₂[i] + D₂₂*B₂₂[i]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                θⱼ = B₁[j]*n₁ + B₂[j]*n₂
                Mⱼ = D₁₁*B₁₁[j] + D₁₂*B₁₂[j] + D₂₂*B₂₂[j]
                k[I,J] += (Mᵢ*θⱼ+θᵢ*Mⱼ+α*θᵢ*θⱼ)*𝑤
            end
            f[I] += (Mᵢ+α*θᵢ)*θ*𝑤
        end
    end
end

"""
    ∫VgdΓ(ap, k, f; D, ν, α)

Kirchhoff plate Kirchhoff-shear boundary term, `∫ (Vₙ(w) + α w) ḡ dΓ`.

    Vₙ = -D[∂ₙ∇²w + (1-ν)∂ₛ(w,ₙₛ)]

The effective shear force is the combination of the twisting-moment derivative
and the normal derivative of the Laplacian; the `α` term is the penalty that
closes the weak form.

# Keyword arguments
- `D` : bending stiffness
- `ν` : Poisson ratio
- `α` : Nitsche/penalty parameter

# Element requirements
`set∇²𝝭!` — reads the second derivatives, `ξ.n₁`, `ξ.n₂`, `ξ.s₁`, `ξ.s₂`,
`ξ.g`, `ξ.𝑤`.
"""
function ∫VgdΓ(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}; ν::Real, D::Real, α::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        s₁ = -ξ.n₂
        s₂ = ξ.n₁
        D₁₁₁ = -D*(n₁ + n₁*s₁*s₁ + ν*n₂*s₁*s₂)
        D₁₁₂ = -D*(n₂ + n₂*s₁*s₁ + 2*n₁*s₁*s₂ + (n₂*s₂*s₂ - n₂*s₁*s₁ - n₁*s₁*s₂)*ν)
        D₁₂₂ = -D*(n₁ + n₁*s₂*s₂ + 2*n₂*s₁*s₂ + (n₁*s₁*s₁ - n₁*s₂*s₂ - n₂*s₁*s₂)*ν)
        D₂₂₂ = -D*(n₂ + n₂*s₂*s₂ + ν*n₁*s₁*s₂)
        N = ξ[:𝝭]
        B₁₁₁ = ξ[:∂³𝝭∂x³]
        B₁₁₂ = ξ[:∂³𝝭∂x²∂y]
        B₁₂₂ = ξ[:∂³𝝭∂x∂y²]
        B₂₂₂ = ξ[:∂³𝝭∂y³]
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            Vᵢ = D₁₁₁*B₁₁₁[i] + D₁₁₂*B₁₁₂[i] + D₁₂₂*B₁₂₂[i] + D₂₂₂*B₂₂₂[i]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                Vⱼ = D₁₁₁*B₁₁₁[j] + D₁₁₂*B₁₁₂[j] + D₁₂₂*B₁₂₂[j] + D₂₂₂*B₂₂₂[j]
                k[I,J] += (-Vᵢ*N[j]-N[i]*Vⱼ+α*N[i]*N[j])*𝑤
            end
            f[I] += (-Vᵢ+α*N[i])*g*𝑤
        end
    end
end

"""
    ∫M̃ₙₙθdΓ(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64})

Element-set overload of `∫M̃ₙₙθdΓ`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function ∫M̃ₙₙθdΓ(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    n₁ = 𝓖[1].n₁
    n₂ = 𝓖[1].n₂
    s₁ = 𝓖[1].s₁
    s₂ = 𝓖[1].s₂
    D = 𝓖[1].D
    ν = 𝓖[1].ν
    D₁₁ = -D*(n₁^2+ν*n₂^2)
    D₁₂ = -2*D*n₁*n₂*(1-ν)
    D₂₂ = -D*(ν*n₁^2+n₂^2)
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₁₁ = ξ[:∂²𝝭∂x²]
        B₁₂ = ξ[:∂²𝝭∂x∂y]
        B₂₂ = ξ[:∂²𝝭∂y²]
        B̄₁₁ = ξ[:∂²𝝭∂x²_]
        B̄₁₂ = ξ[:∂²𝝭∂x∂y_]
        B̄₂₂ = ξ[:∂²𝝭∂y²_]
        θ = ξ.θ
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            θᵢ = B₁[i]*n₁ + B₂[i]*n₂
            Mᵢ = D₁₁*B₁₁[i] + D₁₂*B₁₂[i] + D₂₂*B₂₂[i]
            M̄ᵢ = D₁₁*B̄₁₁[i] + D₁₂*B̄₁₂[i] + D₂₂*B̄₂₂[i]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                θⱼ = B₁[j]*n₁ + B₂[j]*n₂
                Mⱼ = D₁₁*B₁₁[j] + D₁₂*B₁₂[j] + D₂₂*B₂₂[j]
                k[I,J] += (Mᵢ*θⱼ+θᵢ*Mⱼ-M̄ᵢ*θⱼ)*𝑤
            end
            f[I] += (Mᵢ-M̄ᵢ)*θ*𝑤
        end
    end
end

"""
    ∫ṼgdΓ(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64})

Element-set overload of `∫ṼgdΓ`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function ∫ṼgdΓ(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    n₁ = 𝓖[1].n₁
    n₂ = 𝓖[1].n₂
    s₁ = 𝓖[1].s₁
    s₂ = 𝓖[1].s₂
    D = 𝓖[1].D
    ν = 𝓖[1].ν
    D₁₁₁ = -D*(n₁ + n₁*s₁*s₁ + ν*n₂*s₁*s₂)
    D₁₁₂ = -D*(n₁*s₁*s₂ + (n₂*s₂*s₂ + n₂)*ν)
    D₁₂₁ = -D*(n₂ + n₂*s₁*s₁ + n₁*s₁*s₂ + (-n₂ - n₂*s₁*s₁ - n₁*s₁*s₂)*ν)
    D₁₂₂ = -D*(n₁ + n₁*s₂*s₂ + n₂*s₁*s₂ + (-n₁ - n₁*s₂*s₂ - n₂*s₁*s₂)*ν)
    D₂₂₁ = -D*(n₂*s₁*s₂ + (n₁*s₁*s₁ + n₁)*ν)
    D₂₂₂ = -D*(n₂ + n₂*s₂*s₂ + ν*n₁*s₁*s₂)
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁₁₁ = ξ[:∂∂²𝝭∂x²∂x]
        B₁₁₂ = ξ[:∂∂²𝝭∂x²∂y]
        B₁₂₁ = ξ[:∂∂²𝝭∂x∂y∂x]
        B₁₂₂ = ξ[:∂∂²𝝭∂x∂y∂y]
        B₂₂₁ = ξ[:∂∂²𝝭∂y²∂x]
        B₂₂₂ = ξ[:∂∂²𝝭∂y²∂y]
        B̄₁₁₁ = ξ[:∂∂²𝝭∂x²∂x_]
        B̄₁₁₂ = ξ[:∂∂²𝝭∂x²∂y_]
        B̄₁₂₁ = ξ[:∂∂²𝝭∂x∂y∂x_]
        B̄₁₂₂ = ξ[:∂∂²𝝭∂x∂y∂y_]
        B̄₂₂₁ = ξ[:∂∂²𝝭∂y²∂x_]
        B̄₂₂₂ = ξ[:∂∂²𝝭∂y²∂y_]
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            Vᵢ = D₁₁₁*B₁₁₁[i] + D₁₁₂*B₁₁₂[i] + D₁₂₁*B₁₂₁[i] + D₁₂₂*B₁₂₂[i] + D₂₂₁*B₂₂₁[i] + D₂₂₂*B₂₂₂[i]
            V̄ᵢ = D₁₁₁*B̄₁₁₁[i] + D₁₁₂*B̄₁₁₂[i] + D₁₂₁*B̄₁₂₁[i] + D₁₂₂*B̄₁₂₂[i] + D₂₂₁*B̄₂₂₁[i] + D₂₂₂*B̄₂₂₂[i]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                Vⱼ = D₁₁₁*B₁₁₁[j] + D₁₁₂*B₁₁₂[j] + D₁₂₁*B₁₂₁[j] + D₁₂₂*B₁₂₂[j] + D₂₂₁*B₂₂₁[j] + D₂₂₂*B₂₂₂[j]
                k[I,J] -= (Vᵢ*N[j]+N[i]*Vⱼ-V̄ᵢ*N[j])*𝑤
            end
            f[I] -= (Vᵢ-V̄ᵢ)*g*𝑤
        end
    end
end

"""
    wΔMₙₛ(ap::T,f::AbstractVector{Float64})

Corner-force load of the Kirchhoff plate, `∫ w ΔMₙₛ`.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.Δ`, `ξ.𝑤`.
"""
function wΔMₙₛ(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        ΔM = ξ.ΔM
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] -= N[i]*ΔM
        end
    end
end

"""
    ΔMₙₛg(ap, k, f; D, ν, α)
    ΔM̃ₙₛg(ap, k, f; D, ν)

Corner-force terms of the Kirchhoff plate, `∫ ΔMₙₛ ḡ`.

At a boundary corner the twisting moment `Mₙₛ` jumps; the jump `ΔMₙₛ` acts as a
concentrated force. `ξ.Δ` carries the corner sign (+1, -1 or 0 at a smooth
point), which is why the quadrature rule includes the boundary endpoints.

# Keyword arguments
- `D` : bending stiffness
- `ν` : Poisson ratio
- `α` : Nitsche/penalty parameter (used by `ΔMₙₛg`)

# Element requirements
`set∇²𝝭!` — reads the second derivatives, `ξ.n₁`, `ξ.n₂`, `ξ.s₁`, `ξ.s₂`,
`ξ.Δ`, `ξ.g`, `ξ.𝑤`.
"""
function ΔMₙₛg(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}; ν::Real, D::Real, α::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        Δn₁s₁ = ξ.Δn₁s₁
        Δn₁s₂n₂s₁ = ξ.Δn₁s₂n₂s₁
        Δn₂s₂ = ξ.Δn₂s₂
        D₁₁ = - D*(Δn₁s₁+Δn₂s₂*ν)
        D₁₂ = - D*(1-ν)*Δn₁s₂n₂s₁
        D₂₂ = - D*(Δn₁s₁*ν+Δn₂s₂)
        N = ξ[:𝝭]
        B₁₁ = ξ[:∂²𝝭∂x²]
        B₁₂ = ξ[:∂²𝝭∂x∂y]
        B₂₂ = ξ[:∂²𝝭∂y²]
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            ΔMₙₛᵢ = D₁₁*B₁₁[i] + D₁₂*B₁₂[i] + D₂₂*B₂₂[i]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                ΔMₙₛⱼ = D₁₁*B₁₁[j] + D₁₂*B₁₂[j] + D₂₂*B₂₂[j]
                k[I,J] += ΔMₙₛᵢ*N[j] + N[i]*ΔMₙₛⱼ + α*N[i]*N[j]
            end
            f[I] += (ΔMₙₛᵢ + α*N[i])*g
        end
    end
end

"""
    ΔM̃ₙₛg(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}; ν, D)

Element-set overload of `ΔM̃ₙₛg`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.


    ∫ ΔM̃ₙₛ ḡ dΓ,  the symmetric-mixed-derivative twin of `ΔMₙₛg`.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function ΔM̃ₙₛg(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}; ν::Real, D::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        Δn₁s₁ = ξ.Δn₁s₁
        Δn₁s₂n₂s₁ = ξ.Δn₁s₂n₂s₁
        Δn₂s₂ = ξ.Δn₂s₂
        D₁₁ = - D*(Δn₁s₁+Δn₂s₂*ν)
        D₁₂ = - D*(1-ν)*Δn₁s₂n₂s₁
        D₂₂ = - D*(Δn₁s₁*ν+Δn₂s₂)
        N = ξ[:𝝭]
        B₁₁ = ξ[:∂²𝝭∂x²]
        B₁₂ = ξ[:∂²𝝭∂x∂y]
        B₂₂ = ξ[:∂²𝝭∂y²]
        B̄₁₁ = ξ[:∂²𝝭∂x²_]
        B̄₁₂ = ξ[:∂²𝝭∂x∂y_]
        B̄₂₂ = ξ[:∂²𝝭∂y²_]
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            ΔMₙₛᵢ = D₁₁*B₁₁[i] + D₁₂*B₁₂[i] + D₂₂*B₂₂[i]
            ΔM̄ₙₛᵢ = D₁₁*B̄₁₁[i] + D₁₂*B̄₁₂[i] + D₂₂*B̄₂₂[i]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                ΔMₙₛⱼ = D₁₁*B₁₁[j] + D₁₂*B₁₂[j] + D₂₂*B₂₂[j]
                k[I,J] += ΔMₙₛᵢ*N[j] + N[i]*ΔMₙₛⱼ - ΔM̄ₙₛᵢ*N[j]
            end
            f[I] += (ΔMₙₛᵢ - ΔM̄ₙₛᵢ)*g
        end
    end
end

"""
    L₂(ap::T)

Relative L₂ error of the deflection against an exact field,
`‖wʰ - w̄‖_L₂ / ‖w̄‖_L₂` with `wʰ = Σᵢ Nᵢ wᵢ`. Returns the two squared norms
`(Δw², w̄²)`; use [`L₂`](@ref)`(aps)` for the relative error over a set.

The exact field is supplied with

    prescribe!(aps, :w => f)

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.w`, `ξ.𝑤`.
"""
function L₂(ap::T) where T<:AbstractElement
    Δw² = 0.0
    w̄² = 0.0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        w̄ = ξ.w
        w = 0.0
        for (i,xᵢ) in enumerate(ap.𝓒)
            w += N[i]*xᵢ.d
        end
        Δw² += (w - w̄)^2*𝑤
        w̄²  += w̄^2*𝑤
    end
    return Δw², w̄²
end

"""
    L₂(aps::Vector{T})

Relative L₂ error of the deflection, `‖wʰ - w̄‖_L₂ / ‖w̄‖_L₂`, summed over a
set of elements. See the element-level method for the exact-field convention.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.w`, `ξ.𝑤`.
"""
function L₂(aps::Vector{T}) where T<:AbstractElement
    L₂Norm_Δw² = 0.0
    L₂Norm_w̄² = 0.0
    for ap in aps
        Δw², w̄² = L₂(ap)
        L₂Norm_Δw² += Δw²
        L₂Norm_w̄²  += w̄²
    end
    return (L₂Norm_Δw²/L₂Norm_w̄²)^0.5
end

"""
    H₁(ap::T)

Squared L₂ and H¹-seminorm errors of the deflection against an exact field,
returned as `(Δ∇w², ∇w̄², Δw², w̄²)`. Use [`H₁`](@ref)`(aps)` for the relative
H¹ and L₂ errors over a set.

The exact field is supplied with

    prescribe!(aps, :w => f, :∂w∂x => g₁, :∂w∂y => g₂)

# Element requirements
`set∇𝝭!` — reads `ξ[:𝝭]`, `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.w`, `ξ.∂w∂x`,
`ξ.∂w∂y`, `ξ.𝑤`.
"""
function H₁(ap::T) where T<:AbstractElement
    Δ∇w² = 0.0
    ∇w̄² = 0.0
    Δw² = 0.0
    w̄² = 0.0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        w̄ = ξ.w
        ∂w̄∂x = ξ.∂w∂x
        ∂w̄∂y = ξ.∂w∂y
        w = 0.0
        ∂w∂x = 0.0
        ∂w∂y = 0.0
        for (i,xᵢ) in enumerate(ap.𝓒)
            w += N[i]*xᵢ.d
            ∂w∂x += B₁[i]*xᵢ.d
            ∂w∂y += B₂[i]*xᵢ.d
        end
        Δw² += (w - w̄)^2*𝑤
        w̄²  += w̄^2*𝑤
        Δ∇w² += ((∂w∂x - ∂w̄∂x)^2 + (∂w∂y - ∂w̄∂y)^2)*𝑤
        ∇w̄² += (∂w̄∂x^2 + ∂w̄∂y^2)*𝑤
    end
    return Δ∇w², ∇w̄², Δw², w̄²
end

"""
    H₁(aps::Vector{T})

Relative H¹ and L₂ errors of the deflection over a set of elements,
returned as `(H₁, L₂)`. See the element-level method for the field convention.

# Element requirements
`set∇𝝭!` — reads `ξ[:𝝭]`, `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.w`, `ξ.∂w∂x`,
`ξ.∂w∂y`, `ξ.𝑤`.
"""
function H₁(aps::Vector{T}) where T<:AbstractElement
    H₁Norm_Δw² = 0.0
    H₁Norm_w̄² = 0.0
    L₂Norm_Δw² = 0.0
    L₂Norm_w̄² = 0.0
    for ap in aps
        Δ∇w², ∇w̄², Δw², w̄² = H₁(ap)
        H₁Norm_Δw² += Δw² + Δ∇w²
        H₁Norm_w̄²  += w̄² + ∇w̄²
        L₂Norm_Δw² += Δw²
        L₂Norm_w̄²  += w̄²
    end
    return (H₁Norm_Δw²/H₁Norm_w̄²)^0.5, (L₂Norm_Δw²/L₂Norm_w̄²)^0.5
end

"""
    Hₑ(ap::T; D, ν)

Squared energy-norm error of the deflection against an exact field,
returned as `(ΔW², W̄²)`. Use [`Hₑ`](@ref)`(aps; D, ν)` for the relative
error over a set.

    ‖v‖_E² = ∫_Ω D[v,ₓₓ² + 2ν v,ₓₓ v,ᵧᵧ + v,ᵧᵧ² + 2(1-ν) v,ₓᵧ²] dΩ

The exact field is supplied with

    prescribe!(aps, :∂²w∂x² => f₁₁, :∂²w∂x∂y => f₁₂, :∂²w∂y² => f₂₂)

# Keyword arguments
- `D` : bending stiffness `E h³ / (12(1-ν²))`
- `ν` : Poisson ratio

# Element requirements
`set∇²𝝭!` — reads `ξ[:∂²𝝭∂x²]`, `ξ[:∂²𝝭∂x∂y]`, `ξ[:∂²𝝭∂y²]`,
`ξ.∂²w∂x²`, `ξ.∂²w∂x∂y`, `ξ.∂²w∂y²`, `ξ.𝑤`.
"""
function Hₑ(ap::T; D::Real, ν::Real) where T<:AbstractElement
    ΔW² = 0.0
    W̄² = 0.0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        B₁₁ = ξ[:∂²𝝭∂x²]
        B₁₂ = ξ[:∂²𝝭∂x∂y]
        B₂₂ = ξ[:∂²𝝭∂y²]
        w̄₁₁ = ξ.∂²w∂x²
        w̄₁₂ = ξ.∂²w∂x∂y
        w̄₂₂ = ξ.∂²w∂y²
        w₁₁ = 0.0
        w₁₂ = 0.0
        w₂₂ = 0.0
        for (i,xᵢ) in enumerate(ap.𝓒)
            w₁₁ += B₁₁[i]*xᵢ.d
            w₁₂ += B₁₂[i]*xᵢ.d
            w₂₂ += B₂₂[i]*xᵢ.d
        end
        Δw₁₁ = w₁₁ - w̄₁₁
        Δw₁₂ = w₁₂ - w̄₁₂
        Δw₂₂ = w₂₂ - w̄₂₂
        ΔW² += D*(Δw₁₁^2 + 2*ν*Δw₁₁*Δw₂₂ + Δw₂₂^2 + 2*(1-ν)*Δw₁₂^2)*𝑤
        W̄² += D*(w̄₁₁^2 + 2*ν*w̄₁₁*w̄₂₂ + w̄₂₂^2 + 2*(1-ν)*w̄₁₂^2)*𝑤
    end
    return ΔW², W̄²
end

"""
    Hₑ(aps::Vector{T}; D, ν)

Relative energy-norm error of the deflection, `‖e‖_E / ‖w̄‖_E`, summed over a
set of elements. See the element-level method for the formula and the
exact-field convention.

# Keyword arguments
- `D` : bending stiffness `E h³ / (12(1-ν²))`
- `ν` : Poisson ratio

# Element requirements
`set∇²𝝭!` — reads the second derivatives, `ξ.∂²w∂x²`, `ξ.∂²w∂x∂y`,
`ξ.∂²w∂y²`, `ξ.𝑤`.
"""
function Hₑ(aps::Vector{T}; D::Real, ν::Real) where T<:AbstractElement
    HₑNorm_ΔW² = 0.0
    HₑNorm_W̄² = 0.0
    for ap in aps
        ΔW², W̄² = Hₑ(ap; D = D, ν = ν)
        HₑNorm_ΔW² += ΔW²
        HₑNorm_W̄²  += W̄²
    end
    return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5
end

# ---------------------------------------------------------------------------
# Mixed Hellinger–Reissner blocks for the Kirchhoff plate
#
# Unknowns
#   w  : RK scalar, one DOF per node,        index I
#   M  : piecewise polynomial, three DOFs per M node:
#          M11 -> 3I-2,  M22 -> 3I-1,  M12 -> 3I
#
# The operators below fill *block* matrices with local DOF numbering, as in
# the mixed Mindlin plate module.  The caller assembles the global saddle-point
# matrix with offsets; no Schur complement is formed.
#
# The domain coupling uses the second-derivative form of the manuscript:
#
#     ∫Ω w δMᵢⱼ,ᵢⱼ dΩ
#
# For a piecewise-linear M field this term is zero; it becomes active for
# quadratic M (q = 4 in the manuscript).  The boundary terms use only first
# derivatives of M.
# ---------------------------------------------------------------------------

"""
    ∫δMC⁻¹MdΩ(ap, k; D, ν)

Mixed Hellinger–Reissner moment compliance block,

    ∫Ω δMᵢⱼ C⁻¹ᵢⱼₖₗ Mₖₗ dΩ

with the Kirchhoff plate compliance

    C⁻¹₁₁₁₁ = 1/(D(1-ν²)),   C⁻¹₁₁₂₂ = -ν/(D(1-ν²)),
    C⁻¹₁₂₁₂ = 2/(D(1-ν)).

The M element uses three DOFs per node: `3I-2` (M11), `3I-1` (M22),
`3I` (M12).

# Keyword arguments
- `D` : bending stiffness
- `ν` : Poisson ratio

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.𝑤`.
"""
function ∫MC⁻¹MdΩ(ap::T, k::AbstractMatrix{Float64};
                   D::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    C₁₁ = 1.0/(D*(1.0 - ν^2))
    C₁₂ = -ν/(D*(1.0 - ν^2))
    C₃₃ = 2.0/(D*(1.0 - ν))
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3I-2,3J-2] -= N[i]*C₁₁*N[j]*𝑤
                k[3I-2,3J-1] -= N[i]*C₁₂*N[j]*𝑤
                k[3I-1,3J-2] -= N[i]*C₁₂*N[j]*𝑤
                k[3I-1,3J-1] -= N[i]*C₁₁*N[j]*𝑤
                k[3I,3J]     -= N[i]*C₃₃*N[j]*𝑤
            end
        end
    end
    return k
end

"""
    ∫∇²MwdΩ(aᵐ, aʷ, k)

Row `M`, column `w` domain coupling,

    -∫Ω Mᵢⱼ,ᵢⱼ w dΩ

The moment field is the independent piecewise field: its second derivatives
enter the operator, while `w` enters through its shape-function values.

`k` is the independently numbered `M`-`w` block of size `(3 n_M, n_w)`.
The transpose block is `k'`.

# Element requirements
`set∇²𝝭!` for `aᵐ`; `set𝝭!` for `aʷ`.
"""
function ∫∇²MwdΩ(aᵐ::T, aʷ::S,
                  k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    for (ξᵐ,ξʷ) in zip(aᵐ.𝓖, aʷ.𝓖)
        B₁₁ = ξᵐ[:∂²𝝭∂x²]
        B₁₂ = ξᵐ[:∂²𝝭∂x∂y]
        B₂₂ = ξᵐ[:∂²𝝭∂y²]
        Nʷ = ξʷ[:𝝭]
        𝑤 = ξʷ.𝑤
        for (i,xᵐ) in enumerate(aᵐ.𝓒)
            I = xᵐ.𝐼
            for (j,xʷ) in enumerate(aʷ.𝓒)
                J = xʷ.𝐼
                k[3I-2,J] -= B₁₁[i]*Nʷ[j]*𝑤
                k[3I-1,J] -= B₂₂[i]*Nʷ[j]*𝑤
                k[3I,J]   -= 2.0*B₁₂[i]*Nʷ[j]*𝑤
            end
        end
    end
    return k
end

"""
    ∫VwdΓ(aᵐ, aʷ, k)

Row `M`, column `w` boundary term on `Γw`,

    -∫Γ V(M) δw dΓ

This `k`-only overload serves the cell boundary `Γ_C`; the `(aᵐ, aʷ, k, f)`
overload serves the essential boundary `Γ_w`.

with `V(M) = Mᵢⱼ,ⱼ nᵢ + Mₙₛ,ₛ`.  The essential-boundary field `w` is the last
argument and occupies the column index `J` of `k`.

`k` has size `(3 n_M, n_w)`; the transpose block is `k'`.

# Element requirements
`set∇𝝭!` for `aᵐ`; `set𝝭!` for `aʷ`; boundary `ξ.n₁`, `ξ.n₂`.
"""
function ∫VwdΓ(aᵐ::T, aʷ::S,
                k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    for (ξᵐ,ξʷ) in zip(aᵐ.𝓖, aʷ.𝓖)
        Nʷ = ξʷ[:𝝭]
        B₁ = ξᵐ[:∂𝝭∂x]
        B₂ = ξᵐ[:∂𝝭∂y]
        n₁ = ξᵐ.n₁; n₂ = ξᵐ.n₂
        s₁ = -n₂; s₂ = n₁
        𝑤 = ξʷ.𝑤
        for (i,xᵐ) in enumerate(aᵐ.𝓒)
            I = xᵐ.𝐼
            V₁₁ = B₁[i]*n₁ + B₁[i]*n₁*s₁*s₁ + B₂[i]*n₁*s₁*s₂
            V₂₂ = B₂[i]*n₂ + B₁[i]*n₂*s₂*s₁ + B₂[i]*n₂*s₂*s₂
            V₁₂ = B₁[i]*n₂ + B₂[i]*n₁ + B₁[i]*(n₁*s₂+n₂*s₁)*s₁ + B₂[i]*(n₁*s₂+n₂*s₁)*s₂
            for (j,xʷ) in enumerate(aʷ.𝓒)
                J = xʷ.𝐼
                k[3I-2,J] += V₁₁*Nʷ[j]*𝑤
                k[3I-1,J] += V₂₂*Nʷ[j]*𝑤
                k[3I,J]   += V₁₂*Nʷ[j]*𝑤
            end
        end
    end
    return k
end

"""
    ∫VwdΓ(aᵐ, aʷ, k, f)

Essential-boundary `Γ_w` overload, assembling the matrix contribution

    +∫Γw V(M) δw dΓ

and the right-hand side

    +∫Γw V(δM) w̄ dΓ

The corner/essential field `w` is the last argument and occupies column `J`.
`k` has size `(3 n_M, n_w)`; `f` has length `3 n_M`.
"""
function ∫VwdΓ(aᵐ::T, aʷ::S,
                k::AbstractMatrix{Float64},
                f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    for (ξᵐ,ξʷ) in zip(aᵐ.𝓖, aʷ.𝓖)
        Nʷ = ξʷ[:𝝭]
        B₁ = ξᵐ[:∂𝝭∂x]
        B₂ = ξᵐ[:∂𝝭∂y]
        n₁ = ξᵐ.n₁; n₂ = ξᵐ.n₂
        s₁ = -n₂; s₂ = n₁
        w̄ = ξʷ.w
        𝑤 = ξʷ.𝑤
        for (i,xᵐ) in enumerate(aᵐ.𝓒)
            I = xᵐ.𝐼
            V₁₁ = B₁[i]*n₁ + B₁[i]*n₁*s₁*s₁ + B₂[i]*n₁*s₁*s₂
            V₂₂ = B₂[i]*n₂ + B₁[i]*n₂*s₂*s₁ + B₂[i]*n₂*s₂*s₂
            V₁₂ = B₁[i]*n₂ + B₂[i]*n₁ + B₁[i]*(n₁*s₂+n₂*s₁)*s₁ + B₂[i]*(n₁*s₂+n₂*s₁)*s₂
            for (j,xʷ) in enumerate(aʷ.𝓒)
                J = xʷ.𝐼
                k[3I-2,J] -= V₁₁*Nʷ[j]*𝑤
                k[3I-1,J] -= V₂₂*Nʷ[j]*𝑤
                k[3I,J]   -= V₁₂*Nʷ[j]*𝑤
            end
            f[3I-2] -= V₁₁*w̄*𝑤
            f[3I-1] -= V₂₂*w̄*𝑤
            f[3I]   -= V₁₂*w̄*𝑤
        end
    end
    return k,f
end

"""
    ∫θwMₙₙδMdΓ(aᵐ, aʷ, k)

Row `M`, column `w` boundary term on `Γθ`,

    +∫Γ θ(w) Mₙₙ(δM) dΓ

This `k`-only overload serves the cell boundary `Γ_C`; the
`(aᵐ, aʷ, k, f)` overload serves the essential boundary `Γ_θ`.

with `θ(w) = w,ᵢ nᵢ` and `Mₙₙ = M₁₁n₁² + M₂₂n₂² + 2M₁₂n₁n₂`.
The essential-boundary field `w` is the last argument and occupies the column
index `J` of `k`.

`k` has size `(3 n_M, n_w)`; the transpose block is `k'`.

# Element requirements
`set𝝭!` for `aᵐ`; `set∇𝝭!` for `aʷ`; boundary `ξ.n₁`, `ξ.n₂`.
"""
function ∫MθdΓ(aᵐ::T, aʷ::S,
                    k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    for (ξᵐ,ξʷ) in zip(aᵐ.𝓖, aʷ.𝓖)
        Nᵐ = ξᵐ[:𝝭]
        Bʷ₁ = ξʷ[:∂𝝭∂x]
        Bʷ₂ = ξʷ[:∂𝝭∂y]
        n₁ = ξᵐ.n₁; n₂ = ξᵐ.n₂
        𝑤 = ξʷ.𝑤
        for (i,xᵐ) in enumerate(aᵐ.𝓒)
            I = xᵐ.𝐼
            for (j,xʷ) in enumerate(aʷ.𝓒)
                J = xʷ.𝐼
                θ = Bʷ₁[j]*n₁ + Bʷ₂[j]*n₂
                k[3I-2,J] -= θ*Nᵐ[i]*n₁*n₁*𝑤
                k[3I-1,J] -= θ*Nᵐ[i]*n₂*n₂*𝑤
                k[3I,J]   -= 2.0*θ*Nᵐ[i]*n₁*n₂*𝑤
            end
        end
    end
    return k
end

"""
    ∫MθdΓ(aᵐ, aʷ, k, f)

Essential-boundary `Γ_θ` overload, assembling the matrix contribution

    -∫Γθ θ(w) Mₙₙ(δM) dΓ

and the right-hand side

    -∫Γθ θ̄ Mₙₙ(δM) dΓ

`k` has size `(3 n_M, n_w)`; `f` has length `3 n_M`.
"""
function ∫MθdΓ(aᵐ::T, aʷ::S,
                    k::AbstractMatrix{Float64},
                    f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    for (ξᵐ,ξʷ) in zip(aᵐ.𝓖, aʷ.𝓖)
        Nᵐ = ξᵐ[:𝝭]
        Bʷ₁ = ξʷ[:∂𝝭∂x]
        Bʷ₂ = ξʷ[:∂𝝭∂y]
        n₁ = ξᵐ.n₁; n₂ = ξᵐ.n₂
        θ̄ = ξʷ.θ
        𝑤 = ξʷ.𝑤
        for (i,xᵐ) in enumerate(aᵐ.𝓒)
            I = xᵐ.𝐼
            for (j,xʷ) in enumerate(aʷ.𝓒)
                J = xʷ.𝐼
                θ = Bʷ₁[j]*n₁ + Bʷ₂[j]*n₂
                k[3I-2,J] += θ*Nᵐ[i]*n₁*n₁*𝑤
                k[3I-1,J] += θ*Nᵐ[i]*n₂*n₂*𝑤
                k[3I,J]   += 2.0*θ*Nᵐ[i]*n₁*n₂*𝑤
            end
            f[3I-2] += θ̄*Nᵐ[i]*n₁*n₁*𝑤
            f[3I-1] += θ̄*Nᵐ[i]*n₂*n₂*𝑤
            f[3I]   += 2.0*θ̄*Nᵐ[i]*n₁*n₂*𝑤
        end
    end
    return k,f
end

"""
    wP(aᵐ, aʷ, k)

Row `M`, column `w` corner term,

    +∫ w [Mₙₛ(δM)] dΓ

This `k`-only overload serves the cell boundary `Γ_C`; the
`(aᵐ, aʷ, k, f)` overload serves the essential boundary `c_w`.

The corner contribution is a **point contribution**, so no quadrature weight
is used.  It is assembled like [`ΔMₙₛg`](@ref), with the prescribed corner
jump coefficients

    ξ.Δn₁s₁,   ξ.Δn₁s₂n₂s₁,   ξ.Δn₂s₂

instead of the smooth-boundary normals and tangents:

    [Mₙₛ(δM)] = Δn₁s₁ M₁₁(δM) + Δn₂s₂ M₂₂(δM)
              + Δn₁s₂n₂s₁ M₁₂(δM)

The essential-boundary field `w` is the last argument and occupies the column
index `J` of `k`.

`k` has size `(3 n_M, n_w)`; the transpose block is `k'`.

# Element requirements
`set𝝭!` for both `aᵐ` and `aʷ`; corner fields `ξ.Δn₁s₁`,
`ξ.Δn₁s₂n₂s₁`, `ξ.Δn₂s₂`.
"""
function wP(aᵐ::T, aʷ::S,
            k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    for (ξᵐ,ξʷ) in zip(aᵐ.𝓖, aʷ.𝓖)
        data = getfield(ξᵐ,:data)
        haskey(data, :Δn₁s₁) || continue
        Δn₁s₁ = ξᵐ.Δn₁s₁
        Δn₁s₂n₂s₁ = ξᵐ.Δn₁s₂n₂s₁
        Δn₂s₂ = ξᵐ.Δn₂s₂
        Nᵐ = ξᵐ[:𝝭]
        Nʷ = ξʷ[:𝝭]
        for (i,xᵐ) in enumerate(aᵐ.𝓒)
            I = xᵐ.𝐼
            for (j,xʷ) in enumerate(aʷ.𝓒)
                J = xʷ.𝐼
                k[3I-2,J] -= Nʷ[j]*Δn₁s₁*Nᵐ[i]
                k[3I-1,J] -= Nʷ[j]*Δn₂s₂*Nᵐ[i]
                k[3I,J]   -= Nʷ[j]*Δn₁s₂n₂s₁*Nᵐ[i]
            end
        end
    end
    return k
end

"""
    wP(aᵐ, aʷ, k, f)

Same corner term as [`wP`](@ref)`(aᵐ, aʷ, k)`, with the prescribed boundary
deflection `ξʷ.w` added to the right-hand side.

The corner contribution is a point contribution, so no quadrature weight is
used.

`k` has size `(3 n_M, n_w)`; `f` has length `3 n_M`.

# Element requirements
`set𝝭!` for both `aᵐ` and `aʷ`; corner fields `ξ.Δn₁s₁`,
`ξ.Δn₁s₂n₂s₁`, `ξ.Δn₂s₂`; prescribed `ξʷ.w`.
"""
function wP(aᵐ::T, aʷ::S,
            k::AbstractMatrix{Float64},
            f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    for (ξᵐ,ξʷ) in zip(aᵐ.𝓖, aʷ.𝓖)
        data = getfield(ξᵐ,:data)
        haskey(data, :Δn₁s₁) || continue
        Δn₁s₁ = ξᵐ.Δn₁s₁
        Δn₁s₂n₂s₁ = ξᵐ.Δn₁s₂n₂s₁
        Δn₂s₂ = ξᵐ.Δn₂s₂
        Nᵐ = ξᵐ[:𝝭]
        Nʷ = ξʷ[:𝝭]
        w̄ = ξʷ.w
        for (i,xᵐ) in enumerate(aᵐ.𝓒)
            I = xᵐ.𝐼
            for (j,xʷ) in enumerate(aʷ.𝓒)
                J = xʷ.𝐼
                k[3I-2,J] += Nʷ[j]*Δn₁s₁*Nᵐ[i]
                k[3I-1,J] += Nʷ[j]*Δn₂s₂*Nᵐ[i]
                k[3I,J]   += Nʷ[j]*Δn₁s₂n₂s₁*Nᵐ[i]
            end
            f[3I-2] += w̄*Δn₁s₁*Nᵐ[i]
            f[3I-1] += w̄*Δn₂s₂*Nᵐ[i]
            f[3I]   += w̄*Δn₁s₂n₂s₁*Nᵐ[i]
        end
    end
    return k, f
end


end
