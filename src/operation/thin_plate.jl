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
    n₁,n₂ = get𝒏(ap)
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        M₁ = ξ.M₁
        M₂ = ξ.M₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] -= (B₁[i]*M₁+B₂[i]*M₂)*𝑤
        end
    end
end

"""
    ∫∇𝑛vθdΓ(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64})

Boundary penalty on the normal derivative, `α ∫ (∂ₙw) θ dΓ`.

# Element requirements
`set∇𝝭!` — reads gradients, `ξ.n₁`, `ξ.n₂`, `ξ.θ`, `ξ.𝑤`.
"""
function ∫∇𝑛vθdΓ(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
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
    ∫MₙₙθdΓ(ap, k, f; D, ν)

Kirchhoff plate boundary term `∫ (Mₙₙ(w) + α θ) θ̄ dΓ` — normal bending moment
plus penalty, as used in the weak (Nitsche/penalty) free-edge condition.

    Mₙₙ = -D[(n₁² + νn₂²) w,ₓₓ + 2(1-ν)n₁n₂ w,ₓᵧ + (νn₁² + n₂²) w,ᵧᵧ]

# Keyword arguments
- `D` : bending stiffness
- `ν` : Poisson ratio

# Element requirements
`set∇²𝝭!` — reads the second derivatives, `ξ.n₁`, `ξ.n₂`, `ξ.θ`, `ξ.𝑤`.
"""
function ∫MₙₙθdΓ(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}; ν::Real, D::Real) where T<:AbstractElement
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
    ∫VgdΓ(ap, k, f; D, ν)

Kirchhoff plate Kirchhoff-shear boundary term, `∫ (Vₙ(w) + α w) ḡ dΓ`.

    Vₙ = -D[∂ₙ∇²w + (1-ν)∂ₛ(w,ₙₛ)]

The effective shear force is the combination of the twisting-moment derivative
and the normal derivative of the Laplacian; the `α` term is the penalty that
closes the weak form.

# Keyword arguments
- `D` : bending stiffness
- `ν` : Poisson ratio

# Element requirements
`set∇²𝝭!` — reads the second derivatives, `ξ.n₁`, `ξ.n₂`, `ξ.s₁`, `ξ.s₂`,
`ξ.g`, `ξ.𝑤`.
"""
function ∫VgdΓ(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}; ν::Real, D::Real) where T<:AbstractElement
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
    ΔMₙₛg(ap, k, f; D, ν)
    ΔM̃ₙₛg(ap, k, f; D, ν)

Corner-force terms of the Kirchhoff plate, `∫ ΔMₙₛ ḡ`.

At a boundary corner the twisting moment `Mₙₛ` jumps; the jump `ΔMₙₛ` acts as a
concentrated force. `ξ.Δ` carries the corner sign (+1, -1 or 0 at a smooth
point), which is why the quadrature rule includes the boundary endpoints.

# Keyword arguments
- `D` : bending stiffness
- `ν` : Poisson ratio

# Element requirements
`set∇²𝝭!` — reads the second derivatives, `ξ.n₁`, `ξ.n₂`, `ξ.s₁`, `ξ.s₂`,
`ξ.Δ`, `ξ.g`, `ξ.𝑤`.
"""
function ΔMₙₛg(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}; ν::Real, D::Real) where T<:AbstractElement
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

end