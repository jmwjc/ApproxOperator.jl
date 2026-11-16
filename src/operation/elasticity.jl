module Elasticity
    
using ..ApproxOperator: AbstractElement

"""
    ∫εᵢⱼσᵢⱼdΩ(ap, k; E, ν)

Three-dimensional linear-elastic stiffness, `∫ ε(u):σ(v) dΩ`.

    a(u,v) = ∫_Ω ε(u):C:ε(v) dΩ,        C = λ I⊗I + 2μ I

expressed through the Lamé constants

    λ = Eν / ((1+ν)(1-2ν)),    μ = E / (2(1+ν))

The 3×3 block is accumulated at `k[3I-a, 3J-b]`, `a,b ∈ {1,2,3}`.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ[:∂𝝭∂z]`, `ξ.𝑤`.
"""
function ∫εᵢⱼσᵢⱼdΩ(ap::T,k::AbstractMatrix{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z]
        𝑤 = ξ.𝑤
        Cᵢᵢᵢᵢ = E*(1-ν)/(1+ν)/(1-2ν)
        Cᵢᵢⱼⱼ = E*ν/(1+ν)/(1-2ν)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += (Cᵢᵢᵢᵢ*B₁[i]*B₁[j] + Cᵢⱼᵢⱼ*B₂[i]*B₂[j] + Cᵢⱼᵢⱼ*B₃[i]*B₃[j])*𝑤
                k[3*I-2,3*J-1] += (Cᵢᵢⱼⱼ*B₁[i]*B₂[j] + Cᵢⱼᵢⱼ*B₂[i]*B₁[j])*𝑤
                k[3*I-2,3*J]   += (Cᵢᵢⱼⱼ*B₁[i]*B₃[j] + Cᵢⱼᵢⱼ*B₃[i]*B₁[j])*𝑤
                k[3*I-1,3*J-2] += (Cᵢᵢⱼⱼ*B₂[i]*B₁[j] + Cᵢⱼᵢⱼ*B₁[i]*B₂[j])*𝑤
                k[3*I-1,3*J-1] += (Cᵢⱼᵢⱼ*B₁[i]*B₁[j] + Cᵢᵢᵢᵢ*B₂[i]*B₂[j] + Cᵢⱼᵢⱼ*B₃[i]*B₃[j])*𝑤
                k[3*I-1,3*J]   += (Cᵢᵢⱼⱼ*B₂[i]*B₃[j] + Cᵢⱼᵢⱼ*B₃[i]*B₂[j])*𝑤
                k[3*I,3*J-2]   += (Cᵢᵢⱼⱼ*B₃[i]*B₁[j] + Cᵢⱼᵢⱼ*B₁[i]*B₃[j])*𝑤
                k[3*I,3*J-1]   += (Cᵢᵢⱼⱼ*B₃[i]*B₂[j] + Cᵢⱼᵢⱼ*B₂[i]*B₃[j])*𝑤
                k[3*I,3*J]     += (Cᵢⱼᵢⱼ*B₁[i]*B₁[j] + Cᵢⱼᵢⱼ*B₂[i]*B₂[j] + Cᵢᵢᵢᵢ*B₃[i]*B₃[j])*𝑤
            end
        end
    end
end

"""
    ∫εᵈᵢⱼσᵈᵢⱼdΩ(ap, k; E, ν)

Three-dimensional **deviatoric** contribution to the stiffness,
`∫ εᵈ(u):σᵈ(v) dΩ`.

    a(u,v) = ∫_Ω 2μ εᵈ(u):εᵈ(v) dΩ

with the deviatoric projector `Pᵈ = I - ⅓ I⊗I` folded into the coefficients,
so only the shear modulus `μ = E/(2(1+ν))` appears.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ[:∂𝝭∂z]`, `ξ.𝑤`.
"""
function ∫εᵈᵢⱼσᵈᵢⱼdΩ(ap::T,k::AbstractMatrix{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z]
        𝑤 = ξ.𝑤
        Cᵈ = E/(1+ν)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += Cᵈ*( 2/3*B₁[i]*B₁[j]+1/2*B₂[i]*B₂[j]+1/2*B₃[i]*B₃[j])*𝑤
                k[3*I-2,3*J-1] += Cᵈ*(-1/3*B₁[i]*B₂[j]+1/2*B₂[i]*B₁[j])*𝑤
                k[3*I-2,3*J]   += Cᵈ*(-1/3*B₁[i]*B₃[j]+1/2*B₃[i]*B₁[j])*𝑤
                k[3*I-1,3*J-2] += Cᵈ*(-1/3*B₂[i]*B₁[j]+1/2*B₁[i]*B₂[j])*𝑤
                k[3*I-1,3*J-1] += Cᵈ*( 2/3*B₂[i]*B₂[j]+1/2*B₃[i]*B₃[j]+1/2*B₁[i]*B₁[j])*𝑤
                k[3*I-1,3*J]   += Cᵈ*(-1/3*B₂[i]*B₃[j]+1/2*B₃[i]*B₂[j])*𝑤
                k[3*I,3*J-2]   += Cᵈ*(-1/3*B₃[i]*B₁[j]+1/2*B₁[i]*B₃[j])*𝑤
                k[3*I,3*J-1]   += Cᵈ*(-1/3*B₃[i]*B₂[j]+1/2*B₂[i]*B₃[j])*𝑤
                k[3*I,3*J]     += Cᵈ*( 2/3*B₃[i]*B₃[j]+1/2*B₁[i]*B₁[j]+1/2*B₂[i]*B₂[j])*𝑤
            end
        end
    end
end

"""
    ∫∫ρvᵢuᵢdxdy(ap, k; ρ)

Consistent mass matrix, `∫ ρ vᵢ uᵢ dΩ`.

    m[2I-1, 2J-1] += ρ Nᵢ Nⱼ w
    m[2I,   2J  ] += ρ Nᵢ Nⱼ w

Two degrees of freedom per node, both sharing the same scalar integrand.

# Keyword arguments
- `ρ` : mass density

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.𝑤`.
"""
function ∫∫ρvᵢuᵢdxdy(ap::T,k::AbstractMatrix{Float64}; ρ::Real) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += ρ*N[i]*N[j]*𝑤
                k[2*I,2*J]     += ρ*N[i]*N[j]*𝑤
            end
        end
    end
end

"""
    ∫∫εᵢⱼσᵢⱼdxdy(ap, k; E, ν)                         # plane stress, one element set
    ∫∫εᵢⱼσᵢⱼdxdy(aᵤ, aₛ, k; E, ν)                      # test/trial on different sets
    ∫∫εᵢⱼσᵢⱼdxdy_PlaneStrian(aᵤ, aₛ, k; E, ν)           # plane **strain**

Two-dimensional linear-elastic stiffness, `∫ ε(u):σ(v) dΩ`.

Plane stress uses the reduced moduli

    C₁₁₁₁ = E/(1-ν²),   C₁₁₂₂ = Eν/(1-ν²),   C₁₂₁₂ = E/(2(1+ν))

The two-element-set overloads are for mixed or Petrov–Galerkin formulations,
where test and trial functions live on different approximations; the quadrature
points of the two sets are zipped.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫∫εᵢⱼσᵢⱼdxdy(ap::T,k::AbstractMatrix{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        Cᵢᵢᵢᵢ = E/(1-ν^2)
        Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += (Cᵢᵢᵢᵢ*B₁[i]*B₁[j] + Cᵢⱼᵢⱼ*B₂[i]*B₂[j])*𝑤
                k[2*I-1,2*J]   += (Cᵢᵢⱼⱼ*B₁[i]*B₂[j] + Cᵢⱼᵢⱼ*B₂[i]*B₁[j])*𝑤
                k[2*I,2*J-1]   += (Cᵢᵢⱼⱼ*B₂[i]*B₁[j] + Cᵢⱼᵢⱼ*B₁[i]*B₂[j])*𝑤
                k[2*I,2*J]     += (Cᵢᵢᵢᵢ*B₂[i]*B₂[j] + Cᵢⱼᵢⱼ*B₁[i]*B₁[j])*𝑤
            end
        end
    end
end

"""
    ∫∫εᵢⱼσᵢⱼdxdy(aᵤ::T,aₛ::S,k::AbstractMatrix{Float64}; E, ν)

Element-set overload of `∫∫εᵢⱼσᵢⱼdxdy`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function ∫∫εᵢⱼσᵢⱼdxdy(aᵤ::T,aₛ::S,k::AbstractMatrix{Float64}; E::Real, ν::Real) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    for (ξᵤ,ξₛ) in zip(𝓖ᵤ,𝓖ₛ)
        Bᵤ₁ = ξᵤ[:∂𝝭∂x]
        Bᵤ₂ = ξᵤ[:∂𝝭∂y]
        Bₛ₁ = ξₛ[:∂𝝭∂x]
        Bₛ₂ = ξₛ[:∂𝝭∂y]
        𝑤 = ξᵤ.𝑤
        Cᵢᵢᵢᵢ = E/(1-ν^2)
        Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += (Cᵢᵢᵢᵢ*Bₛ₁[i]*Bᵤ₁[j]+Cᵢⱼᵢⱼ*Bₛ₂[i]*Bᵤ₂[j])*𝑤
                k[2*I-1,2*J]   += (Cᵢᵢⱼⱼ*Bₛ₁[i]*Bᵤ₂[j]+Cᵢⱼᵢⱼ*Bₛ₂[i]*Bᵤ₁[j])*𝑤
                k[2*I,2*J-1]   += (Cᵢᵢⱼⱼ*Bₛ₂[i]*Bᵤ₁[j]+Cᵢⱼᵢⱼ*Bₛ₁[i]*Bᵤ₂[j])*𝑤
                k[2*I,2*J]     += (Cᵢᵢᵢᵢ*Bₛ₂[i]*Bᵤ₂[j]+Cᵢⱼᵢⱼ*Bₛ₁[i]*Bᵤ₁[j])*𝑤
            end
        end
    end
end

"""
    ∫∫εᵢⱼσᵢⱼdxdy_PlaneStrian(aᵤ::T,aₛ::S,k::AbstractMatrix{Float64}; E, ν)

Element-set overload of `∫∫εᵢⱼσᵢⱼdxdy_PlaneStrian`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function ∫∫εᵢⱼσᵢⱼdxdy_PlaneStrian(aᵤ::T,aₛ::S,k::AbstractMatrix{Float64}; E::Real, ν::Real) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    for (ξᵤ,ξₛ) in zip(𝓖ᵤ,𝓖ₛ)
        Bᵤ₁ = ξᵤ[:∂𝝭∂x]
        Bᵤ₂ = ξᵤ[:∂𝝭∂y]
        Bₛ₁ = ξₛ[:∂𝝭∂x]
        Bₛ₂ = ξₛ[:∂𝝭∂y]
        𝑤 = ξᵤ.𝑤
        Cᵢᵢᵢᵢ = E*(1-ν)/(1-2*ν)/(1+ν)
        Cᵢᵢⱼⱼ = E*ν/(1-2*ν)/(1+ν)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += (Cᵢᵢᵢᵢ*Bₛ₁[i]*Bᵤ₁[j]+Cᵢⱼᵢⱼ*Bₛ₂[i]*Bᵤ₂[j])*𝑤
                k[2*I-1,2*J]   += (Cᵢᵢⱼⱼ*Bₛ₁[i]*Bᵤ₂[j]+Cᵢⱼᵢⱼ*Bₛ₂[i]*Bᵤ₁[j])*𝑤
                k[2*I,2*J-1]   += (Cᵢᵢⱼⱼ*Bₛ₂[i]*Bᵤ₁[j]+Cᵢⱼᵢⱼ*Bₛ₁[i]*Bᵤ₂[j])*𝑤
                k[2*I,2*J]     += (Cᵢᵢᵢᵢ*Bₛ₂[i]*Bᵤ₂[j]+Cᵢⱼᵢⱼ*Bₛ₁[i]*Bᵤ₁[j])*𝑤
            end
        end
    end
end

"""
    ∫∫εᵛᵢⱼσᵛᵢⱼdxdy(ap, k; E, ν)

Volumetric (dilatational) contribution to the 2-D stiffness,
`∫ εᵛ(u):σᵛ(v) dΩ`.

    a(u,v) = ∫_Ω (E/(1-2ν)) (∇·u)(∇·v) dΩ   up to the 2-D/3-D factor

Only the divergence–divergence term is kept; combined with
`∫∫εᵈᵢⱼσᵈᵢⱼdxdy` it reproduces the full elastic bilinear form.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫∫εᵛᵢⱼσᵛᵢⱼdxdy(ap::T,k::AbstractMatrix{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        Cᵛ = E/(1-2*ν)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += Cᵛ/3*B₁[i]*B₁[j]*𝑤
                k[2*I-1,2*J]   += Cᵛ/3*B₁[i]*B₂[j]*𝑤
                k[2*I,2*J-1]   += Cᵛ/3*B₂[i]*B₁[j]*𝑤
                k[2*I,2*J]     += Cᵛ/3*B₂[i]*B₂[j]*𝑤
            end
        end
    end
end

"""
    ∫∫εᵈᵢⱼσᵈᵢⱼdxdy(ap, k; E, ν)                       # single element set
    ∫∫εᵈᵢⱼσᵈᵢⱼdxdy(aₛ, aᵤ, k; E, ν)                    # test/trial on different sets

Deviatoric contribution to the 2-D stiffness, `∫ εᵈ(u):σᵈ(v) dΩ`.

    Pᵈ = I - ⅓ I⊗I,      a(u,v) = ∫_Ω (E/(1+ν)) εᵈ(u):εᵈ(v) dΩ

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫∫εᵈᵢⱼσᵈᵢⱼdxdy(ap::T,k::AbstractMatrix{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        Cᵈ = E/(1+ν)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += Cᵈ*( 2/3*B₁[i]*B₁[j]+1/2*B₂[i]*B₂[j])*𝑤
                k[2*I-1,2*J]   += Cᵈ*(-1/3*B₁[i]*B₂[j]+1/2*B₂[i]*B₁[j])*𝑤
                k[2*I,2*J-1]   += Cᵈ*(-1/3*B₂[i]*B₁[j]+1/2*B₁[i]*B₂[j])*𝑤
                k[2*I,2*J]     += Cᵈ*( 2/3*B₂[i]*B₂[j]+1/2*B₁[i]*B₁[j])*𝑤
            end
        end
    end
end

"""
    ∫∫εᵈᵢⱼσᵈᵢⱼdxdy(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}; E, ν)

Element-set overload of `∫∫εᵈᵢⱼσᵈᵢⱼdxdy`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function ∫∫εᵈᵢⱼσᵈᵢⱼdxdy(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}; E::Real, ν::Real) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    for (ξᵤ,ξₛ) in zip(𝓖ᵤ,𝓖ₛ)
        Bᵤ₁ = ξᵤ[:∂𝝭∂x]
        Bᵤ₂ = ξᵤ[:∂𝝭∂y]
        Bₛ₁ = ξₛ[:∂𝝭∂x]
        Bₛ₂ = ξₛ[:∂𝝭∂y]
        𝑤 = ξᵤ.𝑤
        Cᵈ = E/(1+ν)
        for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += Cᵈ*( 2/3*Bₛ₁[i]*Bᵤ₁[j]+1/2*Bₛ₂[i]*Bᵤ₂[j])*𝑤
                k[2*I-1,2*J]   += Cᵈ*(-1/3*Bₛ₁[i]*Bᵤ₂[j]+1/2*Bₛ₂[i]*Bᵤ₁[j])*𝑤
                k[2*I,2*J-1]   += Cᵈ*(-1/3*Bₛ₂[i]*Bᵤ₁[j]+1/2*Bₛ₁[i]*Bᵤ₂[j])*𝑤
                k[2*I,2*J]     += Cᵈ*( 2/3*Bₛ₂[i]*Bᵤ₂[j]+1/2*Bₛ₁[i]*Bᵤ₁[j])*𝑤
            end
        end
    end
end

"""
    ∫∫σᵢⱼσₖₗdxdy(ap, k; E, ν)                          # plane stress
    ∫∫σᵢⱼσₖₗdxdy_PlaneStrian(ap, k; E, ν)               # plane strain

Compliance form, `∫ σ(u):C⁻¹:σ(v) dΩ`.

The dual (Hellinger–Reissner / mixed) formulations need the inverse
constitutive law rather than the stiffness:

    plane stress:  C⁻¹₁₁₁₁ = 1/E,          C⁻¹₁₁₂₂ = -ν/E,      C⁻¹₁₂₁₂ = 2(1+ν)/E
    plane strain:  C⁻¹₁₁₁₁ = (1-ν²)/E,     C⁻¹₁₁₂₂ = -(ν+ν²)/E

Requires **one** shape function per stress component, so it uses
`ξ[:𝝭]` rather than a gradient.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.𝑤`.
"""
function ∫∫σᵢⱼσₖₗdxdy(ap::T,k::AbstractMatrix{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        C⁻¹ᵢᵢᵢᵢ = 1/E
        C⁻¹ᵢᵢⱼⱼ = -ν/E
        C⁻¹ᵢⱼᵢⱼ = 2*(1+ν)/E
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += N[i]*C⁻¹ᵢᵢᵢᵢ*N[j]*𝑤
                k[3*I-2,3*J-1] += N[i]*C⁻¹ᵢᵢⱼⱼ*N[j]*𝑤
                k[3*I-1,3*J-2] += N[i]*C⁻¹ᵢᵢⱼⱼ*N[j]*𝑤
                k[3*I-1,3*J-1] += N[i]*C⁻¹ᵢᵢᵢᵢ*N[j]*𝑤
                k[3*I,3*J]     += N[i]*C⁻¹ᵢⱼᵢⱼ*N[j]*𝑤
            end
        end
    end
end

"""
    ∫∫σᵢⱼσₖₗdxdy_PlaneStrian(ap::T,k::AbstractMatrix{Float64}; E, ν)

Element-set overload of `∫∫σᵢⱼσₖₗdxdy_PlaneStrian`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function ∫∫σᵢⱼσₖₗdxdy_PlaneStrian(ap::T,k::AbstractMatrix{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        C⁻¹ᵢᵢᵢᵢ = (1-ν^2)/E
        C⁻¹ᵢᵢⱼⱼ = -(ν+ν^2)/E
        C⁻¹ᵢⱼᵢⱼ = 2*(1+ν)/E
        
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += N[i]*C⁻¹ᵢᵢᵢᵢ*N[j]*𝑤
                k[3*I-2,3*J-1] += N[i]*C⁻¹ᵢᵢⱼⱼ*N[j]*𝑤
                k[3*I-1,3*J-2] += N[i]*C⁻¹ᵢᵢⱼⱼ*N[j]*𝑤
                k[3*I-1,3*J-1] += N[i]*C⁻¹ᵢᵢᵢᵢ*N[j]*𝑤
                k[3*I,3*J]     += N[i]*C⁻¹ᵢⱼᵢⱼ*N[j]*𝑤
            end
        end
    end
end

# function ∫∫σᵢⱼσₖₗdxdy_PlaneStrian(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
#     𝓒 = ap.𝓒;𝓖 = ap.𝓖
#     for ξ in 𝓖
#         N = ξ[:𝝭]
#         𝑤 = ξ.𝑤
#         E = ξ.E
#         ν = ξ.ν
#         C⁻¹ᵢᵢᵢᵢ = (1+ν)*(1-2*ν)/E/(1-ν)
#         C⁻¹ᵢᵢⱼⱼ = (1+ν)*(1-2*ν)/E/ν
#         C⁻¹ᵢⱼᵢⱼ = 2*(1+ν)/E
#         for (i,xᵢ) in enumerate(𝓒)
#             I = xᵢ.𝐼
#             for (j,xⱼ) in enumerate(𝓒)
#                 J = xⱼ.𝐼
#                 k[4*I-3,4*J-3] += N[i]*C⁻¹ᵢᵢᵢᵢ*N[j]*𝑤
#                 k[4*I-3,4*J-2] += N[i]*C⁻¹ᵢᵢⱼⱼ*N[j]*𝑤
#                 k[4*I-2,4*J-3] += N[i]*C⁻¹ᵢᵢⱼⱼ*N[j]*𝑤
#                 k[4*I-2,4*J-2] += N[i]*C⁻¹ᵢᵢᵢᵢ*N[j]*𝑤
               
#                 k[4*I-1,4*J-3] += N[i]*C⁻¹ᵢᵢⱼⱼ*N[j]*𝑤
#                 k[4*I-1,4*J-2] += N[i]*C⁻¹ᵢᵢⱼⱼ*N[j]*𝑤

#                 k[4*I,4*J]     += N[i]*C⁻¹ᵢⱼᵢⱼ*N[j]*𝑤
                
#             end
#         end
#     end
# end

"""
    ∫σᵢⱼnⱼgᵢds(ap, k, f; E, ν)

Traction boundary term, `∫ (σ(u)·n)·v ds`, for a prescribed traction `g`.

    f[I] += Nᵢ (σ(u)·n)·g w

`g` is supplied with `prescribe!(el, :g => f)`.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set∇𝝭!` and `set𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ[:𝝭]`, `ξ.n₁`, `ξ.n₂`, `ξ.𝑤`.
"""
function ∫σᵢⱼnⱼgᵢds(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        𝑤 = ξₛ.𝑤
        # 𝑤 = ξᵤ.𝑤
        N = ξₛ[:𝝭]
        N̄ = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        n₁₁ = ξᵤ.n₁₁
        n₁₂ = ξᵤ.n₁₂
        n₂₂ = ξᵤ.n₂₂
        g₁ = ξᵤ.g₁
        g₂ = ξᵤ.g₂
        for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[3*I-2,2*J-1] += N[i]*n₁*n₁₁*N̄[j]*𝑤
                k[3*I-2,2*J]   += N[i]*n₁*n₁₂*N̄[j]*𝑤
                k[3*I-1,2*J-1] += N[i]*n₂*n₁₂*N̄[j]*𝑤
                k[3*I-1,2*J]   += N[i]*n₂*n₂₂*N̄[j]*𝑤
                k[3*I,2*J-1]   += N[i]*(n₁*n₁₂ + n₂*n₁₁)*N̄[j]*𝑤
                k[3*I,2*J]     += N[i]*(n₁*n₂₂ + n₂*n₁₂)*N̄[j]*𝑤
            end
            f[3*I-2] += N[i]*(n₁*n₁₁*g₁ + n₁*n₁₂*g₂)*𝑤
            f[3*I-1] += N[i]*(n₂*n₁₂*g₁ + n₂*n₂₂*g₂)*𝑤
            f[3*I]   += N[i]*((n₁*n₁₂+n₂*n₁₁)*g₁ + (n₁*n₂₂+n₂*n₁₂)*g₂)*𝑤 
        end
    end
end

"""
    ∫σᵢⱼnⱼuᵢds(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64})

Traction–displacement coupling on the boundary, `∫ (σ(u)·n)·v ds`.

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.n₁`, `ξ.n₂`, `ξ.𝑤`.
"""
function ∫σᵢⱼnⱼuᵢds(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        𝑤 = ξₛ.𝑤
        # 𝑤 = ξᵤ.𝑤

        N = ξₛ[:𝝭]
        N̄ = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[3*I-2,2*J-1] -= N[i]*n₁*N̄[j]*𝑤
                k[3*I-1,2*J]   -= N[i]*n₂*N̄[j]*𝑤
                k[3*I,2*J-1]   -= N[i]*n₂*N̄[j]*𝑤
                k[3*I,2*J]     -= N[i]*n₁*N̄[j]*𝑤
            end
        end
    end
end


"""
    ∫∫∇σᵢⱼuᵢdxdy(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64})

Divergence of the stress against the displacement, `∫ (∇·σ)·v dΩ`.

# Element requirements
`set∇𝝭!` — reads the gradients, `ξ.𝑤`.
"""
function ∫∫∇σᵢⱼuᵢdxdy(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        𝑤 = ξₛ.𝑤
        # 𝑤 = ξᵤ.𝑤
        B₁ = ξₛ[:∂𝝭∂x]
        B₂ = ξₛ[:∂𝝭∂y]
        N = ξᵤ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[3*I-2,2*J-1] += B₁[i]*N[j]*𝑤
                k[3*I-1,2*J]   += B₂[i]*N[j]*𝑤
                k[3*I,2*J-1]   += B₂[i]*N[j]*𝑤
                k[3*I,2*J]     += B₁[i]*N[j]*𝑤
            end
        end
    end
end

"""
    ∫∫vᵢbᵢdxdy(ap::T,f::AbstractVector{Float64})

Body force load in 2-D, `f[I] += ∫ b·v dΩ`. `b` is supplied with `prescribe!(el, :bᵢ => f)`.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.b`, `ξ.𝑤`.
"""
function ∫∫vᵢbᵢdxdy(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        b₁ = ξ.b₁
        b₂ = ξ.b₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += N[i]*b₁*𝑤
            f[2*I]   += N[i]*b₂*𝑤
        end
    end
end

"""
    ∫vᵢbᵢdΩ(ap::T,f::AbstractVector{Float64})

Body force load in 3-D, `f[I] += ∫ b·v dΩ`.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.b`, `ξ.𝑤`.
"""
function ∫vᵢbᵢdΩ(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        b₁ = ξ.b₁
        b₂ = ξ.b₂
        b₃ = ξ.b₃
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[3*I-2] += N[i]*b₁*𝑤
            f[3*I-1] += N[i]*b₂*𝑤
            f[3*I]   += N[i]*b₃*𝑤
        end
    end
end

"""
    ∫vᵢtᵢds(ap::T,f::AbstractVector{Float64})

Surface traction load in 2-D, `f[I] += ∫ t·v ds`.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.t`, `ξ.𝑤`.
"""
function ∫vᵢtᵢds(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        t₁ = ξ.t₁
        t₂ = ξ.t₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += N[i]*t₁*𝑤
            f[2*I]   += N[i]*t₂*𝑤
        end
    end
end

"""
    ∫vᵢtᵢdΓ(ap::T,f::AbstractVector{Float64})

Surface traction load in 3-D, `f[I] += ∫ t·v dΓ`.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.t`, `ξ.𝑤`.
"""
function ∫vᵢtᵢdΓ(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        t₁ = ξ.t₁
        t₂ = ξ.t₂
        t₃ = ξ.t₃
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[3*I-2] += N[i]*t₁*𝑤
            f[3*I-1] += N[i]*t₂*𝑤
            f[3*I]   += N[i]*t₃*𝑤
        end
    end
end

"""
    ∫qpdΩ(ap, k; E, ν)

Divergence–pressure coupling block, `∫ q (∇·u) dΩ` (or its transpose).

Used by mixed formulations; `q` is the pressure test function and the
divergence comes from the displacement trial function.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.𝑤`.
"""
function ∫qpdΩ(ap::T,k::AbstractMatrix{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        K = E/3/(1-2*ν)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += N[i]*N[j]/K*𝑤
            end
        end
    end
end

"""
    ∫p∇udΩ(aₚ::T,aᵤ::S,k::AbstractMatrix{Float64})

Pressure–displacement coupling in 3-D, `∫ p (∇·u) dΩ`.

# Element requirements
`set∇𝝭!` — reads `ξ[:𝝭]`, gradients, `ξ.𝑤`.
"""
function ∫p∇udΩ(aₚ::T,aᵤ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    𝓒ₚ = aₚ.𝓒;𝓖ₚ = aₚ.𝓖
    for (ξᵤ,ξₚ) in zip(𝓖ᵤ,𝓖ₚ)
        N = ξₚ[:𝝭]
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        B₃ = ξᵤ[:∂𝝭∂z]
        𝑤 = ξᵤ.𝑤
        for (i,xᵢ) in enumerate(𝓒ₚ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[I,3*J-2] -= N[i]*B₁[j]*𝑤
                k[I,3*J-1] -= N[i]*B₂[j]*𝑤
                k[I,3*J]   -= N[i]*B₃[j]*𝑤
            end
        end
    end
end

"""
    ∫∫p∇udxdy(aₚ::T,aᵤ::S,k::AbstractMatrix{Float64})

Pressure–displacement coupling in 2-D, `∫ p (∇·u) dΩ`.

# Element requirements
`set∇𝝭!` — reads `ξ[:𝝭]`, gradients, `ξ.𝑤`.
"""
function ∫∫p∇udxdy(aₚ::T,aᵤ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    𝓒ₚ = aₚ.𝓒;𝓖ₚ = aₚ.𝓖
    for (ξᵤ,ξₚ) in zip(𝓖ᵤ,𝓖ₚ)
        N = ξₚ[:𝝭]
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξᵤ.𝑤
        # 𝑤 = ξₚ.𝑤
        for (i,xᵢ) in enumerate(𝓒ₚ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[I,2*J-1] -= N[i]*B₁[j]*𝑤
                k[I,2*J]   -= N[i]*B₂[j]*𝑤
            end
        end
    end
end

"""
    ∫∫∇puᵢdxdy(aₚ::T,aᵤ::S,k::AbstractMatrix{Float64})

Transpose pressure coupling `∫ ∇p·u dΩ`.

# Element requirements
`set∇𝝭!` — reads `ξ[:𝝭]`, gradients, `ξ.𝑤`.
"""
function ∫∫∇puᵢdxdy(aₚ::T,aᵤ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    𝓒ₚ = aₚ.𝓒;𝓖ₚ = aₚ.𝓖
    for (ξᵤ,ξₚ) in zip(𝓖ᵤ,𝓖ₚ)
        B₁ = ξₚ[:∂𝝭∂x]
        B₂ = ξₚ[:∂𝝭∂y]
        N = ξᵤ[:𝝭]
        𝑤 = ξᵤ.𝑤
        for (i,xᵢ) in enumerate(𝓒ₚ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[I,2*J-1] += B₁[i]*N[j]*𝑤
                k[I,2*J]   += B₂[i]*N[j]*𝑤
            end
        end
    end
end

"""
    ∫pnᵢuᵢds(aₚ::T,aᵤ::S,k::AbstractMatrix{Float64})

Boundary pressure term `∫ p (n·u) ds`.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.n₁`, `ξ.n₂`, `ξ.𝑤`.
"""
function ∫pnᵢuᵢds(aₚ::T,aᵤ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    𝓒ₚ = aₚ.𝓒;𝓖ₚ = aₚ.𝓖
    for (ξᵤ,ξₚ) in zip(𝓖ᵤ,𝓖ₚ)
        Nₚ = ξₚ[:𝝭]
        Nᵤ = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        𝑤 = ξᵤ.𝑤
        # 𝑤 = ξₚ.𝑤
        for (i,xᵢ) in enumerate(𝓒ₚ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[I,2*J-1] -= Nₚ[i]*Nᵤ[j]*n₁*𝑤
                k[I,2*J]   -= Nₚ[i]*Nᵤ[j]*n₂*𝑤
            end
        end
    end
end

"""
    ∫pnᵢgᵢds(aₚ::T,aᵤ::S,k::AbstractMatrix{Float64},f::AbstractVector{Float64})

Boundary pressure load `∫ p (n·g) ds`.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.n₁`, `ξ.n₂`, `ξ.g`, `ξ.𝑤`.
"""
function ∫pnᵢgᵢds(aₚ::T,aᵤ::S,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    𝓒ₚ = aₚ.𝓒;𝓖ₚ = aₚ.𝓖
    for (ξᵤ,ξₚ) in zip(𝓖ᵤ,𝓖ₚ)
        Nₚ = ξₚ[:𝝭]
        Nᵤ = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        g₁ = ξᵤ.g₁
        g₂ = ξᵤ.g₂
        n₁₁ = ξᵤ.n₁₁
        n₁₂ = ξᵤ.n₁₂
        n₂₂ = ξᵤ.n₂₂
        𝑤 = ξᵤ.𝑤
        # 𝑤 = ξₚ.𝑤
        for (i,xᵢ) in enumerate(𝓒ₚ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[I,2*J-1] += Nₚ[i]*Nᵤ[j]*(n₁*n₁₁+n₂*n₁₂)*𝑤
                k[I,2*J]   += Nₚ[i]*Nᵤ[j]*(n₁*n₁₂+n₂*n₂₂)*𝑤
            end
            f[I] += Nₚ[i]*(n₁*n₁₁*g₁+n₁*n₁₂*g₂+n₂*n₁₂*g₁+n₂*n₂₂*g₂)*𝑤
        end
    end
end

"""
    ∫∫sᵢⱼsᵢⱼdxdy(ap, k; E, ν)

Deviation norm of the stress, `∫ sᵢⱼ sᵢⱼ dΩ` — a post-processing operator, not
a stiffness block.

Reports how much deviatoric stress a solution carries; useful as an error
indicator for nearly-incompressible problems where the pressure dominates.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.𝑤`.
"""
function ∫∫sᵢⱼsᵢⱼdxdy(ap::T,k::AbstractMatrix{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        G = E/(1+ν)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[4*I-3,4*J-3] +=   N[i]*N[j]/G*𝑤
                k[4*I-2,4*J-2] +=   N[i]*N[j]/G*𝑤
                k[4*I-1,4*J-1] +=   N[i]*N[j]/G*𝑤
                k[4*I,4*J]     += 2*N[i]*N[j]/G*𝑤
            end
        end
    end
end

"""
    ∫∫sᵢⱼεᵢⱼdxdy(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64})

Deviatoric stress–strain coupling `∫ s : ε dΩ`.

# Element requirements
`set∇𝝭!` — reads gradients, `ξ.𝑤`.
"""
function ∫∫sᵢⱼεᵢⱼdxdy(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    for (ξᵤ,ξₛ) in zip(𝓖ᵤ,𝓖ₛ)
        N = ξₛ[:𝝭]
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξᵤ.𝑤
        for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[4*I-3,2*J-1] -= 2/3*N[i]*B₁[j]*𝑤
                k[4*I-3,2*J]   += 1/3*N[i]*B₂[j]*𝑤
                k[4*I-2,2*J-1] += 1/3*N[i]*B₁[j]*𝑤
                k[4*I-2,2*J]   -= 2/3*N[i]*B₂[j]*𝑤
                k[4*I-1,2*J-1] += 1/3*N[i]*B₁[j]*𝑤
                k[4*I-1,2*J]   += 1/3*N[i]*B₂[j]*𝑤
                k[4*I,2*J-1]   -=     N[i]*B₂[j]*𝑤
                k[4*I,2*J]     -=     N[i]*B₁[j]*𝑤
            end
        end
    end
end

"""
    ∫∫∇sᵢⱼuᵢdxdy(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64})

Divergence of the deviatoric stress `∫ (∇·s)·u dΩ`.

# Element requirements
`set∇𝝭!` — reads gradients, `ξ.𝑤`.
"""
function ∫∫∇sᵢⱼuᵢdxdy(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    for (ξᵤ,ξₛ) in zip(𝓖ᵤ,𝓖ₛ)
        B₁ = ξₛ[:∂𝝭∂x]
        B₂ = ξₛ[:∂𝝭∂y]
        Nᵤ = ξᵤ[:𝝭]
        𝑤 = ξᵤ.𝑤
        for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[4*I-3,2*J-1] += 2/3*B₁[i]*Nᵤ[j]*𝑤
                k[4*I-3,2*J]   -= 1/3*B₂[i]*Nᵤ[j]*𝑤
                k[4*I-2,2*J-1] -= 1/3*B₁[i]*Nᵤ[j]*𝑤
                k[4*I-2,2*J]   += 2/3*B₂[i]*Nᵤ[j]*𝑤
                k[4*I-1,2*J-1] -= 1/3*B₁[i]*Nᵤ[j]*𝑤
                k[4*I-1,2*J]   -= 1/3*B₂[i]*Nᵤ[j]*𝑤
                k[4*I,2*J-1]   += B₂[i]*Nᵤ[j]*𝑤
                k[4*I,2*J]     += B₁[i]*Nᵤ[j]*𝑤
            end
        end
    end
end

"""
    ∫sᵢⱼnⱼuᵢds(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64})

Deviatoric traction coupling `∫ (s·n)·u ds`.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.n₁`, `ξ.n₂`, `ξ.𝑤`.
"""
function ∫sᵢⱼnⱼuᵢds(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    for (ξᵤ,ξₛ) in zip(𝓖ᵤ,𝓖ₛ)
        Nₛ = ξₛ[:𝝭]
        Nᵤ = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        𝑤 = ξᵤ.𝑤
        for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[4*I-3,2*J-1] -= 2/3*Nₛ[i]*Nᵤ[j]*n₁*𝑤
                k[4*I-3,2*J]   += 1/3*Nₛ[i]*Nᵤ[j]*n₂*𝑤
                k[4*I-2,2*J-1] += 1/3*Nₛ[i]*Nᵤ[j]*n₁*𝑤
                k[4*I-2,2*J]   -= 2/3*Nₛ[i]*Nᵤ[j]*n₂*𝑤
                k[4*I-1,2*J-1] += 1/3*Nₛ[i]*Nᵤ[j]*n₁*𝑤
                k[4*I-1,2*J]   += 1/3*Nₛ[i]*Nᵤ[j]*n₂*𝑤
                k[4*I,2*J-1]   -= Nₛ[i]*Nᵤ[j]*n₂*𝑤
                k[4*I,2*J]     -= Nₛ[i]*Nᵤ[j]*n₁*𝑤
            end
        end
    end
end

"""
    ∫sᵢⱼnⱼgᵢds(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64},f::AbstractVector{Float64})

Deviatoric traction load `∫ (s·n)·g ds`.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.n₁`, `ξ.n₂`, `ξ.g`, `ξ.𝑤`.
"""
function ∫sᵢⱼnⱼgᵢds(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    for (ξᵤ,ξₛ) in zip(𝓖ᵤ,𝓖ₛ)
        Nₛ = ξₛ[:𝝭]
        Nᵤ = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        g₁ = ξᵤ.g₁
        g₂ = ξᵤ.g₂
        n₁₁ = ξᵤ.n₁₁
        n₁₂ = ξᵤ.n₁₂
        n₂₂ = ξᵤ.n₂₂
        𝑤 = ξᵤ.𝑤
        for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[4*I-3,2*J-1] += Nₛ[i]*Nᵤ[j]*( 2/3*n₁*n₁₁ - 1/3*n₂*n₁₂)*𝑤
                k[4*I-3,2*J]   += Nₛ[i]*Nᵤ[j]*( 2/3*n₁*n₁₂ - 1/3*n₂*n₂₂)*𝑤
                k[4*I-2,2*J-1] += Nₛ[i]*Nᵤ[j]*(-1/3*n₁*n₁₁ + 2/3*n₂*n₁₂)*𝑤
                k[4*I-2,2*J]   += Nₛ[i]*Nᵤ[j]*(-1/3*n₁*n₁₂ + 2/3*n₂*n₂₂)*𝑤
                k[4*I-1,2*J-1] += Nₛ[i]*Nᵤ[j]*(-1/3*n₁*n₁₁ - 1/3*n₂*n₁₂)*𝑤
                k[4*I-1,2*J]   += Nₛ[i]*Nᵤ[j]*(-1/3*n₁*n₁₂ - 1/3*n₂*n₂₂)*𝑤
                k[4*I,2*J-1]   += Nₛ[i]*Nᵤ[j]*(n₁*n₁₂ + n₂*n₁₁)*𝑤
                k[4*I,2*J]     += Nₛ[i]*Nᵤ[j]*(n₁*n₂₂ + n₂*n₁₂)*𝑤
            end
            f[4*I-3] += Nₛ[i]*(( 2/3*n₁*n₁₁-1/3*n₂*n₁₂)*g₁+( 2/3*n₁*n₁₂-1/3*n₂*n₂₂)*g₂)*𝑤
            f[4*I-2] += Nₛ[i]*((-1/3*n₁*n₁₁+2/3*n₂*n₁₂)*g₁+(-1/3*n₁*n₁₂+2/3*n₂*n₂₂)*g₂)*𝑤
            f[4*I-1] += Nₛ[i]*((-1/3*n₁*n₁₁-1/3*n₂*n₁₂)*g₁+(-1/3*n₁*n₁₂-1/3*n₂*n₂₂)*g₂)*𝑤
            f[4*I]   += Nₛ[i]*((n₁*n₁₂+n₂*n₁₁)*g₁+(n₁*n₂₂+n₂*n₁₂)*g₂)*𝑤
        end
    end
end


"""
    g₂(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64},dof::Symbol)

Element-set overload of `g₂`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.


    Assemble the prescribed value into the right-hand side and zero the
    corresponding row and column of the stiffness matrix — the "row/column
    elimination" strong-imposition of a Dirichlet condition:

        f[I] -= K[I,J] g        for the constrained DOF J
        K[I,J] = K[J,I] = 0,    K[J,J] = 1,    f[J] = g

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function g₂(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64},dof::Symbol) where T<:AbstractElement
    x, = ap.𝓒
    if dof == :d₁
        j = 2*x.𝐼-1
    else
        j = 2*x.𝐼
    end
    g = getproperty(x,dof)
    for i in eachindex(f)
        f[i] -= k[i,j]*g
    end
    k[j,:] .= 0.
    k[:,j] .= 0.
    k[j,j] = 1.
    f[j] = g
end

"""
    ∫λᵢgᵢds(aₗ::T,aᵤ::S,k::AbstractMatrix{Float64},f::AbstractVector{Float64})

Lagrange-multiplier boundary term `∫ λ·g ds` for a mixed essential condition.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.g`, `ξ.𝑤`.
"""
function ∫λᵢgᵢds(aₗ::T,aᵤ::S,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₗ = aₗ.𝓒;𝓖ₗ = aₗ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    for (ξₗ,ξᵤ) in zip(𝓖ₗ,𝓖ᵤ)
        𝑤 = ξₗ.𝑤
        N̄ = ξₗ[:𝝭]
        N = ξᵤ[:𝝭]
        g₁ = ξᵤ.g₁
        g₂ = ξᵤ.g₂
        n₁₁ = ξᵤ.n₁₁
        n₂₂ = ξᵤ.n₂₂
        n₁₂ = ξᵤ.n₁₂
        for (i,xᵢ) in enumerate(𝓒ₗ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] -= n₁₁*N̄[i]*N[j]*𝑤
                k[2*I-1,2*J]   -= n₁₂*N̄[i]*N[j]*𝑤
                k[2*I,2*J-1]   -= n₁₂*N̄[i]*N[j]*𝑤
                k[2*I,2*J]     -= n₂₂*N̄[i]*N[j]*𝑤
            end
            f[2*I-1] -= N̄[i]*(n₁₁*g₁+n₁₂*g₂)*𝑤
            f[2*I]   -= N̄[i]*(n₁₂*g₁+n₂₂*g₂)*𝑤
        end
    end
end

"""
    ∫σᵢⱼnⱼgᵢds(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}; E, ν)

Element-set overload of `∫σᵢⱼnⱼgᵢds`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function ∫σᵢⱼnⱼgᵢds(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        n₁₁ = ξ.n₁₁
        n₁₂ = ξ.n₁₂
        n₂₂ = ξ.n₂₂
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        Cᵢᵢᵢᵢ = E/(1-ν^2)
        Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        C₁₁₁ = Cᵢᵢᵢᵢ*n₁*n₁₁+Cᵢᵢⱼⱼ*n₂*n₁₂
        C₁₁₂ = Cᵢᵢᵢᵢ*n₁*n₁₂+Cᵢᵢⱼⱼ*n₂*n₂₂
        C₂₂₁ = Cᵢᵢⱼⱼ*n₁*n₁₁+Cᵢᵢᵢᵢ*n₂*n₁₂
        C₂₂₂ = Cᵢᵢⱼⱼ*n₁*n₁₂+Cᵢᵢᵢᵢ*n₂*n₂₂
        C₁₂₁ = Cᵢⱼᵢⱼ*(n₁*n₁₂+n₂*n₁₁)
        C₁₂₂ = Cᵢⱼᵢⱼ*(n₂*n₁₂+n₁*n₂₂)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] -= (C₁₁₁*(N[i]*B₁[j]+B₁[i]*N[j]) + C₁₂₁*(N[i]*B₂[j]+B₂[i]*N[j]))*𝑤
                k[2*I-1,2*J]   -= (C₁₂₁*N[i]*B₁[j] + C₁₁₂*B₁[i]*N[j] + C₂₂₁*N[i]*B₂[j] + C₁₂₂*B₂[i]*N[j])*𝑤
                k[2*I,2*J-1]   -= (C₁₁₂*N[i]*B₁[j] + C₁₂₁*B₁[i]*N[j] + C₁₂₂*N[i]*B₂[j] + C₂₂₁*B₂[i]*N[j])*𝑤
                k[2*I,2*J]     -= (C₁₂₂*(N[i]*B₁[j]+B₁[i]*N[j]) + C₂₂₂*(N[i]*B₂[j]+B₂[i]*N[j]))*𝑤
            end
            f[2*I-1] -= ((C₁₁₁*B₁[i]+C₁₂₁*B₂[i])*g₁ + (C₁₁₂*B₁[i]+C₁₂₂*B₂[i])*g₂)*𝑤
            f[2*I]   -= ((C₁₂₁*B₁[i]+C₂₂₁*B₂[i])*g₁ + (C₁₂₂*B₁[i]+C₂₂₂*B₂[i])*g₂)*𝑤
        end
    end
end

"""
    ∫σᵢⱼnⱼgᵢvᵢgᵢds(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}; E, ν, α)

Element-set overload of `∫σᵢⱼnⱼgᵢvᵢgᵢds`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function ∫σᵢⱼnⱼgᵢvᵢgᵢds(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}; E::Real, ν::Real, α::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        n₁₁ = ξ.n₁₁
        n₁₂ = ξ.n₁₂
        n₂₂ = ξ.n₂₂
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        Cᵢᵢᵢᵢ = E/(1-ν^2)
        Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        C₁₁₁ = Cᵢᵢᵢᵢ*n₁*n₁₁+Cᵢᵢⱼⱼ*n₂*n₁₂
        C₁₁₂ = Cᵢᵢᵢᵢ*n₁*n₁₂+Cᵢᵢⱼⱼ*n₂*n₂₂
        C₂₂₁ = Cᵢᵢⱼⱼ*n₁*n₁₁+Cᵢᵢᵢᵢ*n₂*n₁₂
        C₂₂₂ = Cᵢᵢⱼⱼ*n₁*n₁₂+Cᵢᵢᵢᵢ*n₂*n₂₂
        C₁₂₁ = Cᵢⱼᵢⱼ*(n₁*n₁₂+n₂*n₁₁)
        C₁₂₂ = Cᵢⱼᵢⱼ*(n₂*n₁₂+n₁*n₂₂)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] -= (C₁₁₁*(N[i]*B₁[j]+B₁[i]*N[j]) + C₁₂₁*(N[i]*B₂[j]+B₂[i]*N[j]) - α*N[i]*n₁₁*N[j])*𝑤
                k[2*I-1,2*J]   -= (C₁₂₁*N[i]*B₁[j] + C₁₁₂*B₁[i]*N[j] + C₂₂₁*N[i]*B₂[j] + C₁₂₂*B₂[i]*N[j] -α*N[i]*n₁₂*N[j])*𝑤
                k[2*I,2*J-1]   -= (C₁₁₂*N[i]*B₁[j] + C₁₂₁*B₁[i]*N[j] + C₁₂₂*N[i]*B₂[j] + C₂₂₁*B₂[i]*N[j] -α*N[i]*n₁₂*N[j])*𝑤
                k[2*I,2*J]     -= (C₁₂₂*(N[i]*B₁[j]+B₁[i]*N[j]) + C₂₂₂*(N[i]*B₂[j]+B₂[i]*N[j]) -α*N[i]*n₂₂*N[j])*𝑤
            end
            f[2*I-1] -= ((C₁₁₁*B₁[i]+C₁₂₁*B₂[i]-α*N[i]*n₁₁)*g₁ + (C₁₁₂*B₁[i]+C₁₂₂*B₂[i]-α*N[i]*n₁₂)*g₂)*𝑤
            f[2*I]   -= ((C₁₂₁*B₁[i]+C₂₂₁*B₂[i]-α*N[i]*n₁₂)*g₁ + (C₁₂₂*B₁[i]+C₂₂₂*B₂[i]-α*N[i]*n₂₂)*g₂)*𝑤
        end
    end
end

"""
    ∫vᵢgᵢds(ap, k, f; α)
    ∫vᵢgᵢds(ap, k; α)
    ∫vᵢgᵢdΓ(ap, k, f; α)

Penalty enforcement of an essential boundary condition, `α ∫ (u-g)·v ds`.

    k[2I-1, 2J-1] += α Nᵢ Nⱼ w        (and the 2-2 block)
    f[2I-1]       += α Nᵢ g₁ w

`g` is the prescribed boundary value, supplied with `prescribe!(el, :g => f)`.

# Keyword arguments
- `α` : penalty coefficient. Large (≳ 10⁶ × the diagonal scale) enforces the
  constraint more strongly at the cost of conditioning.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ[:g]`, `ξ.𝑤`.
"""
function ∫vᵢgᵢds(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}; α::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        n₁₁ = ξ.n₁₁
        n₂₂ = ξ.n₂₂
        n₁₂ = ξ.n₁₂
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += α*N[i]*n₁₁*N[j]*𝑤
                k[2*I,2*J-1]   += α*N[i]*n₁₂*N[j]*𝑤
                k[2*I-1,2*J]   += α*N[i]*n₁₂*N[j]*𝑤
                k[2*I,2*J]     += α*N[i]*n₂₂*N[j]*𝑤
            end
            f[2*I-1] += α*N[i]*(n₁₁*g₁+n₁₂*g₂)*𝑤
            f[2*I]   += α*N[i]*(n₁₂*g₁+n₂₂*g₂)*𝑤
        end
    end
end

"""
    ∫vᵢgᵢds(ap::T,k::AbstractMatrix{Float64}; α)

Element-set overload of `∫vᵢgᵢds`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function ∫vᵢgᵢds(ap::T,k::AbstractMatrix{Float64}; α::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        n₁₁ = ξ.n₁₁
        n₂₂ = ξ.n₂₂
        n₁₂ = ξ.n₁₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += α*N[i]*n₁₁*N[j]*𝑤
                k[2*I,2*J-1]   += α*N[i]*n₁₂*N[j]*𝑤
                k[2*I-1,2*J]   += α*N[i]*n₁₂*N[j]*𝑤
                k[2*I,2*J]     += α*N[i]*n₂₂*N[j]*𝑤
            end
        end
    end
end


"""
    ∫vᵢgᵢdΓ(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}; α)

Element-set overload of `∫vᵢgᵢdΓ`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function ∫vᵢgᵢdΓ(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}; α::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        n₁₁ = ξ.n₁₁
        n₁₂ = ξ.n₁₂
        n₁₃ = ξ.n₁₃
        n₂₂ = ξ.n₂₂
        n₂₃ = ξ.n₂₃
        n₃₃ = ξ.n₃₃
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        g₃ = ξ.g₃
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += α*N[i]*n₁₁*N[j]*𝑤
                k[3*I-2,3*J-1] += α*N[i]*n₁₂*N[j]*𝑤
                k[3*I-2,3*J]   += α*N[i]*n₁₃*N[j]*𝑤
                k[3*I-1,3*J-2] += α*N[i]*n₁₂*N[j]*𝑤
                k[3*I-1,3*J-1] += α*N[i]*n₂₂*N[j]*𝑤
                k[3*I-1,3*J]   += α*N[i]*n₂₃*N[j]*𝑤
                k[3*I,3*J-2]   += α*N[i]*n₁₃*N[j]*𝑤
                k[3*I,3*J-1]   += α*N[i]*n₂₃*N[j]*𝑤
                k[3*I,3*J]     += α*N[i]*n₃₃*N[j]*𝑤
            end
            f[3*I-2] += α*N[i]*(n₁₁*g₁+n₁₂*g₂+n₁₃*g₃)*𝑤
            f[3*I-1] += α*N[i]*(n₁₂*g₁+n₂₂*g₂+n₂₃*g₃)*𝑤
            f[3*I]   += α*N[i]*(n₁₃*g₁+n₂₃*g₂+n₃₃*g₃)*𝑤
        end
    end
end

"""
    ∫∫τ∇q∇pdxdy(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64})

Stabilisation `τ ∫ ∇q·∇p dΩ` with the local time/streamline scale `τ`.

# Element requirements
`set∇𝝭!` — reads gradients, `ξ.τ`, `ξ.𝑤`.
"""
function ∫∫τ∇q∇pdxdy(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        τ = ξ.τ
        b₁ = ξ.b₁
        b₂ = ξ.b₂
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += τ*(B₁[i]*B₁[j] + B₂[i]*B₂[j])*𝑤
            end
            f[I] += τ*(B₁[i]*b₁ + B₂[i]*b₂)*𝑤
        end
    end
end



"""
    ∫∫τ∇sᵢⱼ∇sᵢₖdxdy(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64})

Stabilisation of the deviatoric stress `τ ∫ ∇s:∇s dΩ`.

# Element requirements
`set∇𝝭!` — reads gradients, `ξ.τ`, `ξ.𝑤`.
"""
function ∫∫τ∇sᵢⱼ∇sᵢₖdxdy(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        τ = ξ.τ
        b₁ = ξ.b₁
        b₂ = ξ.b₂
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[4*I-3,4*J-3] += τ*B₁[i]*B₁[j]*𝑤
                k[4*I-3,4*J]   += τ*B₁[i]*B₂[j]*𝑤
                k[4*I-2,4*J-2] += τ*B₂[i]*B₂[j]*𝑤
                k[4*I-2,4*J]   += τ*B₂[i]*B₁[j]*𝑤
                k[4*I,4*J-3]   += τ*B₂[i]*B₁[j]*𝑤
                k[4*I,4*J-2]   += τ*B₁[i]*B₂[j]*𝑤
                k[4*I,4*J]     += τ*(B₁[i]*B₁[j] + B₂[i]*B₂[j])*𝑤
            end
            f[4*I-3] += τ*B₁[i]*b₁*𝑤
            f[4*I-2] += τ*B₂[i]*b₂*𝑤
            f[4*I]   += τ*(B₁[i]*b₂ + B₂[i]*b₁)*𝑤
        end
    end
end
"""
    ∫∫τ∇sᵢⱼ∇pdxdy(aₛ::T,aₚ::S,k::AbstractMatrix{Float64})

Stress–pressure stabilisation coupling `τ ∫ ∇s:∇p dΩ`.

# Element requirements
`set∇𝝭!` — reads gradients, `ξ.τ`, `ξ.𝑤`.
"""
function ∫∫τ∇sᵢⱼ∇pdxdy(aₛ::T,aₚ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₚ = aₚ.𝓒; 𝓖ₚ = aₚ.𝓖
    𝓒ₛ = aₛ.𝓒; 𝓖ₛ = aₛ.𝓖
    for (ξₚ,ξₛ) in zip(𝓖ₚ,𝓖ₛ)
        𝑤 = ξₚ.𝑤
        τ = ξₚ.τ
        B₁ = ξₛ[:∂𝝭∂x]
        B₂ = ξₛ[:∂𝝭∂y]
        B̄₁ = ξₚ[:∂𝝭∂x]
        B̄₂ = ξₚ[:∂𝝭∂y]
        for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ₚ)
                J = xⱼ.𝐼
                k[4*I-3,J] += τ*B₁[i]*B̄₁[j]*𝑤
                k[4*I-2,J] += τ*B₂[i]*B̄₂[j]*𝑤
                k[4*I,J]   += τ*(B₁[i]*B̄₂[j] + B₂[i]*B̄₁[j])*𝑤
            end
        end
    end
end
"""
    ∫∫τ∇σᵢⱼ∇σᵢₖdxdy(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64})

Stabilisation of the full stress `τ ∫ ∇σ:∇σ dΩ`.

# Element requirements
`set∇𝝭!` — reads gradients, `ξ.τ`, `ξ.𝑤`.
"""
function ∫∫τ∇σᵢⱼ∇σᵢₖdxdy(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        τ = ξ.τ
        b₁ = ξ.b₁
        b₂ = ξ.b₂
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += τ*B₁[i]*B₁[j]*𝑤
                k[3*I-2,3*J]   += τ*B₁[i]*B₂[j]*𝑤
                k[3*I-1,3*J-1] += τ*B₂[i]*B₂[j]*𝑤
                k[3*I-1,3*J]   += τ*B₂[i]*B₁[j]*𝑤
                k[3*I,3*J-2]   += τ*B₂[i]*B₁[j]*𝑤
                k[3*I,3*J-1]   += τ*B₁[i]*B₂[j]*𝑤
                k[3*I,3*J]     += τ*(B₁[i]*B₁[j] + B₂[i]*B₂[j])*𝑤
            end
            f[3*I-2] += τ*B₁[i]*b₁*𝑤
            f[3*I-1] += τ*B₂[i]*b₂*𝑤
            f[3*I]   += τ*(B₁[i]*b₂ + B₂[i]*b₁)*𝑤
        end
    end
end
"""
    getσₙ(σ₁₁, σ₂₂, σ₁₂) -> (σ₁, σ₂, 𝒏₁, 𝒏₂)

Principal stresses and their unit directions of a 2-D stress state.

    σ₁,₂ = (σ₁₁ + σ₂₂)/2 ± √( ((σ₁₁ - σ₂₂)/2)² + σ₁₂² )

`𝒏₁` is the eigenvector belonging to `σ₁`, `𝒏₂` that of `σ₂`. When `σ₁₂ ≈ 0` the
state is already diagonal and the axes are returned directly.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function getσₙ(σ₁₁::Float64,σ₂₂::Float64,σ₁₂::Float64)
    # trace
    t = σ₁₁ + σ₂₂
    # determinant
    d = σ₁₁*σ₂₂ - σ₁₂^2

    σ₁ = t/2 + (t^2/4-d)^0.5
    σ₂ = t/2 - (t^2/4-d)^0.5
    if σ₁₂ ≈ 0.0
        (n₁,n₂) = σ₁₁ > σ₂₂ ? ((1.0,0.0),(0.0,1.0)) : ((0.0,1.0),(1.0,0.0))
    else
        l₁ = ((σ₁-σ₂₂)^2+σ₁₂^2)^0.5
        l₂ = ((σ₂-σ₂₂)^2+σ₁₂^2)^0.5
        n₁ = ((σ₁-σ₂₂)/l₁,σ₁₂/l₁)
        n₂ = ((σ₂-σ₂₂)/l₁,σ₁₂/l₁)
    end
    return σ₁,σ₂,n₁,n₂
end

"""
    getσₙ(σ₁₁, σ₂₂, σ₃₃, σ₁₂, σ₁₃, σ₂₃) -> (σ₁, σ₂, σ₃, 𝒏₁, 𝒏₂, 𝒏₃)

Principal stresses and their unit directions of a 3-D stress state, ordered
`σ₁ ≥ σ₂ ≥ σ₃`. The eigenvalues follow from the closed-form solution of

    det(σ - σᵢ I) = 0

with `p = tr(σ)/3` and `s` the second invariant of the deviator, as used by the
Mohr–Coulomb return mapping.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function getσₙ(σ₁₁::Float64,σ₂₂::Float64,σ₃₃::Float64,σ₁₂::Float64,σ₁₃::Float64,σ₂₃::Float64)
    p₁ = σ₁₂^2+σ₁₃^2+σ₂₃^2
    if p₁ == 0.0
        σ₁ = σ₁₁
        σ₂ = σ₂₂
        σ₃ = σ₃₃
    else
        p = (σ₁₁+σ₂₂+σ₃₃)/3.0
        s = (((σ₁₁-p)^2+(σ₂₂-p)^2+(σ₃₃-p)^2 + 2*p₁)/6)^0.5
        (B₁₁,B₂₂,B₃₃,B₁₂,B₁₃,B₂₃) = 1.0/s .* (σ₁₁-p,σ₂₂-p,σ₃₃-p,σ₁₂,σ₁₃,σ₂₃)
        r = (B₁₁*B₂₂*B₃₃+2.0*B₁₂*B₁₃*B₂₃-B₁₁*B₂₃^2-B₂₂*B₁₃^2-B₃₃*B₁₂^2)/2.0

        if r≤-1.0
            θ = π/3.0
        elseif r≥1.0
            θ = 0.0
        else
            θ = acos(r)/3.0
        end

        σ₁ = p + 2.0*s*cos(θ)
        σ₃ = p + 2.0*s*cos(θ+2.0/3.0*π)
        σ₂ = 3.0*p -σ₁-σ₃
    end

    N₁ = (
        (σ₁₁-σ₂)*(σ₁₁-σ₃)+σ₁₂*σ₁₂     +σ₁₃*σ₁₃
       +(σ₁₁-σ₂)*σ₁₂     +σ₁₂*(σ₂₂-σ₃)+σ₁₃*σ₂₃
       +(σ₁₁-σ₂)*σ₁₃     +σ₁₂*σ₂₃     +σ₁₃*(σ₃₃-σ₃),
        σ₁₂*(σ₁₁-σ₃)+(σ₂₂-σ₂)*σ₁₂     +σ₂₃*σ₁₃
       +σ₁₂*σ₁₂     +(σ₂₂-σ₂)*(σ₂₂-σ₃)+σ₂₃*σ₂₃
       +σ₁₂*σ₁₃     +(σ₂₂-σ₂)*σ₂₃     +σ₂₃*(σ₃₃-σ₃),
        σ₁₃*(σ₁₁-σ₃)+σ₂₃*σ₁₂     +(σ₃₃-σ₂)*σ₁₃
       +σ₁₃*σ₁₂     +σ₂₃*(σ₂₂-σ₃)+(σ₃₃-σ₂)*σ₂₃
       +σ₁₃*σ₁₃     +σ₂₃*σ₂₃     +(σ₃₃-σ₂)*(σ₃₃-σ₃)
    )
    N₂ = (
        (σ₁₁-σ₁)*(σ₁₁-σ₃)+σ₁₂*σ₁₂     +σ₁₃*σ₁₃
       +(σ₁₁-σ₁)*σ₁₂     +σ₁₂*(σ₂₂-σ₃)+σ₁₃*σ₂₃
       +(σ₁₁-σ₁)*σ₁₃     +σ₁₂*σ₂₃     +σ₁₃*(σ₃₃-σ₃),
        σ₁₂*(σ₁₁-σ₃)+(σ₂₂-σ₁)*σ₁₂     +σ₂₃*σ₁₃
       +σ₁₂*σ₁₂     +(σ₂₂-σ₁)*(σ₂₂-σ₃)+σ₂₃*σ₂₃
       +σ₁₂*σ₁₃     +(σ₂₂-σ₁)*σ₂₃     +σ₂₃*(σ₃₃-σ₃),
        σ₁₃*(σ₁₁-σ₃)+σ₂₃*σ₁₂     +(σ₃₃-σ₁)*σ₁₃
       +σ₁₃*σ₁₂     +σ₂₃*(σ₂₂-σ₃)+(σ₃₃-σ₁)*σ₂₃
       +σ₁₃*σ₁₃     +σ₂₃*σ₂₃     +(σ₃₃-σ₁)*(σ₃₃-σ₃)
    )
    N₃ = (
        (σ₁₁-σ₁)*(σ₁₁-σ₂)+σ₁₂*σ₁₂     +σ₁₃*σ₁₃
       +(σ₁₁-σ₁)*σ₁₂     +σ₁₂*(σ₂₂-σ₂)+σ₁₃*σ₂₃
       +(σ₁₁-σ₁)*σ₁₃     +σ₁₂*σ₂₃     +σ₁₃*(σ₃₃-σ₂),
        σ₁₂*(σ₁₁-σ₂)+(σ₂₂-σ₁)*σ₁₂     +σ₂₃*σ₁₃
       +σ₁₂*σ₁₂     +(σ₂₂-σ₁)*(σ₂₂-σ₂)+σ₂₃*σ₂₃
       +σ₁₂*σ₁₃     +(σ₂₂-σ₁)*σ₂₃     +σ₂₃*(σ₃₃-σ₂),
        σ₁₃*(σ₁₁-σ₂)+σ₂₃*σ₁₂     +(σ₃₃-σ₁)*σ₁₃
       +σ₁₃*σ₁₂     +σ₂₃*(σ₂₂-σ₂)+(σ₃₃-σ₁)*σ₂₃
       +σ₁₃*σ₁₃     +σ₂₃*σ₂₃     +(σ₃₃-σ₁)*(σ₃₃-σ₂)
    )
    normN₁ = (N₁[1]^2+N₁[2]^2+N₁[3]^2)^0.5
    normN₂ = (N₂[1]^2+N₂[2]^2+N₂[3]^2)^0.5
    normN₃ = (N₃[1]^2+N₃[2]^2+N₃[3]^2)^0.5
    n₁ = (N₁[1]/normN₁,N₁[2]/normN₁,N₁[3]/normN₁)
    n₂ = (N₂[1]/normN₂,N₂[2]/normN₂,N₂[3]/normN₂)
    n₃ = (N₃[1]/normN₃,N₃[2]/normN₃,N₃[3]/normN₃)
    return σ₁,σ₂,σ₃,n₁,n₂,n₃
end

"""
    L₂(ap::T)

Relative L₂ error of the corresponding field against an exact
solution, `‖e‖/‖u‖`. The exact field is supplied with `prescribe!`.

# Element requirements
`set𝝭!` (or `set∇𝝭!` for an H¹ norm) — reads `ξ[:𝝭]`, `ξ.𝑤` and the exact field.
"""
function L₂(ap::T) where T<:AbstractElement
    Δu²= 0
    ū² = 0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        ū₁ = ξ.u₁
        ū₂ = ξ.u₂
        ū₃ = ξ.u₃
        u₁ = 0.
        u₂ = 0.
        u₃ = 0.
        for (i,xᵢ) in enumerate(ap.𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            u₃ += N[i]*xᵢ.d₃
        end
        Δu² += ((u₁ - ū₁)^2 + (u₂ - ū₂)^2 + (u₃ - ū₃)^2)*𝑤
        ū² += (ū₁^2 + ū₂^2 + ū₃^2)*𝑤
    end
    return Δu², ū²
end

"""
    L₂(aps::Vector{T})

Relative L₂ error of the corresponding field against an exact
solution, `‖e‖/‖u‖`. The exact field is supplied with `prescribe!`.

# Element requirements
`set𝝭!` (or `set∇𝝭!` for an H¹ norm) — reads `ξ[:𝝭]`, `ξ.𝑤` and the exact field.
"""
function L₂(aps::Vector{T}) where T<:AbstractElement
    L₂Norm_Δu²= 0.0
    L₂Norm_ū² = 0.0
    for ap in aps
        Δu², ū² = L₂(ap)
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end

"""
    L₂𝑝(ap::T)

Relative L₂ error of the corresponding field against an exact
solution, `‖e‖/‖u‖`. The exact field is supplied with `prescribe!`.

# Element requirements
`set𝝭!` (or `set∇𝝭!` for an H¹ norm) — reads `ξ[:𝝭]`, `ξ.𝑤` and the exact field.
"""
function L₂𝑝(ap::T) where T<:AbstractElement
    Δp²= 0
    p̄² = 0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        p̄ = ξ.p
        p = 0.
        for (i,xᵢ) in enumerate(ap.𝓒)
            p += N[i]*xᵢ.p
        end
        Δp² += (p - p̄)^2*𝑤
        p̄² += p̄^2*𝑤
    end
    return Δp², p̄²
end

"""
    L₂𝑝(aps::Vector{T})

Relative L₂ error of the corresponding field against an exact
solution, `‖e‖/‖u‖`. The exact field is supplied with `prescribe!`.

# Element requirements
`set𝝭!` (or `set∇𝝭!` for an H¹ norm) — reads `ξ[:𝝭]`, `ξ.𝑤` and the exact field.
"""
function L₂𝑝(aps::Vector{T}) where T<:AbstractElement
    L₂Norm_Δp²= 0.0
    L₂Norm_p̄² = 0.0
    for ap in aps
        Δp², p̄² = L₂𝑝(ap)
        L₂Norm_Δp² += Δp²
        L₂Norm_p̄²  += p̄²
    end
    return (L₂Norm_Δp²/L₂Norm_p̄²)^0.5
end

"""
    Hₑ(ap; E, ν)        -> (Hₑ, L₂)``        # per element
    Hₑ(aps; E, ν)       -> (Hₑ, L₂)          # over a set of elements

Relative energy-norm error against an exact displacement field.

    ‖e‖_E / ‖u‖_E,   ‖v‖_E² = ∫_Ω ε(v):C:ε(v) dΩ

and the relative L₂ error of the displacement. The exact field is supplied with

    prescribe!(el, :u₁ => f₁, :u₂ => f₂, :u₃ => f₃)      # values
    prescribe!(el, :∂u₁∂x => g, ...)                     # gradients

Both `Hₑ` and `L₂` are returned as relative errors, i.e. already divided by the
exact norm, so a converging discretisation approaches 0.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set∇𝝭!` and `set𝝭!` — reads `ξ[:𝝭]`, `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ[:∂𝝭∂z]`,
`ξ.𝑤`, plus the exact field entries `ξ.uᵢ`, `ξ.∂uᵢ∂x`, ...
"""
function Hₑ(ap::T; E::Real, ν::Real) where T<:AbstractElement
    ΔW²= 0
    W̄² = 0
    Δu²= 0
    ū² = 0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z]
        Cᵢᵢᵢᵢ = E*(1-ν)/(1-2*ν)/(1+ν)
        Cᵢᵢⱼⱼ = E*ν/(1-2*ν)/(1+ν)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        ū₁ = ξ.u₁
        ū₂ = ξ.u₂
        ū₃ = ξ.u₃
        ∂ū₁∂x = ξ.∂u₁∂x
        ∂ū₁∂y = ξ.∂u₁∂y
        ∂ū₁∂z = ξ.∂u₁∂z
        ∂ū₂∂x = ξ.∂u₂∂x
        ∂ū₂∂y = ξ.∂u₂∂y
        ∂ū₂∂z = ξ.∂u₂∂z
        ∂ū₃∂x = ξ.∂u₃∂x
        ∂ū₃∂y = ξ.∂u₃∂y
        ∂ū₃∂z = ξ.∂u₃∂z
        ε̄₁₁ = ∂ū₁∂x
        ε̄₂₂ = ∂ū₂∂y
        ε̄₃₃ = ∂ū₃∂z
        ε̄₁₂ = ∂ū₁∂y + ∂ū₂∂x
        ε̄₁₃ = ∂ū₁∂z + ∂ū₃∂x
        ε̄₂₃ = ∂ū₂∂z + ∂ū₃∂y
        σ̄₁₁ = Cᵢᵢᵢᵢ*ε̄₁₁ + Cᵢᵢⱼⱼ*ε̄₂₂ + Cᵢᵢⱼⱼ*ε̄₃₃
        σ̄₂₂ = Cᵢᵢⱼⱼ*ε̄₁₁ + Cᵢᵢᵢᵢ*ε̄₂₂ + Cᵢᵢⱼⱼ*ε̄₃₃ 
        σ̄₃₃ = Cᵢᵢⱼⱼ*ε̄₁₁ + Cᵢᵢⱼⱼ*ε̄₂₂ + Cᵢᵢᵢᵢ*ε̄₃₃ 
        σ̄₁₂ = Cᵢⱼᵢⱼ*ε̄₁₂
        σ̄₁₃ = Cᵢⱼᵢⱼ*ε̄₁₃
        σ̄₂₃ = Cᵢⱼᵢⱼ*ε̄₂₃
        u₁ = 0.
        u₂ = 0.
        u₃ = 0.
        ε₁₁ = 0.
        ε₂₂ = 0.
        ε₃₃ = 0.
        ε₁₂ = 0.
        ε₁₃ = 0.
        ε₂₃ = 0.
        for (i,xᵢ) in enumerate(ap.𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            u₃ += N[i]*xᵢ.d₃   
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₃₃ += B₃[i]*xᵢ.d₃
            ε₁₂ += B₂[i]*xᵢ.d₁ + B₁[i]*xᵢ.d₂
            ε₁₃ += B₃[i]*xᵢ.d₁ + B₁[i]*xᵢ.d₃
            ε₂₃ += B₃[i]*xᵢ.d₂ + B₂[i]*xᵢ.d₃
        end
        σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂ + Cᵢᵢⱼⱼ*ε₃₃
        σ₂₂ = Cᵢᵢⱼⱼ*ε₁₁ + Cᵢᵢᵢᵢ*ε₂₂ + Cᵢᵢⱼⱼ*ε₃₃ 
        σ₃₃ = Cᵢᵢⱼⱼ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂ + Cᵢᵢᵢᵢ*ε₃₃ 
        σ₁₂ = Cᵢⱼᵢⱼ*ε₁₂
        σ₁₃ = Cᵢⱼᵢⱼ*ε₁₃
        σ₂₃ = Cᵢⱼᵢⱼ*ε₂₃
        ΔW² += 0.5*((σ₁₁-σ̄₁₁)*(ε₁₁-ε̄₁₁) + (σ₂₂-σ̄₂₂)*(ε₂₂-ε̄₂₂) + (σ₃₃-σ̄₃₃)*(ε₃₃-ε̄₃₃) + (σ₁₂-σ̄₁₂)*(ε₁₂-ε̄₁₂) + (σ₁₃-σ̄₁₃)*(ε₁₃-ε̄₁₃) + (σ₂₃-σ̄₂₃)*(ε₂₃-ε̄₂₃))*𝑤
        W̄² += 0.5*(σ̄₁₁*ε̄₁₁ + σ̄₂₂*ε̄₂₂ + σ̄₃₃*ε̄₃₃ + σ̄₁₂*ε̄₁₂ + σ̄₁₃*ε̄₁₃ + σ̄₂₃*ε̄₂₃)*𝑤
        Δu² += ((u₁ - ū₁)^2 + (u₂ - ū₂)^2 + (u₃ - ū₃)^2)*𝑤
        ū² += (ū₁^2 + ū₂^2 + ū₃^2)*𝑤
    end
    return ΔW², W̄², Δu², ū²
end

"""
    Hₑ(aps::Vector{T})

Element-set overload of `Hₑ`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.


    ‖e‖_E / ‖u‖_E,   ‖v‖_E² = Σ_e ∫_Ω ε(v):C:ε(v) dΩ

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function Hₑ(aps::Vector{T}) where T<:AbstractElement
    HₑNorm_ΔW²= 0.0
    HₑNorm_W̄² = 0.0
    L₂Norm_Δu²= 0.0
    L₂Norm_ū² = 0.0
    for ap in aps
        ΔW², W̄², Δu², ū² = Hₑ(ap)
        HₑNorm_ΔW² += ΔW²
        HₑNorm_W̄²  += W̄²
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end

"""
    Hₑ_PlaneStress(ap; E, ν)       # per element
    Hₑ_PlaneStress(aps; E, ν)      # over a set of elements

Relative energy-norm error for 2-D **plane stress** problems.

    ‖e‖_E / ‖u‖_E,   C₁₁₁₁ = E/(1-ν²),  C₁₁₂₂ = Eν/(1-ν²),  C₁₂₁₂ = E/(2(1+ν))

See [`Hₑ`](@ref) for the element requirements and how the exact field is
supplied.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio
"""
function Hₑ_PlaneStress(ap::T; E::Real, ν::Real) where T<:AbstractElement
    ΔW²= 0
    W̄² = 0
    Δu²= 0
    ū² = 0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        Cᵢᵢᵢᵢ = E/(1-ν^2)
        Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        ū₁ = ξ.u
        ū₂ = ξ.v
        ∂ū₁∂x = ξ.∂u∂x
        ∂ū₁∂y = ξ.∂u∂y
        ∂ū₂∂x = ξ.∂v∂x
        ∂ū₂∂y = ξ.∂v∂y
        ε̄₁₁ = ∂ū₁∂x
        ε̄₂₂ = ∂ū₂∂y
        ε̄₁₂ = ∂ū₁∂y + ∂ū₂∂x
        σ̄₁₁ = Cᵢᵢᵢᵢ*ε̄₁₁ + Cᵢᵢⱼⱼ*ε̄₂₂
        σ̄₂₂ = Cᵢᵢⱼⱼ*ε̄₁₁ + Cᵢᵢᵢᵢ*ε̄₂₂ 
        σ̄₁₂ = Cᵢⱼᵢⱼ*ε̄₁₂
        u₁ = 0.
        u₂ = 0.
        ε₁₁ = 0.
        ε₂₂ = 0.
        ε₁₂ = 0.
        for (i,xᵢ) in enumerate(ap.𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += B₂[i]*xᵢ.d₁ + B₁[i]*xᵢ.d₂
        end
        σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂
        σ₂₂ = Cᵢᵢⱼⱼ*ε₁₁ + Cᵢᵢᵢᵢ*ε₂₂ 
        σ₁₂ = Cᵢⱼᵢⱼ*ε₁₂
        ΔW² += 0.5*((σ₁₁-σ̄₁₁)*(ε₁₁-ε̄₁₁) + (σ₂₂-σ̄₂₂)*(ε₂₂-ε̄₂₂) + (σ₁₂-σ̄₁₂)*(ε₁₂-ε̄₁₂))*𝑤
        W̄² += 0.5*(σ̄₁₁*ε̄₁₁ + σ̄₂₂*ε̄₂₂ + σ̄₁₂*ε̄₁₂)*𝑤
        Δu² += ((u₁ - ū₁)^2 + (u₂ - ū₂)^2)*𝑤
        ū² += (ū₁^2 + ū₂^2)*𝑤
    end
    return ΔW², W̄², Δu², ū²
end

"""
    Hₑ_PlaneStress(aps::Vector{T})

Element-set overload of `Hₑ_PlaneStress`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.


    ‖e‖_E / ‖u‖_E,   C₁₁₁₁ = E/(1-ν²),   C₁₁₂₂ = Eν/(1-ν²),   C₁₂₁₂ = E/(2(1+ν))

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function Hₑ_PlaneStress(aps::Vector{T}) where T<:AbstractElement
    HₑNorm_ΔW²= 0.0
    HₑNorm_W̄² = 0.0
    L₂Norm_Δu²= 0.0
    L₂Norm_ū² = 0.0
    for ap in aps
        ΔW², W̄², Δu², ū² = Hₑ_PlaneStress(ap)
        HₑNorm_ΔW² += ΔW²
        HₑNorm_W̄²  += W̄²
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end

"""
    Hₑ_PlaneStrain_Deviatoric(ap; E, ν)      # per element
    Hₑ_PlaneStrain_Deviatoric(aps; E, ν)     # over a set of elements

Relative error of the **deviatoric** part of the energy norm.

    ΔW² = ½ (σ-σ̄):(ε-ε̄),   σ = Cᵈ ε,   Cᵈ = E/(1+ν),   Pᵈ = I - ⅓ I⊗I

Used to separate discretisation error in the shear/volumetric split, e.g. for
nearly-incompressible plane-strain problems. For the full energy norm including
the volumetric part use [`Hₑ`](@ref).

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]` and the exact gradients `ξ.∂u∂x`, ...
"""
function Hₑ_PlaneStrain_Deviatoric(ap::T; E::Real, ν::Real) where T<:AbstractElement
    ΔW²= 0
    W̄² = 0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        Cᵈ = E/(1+ν)
        ∂ū₁∂x = ξ.∂u∂x
        ∂ū₁∂y = ξ.∂u∂y
        ∂ū₂∂x = ξ.∂v∂x
        ∂ū₂∂y = ξ.∂v∂y
        ε̄₁₁ = ∂ū₁∂x
        ε̄₂₂ = ∂ū₂∂y
        ε̄₁₂ = ∂ū₁∂y + ∂ū₂∂x
        σ̄₁₁ = Cᵈ*( 2/3*ε̄₁₁ - 1/3*ε̄₂₂)
        σ̄₂₂ = Cᵈ*(-1/3*ε̄₁₁ + 2/3*ε̄₂₂)
        σ̄₁₂ = Cᵈ*ε̄₁₂/2
        ε₁₁ = 0.
        ε₂₂ = 0.
        ε₁₂ = 0.
        for (i,xᵢ) in enumerate(ap.𝓒)
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += B₂[i]*xᵢ.d₁ + B₁[i]*xᵢ.d₂
        end
        σ₁₁ = Cᵈ*( 2/3*ε₁₁ - 1/3*ε₂₂)
        σ₂₂ = Cᵈ*(-1/3*ε₁₁ + 2/3*ε₂₂)
        σ₁₂ = Cᵈ*ε₁₂/2
        ΔW² += 0.5*((σ₁₁-σ̄₁₁)*(ε₁₁-ε̄₁₁) + (σ₂₂-σ̄₂₂)*(ε₂₂-ε̄₂₂) + (σ₁₂-σ̄₁₂)*(ε₁₂-ε̄₁₂))*𝑤
        W̄² += 0.5*(σ̄₁₁*ε̄₁₁ + σ̄₂₂*ε̄₂₂ + σ̄₁₂*ε̄₁₂)*𝑤
    end
    return ΔW², W̄²
end

"""
    Hₑ_PlaneStrain_Deviatoric(aps::Vector{T})

Element-set overload of `Hₑ_PlaneStrain_Deviatoric`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.


    ‖e‖_E / ‖u‖_E with the deviatoric energy above, summed over the elements.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function Hₑ_PlaneStrain_Deviatoric(aps::Vector{T}) where T<:AbstractElement
    HₑNorm_ΔW²= 0.0
    HₑNorm_W̄² = 0.0
    for ap in aps
        ΔW², W̄² = Hₑ_PlaneStrain_Deviatoric(ap)
        HₑNorm_ΔW² += ΔW²
        HₑNorm_W̄²  += W̄²
    end
    return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5
end


"""
    Hₑ_up_fem(apu; E, ν)
    Hₑ_up_fem(apsu; E, ν)

Relative errors of a displacement–pressure formulation with the **FEM**
(continuous) pressure approximation.

Returns six relative errors for a set of elements:

    (Hₑ, L₂(u), Hₑ_dil, Hₑ_dev, L₂(p))

`Hₑ_dil` isolates the volumetric term `3(p-p̄)²/(2K)` and `Hₑ_dev` the
deviatoric remainder, so a locking or pressure-oscillation problem can be
attributed to one of them. The pressure is recovered elementwise from the
displacement through `p = K(∇·u)`, so no pressure element set is needed.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set∇𝝭!` and `set𝝭!` — reads `ξ[:𝝭]`, `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤` and the
exact field entries `ξ.u`, `ξ.v`, `ξ.∂u∂x`, ...
"""
function Hₑ_up_fem(apu::T; E::Real, ν::Real) where {T<:AbstractElement,S<:AbstractElement}
    ΔW²= 0
    W̄² = 0
    ΔW²_dil= 0.0
    ΔW²_dev= 0.0
    W̄²_dil = 0.0
    W̄²_dev = 0.0
    Δu²= 0
    Δp²= 0
    ū² = 0
    p̄² =0
    p² =0
    K=E/3/(1-2ν)
    G=E/2/(1+ν)
    𝓖u = apu.𝓖
    for ξu in 𝓖u
    𝑤 = ξu.𝑤
        Nᵘ= ξu[:𝝭]
        B₁ = ξu[:∂𝝭∂x]
        B₂ = ξu[:∂𝝭∂y] 
        ū₁ = ξu.u
        ū₂ = ξu.v
        ∂ū₁∂x = ξu.∂u∂x
        ∂ū₁∂y = ξu.∂u∂y
        ∂ū₂∂x = ξu.∂v∂x
        ∂ū₂∂y = ξu.∂v∂y
        
        ε̄₁₁ = ∂ū₁∂x
        ε̄₂₂ = ∂ū₂∂y
        ε̄₁₂ = ∂ū₁∂y + ∂ū₂∂x
        σ̄₁₁ = K*(ε̄₁₁+ε̄₂₂)+2G*(2/3*ε̄₁₁-1/3*ε̄₂₂)
        σ̄₂₂ = K*(ε̄₁₁+ε̄₂₂)+2G*(-1/3*ε̄₁₁+2/3*ε̄₂₂)
        σ̄₁₂ = 2G*ε̄₁₂
        σ̄₃₃ = K*(ε̄₁₁+ε̄₂₂)+2G*(-1/3*ε̄₁₁-1/3*ε̄₂₂)

        p̄   = K*(ε̄₁₁+ε̄₂₂)
        s̄₁₁ = 2G*(2/3*ε̄₁₁-1/3*ε̄₂₂)
        s̄₂₂ = 2G*(-1/3*ε̄₁₁+2/3*ε̄₂₂)
        s̄₁₂ = 2G*ε̄₁₂
        s̄₃₃ = 2G*(-1/3*ε̄₁₁-1/3*ε̄₂₂)
        u₁ = 0.
        u₂ = 0.
        ε₁₁ = 0.
        ε₂₂ = 0.
        ε₁₂ = 0.
        for (i,xᵢ) in enumerate(apu.𝓒)
            
            u₁ += Nᵘ[i]*xᵢ.d₁
            u₂ += Nᵘ[i]*xᵢ.d₂
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += B₂[i]*xᵢ.d₁ + B₁[i]*xᵢ.d₂
        end
        p = K*(ε₁₁+ε₂₂)
        s₁₁ = 2G*(2/3*ε₁₁-1/3*ε₂₂)
        s₂₂ = 2G*(-1/3*ε₁₁+2/3*ε₂₂)
        s₁₂ = 2G*ε₁₂
        s₃₃ = 2G*(-1/3*ε₁₁-1/3*ε₂₂)
        ΔW² += (3*(p-p̄)^2/2/K+(s₁₁-s̄₁₁)^2/4/G+(s₂₂-s̄₂₂)^2/4/G+(s₃₃-s̄₃₃)^2/4/G+(s₁₂-s̄₁₂)^2/4/G)*𝑤
        W̄² += (3*p̄^2/2/K + s̄₁₁^2/4/G + s̄₂₂^2/4/G + s̄₃₃^2/4/G + s̄₁₂^2/4/G)*𝑤
        # ΔW² += (3*(p-p̄)*(ε₁₁-ε̄₁₁+ε₂₂-ε̄₂₂)/2+(s₁₁-s̄₁₁)*(2/3*(ε₁₁-ε̄₁₁)-1/3*(ε₂₂-ε̄₂₂))/2+(s₂₂-s̄₂₂)*(-1/3*(ε₁₁-ε̄₁₁)+2/3*(ε₂₂-ε̄₂₂))/2+(s₃₃-s̄₃₃)*(-1/3*(ε₁₁-ε̄₁₁)-1/3*(ε₂₂-ε̄₂₂))/2+(s₁₂-s̄₁₂)*(ε₁₂-ε̄₁₂)/2)*𝑤
        # W̄² += (3*p̄*(ε̄₁₁+ε̄₂₂)/2+s̄₁₁*(2/3*ε̄₁₁-1/3*ε̄₂₂)/2+s̄₂₂*(-1/3*ε̄₁₁+2/3*ε̄₂₂)/2+s̄₃₃*(-1/3*ε̄₁₁-1/3*ε̄₂₂)/2+s̄₁₂*ε̄₁₂/2)*𝑤
        Δu² += ((u₁ - ū₁)^2 + (u₂ - ū₂)^2)*𝑤
        ū² += (ū₁^2 + ū₂^2)*𝑤
        Δp² += ((p - p̄)^2 )*𝑤
        p̄² += (p̄^2)*𝑤
        # ΔW²_dil += (3*(p-p̄)*(ε₁₁-ε̄₁₁+ε₂₂-ε̄₂₂)/2)*𝑤
        ΔW²_dil +=(3*(p-p̄)^2/2/K)*𝑤
        # ΔW²_dev += ((s₁₁-s̄₁₁)*(2/3*(ε₁₁-ε̄₁₁)-1/3*(ε₂₂-ε̄₂₂))/2+(s₂₂-s̄₂₂)*(-1/3*(ε₁₁-ε̄₁₁)+2/3*(ε₂₂-ε̄₂₂))/2+(s₃₃-s̄₃₃)*(-1/3*(ε₁₁-ε̄₁₁)-1/3*(ε₂₂-ε̄₂₂))/2+(s₁₂-s̄₁₂)*(ε₁₂-ε̄₁₂)/2)*𝑤
        ΔW²_dev +=((s₁₁-s̄₁₁)^2/4/G+(s₂₂-s̄₂₂)^2/4/G+(s₃₃-s̄₃₃)^2/4/G+(s₁₂-s̄₁₂)^2/4/G)*𝑤
        W̄²_dil +=(3*p̄^2/2/K)*𝑤
        W̄²_dev +=(s̄₁₁^2/4/G + s̄₂₂^2/4/G + s̄₃₃^2/4/G + s̄₁₂^2/4/G)*𝑤
        # p²  +=p^2*𝑤
        # println( p̄²)
        # println( p²)
        

    end
    return ΔW², W̄², Δu², ū²,ΔW²_dil, ΔW²_dev, W̄²_dil,W̄²_dev,Δp²,p̄²
end

"""
    Hₑ_up_mix(apu, app; E, ν)
    Hₑ_up_mix(apsu, apsp; E, ν)

Relative errors of a displacement–pressure formulation with an independent
**mixed** pressure approximation.

Same set of returns as [`Hₑ_up_fem`](@ref), but the pressure comes from its own
element set `app` (one pressure element zipped with each displacement element),
so the discrete pressure is not tied to `∇·u`.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio
"""
function Hₑ_up_mix(apu::T,app::S; E::Real, ν::Real) where {T<:AbstractElement,S<:AbstractElement}
    ΔW²= 0
    W̄² = 0
    ΔW²_dil= 0.0
    ΔW²_dev= 0.0
    W̄²_dil = 0.0
    W̄²_dev = 0.0
    Δu²= 0
    Δp²= 0
    ū² = 0
    p̄² =0
    p² =0
    K=E/3/(1-2ν)
    G=E/2/(1+ν)
    𝓖u = apu.𝓖
    𝓖p = app.𝓖
    for (ξu,ξp) in zip(𝓖u,𝓖p)
    𝑤 = ξu.𝑤
        Nᵖ= ξp[:𝝭]
        Nᵘ= ξu[:𝝭]
        B₁ = ξu[:∂𝝭∂x]
        B₂ = ξu[:∂𝝭∂y] 
        ū₁ = ξu.u
        ū₂ = ξu.v
        ∂ū₁∂x = ξu.∂u∂x
        ∂ū₁∂y = ξu.∂u∂y
        ∂ū₂∂x = ξu.∂v∂x
        ∂ū₂∂y = ξu.∂v∂y
        
        ε̄₁₁ = ∂ū₁∂x
        ε̄₂₂ = ∂ū₂∂y
        ε̄₁₂ = ∂ū₁∂y + ∂ū₂∂x
        σ̄₁₁ = K*(ε̄₁₁+ε̄₂₂)+2G*(2/3*ε̄₁₁-1/3*ε̄₂₂)
        σ̄₂₂ = K*(ε̄₁₁+ε̄₂₂)+2G*(-1/3*ε̄₁₁+2/3*ε̄₂₂)
        σ̄₁₂ = 2G*ε̄₁₂
        σ̄₃₃ = K*(ε̄₁₁+ε̄₂₂)+2G*(-1/3*ε̄₁₁-1/3*ε̄₂₂)

        p̄   = K*(ε̄₁₁+ε̄₂₂)
        s̄₁₁ = 2G*(2/3*ε̄₁₁-1/3*ε̄₂₂)
        s̄₂₂ = 2G*(-1/3*ε̄₁₁+2/3*ε̄₂₂)
        s̄₁₂ = 2G*ε̄₁₂
        s̄₃₃ = 2G*(-1/3*ε̄₁₁-1/3*ε̄₂₂)
        p  = 0.
        u₁ = 0.
        u₂ = 0.
        ε₁₁ = 0.
        ε₂₂ = 0.
        ε₁₂ = 0.
        for (i,xᵢ) in enumerate(apu.𝓒)
            
            u₁ += Nᵘ[i]*xᵢ.d₁
            u₂ += Nᵘ[i]*xᵢ.d₂
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += B₂[i]*xᵢ.d₁ + B₁[i]*xᵢ.d₂
        end
        for (j,xⱼ) in enumerate(app.𝓒)
            p  += Nᵖ[j]*xⱼ.q
        end
        s₁₁ = 2G*(2/3*ε₁₁-1/3*ε₂₂)
        s₂₂ = 2G*(-1/3*ε₁₁+2/3*ε₂₂)
        s₁₂ = 2G*ε₁₂
        s₃₃ = 2G*(-1/3*ε₁₁-1/3*ε₂₂)
        ΔW² += (3*(p-p̄)^2/2/K+(s₁₁-s̄₁₁)^2/4/G+(s₂₂-s̄₂₂)^2/4/G+(s₃₃-s̄₃₃)^2/4/G+(s₁₂-s̄₁₂)^2/4/G)*𝑤
        W̄² += (3*p̄^2/2/K + s̄₁₁^2/4/G + s̄₂₂^2/4/G + s̄₃₃^2/4/G + s̄₁₂^2/4/G)*𝑤
        # ΔW² += (3*(p-p̄)*(ε₁₁-ε̄₁₁+ε₂₂-ε̄₂₂)/2+(s₁₁-s̄₁₁)*(2/3*(ε₁₁-ε̄₁₁)-1/3*(ε₂₂-ε̄₂₂))/2+(s₂₂-s̄₂₂)*(-1/3*(ε₁₁-ε̄₁₁)+2/3*(ε₂₂-ε̄₂₂))/2+(s₃₃-s̄₃₃)*(-1/3*(ε₁₁-ε̄₁₁)-1/3*(ε₂₂-ε̄₂₂))/2+(s₁₂-s̄₁₂)*(ε₁₂-ε̄₁₂)/2)*𝑤
        # W̄² += (3*p̄*(ε̄₁₁+ε̄₂₂)/2+s̄₁₁*(2/3*ε̄₁₁-1/3*ε̄₂₂)/2+s̄₂₂*(-1/3*ε̄₁₁+2/3*ε̄₂₂)/2+s̄₃₃*(-1/3*ε̄₁₁-1/3*ε̄₂₂)/2+s̄₁₂*ε̄₁₂/2)*𝑤
        Δu² += ((u₁ - ū₁)^2 + (u₂ - ū₂)^2)*𝑤
        ū² += (ū₁^2 + ū₂^2)*𝑤
        Δp² += ((p - p̄)^2 )*𝑤
        p̄² += (p̄^2)*𝑤
        # ΔW²_dil += (3*(p-p̄)*(ε₁₁-ε̄₁₁+ε₂₂-ε̄₂₂)/2)*𝑤
        ΔW²_dil +=(3*(p-p̄)^2/2/K)*𝑤
        # ΔW²_dev += ((s₁₁-s̄₁₁)*(2/3*(ε₁₁-ε̄₁₁)-1/3*(ε₂₂-ε̄₂₂))/2+(s₂₂-s̄₂₂)*(-1/3*(ε₁₁-ε̄₁₁)+2/3*(ε₂₂-ε̄₂₂))/2+(s₃₃-s̄₃₃)*(-1/3*(ε₁₁-ε̄₁₁)-1/3*(ε₂₂-ε̄₂₂))/2+(s₁₂-s̄₁₂)*(ε₁₂-ε̄₁₂)/2)*𝑤
        ΔW²_dev +=((s₁₁-s̄₁₁)^2/4/G+(s₂₂-s̄₂₂)^2/4/G+(s₃₃-s̄₃₃)^2/4/G+(s₁₂-s̄₁₂)^2/4/G)*𝑤
        W̄²_dil +=(3*p̄^2/2/K)*𝑤
        W̄²_dev +=(s̄₁₁^2/4/G + s̄₂₂^2/4/G + s̄₃₃^2/4/G + s̄₁₂^2/4/G)*𝑤
        # p²  +=p^2*𝑤
        # println( p̄²)
        # println( p²)
        

    end
    return ΔW², W̄², Δu², ū²,ΔW²_dil, ΔW²_dev, W̄²_dil,W̄²_dev,Δp²,p̄²
end

"""
    Hₑ_up_fem(apsu::Vector{T})

Element-set overload of `Hₑ_up_fem`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.


    Element-set overload of `Hₑ_up_fem(apu)`, summed over the elements.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function Hₑ_up_fem(apsu::Vector{T}) where {T<:AbstractElement,S<:AbstractElement}
    HₑNorm_ΔW²= 0.0
    HₑNorm_W̄² = 0.0
    L₂Norm_Δp²= 0.0
    L₂Norm_p̄² = 0.0
    L₂Norm_Δu²= 0.0
    L₂Norm_ū² = 0.0
    HₑNorm_ΔW²_dil= 0.0
    HₑNorm_ΔW²_dev= 0.0
    HₑNorm_W̄²_dil = 0.0
    HₑNorm_W̄²_dev = 0.0
    for apu in apsu
    # for apu in apsu
        # for app in apsp
            ΔW², W̄², Δu², ū²,ΔW²_dil, ΔW²_dev, W̄²_dil,W̄²_dev, Δp²,p̄²  = Hₑ_up_fem(apu)
            HₑNorm_ΔW² += ΔW²
            HₑNorm_W̄²  += W̄²
            L₂Norm_Δu² += Δu²
            L₂Norm_ū²  += ū²
            L₂Norm_Δp² += Δp²
            L₂Norm_p̄²  += p̄²
            
            HₑNorm_ΔW²_dil += ΔW²_dil
            HₑNorm_ΔW²_dev += ΔW²_dev
            HₑNorm_W̄²_dil  += W̄²_dil
            HₑNorm_W̄²_dev  += W̄²_dev
        # end
    end
    return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5,(HₑNorm_ΔW²_dil/HₑNorm_W̄²_dil)^0.5,(HₑNorm_ΔW²_dev/HₑNorm_W̄²_dev)^0.5,(L₂Norm_Δp²/L₂Norm_p̄²)^0.5
    # return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end    

"""
    Hₑ_up_mix(apsu::Vector{T},apsp::Vector{S})

Element-set overload of `Hₑ_up_mix`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.


    Element-set overload of `Hₑ_up_mix(apu, app)`, summed over the elements.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function Hₑ_up_mix(apsu::Vector{T},apsp::Vector{S}) where {T<:AbstractElement,S<:AbstractElement}
    HₑNorm_ΔW²= 0.0
    HₑNorm_W̄² = 0.0
    L₂Norm_Δp²= 0.0
    L₂Norm_p̄² = 0.0
    L₂Norm_Δu²= 0.0
    L₂Norm_ū² = 0.0
    HₑNorm_ΔW²_dil= 0.0
    HₑNorm_ΔW²_dev= 0.0
    HₑNorm_W̄²_dil = 0.0
    HₑNorm_W̄²_dev = 0.0
    for (apu,app) in zip(apsu,apsp)
    # for apu in apsu
        # for app in apsp
            ΔW², W̄², Δu², ū²,ΔW²_dil, ΔW²_dev, W̄²_dil,W̄²_dev, Δp²,p̄²  = Hₑ_up_mix(apu,app)
            HₑNorm_ΔW² += ΔW²
            HₑNorm_W̄²  += W̄²
            L₂Norm_Δu² += Δu²
            L₂Norm_ū²  += ū²
            L₂Norm_Δp² += Δp²
            L₂Norm_p̄²  += p̄²
            
            HₑNorm_ΔW²_dil += ΔW²_dil
            HₑNorm_ΔW²_dev += ΔW²_dev
            HₑNorm_W̄²_dil  += W̄²_dil
            HₑNorm_W̄²_dev  += W̄²_dev
        # end
    end
    return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5,(HₑNorm_ΔW²_dil/HₑNorm_W̄²_dil)^0.5,(HₑNorm_ΔW²_dev/HₑNorm_W̄²_dev)^0.5,(L₂Norm_Δp²/L₂Norm_p̄²)^0.5
    # return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end    

"""
    Hₑ_up_mix_Q4P1(apu, app; E, ν)
    Hₑ_up_mix_Q4P1(apsu, apsp; E, ν)

Relative errors for the Q4/P1 pair (bilinear displacement, elementwise-constant
pressure).

Identical to [`Hₑ_up_mix`](@ref) except that the discrete pressure is taken as
the value at the first quadrature point only — the constant per element that
P1 provides — instead of being re-interpolated at every point.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio
"""
function Hₑ_up_mix_Q4P1(apu::T,app::S; E::Real, ν::Real) where {T<:AbstractElement,S<:AbstractElement}
    ΔW²= 0
    W̄² = 0
    ΔW²_dil= 0.0
    ΔW²_dev= 0.0
    W̄²_dil = 0.0
    W̄²_dev = 0.0
    Δu²= 0
    Δp²= 0
    ū² = 0
    p̄² =0
    p² =0
    K=E/3/(1-2ν)
    G=E/2/(1+ν)
    𝓖u = apu.𝓖
    𝓖p = app.𝓖
    
    for (ξu,ξp) in zip(𝓖u,𝓖p)
    𝑤 = ξu.𝑤
        Nᵖ= ξp[:𝝭]
        Nᵘ= ξu[:𝝭]
        B₁ = ξu[:∂𝝭∂x]
        B₂ = ξu[:∂𝝭∂y] 
        ū₁ = ξu.u
        ū₂ = ξu.v
        ∂ū₁∂x = ξu.∂u∂x
        ∂ū₁∂y = ξu.∂u∂y
        ∂ū₂∂x = ξu.∂v∂x
        ∂ū₂∂y = ξu.∂v∂y
        
        ε̄₁₁ = ∂ū₁∂x
        ε̄₂₂ = ∂ū₂∂y
        ε̄₁₂ = ∂ū₁∂y + ∂ū₂∂x
        σ̄₁₁ = K*(ε̄₁₁+ε̄₂₂)+2G*(2/3*ε̄₁₁-1/3*ε̄₂₂)
        σ̄₂₂ = K*(ε̄₁₁+ε̄₂₂)+2G*(-1/3*ε̄₁₁+2/3*ε̄₂₂)
        σ̄₁₂ = 2G*ε̄₁₂
        σ̄₃₃ = K*(ε̄₁₁+ε̄₂₂)+2G*(-1/3*ε̄₁₁-1/3*ε̄₂₂)

        p̄   = K*(ε̄₁₁+ε̄₂₂)
        s̄₁₁ = 2G*(2/3*ε̄₁₁-1/3*ε̄₂₂)
        s̄₂₂ = 2G*(-1/3*ε̄₁₁+2/3*ε̄₂₂)
        s̄₁₂ = 2G*ε̄₁₂
        s̄₃₃ = 2G*(-1/3*ε̄₁₁-1/3*ε̄₂₂)
        p  = 0.
        u₁ = 0.
        u₂ = 0.
        ε₁₁ = 0.
        ε₂₂ = 0.
        ε₁₂ = 0.
        for (i,xᵢ) in enumerate(apu.𝓒)
            
            u₁ += Nᵘ[i]*xᵢ.d₁
            u₂ += Nᵘ[i]*xᵢ.d₂
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += B₂[i]*xᵢ.d₁ + B₁[i]*xᵢ.d₂
        end
        
        p  = p
        
        s₁₁ = 2G*(2/3*ε₁₁-1/3*ε₂₂)
        s₂₂ = 2G*(-1/3*ε₁₁+2/3*ε₂₂)
        s₁₂ = 2G*ε₁₂
        s₃₃ = 2G*(-1/3*ε₁₁-1/3*ε₂₂)
        ΔW² += (3*(p-p̄)^2/2/K+(s₁₁-s̄₁₁)^2/4/G+(s₂₂-s̄₂₂)^2/4/G+(s₃₃-s̄₃₃)^2/4/G+(s₁₂-s̄₁₂)^2/4/G)*𝑤
        W̄² += (3*p̄^2/2/K + s̄₁₁^2/4/G + s̄₂₂^2/4/G + s̄₃₃^2/4/G + s̄₁₂^2/4/G)*𝑤
        # ΔW² += (3*(p-p̄)*(ε₁₁-ε̄₁₁+ε₂₂-ε̄₂₂)/2+(s₁₁-s̄₁₁)*(2/3*(ε₁₁-ε̄₁₁)-1/3*(ε₂₂-ε̄₂₂))/2+(s₂₂-s̄₂₂)*(-1/3*(ε₁₁-ε̄₁₁)+2/3*(ε₂₂-ε̄₂₂))/2+(s₃₃-s̄₃₃)*(-1/3*(ε₁₁-ε̄₁₁)-1/3*(ε₂₂-ε̄₂₂))/2+(s₁₂-s̄₁₂)*(ε₁₂-ε̄₁₂)/2)*𝑤
        # W̄² += (3*p̄*(ε̄₁₁+ε̄₂₂)/2+s̄₁₁*(2/3*ε̄₁₁-1/3*ε̄₂₂)/2+s̄₂₂*(-1/3*ε̄₁₁+2/3*ε̄₂₂)/2+s̄₃₃*(-1/3*ε̄₁₁-1/3*ε̄₂₂)/2+s̄₁₂*ε̄₁₂/2)*𝑤
        Δu² += ((u₁ - ū₁)^2 + (u₂ - ū₂)^2)*𝑤
        ū² += (ū₁^2 + ū₂^2)*𝑤
        Δp² += ((p - p̄)^2 )*𝑤
        p̄² += (p̄^2)*𝑤
        # ΔW²_dil += (3*(p-p̄)*(ε₁₁-ε̄₁₁+ε₂₂-ε̄₂₂)/2)*𝑤
        ΔW²_dil +=(3*(p-p̄)^2/2/K)*𝑤
        # ΔW²_dev += ((s₁₁-s̄₁₁)*(2/3*(ε₁₁-ε̄₁₁)-1/3*(ε₂₂-ε̄₂₂))/2+(s₂₂-s̄₂₂)*(-1/3*(ε₁₁-ε̄₁₁)+2/3*(ε₂₂-ε̄₂₂))/2+(s₃₃-s̄₃₃)*(-1/3*(ε₁₁-ε̄₁₁)-1/3*(ε₂₂-ε̄₂₂))/2+(s₁₂-s̄₁₂)*(ε₁₂-ε̄₁₂)/2)*𝑤
        ΔW²_dev +=((s₁₁-s̄₁₁)^2/4/G+(s₂₂-s̄₂₂)^2/4/G+(s₃₃-s̄₃₃)^2/4/G+(s₁₂-s̄₁₂)^2/4/G)*𝑤
        W̄²_dil +=(3*p̄^2/2/K)*𝑤
        W̄²_dev +=(s̄₁₁^2/4/G + s̄₂₂^2/4/G + s̄₃₃^2/4/G + s̄₁₂^2/4/G)*𝑤
        # p²  +=p^2*𝑤
        # println( p̄²)
        # println( p²)
        

    end
    return ΔW², W̄², Δu², ū²,ΔW²_dil, ΔW²_dev, W̄²_dil,W̄²_dev,Δp²,p̄²
end

"""
    Hₑ_up_mix_Q4P1(apsu::Vector{T},apsp::Vector{S})

Element-set overload of `Hₑ_up_mix_Q4P1`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.


    Element-set overload of `Hₑ_up_mix_Q4P1(apu, app)`.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function Hₑ_up_mix_Q4P1(apsu::Vector{T},apsp::Vector{S}) where {T<:AbstractElement,S<:AbstractElement}
    HₑNorm_ΔW²= 0.0
    HₑNorm_W̄² = 0.0
    L₂Norm_Δp²= 0.0
    L₂Norm_p̄² = 0.0
    L₂Norm_Δu²= 0.0
    L₂Norm_ū² = 0.0
    HₑNorm_ΔW²_dil= 0.0
    HₑNorm_ΔW²_dev= 0.0
    HₑNorm_W̄²_dil = 0.0
    HₑNorm_W̄²_dev = 0.0
    for (apu,app) in zip(apsu,apsp)
    # for apu in apsu
        # for app in apsp
            ΔW², W̄², Δu², ū²,ΔW²_dil, ΔW²_dev, W̄²_dil,W̄²_dev, Δp²,p̄²  = Hₑ_up_mix_Q4P1(apu,app)
            HₑNorm_ΔW² += ΔW²
            HₑNorm_W̄²  += W̄²
            L₂Norm_Δu² += Δu²
            L₂Norm_ū²  += ū²
            L₂Norm_Δp² += Δp²
            L₂Norm_p̄²  += p̄²
            
            HₑNorm_ΔW²_dil += ΔW²_dil
            HₑNorm_ΔW²_dev += ΔW²_dev
            HₑNorm_W̄²_dil  += W̄²_dil
            HₑNorm_W̄²_dev  += W̄²_dev
        # end
    end
    return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5,(HₑNorm_ΔW²_dil/HₑNorm_W̄²_dil)^0.5,(HₑNorm_ΔW²_dev/HₑNorm_W̄²_dev)^0.5,(L₂Norm_Δp²/L₂Norm_p̄²)^0.5
    # return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end    


"""
    Hₑ_Incompressible(ap; E, ν)
    Hₑ_Incompressible(aps; E, ν)    -> (Hₑ, L₂(u))

Relative errors for a purely deviatoric (incompressible-limit) material.

    ε₁₂ = ½(u,ᵧ + v,ₓ),   σ = Cᵈ ε,   Cᵈ = E/(1+ν)

Note the **half-shear convention**: this estimator uses `ε₁₂ = ½(u,ᵧ+v,ₓ)` with
`σ₁₂ = Cᵈ ε₁₂`, whereas [`Hₑ_PlaneStrain_Deviatoric`](@ref) uses the engineering
shear `ε₁₂ = u,ᵧ+v,ₓ` with `σ₁₂ = Cᵈ ε₁₂/2`. The two are **not** the same
functional: the `σ₁₁`/`σ₂₂` contributions do not rescale, so the resulting
norms differ. Kept as a separate estimator for that reason.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio
"""
function Hₑ_Incompressible(ap::T; E::Real, ν::Real) where T<:AbstractElement
    ΔW²= 0
    W̄² = 0
    Δu²= 0
    ū² = 0
    Cᵈ = E/(1+ν)
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        ū₁ = ξ.u
        ū₂ = ξ.v
        ∂ū₁∂x = ξ.∂u∂x
        ∂ū₁∂y = ξ.∂u∂y
        ∂ū₂∂x = ξ.∂v∂x
        ∂ū₂∂y = ξ.∂v∂y
        ε̄₁₁ = ∂ū₁∂x
        ε̄₂₂ = ∂ū₂∂y
        ε̄₁₂ = 0.5*(∂ū₁∂y + ∂ū₂∂x)
        σ̄₁₁ = Cᵈ*(2*ε̄₁₁ -   ε̄₂₂)/3
        σ̄₂₂ = Cᵈ*(- ε̄₁₁ + 2*ε̄₂₂)/3
        σ̄₁₂ = Cᵈ*ε̄₁₂
        u₁ = 0.
        u₂ = 0.
        ε₁₁ = 0.
        ε₂₂ = 0.
        ε₁₂ = 0.
        for (i,xᵢ) in enumerate(ap.𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += 0.5*(B₂[i]*xᵢ.d₁ + B₁[i]*xᵢ.d₂)
        end
        σ₁₁ = Cᵈ*(2*ε₁₁ -   ε₂₂)/3
        σ₂₂ = Cᵈ*(- ε₁₁ + 2*ε₂₂)/3
        σ₁₂ = Cᵈ*ε₁₂
        ΔW² += 0.5*((σ₁₁-σ̄₁₁)*(ε₁₁-ε̄₁₁) + (σ₂₂-σ̄₂₂)*(ε₂₂-ε̄₂₂) + (σ₁₂-σ̄₁₂)*(ε₁₂-ε̄₁₂))*𝑤
        W̄² += 0.5*(σ̄₁₁*ε̄₁₁ + σ̄₂₂*ε̄₂₂ + σ̄₁₂*ε̄₁₂)*𝑤
        Δu² += ((u₁ - ū₁)^2 + (u₂ - ū₂)^2)*𝑤
        ū² += (ū₁^2 + ū₂^2)*𝑤
    end
    return ΔW², W̄², Δu², ū²
end

"""
    Hₑ_Incompressible(aps::Vector{T})

Element-set overload of `Hₑ_Incompressible`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.


    Element-set overload of `Hₑ_Incompressible(ap)`.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function Hₑ_Incompressible(aps::Vector{T}) where T<:AbstractElement
    HₑNorm_ΔW²= 0.0
    HₑNorm_W̄² = 0.0
    L₂Norm_Δu²= 0.0
    L₂Norm_ū² = 0.0
    for ap in aps
        ΔW², W̄², Δu², ū² = Hₑ_Incompressible(ap)
        HₑNorm_ΔW² += ΔW²
        HₑNorm_W̄²  += W̄²
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end

end