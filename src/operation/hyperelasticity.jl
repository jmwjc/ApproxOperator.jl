module Hyperelasticity
    
using ..ApproxOperator: AbstractElement

"""
    ∫ESdx(ap, k, f; E)

1-D hyperelastic bar: tangent stiffness and internal force from the
Saint-Venant–Kirchhoff stored energy.

    F = 1 + du/dx,   E_GL = ½(F² - 1),   S = E·E_GL
    k[I,J] += Bᵢ (S + F²E) Bⱼ w
    f[I]   += Bᵢ E F E_GL w

# Keyword arguments
- `E` : Young's modulus

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ.𝑤`, and the nodal displacement `xᵢ.d`.
"""
function ∫ESdx(ap::T, k::AbstractMatrix{Float64}, f::AbstractVector{Float64}; E::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖 
    for ξ in 𝓖
        B = ξ[:∂𝝭∂x]
        𝑤 = ξ.𝑤
        F = 1.0
        for (i,xᵢ) in enumerate(𝓒)
            F += B[i]*xᵢ.d
        end
        E = 0.5*(F^2-1.0)
        S = Eᵉ*E
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += B[i]*(S+F^2*Eᵉ)*B[j]*𝑤
            end
            f[I] += B[i]*Eᵉ*F*E*𝑤
        end
    end
end

"""
    Δ∫∫EᵢⱼSᵢⱼdxdy_SaintVenantKirchhoff(ap, k; E, ν)

Incremental tangent stiffness of a 2-D
Saint-Venant–Kirchhoff material.

    E = ½(FᵀF - I),   S = C:E,   Δf = ∫ Bᵀ (C + S:∂²E/∂u²) B Δu w

The Green–Lagrange strain is used with the second Piola–Kirchhoff stress
so the tangent carries both the material and the geometric term.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`, and the nodal displacements.
"""
function Δ∫∫EᵢⱼSᵢⱼdxdy_SaintVenantKirchhoff(ap::T, k::AbstractMatrix{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        Cᵢᵢᵢᵢ = E*(1-ν)/(1+ν)/(1-2*ν)
        Cᵢᵢⱼⱼ = E*ν/(1+ν)/(1-2*ν)
        Cᵢⱼᵢⱼ = E/(1+ν)/2
       
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        
        E₁₁ = 0.5*(F₁₁*F₁₁+F₂₁*F₂₁-1.0)
        E₁₂ = 0.5*(F₁₁*F₁₂+F₂₁*F₂₂)
        E₂₂ = 0.5*(F₁₂*F₁₂+F₂₂*F₂₂-1.0)
        S₁₁ = Cᵢᵢᵢᵢ*E₁₁+Cᵢᵢⱼⱼ*E₂₂
        S₂₂ = Cᵢᵢⱼⱼ*E₁₁+Cᵢᵢᵢᵢ*E₂₂
        S₁₂ = 2.0*Cᵢⱼᵢⱼ*E₁₂
        
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] +=((B₁[i]*F₁₁*B₁[j]*F₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂)*Cᵢᵢᵢᵢ
                               +  (B₁[i]*F₁₁*B₂[j]*F₁₂ + B₂[i]*F₁₂*B₁[j]*F₁₁)*Cᵢᵢⱼⱼ
                               +  (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*Cᵢⱼᵢⱼ
                               +   B₁[i]*B₁[j]*S₁₁+B₂[i]*B₂[j]*S₂₂+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂)*𝑤
                              
                k[2*I-1,2*J]   += ((B₁[i]*F₁₁*B₁[j]*F₂₁ + B₂[i]*F₁₂*B₂[j]*F₂₂)*Cᵢᵢᵢᵢ
                               +   (B₁[i]*F₁₁*B₂[j]*F₂₂ + B₂[i]*F₁₂*B₁[j]*F₂₁)*Cᵢᵢⱼⱼ
                               +   (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*Cᵢⱼᵢⱼ)*𝑤
                               
                k[2*I,2*J-1]   += ((B₁[i]*F₂₁*B₁[j]*F₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂)*Cᵢᵢᵢᵢ
                               +   (B₁[i]*F₂₁*B₂[j]*F₁₂ + B₂[i]*F₂₂*B₁[j]*F₁₁)*Cᵢᵢⱼⱼ
                               +   (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*Cᵢⱼᵢⱼ)*𝑤

                k[2*I,2*J]     += ((B₁[i]*F₂₁*B₁[j]*F₂₁ + B₂[i]*F₂₂*B₂[j]*F₂₂)*Cᵢᵢᵢᵢ 
                               +   (B₁[i]*F₂₁*B₂[j]*F₂₂ + B₂[i]*F₂₂*B₁[j]*F₂₁)*Cᵢᵢⱼⱼ
                               +   (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*Cᵢⱼᵢⱼ
                               +    B₁[i]*B₁[j]*S₁₁+B₂[i]*B₂[j]*S₂₂+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂)*𝑤

            end
        end
    end
end

"""
    ∫∫EᵢⱼSᵢⱼdxdy_SaintVenantKirchhoff(ap, k; E, ν)

Incremental internal force of a 2-D
Saint-Venant–Kirchhoff material.

    E = ½(FᵀF - I),   S = C:E,   f[I] += ∫ Bᵀ S w

The Green–Lagrange strain is used with the second Piola–Kirchhoff stress
so the force is nonlinear in the displacement.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`, and the nodal displacements.
"""
function ∫∫EᵢⱼSᵢⱼdxdy_SaintVenantKirchhoff(ap::T, f::AbstractVector{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        Cᵢᵢᵢᵢ = E*(1-ν)/(1+ν)/(1-2*ν)
        Cᵢᵢⱼⱼ = E*ν/(1+ν)/(1-2*ν)
        Cᵢⱼᵢⱼ = E/(1+ν)/2
       
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        
        E₁₁ = 0.5*(F₁₁*F₁₁+F₂₁*F₂₁-1.0)
        E₁₂ = 0.5*(F₁₁*F₁₂+F₂₁*F₂₂)
        E₂₂ = 0.5*(F₁₂*F₁₂+F₂₂*F₂₂-1.0)
        S₁₁ = Cᵢᵢᵢᵢ*E₁₁+Cᵢᵢⱼⱼ*E₂₂
        S₂₂ = Cᵢᵢⱼⱼ*E₁₁+Cᵢᵢᵢᵢ*E₂₂
        S₁₂ = 2.0*Cᵢⱼᵢⱼ*E₁₂
        
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += (B₁[i]*F₁₁*S₁₁+B₂[i]*F₁₂*S₂₂+(B₁[i]*F₁₂+B₂[i]*F₁₁)*S₁₂)*𝑤
            f[2*I]   += (B₁[i]*F₂₁*S₁₁+B₂[i]*F₂₂*S₂₂+(B₁[i]*F₂₂+B₂[i]*F₂₁)*S₁₂)*𝑤
        end
    end
end

"""
    Δ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean(ap, k; E, ν)

Full part of the tangent of a 2-D Neo-Hookean material.

    Δf = ∫ Bᵀ (C + S:∂²E/∂u²) B Δu w

Energy:  ψ = μ/2 (I₁ - 3) - μ ln J + λ/2 (ln J)²

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`, and the nodal displacements.
"""
function Δ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean(ap::T, k::AbstractMatrix{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        C₁₁₁₁=λ*J*(2*J-1.0)*C⁻¹₁₁*C⁻¹₁₁+(μ- λ*J*(J-1.0))*2*C⁻¹₁₁*C⁻¹₁₁
        C₂₂₂₂=λ*J*(2*J-1.0)*C⁻¹₂₂*C⁻¹₂₂+(μ- λ*J*(J-1.0))*2*C⁻¹₂₂*C⁻¹₂₂
        C₁₁₂₂=λ*J*(2*J-1.0)*C⁻¹₁₁*C⁻¹₂₂+(μ- λ*J*(J-1.0))*2*C⁻¹₁₂*C⁻¹₁₂
        C₁₁₁₂=λ*J*(2*J-1.0)*C⁻¹₁₁*C⁻¹₁₂+(μ- λ*J*(J-1.0))*2*C⁻¹₁₁*C⁻¹₁₂ 
        C₂₂₁₂=λ*J*(2*J-1.0)*C⁻¹₂₂*C⁻¹₁₂+(μ- λ*J*(J-1.0))*2*C⁻¹₁₂*C⁻¹₂₂     
        C₁₂₁₂=λ*J*(2*J-1.0)*C⁻¹₁₂*C⁻¹₁₂+(μ- λ*J*(J-1.0))*(C⁻¹₁₁*C⁻¹₂₂+C⁻¹₁₂*C⁻¹₁₂)   

        S₁₁ = λ*J*(J-1)*C⁻¹₁₁ + μ*(1-C⁻¹₁₁)
        S₂₂ = λ*J*(J-1)*C⁻¹₂₂ + μ*(1-C⁻¹₂₂)
        S₁₂ = λ*J*(J-1)*C⁻¹₁₂ - μ*C⁻¹₁₂
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] +=(B₁[i]*F₁₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               + (B₁[i]*F₁₁*B₂[j]*F₁₂ + B₂[i]*F₁₂*B₁[j]*F₁₁)*C₁₁₂₂
                               + (B₁[i]*F₁₁*B₂[j]*F₁₁ + B₁[i]*F₁₁*B₁[j]*F₁₂)*C₁₁₁₂
                               + (B₁[i]*F₁₂*B₁[j]*F₁₁ + B₂[i]*F₁₁*B₁[j]*F₁₁)*C₁₁₁₂
                               + (B₂[i]*F₁₂*B₂[j]*F₁₁ + B₂[i]*F₁₂*B₁[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂*B₂[j]*F₁₂ + B₂[i]*F₁₁*B₂[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂
                               +  B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤
                              
                k[2*I-1,2*J]   += (B₁[i]*F₁₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₂₂*C₂₂₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₂ + B₂[i]*F₁₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₁ + B₁[i]*F₁₁*B₁[j]*F₂₂)*C₁₁₁₂ 
                               +  (B₁[i]*F₁₂*B₁[j]*F₂₁ + B₂[i]*F₁₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₁₂*B₂[j]*F₂₁ + B₂[i]*F₁₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂*B₂[j]*F₂₂ + B₂[i]*F₁₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂)*𝑤
                               
                k[2*I,2*J-1]   += (B₁[i]*F₂₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₂ + B₂[i]*F₂₂*B₁[j]*F₁₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₁ + B₁[i]*F₂₁*B₁[j]*F₁₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₁₁ + B₂[i]*F₂₁*B₁[j]*F₁₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₁₁ + B₂[i]*F₂₂*B₁[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₁₂ + B₂[i]*F₂₁*B₂[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂)*𝑤

                k[2*I,2*J]     += (B₁[i]*F₂₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₂₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₂ + B₂[i]*F₂₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₁ + B₁[i]*F₂₁*B₁[j]*F₂₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₂₁ + B₂[i]*F₂₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₂₁ + B₂[i]*F₂₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₂₂ + B₂[i]*F₂₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂
                               +   B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤

            end
        end
    end
end

"""
    Δ∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean(ap, k; E, ν)

Volumetric part of the tangent of a 2-D Neo-Hookean material.

    Δf = ∫ Bᵀ (C + S:∂²E/∂u²) B Δu w
Combined with the deviatoric companion it reproduces the full operator.

Energy:  ψ = μ/2 (I₁ - 3) - μ ln J + λ/2 (ln J)²

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`, and the nodal displacements.
"""
function Δ∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean(ap::T, k::AbstractMatrix{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        C₁₁₁₁=λ*J*(2*J-1.0)*C⁻¹₁₁*C⁻¹₁₁- λ*J*(J-1.0)*2*C⁻¹₁₁*C⁻¹₁₁
        C₂₂₂₂=λ*J*(2*J-1.0)*C⁻¹₂₂*C⁻¹₂₂- λ*J*(J-1.0)*2*C⁻¹₂₂*C⁻¹₂₂
        C₁₁₂₂=λ*J*(2*J-1.0)*C⁻¹₁₁*C⁻¹₂₂- λ*J*(J-1.0)*2*C⁻¹₁₂*C⁻¹₁₂
        C₁₁₁₂=λ*J*(2*J-1.0)*C⁻¹₁₁*C⁻¹₁₂- λ*J*(J-1.0)*2*C⁻¹₁₁*C⁻¹₁₂ 
        C₂₂₁₂=λ*J*(2*J-1.0)*C⁻¹₂₂*C⁻¹₁₂- λ*J*(J-1.0)*2*C⁻¹₁₂*C⁻¹₂₂     
        C₁₂₁₂=λ*J*(2*J-1.0)*C⁻¹₁₂*C⁻¹₁₂- λ*J*(J-1.0)*(C⁻¹₁₁*C⁻¹₂₂+C⁻¹₁₂*C⁻¹₁₂)   

        S₁₁ = λ*J*(J-1)*C⁻¹₁₁
        S₂₂ = λ*J*(J-1)*C⁻¹₂₂
        S₁₂ = λ*J*(J-1)*C⁻¹₁₂
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] +=(B₁[i]*F₁₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               + (B₁[i]*F₁₁*B₂[j]*F₁₂ + B₂[i]*F₁₂*B₁[j]*F₁₁)*C₁₁₂₂
                               + (B₁[i]*F₁₁*B₂[j]*F₁₁ + B₁[i]*F₁₁*B₁[j]*F₁₂)*C₁₁₁₂
                               + (B₁[i]*F₁₂*B₁[j]*F₁₁ + B₂[i]*F₁₁*B₁[j]*F₁₁)*C₁₁₁₂
                               + (B₂[i]*F₁₂*B₂[j]*F₁₁ + B₂[i]*F₁₂*B₁[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂*B₂[j]*F₁₂ + B₂[i]*F₁₁*B₂[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂
                               +  B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤
                              
                k[2*I-1,2*J]   += (B₁[i]*F₁₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₂₂*C₂₂₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₂ + B₂[i]*F₁₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₁ + B₁[i]*F₁₁*B₁[j]*F₂₂)*C₁₁₁₂ 
                               +  (B₁[i]*F₁₂*B₁[j]*F₂₁ + B₂[i]*F₁₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₁₂*B₂[j]*F₂₁ + B₂[i]*F₁₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂*B₂[j]*F₂₂ + B₂[i]*F₁₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂)*𝑤
                               
                k[2*I,2*J-1]   += (B₁[i]*F₂₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₂ + B₂[i]*F₂₂*B₁[j]*F₁₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₁ + B₁[i]*F₂₁*B₁[j]*F₁₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₁₁ + B₂[i]*F₂₁*B₁[j]*F₁₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₁₁ + B₂[i]*F₂₂*B₁[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₁₂ + B₂[i]*F₂₁*B₂[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂)*𝑤

                k[2*I,2*J]     += (B₁[i]*F₂₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₂₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₂ + B₂[i]*F₂₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₁ + B₁[i]*F₂₁*B₁[j]*F₂₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₂₁ + B₂[i]*F₂₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₂₁ + B₂[i]*F₂₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₂₂ + B₂[i]*F₂₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂
                               +   B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤

            end
        end
    end
end

"""
    Δ∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean(ap, k; E, ν)

Deviatoric part of the tangent of a 2-D Neo-Hookean material.

    Δf = ∫ Bᵀ (C + S:∂²E/∂u²) B Δu w
Combined with the volumetric companion it reproduces the full operator.

Energy:  ψ = μ/2 (I₁ - 3) - μ ln J + λ/2 (ln J)²

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`, and the nodal displacements.
"""
function Δ∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean(ap::T, k::AbstractMatrix{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        C₁₁₁₁=μ*2*C⁻¹₁₁*C⁻¹₁₁
        C₂₂₂₂=μ*2*C⁻¹₂₂*C⁻¹₂₂
        C₁₁₂₂=μ*2*C⁻¹₁₂*C⁻¹₁₂
        C₁₁₁₂=μ*2*C⁻¹₁₁*C⁻¹₁₂ 
        C₂₂₁₂=μ*2*C⁻¹₁₂*C⁻¹₂₂     
        C₁₂₁₂=μ*(C⁻¹₁₁*C⁻¹₂₂+C⁻¹₁₂*C⁻¹₁₂)   

        S₁₁ = μ*(1-C⁻¹₁₁)
        S₂₂ = μ*(1-C⁻¹₂₂)
        S₁₂ = - μ*C⁻¹₁₂
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] +=(B₁[i]*F₁₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               + (B₁[i]*F₁₁*B₂[j]*F₁₂ + B₂[i]*F₁₂*B₁[j]*F₁₁)*C₁₁₂₂
                               + (B₁[i]*F₁₁*B₂[j]*F₁₁ + B₁[i]*F₁₁*B₁[j]*F₁₂)*C₁₁₁₂
                               + (B₁[i]*F₁₂*B₁[j]*F₁₁ + B₂[i]*F₁₁*B₁[j]*F₁₁)*C₁₁₁₂
                               + (B₂[i]*F₁₂*B₂[j]*F₁₁ + B₂[i]*F₁₂*B₁[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂*B₂[j]*F₁₂ + B₂[i]*F₁₁*B₂[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂
                               +  B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤
                              
                k[2*I-1,2*J]   += (B₁[i]*F₁₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₂₂*C₂₂₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₂ + B₂[i]*F₁₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₁ + B₁[i]*F₁₁*B₁[j]*F₂₂)*C₁₁₁₂ 
                               +  (B₁[i]*F₁₂*B₁[j]*F₂₁ + B₂[i]*F₁₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₁₂*B₂[j]*F₂₁ + B₂[i]*F₁₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂*B₂[j]*F₂₂ + B₂[i]*F₁₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂)*𝑤
                               
                k[2*I,2*J-1]   += (B₁[i]*F₂₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₂ + B₂[i]*F₂₂*B₁[j]*F₁₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₁ + B₁[i]*F₂₁*B₁[j]*F₁₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₁₁ + B₂[i]*F₂₁*B₁[j]*F₁₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₁₁ + B₂[i]*F₂₂*B₁[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₁₂ + B₂[i]*F₂₁*B₂[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂)*𝑤

                k[2*I,2*J]     += (B₁[i]*F₂₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₂₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₂ + B₂[i]*F₂₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₁ + B₁[i]*F₂₁*B₁[j]*F₂₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₂₁ + B₂[i]*F₂₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₂₁ + B₂[i]*F₂₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₂₂ + B₂[i]*F₂₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂
                               +   B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤

            end
        end
    end
end

"""
    ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean(ap, k; E, ν)

Full part of the internal force of a 2-D Neo-Hookean material.

    f[I] += ∫ Bᵀ S w

Energy:  ψ = μ/2 (I₁ - 3) - μ ln J + λ/2 (ln J)²

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`, and the nodal displacements.
"""
function ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean(ap::T, f::AbstractVector{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        S₁₁ = λ*J*(J-1)*C⁻¹₁₁ + μ*(1-C⁻¹₁₁)
        S₂₂ = λ*J*(J-1)*C⁻¹₂₂ + μ*(1-C⁻¹₂₂)
        S₁₂ = λ*J*(J-1)*C⁻¹₁₂ - μ*C⁻¹₁₂
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += (B₁[i]*F₁₁*S₁₁+B₂[i]*F₁₂*S₂₂+(B₁[i]*F₁₂+B₂[i]*F₁₁)*S₁₂)*𝑤
            f[2*I]   += (B₁[i]*F₂₁*S₁₁+B₂[i]*F₂₂*S₂₂+(B₁[i]*F₂₂+B₂[i]*F₂₁)*S₁₂)*𝑤
        end
    end
end

"""
    ∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean(ap, k; E, ν)

Volumetric part of the internal force of a 2-D Neo-Hookean material.

    f[I] += ∫ Bᵀ S w
Combined with the deviatoric companion it reproduces the full operator.

Energy:  ψ = μ/2 (I₁ - 3) - μ ln J + λ/2 (ln J)²

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`, and the nodal displacements.
"""
function ∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean(ap::T, f::AbstractVector{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        S₁₁ = λ*J*(J-1)*C⁻¹₁₁
        S₂₂ = λ*J*(J-1)*C⁻¹₂₂
        S₁₂ = λ*J*(J-1)*C⁻¹₁₂
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += (B₁[i]*F₁₁*S₁₁+B₂[i]*F₁₂*S₂₂+(B₁[i]*F₁₂+B₂[i]*F₁₁)*S₁₂)*𝑤
            f[2*I]   += (B₁[i]*F₂₁*S₁₁+B₂[i]*F₂₂*S₂₂+(B₁[i]*F₂₂+B₂[i]*F₂₁)*S₁₂)*𝑤
        end
    end
end

"""
    ∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean(ap, k; E, ν)

Deviatoric part of the internal force of a 2-D Neo-Hookean material.

    f[I] += ∫ Bᵀ S w
Combined with the volumetric companion it reproduces the full operator.

Energy:  ψ = μ/2 (I₁ - 3) - μ ln J + λ/2 (ln J)²

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`, and the nodal displacements.
"""
function ∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean(ap::T, f::AbstractVector{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        S₁₁ = μ*(1-C⁻¹₁₁)
        S₂₂ = μ*(1-C⁻¹₂₂)
        S₁₂ = - μ*C⁻¹₁₂
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += (B₁[i]*F₁₁*S₁₁+B₂[i]*F₁₂*S₂₂+(B₁[i]*F₁₂+B₂[i]*F₁₁)*S₁₂)*𝑤
            f[2*I]   += (B₁[i]*F₂₁*S₁₁+B₂[i]*F₂₂*S₂₂+(B₁[i]*F₂₂+B₂[i]*F₂₁)*S₁₂)*𝑤
        end
    end
end

"""
    ∫vᵢuᵢds(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64})

Element-set overload of `∫vᵢuᵢds`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function ∫vᵢuᵢds(ap::T, k::AbstractMatrix{Float64}, f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        n₁₁ = ξ.n₁₁
        n₂₂ = ξ.n₂₂
        n₁₂ = ξ.n₁₂
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
        end
        Δu₁ = g₁-u₁
        Δu₂ = g₂-u₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += α*N[i]*n₁₁*N[j]*𝑤
                k[2*I,2*J-1]   += α*N[i]*n₁₂*N[j]*𝑤
                k[2*I-1,2*J]   += α*N[i]*n₁₂*N[j]*𝑤
                k[2*I,2*J]     += α*N[i]*n₂₂*N[j]*𝑤
            end
            f[2*I-1] += α*N[i]*(n₁₁*Δu₁+n₁₂*Δu₂)*𝑤
            f[2*I]   += α*N[i]*(n₁₂*Δu₁+n₂₂*Δu₂)*𝑤
        end
    end
end
"""
    Δ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean2(ap, k; E, ν)

Full part of the tangent of a 2-D Neo-Hookean material (the `_2` variant).

    Δf = ∫ Bᵀ (C + S:∂²E/∂u²) B Δu w

Energy:  ψ = μ/2 (I₁ - 3) - μ ln J + λ/2 (ln J)²

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`, and the nodal displacements.
"""
function Δ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean2(ap::T, k::AbstractMatrix{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    K=E/(1-2*ν)/3
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        J⁻²=(1.0/J)^2
        # println(J⁻²)
        J⁻²³=cbrt(J⁻²)
        # println(J⁻²³)
        C₁₁₁₁=(-2/3*C⁻¹₁₁-2/3*C⁻¹₁₁)*J⁻²³*μ+(J^2*K+2/9*J⁻²³*I₁*μ)*C⁻¹₁₁*C⁻¹₁₁+(1/3*μ*J⁻²³*I₁-0.5*K*(J^2-1.0))*2*C⁻¹₁₁*C⁻¹₁₁
            #   println(C₁₁₁₁)
        C₂₂₂₂=(-2/3*C⁻¹₂₂-2/3*C⁻¹₂₂)*J⁻²³*μ+(J^2*K+2/9*J⁻²³*I₁*μ)*C⁻¹₂₂*C⁻¹₂₂+(1/3*μ*J⁻²³*I₁-0.5*K*(J^2-1.0))*2*C⁻¹₂₂*C⁻¹₂₂
        C₁₁₂₂=(-2/3*C⁻¹₂₂-2/3*C⁻¹₁₁)*J⁻²³*μ+(J^2*K+2/9*J⁻²³*I₁*μ)*C⁻¹₁₁*C⁻¹₂₂+(1/3*μ*J⁻²³*I₁-0.5*K*(J^2-1.0))*2*C⁻¹₁₂*C⁻¹₁₂
        C₁₁₁₂=(-2/3*C⁻¹₁₂)*J⁻²³*μ+(J^2*K+2/9*J⁻²³*I₁*μ)*C⁻¹₁₁*C⁻¹₁₂+(1/3*μ*J⁻²³*I₁-0.5*K*(J^2-1.0))*2*C⁻¹₁₁*C⁻¹₁₂ 
        C₂₂₁₂=(-2/3*C⁻¹₁₂)*J⁻²³*μ+(J^2*K+2/9*J⁻²³*I₁*μ)*C⁻¹₂₂*C⁻¹₁₂+(1/3*μ*J⁻²³*I₁-0.5*K*(J^2-1.0))*2*C⁻¹₂₂*C⁻¹₁₂     
        C₁₂₁₂=(J^2*K+2/9*J⁻²³*I₁*μ)*C⁻¹₁₂*C⁻¹₁₂+(1/3*μ*J⁻²³*I₁-0.5*K*(J^2-1.0))*(C⁻¹₁₁*C⁻¹₂₂+C⁻¹₁₂*C⁻¹₁₂)   
              
        S₁₁ = μ*J⁻²³*(1.0-1/3*I₁*C⁻¹₁₁)+0.5*K*(J^2-1.0)*C⁻¹₁₁
        S₂₂ = μ*J⁻²³*(1.0-1/3*I₁*C⁻¹₂₂)+0.5*K*(J^2-1.0)*C⁻¹₂₂
        S₁₂ = μ*J⁻²³*(-1/3*I₁*C⁻¹₁₂)+0.5*K*(J^2-1.0)*C⁻¹₁₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] +=(B₁[i]*F₁₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               + (B₁[i]*F₁₁*B₂[j]*F₁₂ + B₂[i]*F₁₂*B₁[j]*F₁₁)*C₁₁₂₂
                               + (B₁[i]*F₁₁*B₂[j]*F₁₁ + B₁[i]*F₁₁*B₁[j]*F₁₂)*C₁₁₁₂
                               + (B₁[i]*F₁₂*B₁[j]*F₁₁ + B₂[i]*F₁₁*B₁[j]*F₁₁)*C₁₁₁₂
                               + (B₂[i]*F₁₂*B₂[j]*F₁₁ + B₂[i]*F₁₂*B₁[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂*B₂[j]*F₁₂ + B₂[i]*F₁₁*B₂[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂
                               +  B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤
                              
                k[2*I-1,2*J]   += (B₁[i]*F₁₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₂₂*C₂₂₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₂ + B₂[i]*F₁₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₁ + B₁[i]*F₁₁*B₁[j]*F₂₂)*C₁₁₁₂ 
                               +  (B₁[i]*F₁₂*B₁[j]*F₂₁ + B₂[i]*F₁₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₁₂*B₂[j]*F₂₁ + B₂[i]*F₁₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂*B₂[j]*F₂₂ + B₂[i]*F₁₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂)*𝑤
                               
                k[2*I,2*J-1]   += (B₁[i]*F₂₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₂ + B₂[i]*F₂₂*B₁[j]*F₁₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₁ + B₁[i]*F₂₁*B₁[j]*F₁₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₁₁ + B₂[i]*F₂₁*B₁[j]*F₁₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₁₁ + B₂[i]*F₂₂*B₁[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₁₂ + B₂[i]*F₂₁*B₂[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂)*𝑤

                k[2*I,2*J]     += (B₁[i]*F₂₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₂₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₂ + B₂[i]*F₂₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₁ + B₁[i]*F₂₁*B₁[j]*F₂₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₂₁ + B₂[i]*F₂₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₂₁ + B₂[i]*F₂₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₂₂ + B₂[i]*F₂₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂
                               +   B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤

            end
        end
    end
end

"""
    Δ∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean2(ap, k; E, ν)

Volumetric part of the tangent of a 2-D Neo-Hookean material (the `_2` variant).

    Δf = ∫ Bᵀ (C + S:∂²E/∂u²) B Δu w
Combined with the deviatoric companion it reproduces the full operator.

Energy:  ψ = μ/2 (I₁ - 3) - μ ln J + λ/2 (ln J)²

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`, and the nodal displacements.
"""
function Δ∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean2(ap::T, k::AbstractMatrix{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    K=E/(1-2*ν)/3
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        J⁻²=(1.0/J)^2
        J⁻²³=cbrt(J⁻²)

        C₁₁₁₁=J^2*K*C⁻¹₁₁*C⁻¹₁₁+0.5*K*(J^2-1.0)*2*C⁻¹₁₁*C⁻¹₁₁
        C₂₂₂₂=J^2*K*C⁻¹₂₂*C⁻¹₂₂+0.5*K*(J^2-1.0)*2*C⁻¹₂₂*C⁻¹₂₂
        C₁₁₂₂=J^2*K*C⁻¹₁₁*C⁻¹₂₂+0.5*K*(J^2-1.0)*2*C⁻¹₁₂*C⁻¹₁₂
        C₁₁₁₂=J^2*K*C⁻¹₁₁*C⁻¹₁₂+0.5*K*(J^2-1.0)*2*C⁻¹₁₁*C⁻¹₁₂ 
        C₂₂₁₂=J^2*K*C⁻¹₂₂*C⁻¹₁₂+0.5*K*(J^2-1.0)*2*C⁻¹₂₂*C⁻¹₁₂     
        C₁₂₁₂=J^2*K*C⁻¹₁₂*C⁻¹₁₂+0.5*K*(J^2-1.0)*(C⁻¹₁₁*C⁻¹₂₂+C⁻¹₁₂*C⁻¹₁₂)   
              
        S₁₁ = 0.5*K*(J^2-1.0)*C⁻¹₁₁
        S₂₂ = 0.5*K*(J^2-1.0)*C⁻¹₂₂
        S₁₂ = 0.5*K*(J^2-1.0)*C⁻¹₁₂
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] +=(B₁[i]*F₁₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               + (B₁[i]*F₁₁*B₂[j]*F₁₂ + B₂[i]*F₁₂*B₁[j]*F₁₁)*C₁₁₂₂
                               + (B₁[i]*F₁₁*B₂[j]*F₁₁ + B₁[i]*F₁₁*B₁[j]*F₁₂)*C₁₁₁₂
                               + (B₁[i]*F₁₂*B₁[j]*F₁₁ + B₂[i]*F₁₁*B₁[j]*F₁₁)*C₁₁₁₂
                               + (B₂[i]*F₁₂*B₂[j]*F₁₁ + B₂[i]*F₁₂*B₁[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂*B₂[j]*F₁₂ + B₂[i]*F₁₁*B₂[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂
                               +  B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤
                              
                k[2*I-1,2*J]   += (B₁[i]*F₁₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₂₂*C₂₂₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₂ + B₂[i]*F₁₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₁ + B₁[i]*F₁₁*B₁[j]*F₂₂)*C₁₁₁₂ 
                               +  (B₁[i]*F₁₂*B₁[j]*F₂₁ + B₂[i]*F₁₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₁₂*B₂[j]*F₂₁ + B₂[i]*F₁₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂*B₂[j]*F₂₂ + B₂[i]*F₁₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂)*𝑤
                               
                k[2*I,2*J-1]   += (B₁[i]*F₂₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₂ + B₂[i]*F₂₂*B₁[j]*F₁₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₁ + B₁[i]*F₂₁*B₁[j]*F₁₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₁₁ + B₂[i]*F₂₁*B₁[j]*F₁₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₁₁ + B₂[i]*F₂₂*B₁[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₁₂ + B₂[i]*F₂₁*B₂[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂)*𝑤

                k[2*I,2*J]     += (B₁[i]*F₂₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₂₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₂ + B₂[i]*F₂₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₁ + B₁[i]*F₂₁*B₁[j]*F₂₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₂₁ + B₂[i]*F₂₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₂₁ + B₂[i]*F₂₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₂₂ + B₂[i]*F₂₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂
                               +   B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤

            end
        end
    end
end
"""
    Δ∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean2(ap, k; E, ν)

Deviatoric part of the tangent of a 2-D Neo-Hookean material (the `_2` variant).

    Δf = ∫ Bᵀ (C + S:∂²E/∂u²) B Δu w
Combined with the volumetric companion it reproduces the full operator.

Energy:  ψ = μ/2 (I₁ - 3) - μ ln J + λ/2 (ln J)²

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`, and the nodal displacements.
"""
function Δ∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean2(ap::T, k::AbstractMatrix{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    K=E/(1-2*ν)/3
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        J⁻²=(1.0/J)^2
        J⁻²³=cbrt(J⁻²)
        C₁₁₁₁=(-2/3*C⁻¹₁₁-2/3*C⁻¹₁₁)*J⁻²³*μ + (2/9*J⁻²³*I₁*μ)*C⁻¹₁₁*C⁻¹₁₁ + (1/3*μ*J⁻²³*I₁)*2*C⁻¹₁₁*C⁻¹₁₁
        C₂₂₂₂=(-2/3*C⁻¹₂₂-2/3*C⁻¹₂₂)*J⁻²³*μ + (2/9*J⁻²³*I₁*μ)*C⁻¹₂₂*C⁻¹₂₂ + (1/3*μ*J⁻²³*I₁)*2*C⁻¹₂₂*C⁻¹₂₂
        C₁₁₂₂=(-2/3*C⁻¹₂₂-2/3*C⁻¹₁₁)*J⁻²³*μ + (2/9*J⁻²³*I₁*μ)*C⁻¹₁₁*C⁻¹₂₂ + (1/3*μ*J⁻²³*I₁)*2*C⁻¹₁₂*C⁻¹₁₂
        C₁₁₁₂=(-2/3*C⁻¹₁₂)*J⁻²³*μ + (2/9*J⁻²³*I₁*μ)*C⁻¹₁₁*C⁻¹₁₂ + (1/3*μ*J⁻²³*I₁)*2*C⁻¹₁₁*C⁻¹₁₂ 
        C₂₂₁₂=(-2/3*C⁻¹₁₂)*J⁻²³*μ + (2/9*J⁻²³*I₁*μ)*C⁻¹₂₂*C⁻¹₁₂ + (1/3*μ*J⁻²³*I₁)*2*C⁻¹₂₂*C⁻¹₁₂     
        C₁₂₁₂=(2/9*J⁻²³*I₁*μ)*C⁻¹₁₂*C⁻¹₁₂ + (1/3*μ*J⁻²³*I₁)*(C⁻¹₁₁*C⁻¹₂₂+C⁻¹₁₂*C⁻¹₁₂)   
              
        S₁₁ = μ*J⁻²³*(1.0-1/3*I₁*C⁻¹₁₁)
        S₂₂ = μ*J⁻²³*(1.0-1/3*I₁*C⁻¹₂₂)
        S₁₂ = μ*J⁻²³*(-1/3*I₁*C⁻¹₁₂)
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] +=(B₁[i]*F₁₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               + (B₁[i]*F₁₁*B₂[j]*F₁₂ + B₂[i]*F₁₂*B₁[j]*F₁₁)*C₁₁₂₂
                               + (B₁[i]*F₁₁*B₂[j]*F₁₁ + B₁[i]*F₁₁*B₁[j]*F₁₂)*C₁₁₁₂
                               + (B₁[i]*F₁₂*B₁[j]*F₁₁ + B₂[i]*F₁₁*B₁[j]*F₁₁)*C₁₁₁₂
                               + (B₂[i]*F₁₂*B₂[j]*F₁₁ + B₂[i]*F₁₂*B₁[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂*B₂[j]*F₁₂ + B₂[i]*F₁₁*B₂[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂
                               +  B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤
                              
                k[2*I-1,2*J]   += (B₁[i]*F₁₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₂₂*C₂₂₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₂ + B₂[i]*F₁₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₁ + B₁[i]*F₁₁*B₁[j]*F₂₂)*C₁₁₁₂ 
                               +  (B₁[i]*F₁₂*B₁[j]*F₂₁ + B₂[i]*F₁₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₁₂*B₂[j]*F₂₁ + B₂[i]*F₁₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂*B₂[j]*F₂₂ + B₂[i]*F₁₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂)*𝑤
                               
                k[2*I,2*J-1]   += (B₁[i]*F₂₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₂ + B₂[i]*F₂₂*B₁[j]*F₁₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₁ + B₁[i]*F₂₁*B₁[j]*F₁₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₁₁ + B₂[i]*F₂₁*B₁[j]*F₁₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₁₁ + B₂[i]*F₂₂*B₁[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₁₂ + B₂[i]*F₂₁*B₂[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂)*𝑤

                k[2*I,2*J]     += (B₁[i]*F₂₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₂₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₂ + B₂[i]*F₂₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₁ + B₁[i]*F₂₁*B₁[j]*F₂₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₂₁ + B₂[i]*F₂₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₂₁ + B₂[i]*F₂₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₂₂ + B₂[i]*F₂₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂
                               +   B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤

            end
        end
    end
end
"""
    ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean2(ap, k; E, ν)

Full part of the internal force of a 2-D Neo-Hookean material (the `_2` variant).

    f[I] += ∫ Bᵀ S w

Energy:  ψ = μ/2 (I₁ - 3) - μ ln J + λ/2 (ln J)²

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`, and the nodal displacements.
"""
function ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean2(ap::T, f::AbstractVector{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    K=E/(1-2*ν)/3
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)
        
        J⁻²=(1.0/J)^2
        J⁻²³=cbrt(J⁻²)
        S₁₁ = μ*J⁻²³*(1.0-1/3*I₁*C⁻¹₁₁) + 0.5*K*(J^2-1.0)*C⁻¹₁₁
        S₂₂ = μ*J⁻²³*(1.0-1/3*I₁*C⁻¹₂₂) + 0.5*K*(J^2-1.0)*C⁻¹₂₂
        S₁₂ = μ*J⁻²³*(-1/3*I₁*C⁻¹₁₂) + 0.5*K*(J^2-1.0)*C⁻¹₁₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += (B₁[i]*F₁₁*S₁₁+B₂[i]*F₁₂*S₂₂+(B₁[i]*F₁₂+B₂[i]*F₁₁)*S₁₂)*𝑤
            f[2*I]   += (B₁[i]*F₂₁*S₁₁+B₂[i]*F₂₂*S₂₂+(B₁[i]*F₂₂+B₂[i]*F₂₁)*S₁₂)*𝑤
        end
    end
end
"""
    ∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean2(ap, k; E, ν)

Volumetric part of the internal force of a 2-D Neo-Hookean material (the `_2` variant).

    f[I] += ∫ Bᵀ S w
Combined with the deviatoric companion it reproduces the full operator.

Energy:  ψ = μ/2 (I₁ - 3) - μ ln J + λ/2 (ln J)²

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`, and the nodal displacements.
"""
function ∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean2(ap::T, f::AbstractVector{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    K=E/(1-2*ν)/3
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        
        S₁₁ = 0.5*K*(J^2-1.0)*C⁻¹₁₁
        S₂₂ = 0.5*K*(J^2-1.0)*C⁻¹₂₂
        S₁₂ = 0.5*K*(J^2-1.0)*C⁻¹₁₂
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += (B₁[i]*F₁₁*S₁₁+B₂[i]*F₁₂*S₂₂+(B₁[i]*F₁₂+B₂[i]*F₁₁)*S₁₂)*𝑤
            f[2*I]   += (B₁[i]*F₂₁*S₁₁+B₂[i]*F₂₂*S₂₂+(B₁[i]*F₂₂+B₂[i]*F₂₁)*S₁₂)*𝑤
        end
    end
end
"""
    ∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean2(ap, k; E, ν)

Deviatoric part of the internal force of a 2-D Neo-Hookean material (the `_2` variant).

    f[I] += ∫ Bᵀ S w
Combined with the volumetric companion it reproduces the full operator.

Energy:  ψ = μ/2 (I₁ - 3) - μ ln J + λ/2 (ln J)²

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`, and the nodal displacements.
"""
function ∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean2(ap::T, f::AbstractVector{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    K=E/(1-2*ν)/3
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        J⁻²=(1.0/J)^2
        J⁻²³=cbrt(J⁻²)
        S₁₁ = μ*J⁻²³*(1.0-1/3*I₁*C⁻¹₁₁)
        S₂₂ = μ*J⁻²³*(1.0-1/3*I₁*C⁻¹₂₂)
        S₁₂ = μ*J⁻²³*(-1/3*I₁*C⁻¹₁₂)
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += (B₁[i]*F₁₁*S₁₁+B₂[i]*F₁₂*S₂₂+(B₁[i]*F₁₂+B₂[i]*F₁₁)*S₁₂)*𝑤
            f[2*I]   += (B₁[i]*F₂₁*S₁₁+B₂[i]*F₂₂*S₂₂+(B₁[i]*F₂₂+B₂[i]*F₂₁)*S₁₂)*𝑤
        end
    end
end

end