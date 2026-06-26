module Hyperelasticity

using ..ApproxOperator: AbstractElement
using LinearAlgebra,Printf
function ∫ESdx(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖 
    Eᵉ = op.E
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

function Δ∫∫δSᵢⱼXᵢⱼdxdy_HR_SaintVenantKirchhoff(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    for ξ in 𝓖
       E = ξ.E
       ν = ξ.ν
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

function Δ∫∫δSᵢⱼSᵢⱼdxdy_HR_SaintVenantKirchhoff(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    for ξ in 𝓖
       E = ξ.E
       ν = ξ.ν
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        C⁻¹ᵢᵢᵢᵢ = 1/E
        C⁻¹ᵢᵢⱼⱼ = -ν/E
        C⁻¹ᵢⱼᵢⱼ = 2*(1+ν)/E

        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼

                k[3*I-2,3*J-2] += N[i]*N[j]*𝑤
                k[3*I-2,3*J-1] += N[i]*N[j]*𝑤
                k[3*I-1,3*J-2] += N[i]*N[j]*𝑤
                k[3*I-1,3*J-1] += N[i]*N[j]*𝑤
                k[3*I,3*J]     += N[i]*N[j]*𝑤

            end
        end
    end
end

function Δ∫∫∫δSᵢⱼXᵢⱼdxdydz_HR_SaintVenantKirchhoff(ap::T, k::AbstractMatrix{Float64}) where {T<:AbstractElement}
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖

    for ξ in 𝓖
        E = ξ.E
        ν = ξ.ν
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤   

        C⁻¹ᵢᵢᵢᵢ = 1/E
        C⁻¹ᵢᵢⱼⱼ = -ν/E
        C⁻¹ᵢⱼᵢⱼ = 2*(1+ν)/E

        for (i, xᵃ) in enumerate(𝓒)
            I = xᵃ.𝐼
            for (j, xᵇ) in enumerate(𝓒)
                J = xᵇ.𝐼

                k[6I-5, 6J-5] += N[i]*N[j]*𝑤 * C⁻¹ᵢᵢᵢᵢ

                k[6I-5, 6J-4] += N[i]*N[j]*𝑤 * C⁻¹ᵢᵢⱼⱼ

                k[6I-5, 6J-3] += N[i]*N[j]*𝑤 * C⁻¹ᵢᵢⱼⱼ

                k[6I-4, 6J-5] += N[i]*N[j]*𝑤 * C⁻¹ᵢᵢⱼⱼ

                k[6I-4, 6J-4] += N[i]*N[j]*𝑤 * C⁻¹ᵢᵢᵢᵢ

                k[6I-4, 6J-3] += N[i]*N[j]*𝑤 * C⁻¹ᵢᵢⱼⱼ

                k[6I-3, 6J-5] += N[i]*N[j]*𝑤 * C⁻¹ᵢᵢⱼⱼ

                k[6I-3, 6J-4] += N[i]*N[j]*𝑤 * C⁻¹ᵢᵢⱼⱼ

                k[6I-3, 6J-3] += N[i]*N[j]*𝑤 * C⁻¹ᵢᵢᵢᵢ

                k[6I-2, 6J-2] += N[i]*N[j]*𝑤 * C⁻¹ᵢⱼᵢⱼ

                k[6I-1, 6J-1] += N[i]*N[j]*𝑤 * C⁻¹ᵢⱼᵢⱼ

                k[6I  , 6J  ] += N[i]*N[j]*𝑤 * C⁻¹ᵢⱼᵢⱼ

            end
        end
    end
end

function ∫∫δSᵢⱼXᵢⱼdxdy_HR_SaintVenantKirchhoff(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    for ξ in 𝓖
        E = ξ.E
        ν = ξ.ν
        𝑤 = ξ.𝑤
        C⁻¹ᵢᵢᵢᵢ = 1/E
        C⁻¹ᵢᵢⱼⱼ = -ν/E
        C⁻¹ᵢⱼᵢⱼ = 2*(1+ν)/E

        N = ξ[:𝝭]

        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒)
           S₁₁ += N[i]*xᵢ.dₛ₁₁
           S₂₂ += N[i]*xᵢ.dₛ₂₂
           S₁₂ += N[i]*xᵢ.dₛ₁₂
        end

        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[3*I-2] += (N[i]*C⁻¹ᵢᵢᵢᵢ*S₁₁+N[i]*C⁻¹ᵢᵢⱼⱼ*S₂₂)*𝑤
            f[3*I-1] += (N[i]*C⁻¹ᵢᵢⱼⱼ*S₁₁+N[i]*C⁻¹ᵢᵢᵢᵢ*S₂₂)*𝑤
            f[3*I]   += (N[i]*C⁻¹ᵢⱼᵢⱼ*S₁₂ )*𝑤
        end
    end
end

function ∫∫∫δSᵢⱼXᵢⱼdxdydz_HR_SaintVenantKirchhoff(ap::T, f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    for ξ in 𝓖
        E = ξ.E
        ν = ξ.ν
        𝑤 = ξ.𝑤

        C⁻¹ᵢᵢᵢᵢ = 1.0 / E                
        C⁻¹ᵢᵢⱼⱼ = -ν / E                 
        C⁻¹ᵢⱼᵢⱼ = 2.0 * (1.0 + ν) / E    

        N = ξ[:𝝭]

        S₁₁ = 0.0; S₂₂ = 0.0; S₃₃ = 0.0
        S₁₂ = 0.0; S₂₃ = 0.0; S₁₃ = 0.0

        for (i, xᵢ) in enumerate(𝓒)
           S₁₁ += N[i]*xᵢ.dₛ₁₁
           S₂₂ += N[i]*xᵢ.dₛ₂₂
           S₃₃ += N[i]*xᵢ.dₛ₃₃
           S₁₂ += N[i]*xᵢ.dₛ₁₂
           S₂₃ += N[i]*xᵢ.dₛ₂₃
           S₁₃ += N[i]*xᵢ.dₛ₁₃
        end

        X₁₁ = C⁻¹ᵢᵢᵢᵢ*S₁₁ + C⁻¹ᵢᵢⱼⱼ*S₂₂ + C⁻¹ᵢᵢⱼⱼ*S₃₃
        X₂₂ = C⁻¹ᵢᵢⱼⱼ*S₁₁ + C⁻¹ᵢᵢᵢᵢ*S₂₂ + C⁻¹ᵢᵢⱼⱼ*S₃₃
        X₃₃ = C⁻¹ᵢᵢⱼⱼ*S₁₁ + C⁻¹ᵢᵢⱼⱼ*S₂₂ + C⁻¹ᵢᵢᵢᵢ*S₃₃
        X₁₂ = C⁻¹ᵢⱼᵢⱼ*S₁₂
        X₂₃ = C⁻¹ᵢⱼᵢⱼ*S₂₃
        X₁₃ = C⁻¹ᵢⱼᵢⱼ*S₁₃

        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            N_w = N[i] * 𝑤  

            f[6*I-5] += N_w * X₁₁
            f[6*I-4] += N_w * X₂₂
            f[6*I-3] += N_w * X₃₃
            f[6*I-2] += N_w * X₁₂
            f[6*I-1] += N_w * X₂₃
            f[6*I]   += N_w * X₁₃
        end
    end
end

function Δ∫∫δSFnudxdy_HR(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖

    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

      for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼

                k[3*I-2,2*J-1] -= Nₛ[i]*n₁*u₁*B₁[j]*𝑤
                k[3*I-2,2*J]   -= Nₛ[i]*n₁*u₂*B₁[j]*𝑤
                k[3*I-1,2*J-1] -= Nₛ[i]*n₂*u₁*B₂[j]*𝑤
                k[3*I-1,2*J]   -= Nₛ[i]*n₂*u₂*B₂[j]*𝑤
                k[3*I,2*J-1]   -= (Nₛ[i]*n₂*u₁*B₁[j]+Nₛ[i]*n₁*u₁*B₂[j])*𝑤
                k[3*I,2*J]     -= (Nₛ[i]*n₁*u₂*B₂[j]+Nₛ[i]*n₂*u₂*B₁[j])*𝑤

                k[3*I-2,2*J-1] -= Nₛ[i]*n₁*F₁₁*N[j]*𝑤
                k[3*I-2,2*J]   -= Nₛ[i]*n₁*F₂₁*N[j]*𝑤
                k[3*I-1,2*J-1] -= Nₛ[i]*n₂*F₁₂*N[j]*𝑤
                k[3*I-1,2*J]   -= Nₛ[i]*n₂*F₂₂*N[j]*𝑤
                k[3*I,2*J-1]   -= (Nₛ[i]*n₂*F₁₁*N[j]+Nₛ[i]*n₁*F₁₂*N[j])*𝑤
                k[3*I,2*J]     -= (Nₛ[i]*n₁*F₂₂*N[j]+Nₛ[i]*n₂*F₂₁*N[j])*𝑤

            end
        end
    end
end

function Δ∫∫δSFnudxdy_HR_uS(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖

    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

      for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼

                k[3*I-2,2*J-1] -= Nₛ[i]*n₁*F₁₁*N[j]*𝑤
                k[3*I-2,2*J]   -= Nₛ[i]*n₁*F₂₁*N[j]*𝑤
                k[3*I-1,2*J-1] -= Nₛ[i]*n₂*F₁₂*N[j]*𝑤
                k[3*I-1,2*J]   -= Nₛ[i]*n₂*F₂₂*N[j]*𝑤
                k[3*I,2*J-1]   -= (Nₛ[i]*n₂*F₁₁*N[j]+Nₛ[i]*n₁*F₁₂*N[j])*𝑤
                k[3*I,2*J]     -= (Nₛ[i]*n₁*F₂₂*N[j]+Nₛ[i]*n₂*F₂₁*N[j])*𝑤

            end
        end
    end
end

function ∫∫∫ΔSFnδudxdydz_HR_uS(aₛ::T, aᵤ::S, k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒; 𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒; 𝓖ᵤ = aᵤ.𝓖

    for (ξₛ, ξᵤ) in zip(𝓖ₛ, 𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        B₃ = ξᵤ[:∂𝝭∂z]
        𝑤  = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N  = ξᵤ[:𝝭]

        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        n₃ = ξᵤ.n₃

        F₁₁ = 1.0; F₁₂ = 0.0; F₁₃ = 0.0
        F₂₁ = 0.0; F₂₂ = 1.0; F₂₃ = 0.0
        F₃₁ = 0.0; F₃₂ = 0.0; F₃₃ = 1.0
        u₁ = 0.0; u₂ = 0.0; u₃ = 0.0

        for (i, xᵢ) in enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁; F₁₂ += B₂[i]*xᵢ.d₁; F₁₃ += B₃[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂; F₂₂ += B₂[i]*xᵢ.d₂; F₂₃ += B₃[i]*xᵢ.d₂
            F₃₁ += B₁[i]*xᵢ.d₃; F₃₂ += B₂[i]*xᵢ.d₃; F₃₃ += B₃[i]*xᵢ.d₃
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
            u₃  += N[i]*xᵢ.d₃
        end

        for (i, xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j, xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼

                k[6*I-5, 3*J-2] -= Nₛ[i] * n₁ * F₁₁ * N[j] * 𝑤
                k[6*I-5, 3*J-1] -= Nₛ[i] * n₁ * F₂₁ * N[j] * 𝑤
                k[6*I-5, 3*J]   -= Nₛ[i] * n₁ * F₃₁ * N[j] * 𝑤

                k[6*I-4, 3*J-2] -= Nₛ[i] * n₂ * F₁₂ * N[j] * 𝑤
                k[6*I-4, 3*J-1] -= Nₛ[i] * n₂ * F₂₂ * N[j] * 𝑤
                k[6*I-4, 3*J]   -= Nₛ[i] * n₂ * F₃₂ * N[j] * 𝑤

                k[6*I-3, 3*J-2] -= Nₛ[i] * n₃ * F₁₃ * N[j] * 𝑤
                k[6*I-3, 3*J-1] -= Nₛ[i] * n₃ * F₂₃ * N[j] * 𝑤
                k[6*I-3, 3*J]   -= Nₛ[i] * n₃ * F₃₃ * N[j] * 𝑤

                k[6*I-2, 3*J-2] -= (Nₛ[i] * n₂ * F₁₁ * N[j] + Nₛ[i] * n₁ * F₁₂ * N[j]) * 𝑤
                k[6*I-2, 3*J-1] -= (Nₛ[i] * n₂ * F₂₁ * N[j] + Nₛ[i] * n₁ * F₂₂ * N[j]) * 𝑤
                k[6*I-2, 3*J]   -= (Nₛ[i] * n₂ * F₃₁ * N[j] + Nₛ[i] * n₁ * F₃₂ * N[j]) * 𝑤

                k[6*I-1, 3*J-2] -= (Nₛ[i] * n₃ * F₁₂ * N[j] + Nₛ[i] * n₂ * F₁₃ * N[j]) * 𝑤
                k[6*I-1, 3*J-1] -= (Nₛ[i] * n₃ * F₂₂ * N[j] + Nₛ[i] * n₂ * F₂₃ * N[j]) * 𝑤
                k[6*I-1, 3*J]   -= (Nₛ[i] * n₃ * F₃₂ * N[j] + Nₛ[i] * n₂ * F₃₃ * N[j]) * 𝑤

                k[6*I,   3*J-2] -= (Nₛ[i] * n₃ * F₁₁ * N[j] + Nₛ[i] * n₁ * F₁₃ * N[j]) * 𝑤
                k[6*I,   3*J-1] -= (Nₛ[i] * n₃ * F₂₁ * N[j] + Nₛ[i] * n₁ * F₂₃ * N[j]) * 𝑤
                k[6*I,   3*J]   -= (Nₛ[i] * n₃ * F₃₁ * N[j] + Nₛ[i] * n₁ * F₃₃ * N[j]) * 𝑤

            end
        end
    end
end

function ∫∫δSFnudxdy_HR(aₛ::T,aᵤ::S,f::AbstractVector{Float64})  where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖

    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)

        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        𝑤 = ξₛ.𝑤

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

        for (i,xᵢ) in enumerate(𝓒ₛ)
              I = xᵢ.𝐼

               f[3*I-2] -= (Nₛ[i]*n₁*F₁₁*u₁ + Nₛ[i]*n₁*F₂₁*u₂)*𝑤
               f[3*I-1] -= (Nₛ[i]*n₂*F₁₂*u₁ + Nₛ[i]*n₂*F₂₂*u₂)*𝑤
               f[3*I]   -= ((Nₛ[i]*n₂*F₁₁+Nₛ[i]*n₁*F₁₂)*u₁ + (Nₛ[i]*n₁*F₂₂+Nₛ[i]*n₂*F₂₁)*u₂)*𝑤

        end

    end
end

function Δ∫∫δ∇SFudxdy_HR(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖

    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        B₁₁ = ξᵤ[:∂²𝝭∂x²]
        B₁₂ = ξᵤ[:∂²𝝭∂x∂y]
        B₂₂ = ξᵤ[:∂²𝝭∂y²]

        Bₛ₁ = ξₛ[:∂𝝭∂x]
        Bₛ₂ = ξₛ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0

        F₁₁_₁ = 0.0
        F₁₁_₂ = 0.0
        F₁₂_₁ = 0.0
        F₁₂_₂ = 0.0
        F₂₁_₁ = 0.0
        F₂₁_₂ = 0.0
        F₂₂_₁ = 0.0
        F₂₂_₂ = 0.0

        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂

            F₁₁_₁  += B₁₁[i]*xᵢ.d₁
            F₁₁_₂  += B₁₂[i]*xᵢ.d₁
            F₁₂_₁  += B₁₂[i]*xᵢ.d₁
            F₁₂_₂  += B₂₂[i]*xᵢ.d₁
            F₂₁_₁  += B₁₁[i]*xᵢ.d₂
            F₂₁_₂  += B₁₂[i]*xᵢ.d₂
            F₂₂_₁  += B₁₂[i]*xᵢ.d₂
            F₂₂_₂  += B₂₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

      for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼

                k[3*I-2,2*J-1] += Bₛ₁[i]*u₁*B₁[j]*𝑤
                k[3*I-2,2*J]   += Bₛ₁[i]*u₂*B₁[j]*𝑤
                k[3*I-1,2*J-1] += Bₛ₂[i]*u₁*B₂[j]*𝑤
                k[3*I-1,2*J]   += Bₛ₂[i]*u₂*B₂[j]*𝑤
                k[3*I,2*J-1]   += (Bₛ₂[i]*u₁*B₁[j]+Bₛ₁[i]*u₁*B₂[j])*𝑤
                k[3*I,2*J]     += (Bₛ₁[i]*u₂*B₂[j]+Bₛ₂[i]*u₂*B₁[j])*𝑤

                k[3*I-2,2*J-1] += Bₛ₁[i]*F₁₁*N[j]*𝑤
                k[3*I-2,2*J]   += Bₛ₁[i]*F₂₁*N[j]*𝑤
                k[3*I-1,2*J-1] += Bₛ₂[i]*F₁₂*N[j]*𝑤
                k[3*I-1,2*J]   += Bₛ₂[i]*F₂₂*N[j]*𝑤
                k[3*I,2*J-1]   += (Bₛ₂[i]*F₁₁*N[j]+Bₛ₁[i]*F₁₂*N[j])*𝑤
                k[3*I,2*J]     += (Bₛ₁[i]*F₂₂*N[j]+Bₛ₂[i]*F₂₁*N[j])*𝑤

                k[3*I-2,2*J-1] += Nₛ[i]*u₁*B₁₁[j]*𝑤
                k[3*I-2,2*J]   += Nₛ[i]*u₂*B₁₁[j]*𝑤
                k[3*I-1,2*J-1] += Nₛ[i]*u₁*B₂₂[j]*𝑤
                k[3*I-1,2*J]   += Nₛ[i]*u₂*B₂₂[j]*𝑤
                k[3*I,2*J-1]   += (Nₛ[i]*u₁*B₁₂[j]+Nₛ[i]*u₁*B₁₂[j])*𝑤
                k[3*I,2*J]     += (Nₛ[i]*u₂*B₁₂[j]+Nₛ[i]*u₂*B₁₂[j])*𝑤

                k[3*I-2,2*J-1] += Nₛ[i]*F₁₁_₁*N[j]*𝑤
                k[3*I-2,2*J]   += Nₛ[i]*F₂₁_₁*N[j]*𝑤
                k[3*I-1,2*J-1] += Nₛ[i]*F₁₂_₂*N[j]*𝑤
                k[3*I-1,2*J]   += Nₛ[i]*F₂₂_₂*N[j]*𝑤
                k[3*I,2*J-1]   += (Nₛ[i]*F₁₁_₂*N[j]+Nₛ[i]*F₁₂_₁*N[j])*𝑤
                k[3*I,2*J]     += (Nₛ[i]*F₂₂_₁*N[j]+Nₛ[i]*F₂₁_₂*N[j])*𝑤

            end
        end
    end
end

function Δ∫∫δ∇SFudxdy_HR_uS(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖

    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        B₁₁ = ξᵤ[:∂²𝝭∂x²]
        B₁₂ = ξᵤ[:∂²𝝭∂x∂y]
        B₂₂ = ξᵤ[:∂²𝝭∂y²]

        Bₛ₁ = ξₛ[:∂𝝭∂x]
        Bₛ₂ = ξₛ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0

        F₁₁_₁ = 0.0
        F₁₁_₂ = 0.0
        F₁₂_₁ = 0.0
        F₁₂_₂ = 0.0
        F₂₁_₁ = 0.0
        F₂₁_₂ = 0.0
        F₂₂_₁ = 0.0
        F₂₂_₂ = 0.0

        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂

            F₁₁_₁  += B₁₁[i]*xᵢ.d₁
            F₁₁_₂  += B₁₂[i]*xᵢ.d₁
            F₁₂_₁  += B₁₂[i]*xᵢ.d₁
            F₁₂_₂  += B₂₂[i]*xᵢ.d₁
            F₂₁_₁  += B₁₁[i]*xᵢ.d₂
            F₂₁_₂  += B₁₂[i]*xᵢ.d₂
            F₂₂_₁  += B₁₂[i]*xᵢ.d₂
            F₂₂_₂  += B₂₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

      for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼

                k[3*I-2,2*J-1] += Bₛ₁[i]*F₁₁*N[j]*𝑤
                k[3*I-2,2*J]   += Bₛ₁[i]*F₂₁*N[j]*𝑤
                k[3*I-1,2*J-1] += Bₛ₂[i]*F₁₂*N[j]*𝑤
                k[3*I-1,2*J]   += Bₛ₂[i]*F₂₂*N[j]*𝑤
                k[3*I,2*J-1]   += (Bₛ₂[i]*F₁₁*N[j]+Bₛ₁[i]*F₁₂*N[j])*𝑤
                k[3*I,2*J]     += (Bₛ₁[i]*F₂₂*N[j]+Bₛ₂[i]*F₂₁*N[j])*𝑤

                k[3*I-2,2*J-1] += Nₛ[i]*F₁₁_₁*N[j]*𝑤
                k[3*I-2,2*J]   += Nₛ[i]*F₂₁_₁*N[j]*𝑤
                k[3*I-1,2*J-1] += Nₛ[i]*F₁₂_₂*N[j]*𝑤
                k[3*I-1,2*J]   += Nₛ[i]*F₂₂_₂*N[j]*𝑤
                k[3*I,2*J-1]   += (Nₛ[i]*F₁₁_₂*N[j]+Nₛ[i]*F₁₂_₁*N[j])*𝑤
                k[3*I,2*J]     += (Nₛ[i]*F₂₂_₁*N[j]+Nₛ[i]*F₂₁_₂*N[j])*𝑤

            end
        end
    end
end

 function Δ∫∫∫δ∇SFudxdydz_HR_uS(aₛ::T, aᵤ::S, k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒; 𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒; 𝓖ᵤ = aᵤ.𝓖

    for (ξₛ, ξᵤ) in zip(𝓖ₛ, 𝓖ᵤ)

        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        B₃ = ξᵤ[:∂𝝭∂z]

        B₁₁ = ξᵤ[:∂²𝝭∂x²]
        B₂₂ = ξᵤ[:∂²𝝭∂y²]
        B₃₃ = ξᵤ[:∂²𝝭∂z²]
        B₁₂ = ξᵤ[:∂²𝝭∂x∂y]
        B₂₃ = ξᵤ[:∂²𝝭∂y∂z]
        B₁₃ = ξᵤ[:∂²𝝭∂x∂z] 

        Bₛ₁ = ξₛ[:∂𝝭∂x]
        Bₛ₂ = ξₛ[:∂𝝭∂y]
        Bₛ₃ = ξₛ[:∂𝝭∂z]

        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]

        F₁₁ = 1.0; F₁₂ = 0.0; F₁₃ = 0.0
        F₂₁ = 0.0; F₂₂ = 1.0; F₂₃ = 0.0
        F₃₁ = 0.0; F₃₂ = 0.0; F₃₃ = 1.0

        F₁₁_₁ = 0.0; F₁₁_₂ = 0.0; F₁₁_₃ = 0.0
        F₁₂_₁ = 0.0; F₁₂_₂ = 0.0; F₁₂_₃ = 0.0
        F₁₃_₁ = 0.0; F₁₃_₂ = 0.0; F₁₃_₃ = 0.0

        F₂₁_₁ = 0.0; F₂₁_₂ = 0.0; F₂₁_₃ = 0.0
        F₂₂_₁ = 0.0; F₂₂_₂ = 0.0; F₂₂_₃ = 0.0
        F₂₃_₁ = 0.0; F₂₃_₂ = 0.0; F₂₃_₃ = 0.0

        F₃₁_₁ = 0.0; F₃₁_₂ = 0.0; F₃₁_₃ = 0.0
        F₃₂_₁ = 0.0; F₃₂_₂ = 0.0; F₃₂_₃ = 0.0
        F₃₃_₁ = 0.0; F₃₃_₂ = 0.0; F₃₃_₃ = 0.0

        u₁ = 0.0; u₂ = 0.0; u₃ = 0.0

        for (i, xᵢ) in enumerate(𝓒ᵤ)

            F₁₁ += B₁[i]*xᵢ.d₁; F₁₂ += B₂[i]*xᵢ.d₁; F₁₃ += B₃[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂; F₂₂ += B₂[i]*xᵢ.d₂; F₂₃ += B₃[i]*xᵢ.d₂
            F₃₁ += B₁[i]*xᵢ.d₃; F₃₂ += B₂[i]*xᵢ.d₃; F₃₃ += B₃[i]*xᵢ.d₃

            F₁₁_₁ += B₁₁[i]*xᵢ.d₁; F₁₁_₂ += B₁₂[i]*xᵢ.d₁; F₁₁_₃ += B₁₃[i]*xᵢ.d₁
            F₁₂_₁ += B₁₂[i]*xᵢ.d₁; F₁₂_₂ += B₂₂[i]*xᵢ.d₁; F₁₂_₃ += B₂₃[i]*xᵢ.d₁
            F₁₃_₁ += B₁₃[i]*xᵢ.d₁; F₁₃_₂ += B₂₃[i]*xᵢ.d₁; F₁₃_₃ += B₃₃[i]*xᵢ.d₁

            F₂₁_₁ += B₁₁[i]*xᵢ.d₂; F₂₁_₂ += B₁₂[i]*xᵢ.d₂; F₂₁_₃ += B₁₃[i]*xᵢ.d₂
            F₂₂_₁ += B₁₂[i]*xᵢ.d₂; F₂₂_₂ += B₂₂[i]*xᵢ.d₂; F₂₂_₃ += B₂₃[i]*xᵢ.d₂
            F₂₃_₁ += B₁₃[i]*xᵢ.d₂; F₂₃_₂ += B₂₃[i]*xᵢ.d₂; F₂₃_₃ += B₃₃[i]*xᵢ.d₂

            F₃₁_₁ += B₁₁[i]*xᵢ.d₃; F₃₁_₂ += B₁₂[i]*xᵢ.d₃; F₃₁_₃ += B₁₃[i]*xᵢ.d₃
            F₃₂_₁ += B₁₂[i]*xᵢ.d₃; F₃₂_₂ += B₂₂[i]*xᵢ.d₃; F₃₂_₃ += B₂₃[i]*xᵢ.d₃
            F₃₃_₁ += B₁₃[i]*xᵢ.d₃; F₃₃_₂ += B₂₃[i]*xᵢ.d₃; F₃₃_₃ += B₃₃[i]*xᵢ.d₃

            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            u₃ += N[i]*xᵢ.d₃
        end

        for (i, xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j, xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼

                k[6*I-5, 3*J-2] += Bₛ₁[i] * F₁₁ * N[j] * 𝑤
                k[6*I-4, 3*J-2] += Bₛ₂[i] * F₁₂ * N[j] * 𝑤
                k[6*I-3, 3*J-2] += Bₛ₃[i] * F₁₃ * N[j] * 𝑤
                k[6*I-2, 3*J-2] += (Bₛ₂[i] * F₁₁ + Bₛ₁[i] * F₁₂) * N[j] * 𝑤
                k[6*I-1, 3*J-2] += (Bₛ₃[i] * F₁₂ + Bₛ₂[i] * F₁₃) * N[j] * 𝑤
                k[6*I,   3*J-2] += (Bₛ₃[i] * F₁₁ + Bₛ₁[i] * F₁₃) * N[j] * 𝑤

                k[6*I-5, 3*J-1] += Bₛ₁[i] * F₂₁ * N[j] * 𝑤
                k[6*I-4, 3*J-1] += Bₛ₂[i] * F₂₂ * N[j] * 𝑤
                k[6*I-3, 3*J-1] += Bₛ₃[i] * F₂₃ * N[j] * 𝑤
                k[6*I-2, 3*J-1] += (Bₛ₂[i] * F₂₁ + Bₛ₁[i] * F₂₂) * N[j] * 𝑤
                k[6*I-1, 3*J-1] += (Bₛ₃[i] * F₂₂ + Bₛ₂[i] * F₂₃) * N[j] * 𝑤
                k[6*I,   3*J-1] += (Bₛ₃[i] * F₂₁ + Bₛ₁[i] * F₂₃) * N[j] * 𝑤

                k[6*I-5, 3*J]   += Bₛ₁[i] * F₃₁ * N[j] * 𝑤
                k[6*I-4, 3*J]   += Bₛ₂[i] * F₃₂ * N[j] * 𝑤
                k[6*I-3, 3*J]   += Bₛ₃[i] * F₃₃ * N[j] * 𝑤
                k[6*I-2, 3*J]   += (Bₛ₂[i] * F₃₁ + Bₛ₁[i] * F₃₂) * N[j] * 𝑤
                k[6*I-1, 3*J]   += (Bₛ₃[i] * F₃₂ + Bₛ₂[i] * F₃₃) * N[j] * 𝑤
                k[6*I,   3*J]   += (Bₛ₃[i] * F₃₁ + Bₛ₁[i] * F₃₃) * N[j] * 𝑤

                k[6*I-5, 3*J-2] += Nₛ[i] * F₁₁_₁ * N[j] * 𝑤
                k[6*I-4, 3*J-2] += Nₛ[i] * F₁₂_₂ * N[j] * 𝑤
                k[6*I-3, 3*J-2] += Nₛ[i] * F₁₃_₃ * N[j] * 𝑤
                k[6*I-2, 3*J-2] += Nₛ[i] * (F₁₁_₂ + F₁₂_₁) * N[j] * 𝑤
                k[6*I-1, 3*J-2] += Nₛ[i] * (F₁₂_₃ + F₁₃_₂) * N[j] * 𝑤
                k[6*I,   3*J-2] += Nₛ[i] * (F₁₁_₃ + F₁₃_₁) * N[j] * 𝑤

                k[6*I-5, 3*J-1] += Nₛ[i] * F₂₁_₁ * N[j] * 𝑤
                k[6*I-4, 3*J-1] += Nₛ[i] * F₂₂_₂ * N[j] * 𝑤
                k[6*I-3, 3*J-1] += Nₛ[i] * F₂₃_₃ * N[j] * 𝑤
                k[6*I-2, 3*J-1] += Nₛ[i] * (F₂₁_₂ + F₂₂_₁) * N[j] * 𝑤
                k[6*I-1, 3*J-1] += Nₛ[i] * (F₂₂_₃ + F₂₃_₂) * N[j] * 𝑤
                k[6*I,   3*J-1] += Nₛ[i] * (F₂₁_₃ + F₂₃_₁) * N[j] * 𝑤

                k[6*I-5, 3*J]   += Nₛ[i] * F₃₁_₁ * N[j] * 𝑤
                k[6*I-4, 3*J]   += Nₛ[i] * F₃₂_₂ * N[j] * 𝑤
                k[6*I-3, 3*J]   += Nₛ[i] * F₃₃_₃ * N[j] * 𝑤
                k[6*I-2, 3*J]   += Nₛ[i] * (F₃₁_₂ + F₃₂_₁) * N[j] * 𝑤
                k[6*I-1, 3*J]   += Nₛ[i] * (F₃₂_₃ + F₃₃_₂) * N[j] * 𝑤
                k[6*I,   3*J]   += Nₛ[i] * (F₃₁_₃ + F₃₃_₁) * N[j] * 𝑤

            end
        end
    end
end

function ∫∫δ∇SFudxdy_HR(aₛ::T,aᵤ::S,f::AbstractVector{Float64})  where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖

    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)

        Bₛ₁ = ξₛ[:∂𝝭∂x]
        Bₛ₂ = ξₛ[:∂𝝭∂y]
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        B₁₁ = ξᵤ[:∂²𝝭∂x²]
        B₁₂ = ξᵤ[:∂²𝝭∂x∂y]
        B₂₂ = ξᵤ[:∂²𝝭∂y²]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]

        𝑤 = ξₛ.𝑤

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0

        F₁₁_₁ = 0.0
        F₁₁_₂ = 0.0
        F₁₂_₁ = 0.0
        F₁₂_₂ = 0.0
        F₂₁_₁ = 0.0
        F₂₁_₂ = 0.0
        F₂₂_₁ = 0.0
        F₂₂_₂ = 0.0

        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂

            F₁₁_₁  += B₁₁[i]*xᵢ.d₁
            F₁₁_₂  += B₁₂[i]*xᵢ.d₁
            F₁₂_₁  += B₁₂[i]*xᵢ.d₁
            F₁₂_₂  += B₂₂[i]*xᵢ.d₁
            F₂₁_₁  += B₁₁[i]*xᵢ.d₂
            F₂₁_₂  += B₁₂[i]*xᵢ.d₂
            F₂₂_₁  += B₁₂[i]*xᵢ.d₂
            F₂₂_₂  += B₂₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

        for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
               f[3*I-2] += (Bₛ₁[i]*F₁₁*u₁ + Bₛ₁[i]*F₂₁*u₂)*𝑤
               f[3*I-1] += (Bₛ₂[i]*F₁₂*u₁ + Bₛ₂[i]*F₂₂*u₂)*𝑤
               f[3*I]   += ((Bₛ₂[i]*F₁₁+Bₛ₁[i]*F₁₂)*u₁ + (Bₛ₁[i]*F₂₂+Bₛ₂[i]*F₂₁)*u₂)*𝑤

               f[3*I-2] += (Nₛ[i]*F₁₁_₁*u₁ + Nₛ[i]*F₂₁_₁*u₂)*𝑤
               f[3*I-1] += (Nₛ[i]*F₁₂_₂*u₁ + Nₛ[i]*F₂₂_₂*u₂)*𝑤
               f[3*I]   += ((Nₛ[i]*F₁₁_₂+Nₛ[i]*F₁₂_₁)*u₁ + (Nₛ[i]*F₂₂_₁+Nₛ[i]*F₂₁_₂)*u₂)*𝑤

        end

    end
end

function ∫∫SδFnΔudxdy_HR(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖

    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end
        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
        end

      for (i,xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼

                k[2*I-1,2*J-1] -= (N[i]*B₁[j]*(S₁₁*n₁+S₁₂*n₂)+N[i]*B₂[j]*(S₁₂*n₁+S₂₂*n₂))*𝑤
                k[2*I,2*J]     -= (N[i]*B₁[j]*(S₁₁*n₁+S₁₂*n₂)+N[i]*B₂[j]*(S₁₂*n₁+S₂₂*n₂))*𝑤

            end
        end
    end
end
function ∫∫∫SδFnΔudxdydz_HR(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖

    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        B₃ = ξᵤ[:∂𝝭∂z]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        n₃ = ξᵤ.n₃

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₁₃ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        F₂₃ = 0.0
        F₃₁ = 0.0
        F₃₂ = 0.0
        F₃₃ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        u₃ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₁₃ += B₃[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            F₂₃ += B₃[i]*xᵢ.d₂
            F₃₁ += B₁[i]*xᵢ.d₃
            F₃₂ += B₂[i]*xᵢ.d₃
            F₃₃ += B₃[i]*xᵢ.d₃
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
            u₃  += N[i]*xᵢ.d₃
        end

        S₁₁ = 0.0
        S₂₂ = 0.0
        S₃₃ = 0.0
        S₁₂ = 0.0
        S₂₃ = 0.0
        S₁₃ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₃₃ += Nₛ[i]*xᵢ.dₛ₃₃
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
           S₂₃ += Nₛ[i]*xᵢ.dₛ₂₃
           S₁₃ += Nₛ[i]*xᵢ.dₛ₁₃
        end

      for (i,xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼

                k[3*I-2,3*J-2] -= (N[i]*B₁[j]*(S₁₁*n₁ + S₁₂*n₂ + S₁₃*n₃) + N[i]*B₂[j]*(S₁₂*n₁ + S₂₂*n₂ + S₂₃*n₃) + N[i]*B₃[j]*(S₁₃*n₁ + S₂₃*n₂ + S₃₃*n₃))*𝑤
                k[3*I-1,3*J-1] -= (N[i]*B₁[j]*(S₁₁*n₁ + S₁₂*n₂ + S₁₃*n₃) + N[i]*B₂[j]*(S₁₂*n₁ + S₂₂*n₂ + S₂₃*n₃) + N[i]*B₃[j]*(S₁₃*n₁ + S₂₃*n₂ + S₃₃*n₃))*𝑤
                k[3*I,3*J]     -= (N[i]*B₁[j]*(S₁₁*n₁ + S₁₂*n₂ + S₁₃*n₃) + N[i]*B₂[j]*(S₁₂*n₁ + S₂₂*n₂ + S₂₃*n₃) + N[i]*B₃[j]*(S₁₃*n₁ + S₂₃*n₂ + S₃₃*n₃))*𝑤

            end
        end
    end
end

function ∫∫SδFnδudxdy_HR(aₛ::T,aᵤ::S,f::AbstractVector{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖

    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end
        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
        end

        for (i,xᵢ) in enumerate(𝓒ᵤ)
                I = xᵢ.𝐼

               f[2*I-1] -= (N[i]*F₁₁*S₁₁*n₁+N[i]*F₁₁*S₁₂*n₂+N[i]*F₁₂*S₁₂*n₁+N[i]*F₁₂*S₂₂*n₂ )*𝑤
               f[2*I]   -= (N[i]*F₂₁*S₁₁*n₁+N[i]*F₂₁*S₁₂*n₂+N[i]*F₂₂*S₁₂*n₁+N[i]*F₂₂*S₂₂*n₂ )*𝑤    

        end
    end
end

function ∫∫∫SδFnδudxdydz_HR(aₛ::T,aᵤ::S,f::AbstractVector{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖

    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        B₃ = ξᵤ[:∂𝝭∂z]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        n₃ = ξᵤ.n₃

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₁₃ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        F₂₃ = 0.0
        F₃₁ = 0.0
        F₃₂ = 0.0
        F₃₃ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        u₃ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₁₃ += B₃[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            F₂₃ += B₃[i]*xᵢ.d₂
            F₃₁ += B₁[i]*xᵢ.d₃
            F₃₂ += B₂[i]*xᵢ.d₃
            F₃₃ += B₃[i]*xᵢ.d₃
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
            u₃  += N[i]*xᵢ.d₃
        end

        S₁₁ = 0.0
        S₂₂ = 0.0
        S₃₃ = 0.0
        S₁₂ = 0.0
        S₂₃ = 0.0
        S₁₃ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₃₃ += Nₛ[i]*xᵢ.dₛ₃₃
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
           S₂₃ += Nₛ[i]*xᵢ.dₛ₂₃
           S₁₃ += Nₛ[i]*xᵢ.dₛ₁₃
        end

        for (i,xᵢ) in enumerate(𝓒ᵤ)
                I = xᵢ.𝐼

               f[3*I-2] -= (N[i]*F₁₁*S₁₁*n₁ + N[i]*F₁₁*S₁₂*n₂ + N[i]*F₁₁*S₁₃*n₃ + N[i]*F₁₂*S₁₂*n₁ + N[i]*F₁₂*S₂₂*n₂ + N[i]*F₁₂*S₂₃*n₃ + N[i]*F₁₃*S₁₃*n₁ + N[i]*F₁₃*S₂₃*n₂ + N[i]*F₁₃*S₃₃*n₃)*𝑤
               f[3*I-1] -= (N[i]*F₂₁*S₁₁*n₁ + N[i]*F₂₁*S₁₂*n₂ + N[i]*F₂₁*S₁₃*n₃ + N[i]*F₂₂*S₁₂*n₁ + N[i]*F₂₂*S₂₂*n₂ + N[i]*F₂₂*S₂₃*n₃ + N[i]*F₂₃*S₁₃*n₁ + N[i]*F₂₃*S₂₃*n₂ + N[i]*F₂₃*S₃₃*n₃)*𝑤    
               f[3*I]   -= (N[i]*F₃₁*S₁₁*n₁ + N[i]*F₃₁*S₁₂*n₂ + N[i]*F₃₁*S₁₃*n₃ + N[i]*F₃₂*S₁₂*n₁ + N[i]*F₃₂*S₂₂*n₂ + N[i]*F₃₂*S₂₃*n₃ + N[i]*F₃₃*S₁₃*n₁ + N[i]*F₃₃*S₂₃*n₂ + N[i]*F₃₃*S₃₃*n₃)*𝑤    

        end
    end
end

function ∫∫∇SδFΔudxdy_HR(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖

    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        B₁₁ = ξᵤ[:∂²𝝭∂x²]
        B₁₂ = ξᵤ[:∂²𝝭∂x∂y]
        B₂₂ = ξᵤ[:∂²𝝭∂y²]
        Bₛ₁ = ξₛ[:∂𝝭∂x]
        Bₛ₂ = ξₛ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]

        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        S₁₁_₁ = 0.0
        S₂₂_₂ = 0.0
        S₁₂_₁ = 0.0
        S₁₂_₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂

           S₁₁_₁ += Bₛ₁[i]*xᵢ.dₛ₁₁
           S₂₂_₂ += Bₛ₂[i]*xᵢ.dₛ₂₂
           S₁₂_₁ += Bₛ₁[i]*xᵢ.dₛ₁₂
           S₁₂_₂ += Bₛ₂[i]*xᵢ.dₛ₁₂

        end

      for (i,xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼

                k[2*I-1,2*J-1] += (N[i]*B₁[j]*(S₁₁_₁+S₁₂_₂)+N[i]*B₂[j]*(S₁₂_₁+S₂₂_₂))*𝑤
                k[2*I,2*J]     += (N[i]*B₁[j]*(S₁₁_₁+S₁₂_₂)+N[i]*B₂[j]*(S₁₂_₁+S₂₂_₂))*𝑤

                k[2*I-1,2*J-1] += (N[i]*(B₁₁[j]*S₁₁+B₁₂[j]*S₁₂)+N[i]*(B₁₂[j]*S₁₂+B₂₂[j]*S₂₂))*𝑤
                k[2*I,2*J]     += (N[i]*(B₁₁[j]*S₁₁+B₁₂[j]*S₁₂)+N[i]*(B₁₂[j]*S₁₂+B₂₂[j]*S₂₂))*𝑤

            end
        end
    end
end

function ∫∫∫∇SδFΔudxdydz_HR(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖

    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        B₃ = ξᵤ[:∂𝝭∂z]
        B₁₁ = ξᵤ[:∂²𝝭∂x²]
        B₂₂ = ξᵤ[:∂²𝝭∂y²]
        B₃₃ = ξᵤ[:∂²𝝭∂z²]
        B₁₂ = ξᵤ[:∂²𝝭∂x∂y]
        B₂₃ = ξᵤ[:∂²𝝭∂y∂z]
        B₁₃ = ξᵤ[:∂²𝝭∂x∂z]
        Bₛ₁ = ξₛ[:∂𝝭∂x]
        Bₛ₂ = ξₛ[:∂𝝭∂y]
        Bₛ₃ = ξₛ[:∂𝝭∂z]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]

        S₁₁ = 0.0
        S₂₂ = 0.0
        S₃₃ = 0.0
        S₁₂ = 0.0
        S₂₃ = 0.0
        S₁₃ = 0.0

        S₁₁_₁ = 0.0
        S₁₂_₂ = 0.0
        S₁₃_₃ = 0.0
        S₁₂_₁ = 0.0
        S₂₂_₂ = 0.0
        S₂₃_₃ = 0.0
        S₁₃_₁ = 0.0
        S₂₃_₂ = 0.0
        S₃₃_₃ = 0.0

        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₃₃ += Nₛ[i]*xᵢ.dₛ₃₃
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
           S₂₃ += Nₛ[i]*xᵢ.dₛ₂₃
           S₁₃ += Nₛ[i]*xᵢ.dₛ₁₃

           S₁₁_₁ += Bₛ₁[i]*xᵢ.dₛ₁₁
           S₁₂_₂ += Bₛ₂[i]*xᵢ.dₛ₁₂
           S₁₃_₃ += Bₛ₃[i]*xᵢ.dₛ₁₃

           S₁₂_₁ += Bₛ₁[i]*xᵢ.dₛ₁₂
           S₂₂_₂ += Bₛ₂[i]*xᵢ.dₛ₂₂
           S₂₃_₃ += Bₛ₃[i]*xᵢ.dₛ₂₃

           S₁₃_₁ += Bₛ₁[i]*xᵢ.dₛ₁₃
           S₂₃_₂ += Bₛ₂[i]*xᵢ.dₛ₂₃
           S₃₃_₃ += Bₛ₃[i]*xᵢ.dₛ₃₃
        end

      for (i,xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼

                k[3*I-2,3*J-2] += (N[i]*B₁[j]*(S₁₁_₁ + S₁₂_₂ + S₁₃_₃) + N[i]*B₂[j]*(S₁₂_₁ + S₂₂_₂ + S₂₃_₃) + N[i]*B₃[j]*(S₁₃_₁ + S₂₃_₂ + S₃₃_₃))*𝑤
                k[3*I-1,3*J-1] += (N[i]*B₁[j]*(S₁₁_₁ + S₁₂_₂ + S₁₃_₃) + N[i]*B₂[j]*(S₁₂_₁ + S₂₂_₂ + S₂₃_₃) + N[i]*B₃[j]*(S₁₃_₁ + S₂₃_₂ + S₃₃_₃))*𝑤
                k[3*I,3*J]     += (N[i]*B₁[j]*(S₁₁_₁ + S₁₂_₂ + S₁₃_₃) + N[i]*B₂[j]*(S₁₂_₁ + S₂₂_₂ + S₂₃_₃) + N[i]*B₃[j]*(S₁₃_₁ + S₂₃_₂ + S₃₃_₃))*𝑤

                k[3*I-2,3*J-2] += (N[i]*(B₁₁[j]*S₁₁ + B₁₂[j]*S₁₂ + B₁₃[j]*S₁₃) + N[i]*(B₁₂[j]*S₁₂ + B₂₂[j]*S₂₂ + B₂₃[j]*S₂₃) + N[i]*(B₁₃[j]*S₁₃ + B₂₃[j]*S₂₃ + B₃₃[j]*S₃₃))*𝑤
                k[3*I-1,3*J-1] += (N[i]*(B₁₁[j]*S₁₁ + B₁₂[j]*S₁₂ + B₁₃[j]*S₁₃) + N[i]*(B₁₂[j]*S₁₂ + B₂₂[j]*S₂₂ + B₂₃[j]*S₂₃) + N[i]*(B₁₃[j]*S₁₃ + B₂₃[j]*S₂₃ + B₃₃[j]*S₃₃))*𝑤
                k[3*I,3*J]     += (N[i]*(B₁₁[j]*S₁₁ + B₁₂[j]*S₁₂ + B₁₃[j]*S₁₃) + N[i]*(B₁₂[j]*S₁₂ + B₂₂[j]*S₂₂ + B₂₃[j]*S₂₃) + N[i]*(B₁₃[j]*S₁₃ + B₂₃[j]*S₂₃ + B₃₃[j]*S₃₃))*𝑤

            end
        end
    end
end

function ∫∫∇SδFδudxdy_HR(aₛ::T,aᵤ::S,f::AbstractVector{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖

    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
         B₁₁ = ξᵤ[:∂²𝝭∂x²]
        B₁₂ = ξᵤ[:∂²𝝭∂x∂y]
        B₂₂ = ξᵤ[:∂²𝝭∂y²]
        Bₛ₁ = ξₛ[:∂𝝭∂x]
        Bₛ₂ = ξₛ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0

        F₁₁_₁ = 0.0
        F₁₁_₂ = 0.0
        F₁₂_₁ = 0.0
        F₁₂_₂ = 0.0
        F₂₁_₁ = 0.0
        F₂₁_₂ = 0.0
        F₂₂_₁ = 0.0
        F₂₂_₂ = 0.0

        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂

            F₁₁_₁  += B₁₁[i]*xᵢ.d₁
            F₁₁_₂  += B₁₂[i]*xᵢ.d₁
            F₁₂_₁  += B₁₂[i]*xᵢ.d₁
            F₁₂_₂  += B₂₂[i]*xᵢ.d₁
            F₂₁_₁  += B₁₁[i]*xᵢ.d₂
            F₂₁_₂  += B₁₂[i]*xᵢ.d₂
            F₂₂_₁  += B₁₂[i]*xᵢ.d₂
            F₂₂_₂  += B₂₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        S₁₁_₁ = 0.0
        S₂₂_₂ = 0.0
        S₁₂_₁ = 0.0
        S₁₂_₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂

           S₁₁_₁ += Bₛ₁[i]*xᵢ.dₛ₁₁
           S₂₂_₂ += Bₛ₂[i]*xᵢ.dₛ₂₂
           S₁₂_₁ += Bₛ₁[i]*xᵢ.dₛ₁₂
           S₁₂_₂ += Bₛ₂[i]*xᵢ.dₛ₁₂

        end

      for (i,xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼

               f[2*I-1] += (N[i]*F₁₁*S₁₁_₁+N[i]*F₁₁*S₁₂_₂+N[i]*F₁₂*S₁₂_₁+N[i]*F₁₂*S₂₂_₂ )*𝑤
               f[2*I]   += (N[i]*F₂₁*S₁₁_₁+N[i]*F₂₁*S₁₂_₂+N[i]*F₂₂*S₁₂_₁+N[i]*F₂₂*S₂₂_₂ )*𝑤

               f[2*I-1] += (N[i]*F₁₁_₁*S₁₁+N[i]*F₁₁_₂*S₁₂+N[i]*F₁₂_₁*S₁₂+N[i]*F₁₂_₂*S₂₂ )*𝑤
               f[2*I]   += (N[i]*F₂₁_₁*S₁₁+N[i]*F₂₁_₂*S₁₂+N[i]*F₂₂_₁*S₁₂+N[i]*F₂₂_₂*S₂₂ )*𝑤

        end
    end
end

function ∫∫∫∇SδFδudxdydz_HR(aₛ::T,aᵤ::S,f::AbstractVector{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖

    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        B₃ = ξᵤ[:∂𝝭∂z]
        B₁₁ = ξᵤ[:∂²𝝭∂x²]
        B₂₂ = ξᵤ[:∂²𝝭∂y²]
        B₃₃ = ξᵤ[:∂²𝝭∂z²]
        B₁₂ = ξᵤ[:∂²𝝭∂x∂y]
        B₂₃ = ξᵤ[:∂²𝝭∂y∂z]
        B₁₃ = ξᵤ[:∂²𝝭∂x∂z]
        Bₛ₁ = ξₛ[:∂𝝭∂x]
        Bₛ₂ = ξₛ[:∂𝝭∂y]
        Bₛ₃ = ξₛ[:∂𝝭∂z]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]

        F₁₁ = 1.0; F₁₂ = 0.0; F₁₃ = 0.0
        F₂₁ = 0.0; F₂₂ = 1.0; F₂₃ = 0.0
        F₃₁ = 0.0; F₃₂ = 0.0; F₃₃ = 1.0

        F₁₁_₁ = 0.0; F₁₁_₂ = 0.0; F₁₁_₃ = 0.0
        F₁₂_₁ = 0.0; F₁₂_₂ = 0.0; F₁₂_₃ = 0.0
        F₁₃_₁ = 0.0; F₁₃_₂ = 0.0; F₁₃_₃ = 0.0

        F₂₁_₁ = 0.0; F₂₁_₂ = 0.0; F₂₁_₃ = 0.0
        F₂₂_₁ = 0.0; F₂₂_₂ = 0.0; F₂₂_₃ = 0.0
        F₂₃_₁ = 0.0; F₂₃_₂ = 0.0; F₂₃_₃ = 0.0

        F₃₁_₁ = 0.0; F₃₁_₂ = 0.0; F₃₁_₃ = 0.0
        F₃₂_₁ = 0.0; F₃₂_₂ = 0.0; F₃₂_₃ = 0.0
        F₃₃_₁ = 0.0; F₃₃_₂ = 0.0; F₃₃_₃ = 0.0

        u₁ = 0.0; u₂ = 0.0; u₃ = 0.0

        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₁₃ += B₃[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            F₂₃ += B₃[i]*xᵢ.d₂
            F₃₁ += B₁[i]*xᵢ.d₃
            F₃₂ += B₂[i]*xᵢ.d₃
            F₃₃ += B₃[i]*xᵢ.d₃

            F₁₁_₁  += B₁₁[i]*xᵢ.d₁
            F₁₁_₂  += B₁₂[i]*xᵢ.d₁
            F₁₁_₃  += B₁₃[i]*xᵢ.d₁
            F₁₂_₁  += B₁₂[i]*xᵢ.d₁
            F₁₂_₂  += B₂₂[i]*xᵢ.d₁
            F₁₂_₃  += B₂₃[i]*xᵢ.d₁
            F₁₃_₁  += B₁₃[i]*xᵢ.d₁
            F₁₃_₂  += B₂₃[i]*xᵢ.d₁
            F₁₃_₃  += B₃₃[i]*xᵢ.d₁

            F₂₁_₁  += B₁₁[i]*xᵢ.d₂
            F₂₁_₂  += B₁₂[i]*xᵢ.d₂
            F₂₁_₃  += B₁₃[i]*xᵢ.d₂
            F₂₂_₁  += B₁₂[i]*xᵢ.d₂
            F₂₂_₂  += B₂₂[i]*xᵢ.d₂
            F₂₂_₃  += B₂₃[i]*xᵢ.d₂
            F₂₃_₁  += B₁₃[i]*xᵢ.d₂
            F₂₃_₂  += B₂₃[i]*xᵢ.d₂
            F₂₃_₃  += B₃₃[i]*xᵢ.d₂

            F₃₁_₁  += B₁₁[i]*xᵢ.d₃
            F₃₁_₂  += B₁₂[i]*xᵢ.d₃
            F₃₁_₃  += B₁₃[i]*xᵢ.d₃
            F₃₂_₁  += B₁₂[i]*xᵢ.d₃
            F₃₂_₂  += B₂₂[i]*xᵢ.d₃
            F₃₂_₃  += B₂₃[i]*xᵢ.d₃
            F₃₃_₁  += B₁₃[i]*xᵢ.d₃
            F₃₃_₂  += B₂₃[i]*xᵢ.d₃
            F₃₃_₃  += B₃₃[i]*xᵢ.d₃

            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
            u₃  += N[i]*xᵢ.d₃
        end

        S₁₁ = 0.0; S₂₂ = 0.0; S₃₃ = 0.0
        S₁₂ = 0.0; S₂₃ = 0.0; S₁₃ = 0.0

        S₁₁_₁ = 0.0; S₁₂_₂ = 0.0; S₁₃_₃ = 0.0
        S₁₂_₁ = 0.0; S₂₂_₂ = 0.0; S₂₃_₃ = 0.0
        S₁₃_₁ = 0.0; S₂₃_₂ = 0.0; S₃₃_₃ = 0.0

        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₃₃ += Nₛ[i]*xᵢ.dₛ₃₃
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
           S₂₃ += Nₛ[i]*xᵢ.dₛ₂₃
           S₁₃ += Nₛ[i]*xᵢ.dₛ₁₃

           S₁₁_₁ += Bₛ₁[i]*xᵢ.dₛ₁₁
           S₁₂_₂ += Bₛ₂[i]*xᵢ.dₛ₁₂
           S₁₃_₃ += Bₛ₃[i]*xᵢ.dₛ₁₃

           S₁₂_₁ += Bₛ₁[i]*xᵢ.dₛ₁₂
           S₂₂_₂ += Bₛ₂[i]*xᵢ.dₛ₂₂
           S₂₃_₃ += Bₛ₃[i]*xᵢ.dₛ₂₃

           S₁₃_₁ += Bₛ₁[i]*xᵢ.dₛ₁₃
           S₂₃_₂ += Bₛ₂[i]*xᵢ.dₛ₂₃
           S₃₃_₃ += Bₛ₃[i]*xᵢ.dₛ₃₃
        end

        for (i,xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼

               f[3*I-2] += (N[i]*F₁₁*S₁₁_₁ + N[i]*F₁₁*S₁₂_₂ + N[i]*F₁₁*S₁₃_₃ + N[i]*F₁₂*S₁₂_₁ + N[i]*F₁₂*S₂₂_₂ + N[i]*F₁₂*S₂₃_₃ + N[i]*F₁₃*S₁₃_₁ + N[i]*F₁₃*S₂₃_₂ + N[i]*F₁₃*S₃₃_₃)*𝑤
               f[3*I-1] += (N[i]*F₂₁*S₁₁_₁ + N[i]*F₂₁*S₁₂_₂ + N[i]*F₂₁*S₁₃_₃ + N[i]*F₂₂*S₁₂_₁ + N[i]*F₂₂*S₂₂_₂ + N[i]*F₂₂*S₂₃_₃ + N[i]*F₂₃*S₁₃_₁ + N[i]*F₂₃*S₂₃_₂ + N[i]*F₂₃*S₃₃_₃)*𝑤
               f[3*I]   += (N[i]*F₃₁*S₁₁_₁ + N[i]*F₃₁*S₁₂_₂ + N[i]*F₃₁*S₁₃_₃ + N[i]*F₃₂*S₁₂_₁ + N[i]*F₃₂*S₂₂_₂ + N[i]*F₃₂*S₂₃_₃ + N[i]*F₃₃*S₁₃_₁ + N[i]*F₃₃*S₂₃_₂ + N[i]*F₃₃*S₃₃_₃)*𝑤

               f[3*I-2] += (N[i]*F₁₁_₁*S₁₁ + N[i]*F₁₁_₂*S₁₂ + N[i]*F₁₁_₃*S₁₃ + N[i]*F₁₂_₁*S₁₂ + N[i]*F₁₂_₂*S₂₂ + N[i]*F₁₂_₃*S₂₃ + N[i]*F₁₃_₁*S₁₃ + N[i]*F₁₃_₂*S₂₃ + N[i]*F₁₃_₃*S₃₃)*𝑤
               f[3*I-1] += (N[i]*F₂₁_₁*S₁₁ + N[i]*F₂₁_₂*S₁₂ + N[i]*F₂₁_₃*S₁₃ + N[i]*F₂₂_₁*S₁₂ + N[i]*F₂₂_₂*S₂₂ + N[i]*F₂₂_₃*S₂₃ + N[i]*F₂₃_₁*S₁₃ + N[i]*F₂₃_₂*S₂₃ + N[i]*F₂₃_₃*S₃₃)*𝑤
               f[3*I]   += (N[i]*F₃₁_₁*S₁₁ + N[i]*F₃₁_₂*S₁₂ + N[i]*F₃₁_₃*S₁₃ + N[i]*F₃₂_₁*S₁₂ + N[i]*F₃₂_₂*S₂₂ + N[i]*F₃₂_₃*S₂₃ + N[i]*F₃₃_₁*S₁₃ + N[i]*F₃₃_₂*S₂₃ + N[i]*F₃₃_₃*S₃₃)*𝑤

        end
    end
end

function ∫∫δSΔFngdxdy_HR(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖

    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        n₁₁ = ξᵤ.n₁₁
        n₁₂ = ξᵤ.n₁₂
        n₂₂ = ξᵤ.n₂₂
        g₁ = ξᵤ.g₁
        g₂ = ξᵤ.g₂

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
        end

        Δu₁ = g₁-u₁
        Δu₂ = g₂-u₂

      for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼

                k[3*I-2,2*J-1] += Nₛ[i]*(n₁*n₁₁*u₁+n₁*n₁₂*u₂)*B₁[j]*𝑤
                k[3*I-2,2*J]   += Nₛ[i]*(n₁*n₁₂*u₁+n₁*n₂₂*u₂)*B₁[j]*𝑤
                k[3*I-1,2*J-1] += Nₛ[i]*(n₂*n₁₁*u₁+n₂*n₁₂*u₂)*B₂[j]*𝑤
                k[3*I-1,2*J]   += Nₛ[i]*(n₂*n₁₂*u₁+n₂*n₂₂*u₂)*B₂[j]*𝑤
                k[3*I,2*J-1]   += (Nₛ[i]*(n₂*n₁₁*u₁+n₂*n₁₂*u₂)*B₁[j]+Nₛ[i]*(n₁*n₁₁*u₁+n₁*n₁₂*u₂)*B₂[j])*𝑤
                k[3*I,2*J]     += (Nₛ[i]*(n₂*n₁₂*u₁+n₂*n₂₂*u₂)*B₁[j]+Nₛ[i]*(n₁*n₁₂*u₁+n₁*n₂₂*u₂)*B₂[j])*𝑤

                k[3*I-2,2*J-1] += Nₛ[i]*(n₁*n₁₁*F₁₁+n₁*n₁₂*F₂₁)*N[j]*𝑤
                k[3*I-2,2*J]   += Nₛ[i]*(n₁*n₁₂*F₁₁+n₁*n₂₂*F₂₁)*N[j]*𝑤
                k[3*I-1,2*J-1] += Nₛ[i]*(n₂*n₁₁*F₁₂+n₂*n₁₂*F₂₂)*N[j]*𝑤
                k[3*I-1,2*J]   += Nₛ[i]*(n₂*n₁₂*F₁₂+n₂*n₂₂*F₂₂)*N[j]*𝑤
                k[3*I,2*J-1]   += (Nₛ[i]*(n₂*n₁₁*F₁₁+n₂*n₁₂*F₂₁)*N[j]+Nₛ[i]*(n₁*n₁₁*F₁₂+n₁*n₁₂*F₂₂)*N[j])*𝑤
                k[3*I,2*J]     += (Nₛ[i]*(n₂*n₁₂*F₁₁+n₂*n₂₂*F₂₁)*N[j]+Nₛ[i]*(n₁*n₁₂*F₁₂+n₁*n₂₂*F₂₂)*N[j])*𝑤

                k[3*I-2,2*J-1] -= Nₛ[i]*(n₁*n₁₁*g₁+n₁*n₁₂*g₂)*B₁[j]*𝑤
                k[3*I-2,2*J]   -= Nₛ[i]*(n₁*n₁₂*g₁+n₁*n₂₂*g₂)*B₁[j]*𝑤
                k[3*I-1,2*J-1] -= Nₛ[i]*(n₂*n₁₁*g₁+n₂*n₁₂*g₂)*B₂[j]*𝑤
                k[3*I-1,2*J]   -= Nₛ[i]*(n₂*n₁₂*g₁+n₂*n₂₂*g₂)*B₂[j]*𝑤
                k[3*I,2*J-1]   -= (Nₛ[i]*(n₂*n₁₁*g₁+n₂*n₁₂*g₂)*B₁[j]+Nₛ[i]*(n₁*n₁₁*g₁+n₁*n₁₂*g₂)*B₂[j])*𝑤
                k[3*I,2*J]     -= (Nₛ[i]*(n₂*n₁₂*g₁+n₂*n₂₂*g₂)*B₁[j]+Nₛ[i]*(n₁*n₁₂*g₁+n₁*n₂₂*g₂)*B₂[j])*𝑤

            end

            f[3*I-2] += Nₛ[i]*((n₁*n₁₁*F₁₁+n₁*n₁₂*F₂₁)*Δu₁ + (n₁*n₁₂*F₁₁+n₁*n₂₂*F₂₁)*Δu₂)*𝑤
            f[3*I-1] += Nₛ[i]*((n₂*n₁₁*F₁₂+n₂*n₁₂*F₂₂)*Δu₁ + (n₂*n₁₂*F₁₂+n₂*n₂₂*F₂₂)*Δu₂)*𝑤
            f[3*I]   += Nₛ[i]*(((n₂*n₁₁*F₁₁+n₂*n₁₂*F₂₁)+(n₁*n₁₁*F₁₂+n₁*n₁₂*F₂₂))*Δu₁ + ((n₂*n₁₂*F₁₁+n₂*n₂₂*F₂₁)+(n₁*n₁₂*F₁₂+n₁*n₂₂*F₂₂))*Δu₂)*𝑤

        end

    end
end

function ∫∫∫δSΔFngdxdydz_HR(aₛ::T, aᵤ::S, k::AbstractMatrix{Float64}, f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒; 𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒; 𝓖ᵤ = aᵤ.𝓖
    Nₛ = zeros(length(𝓒ₛ))
    for (ξₛ, ξᵤ) in zip(𝓖ₛ, 𝓖ᵤ)

        B₁ = ξᵤ[:∂𝝭∂x]; B₂ = ξᵤ[:∂𝝭∂y]; B₃ = ξᵤ[:∂𝝭∂z]
        x = ξᵤ.x
        y = ξᵤ.y
        z = ξᵤ.z
        eval_piecewise_3d!(Nₛ, x, y, z)

        N  = ξᵤ[:𝝭]
        𝑤  = ξₛ.𝑤

        n₁ = ξᵤ.n₁; n₂ = ξᵤ.n₂; n₃ = ξᵤ.n₃

        n₁₁ = ξᵤ.n₁₁; n₂₂ = ξᵤ.n₂₂; n₃₃ = ξᵤ.n₃₃
        n₁₂ = ξᵤ.n₁₂; n₂₃ = ξᵤ.n₂₃; n₁₃ = ξᵤ.n₁₃

        g₁ = ξᵤ.g₁; g₂ = ξᵤ.g₂; g₃ = ξᵤ.g₃

        F₁₁ = 1.0; F₁₂ = 0.0; F₁₃ = 0.0
        F₂₁ = 0.0; F₂₂ = 1.0; F₂₃ = 0.0
        F₃₁ = 0.0; F₃₂ = 0.0; F₃₃ = 1.0
        u₁ = 0.0; u₂ = 0.0; u₃ = 0.0

        S₁₁ = 0.0; S₂₂ = 0.0; S₃₃ = 0.0
        S₁₂ = 0.0; S₂₃ = 0.0; S₁₃ = 0.0

        for (i, xᵢ) in enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁; F₁₂ += B₂[i]*xᵢ.d₁; F₁₃ += B₃[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂; F₂₂ += B₂[i]*xᵢ.d₂; F₂₃ += B₃[i]*xᵢ.d₂
            F₃₁ += B₁[i]*xᵢ.d₃; F₃₂ += B₂[i]*xᵢ.d₃; F₃₃ += B₃[i]*xᵢ.d₃
            u₁  += N[i]*xᵢ.d₁;  u₂  += N[i]*xᵢ.d₂;  u₃  += N[i]*xᵢ.d₃
        end

        for (i, xᵢ) in enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁; S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂; S₃₃ += Nₛ[i]*xᵢ.dₛ₃₃
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂; S₂₃ += Nₛ[i]*xᵢ.dₛ₂₃; S₁₃ += Nₛ[i]*xᵢ.dₛ₁₃
        end

        Δu₁ = g₁ - u₁
        Δu₂ = g₂ - u₂
        Δu₃ = g₃ - u₃

        v_Δu1 = n₁₁*Δu₁ + n₁₂*Δu₂ + n₁₃*Δu₃
        v_Δu2 = n₁₂*Δu₁ + n₂₂*Δu₂ + n₂₃*Δu₃
        v_Δu3 = n₁₃*Δu₁ + n₂₃*Δu₂ + n₃₃*Δu₃

        P₁₁ = F₁₁*n₁₁ + F₂₁*n₁₂ + F₃₁*n₁₃
        P₁₂ = F₁₁*n₁₂ + F₂₁*n₂₂ + F₃₁*n₂₃
        P₁₃ = F₁₁*n₁₃ + F₂₁*n₂₃ + F₃₁*n₃₃

        P₂₁ = F₁₂*n₁₁ + F₂₂*n₁₂ + F₃₂*n₁₃
        P₂₂ = F₁₂*n₁₂ + F₂₂*n₂₂ + F₃₂*n₂₃
        P₂₃ = F₁₂*n₁₃ + F₂₂*n₂₃ + F₃₂*n₃₃

        P₃₁ = F₁₃*n₁₁ + F₂₃*n₁₂ + F₃₃*n₁₃
        P₃₂ = F₁₃*n₁₂ + F₂₃*n₂₂ + F₃₃*n₂₃
        P₃₃ = F₁₃*n₁₃ + F₂₃*n₂₃ + F₃₃*n₃₃

        for (i, xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            Nₛ_i = Nₛ[i]

            w_res1 = F₁₁*v_Δu1 + F₂₁*v_Δu2 + F₃₁*v_Δu3
            w_res2 = F₁₂*v_Δu1 + F₂₂*v_Δu2 + F₃₂*v_Δu3
            w_res3 = F₁₃*v_Δu1 + F₂₃*v_Δu2 + F₃₃*v_Δu3

            f[6*I-5] += Nₛ_i * n₁ * w_res1 * 𝑤
            f[6*I-4] += Nₛ_i * n₂ * w_res2 * 𝑤
            f[6*I-3] += Nₛ_i * n₃ * w_res3 * 𝑤
            f[6*I-2] += Nₛ_i * (n₁ * w_res2 + n₂ * w_res1) * 𝑤
            f[6*I-1] += Nₛ_i * (n₂ * w_res3 + n₃ * w_res2) * 𝑤
            f[6*I]   += Nₛ_i * (n₁ * w_res3 + n₃ * w_res1) * 𝑤

            for (j, xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                N_j = N[j]
                B1_j = B₁[j]; B2_j = B₂[j]; B3_j = B₃[j]

                w_k1_1 = -v_Δu1 * B1_j + P₁₁ * N_j
                w_k2_1 = -v_Δu1 * B2_j + P₂₁ * N_j
                w_k3_1 = -v_Δu1 * B3_j + P₃₁ * N_j

                k[6*I-5, 3*J-2] += Nₛ_i * n₁ * w_k1_1 * 𝑤
                k[6*I-4, 3*J-2] += Nₛ_i * n₂ * w_k2_1 * 𝑤
                k[6*I-3, 3*J-2] += Nₛ_i * n₃ * w_k3_1 * 𝑤
                k[6*I-2, 3*J-2] += Nₛ_i * (n₁ * w_k2_1 + n₂ * w_k1_1) * 𝑤
                k[6*I-1, 3*J-2] += Nₛ_i * (n₂ * w_k3_1 + n₃ * w_k2_1) * 𝑤
                k[6*I,   3*J-2] += Nₛ_i * (n₁ * w_k3_1 + n₃ * w_k1_1) * 𝑤

                w_k1_2 = -v_Δu2 * B1_j + P₁₂ * N_j
                w_k2_2 = -v_Δu2 * B2_j + P₂₂ * N_j
                w_k3_2 = -v_Δu2 * B3_j + P₃₂ * N_j

                k[6*I-5, 3*J-1] += Nₛ_i * n₁ * w_k1_2 * 𝑤
                k[6*I-4, 3*J-1] += Nₛ_i * n₂ * w_k2_2 * 𝑤
                k[6*I-3, 3*J-1] += Nₛ_i * n₃ * w_k3_2 * 𝑤
                k[6*I-2, 3*J-1] += Nₛ_i * (n₁ * w_k2_2 + n₂ * w_k1_2) * 𝑤
                k[6*I-1, 3*J-1] += Nₛ_i * (n₂ * w_k3_2 + n₃ * w_k2_2) * 𝑤
                k[6*I,   3*J-1] += Nₛ_i * (n₁ * w_k3_2 + n₃ * w_k1_2) * 𝑤

                w_k1_3 = -v_Δu3 * B1_j + P₁₃ * N_j
                w_k2_3 = -v_Δu3 * B2_j + P₂₃ * N_j
                w_k3_3 = -v_Δu3 * B3_j + P₃₃ * N_j

                k[6*I-5, 3*J] += Nₛ_i * n₁ * w_k1_3 * 𝑤
                k[6*I-4, 3*J] += Nₛ_i * n₂ * w_k2_3 * 𝑤
                k[6*I-3, 3*J] += Nₛ_i * n₃ * w_k3_3 * 𝑤
                k[6*I-2, 3*J] += Nₛ_i * (n₁ * w_k2_3 + n₂ * w_k1_3) * 𝑤
                k[6*I-1, 3*J] += Nₛ_i * (n₂ * w_k3_3 + n₃ * w_k2_3) * 𝑤
                k[6*I,   3*J] += Nₛ_i * (n₁ * w_k3_3 + n₃ * w_k1_3) * 𝑤
            end
        end
    end
end

function ∫∫δSΔFngdxdy_HR_uS(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖

    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        n₁₁ = ξᵤ.n₁₁
        n₁₂ = ξᵤ.n₁₂
        n₂₂ = ξᵤ.n₂₂
        g₁ = ξᵤ.g₁
        g₂ = ξᵤ.g₂

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
        end

        Δu₁ = g₁-u₁
        Δu₂ = g₂-u₂

      for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼

                k[3*I-2,2*J-1] += Nₛ[i]*(n₁*n₁₁*F₁₁+n₁*n₁₂*F₂₁)*N[j]*𝑤
                k[3*I-2,2*J]   += Nₛ[i]*(n₁*n₁₂*F₁₁+n₁*n₂₂*F₂₁)*N[j]*𝑤
                k[3*I-1,2*J-1] += Nₛ[i]*(n₂*n₁₁*F₁₂+n₂*n₁₂*F₂₂)*N[j]*𝑤
                k[3*I-1,2*J]   += Nₛ[i]*(n₂*n₁₂*F₁₂+n₂*n₂₂*F₂₂)*N[j]*𝑤
                k[3*I,2*J-1]   += (Nₛ[i]*(n₂*n₁₁*F₁₁+n₂*n₁₂*F₂₁)*N[j]+Nₛ[i]*(n₁*n₁₁*F₁₂+n₁*n₁₂*F₂₂)*N[j])*𝑤
                k[3*I,2*J]     += (Nₛ[i]*(n₂*n₁₂*F₁₁+n₂*n₂₂*F₂₁)*N[j]+Nₛ[i]*(n₁*n₁₂*F₁₂+n₁*n₂₂*F₂₂)*N[j])*𝑤

            end
        end

    end
end

function ∫∫∫δSΔFngdxdydz_HR_uS(aₛ::T, aᵤ::S, k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒; 𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒; 𝓖ᵤ = aᵤ.𝓖
   Nₛ = zeros(length(𝓒ₛ))
    for (ξₛ, ξᵤ) in zip(𝓖ₛ, 𝓖ᵤ)

        B₁ = ξᵤ[:∂𝝭∂x]; B₂ = ξᵤ[:∂𝝭∂y]; B₃ = ξᵤ[:∂𝝭∂z]

        x = ξᵤ.x
        y = ξᵤ.y
        z = ξᵤ.z
        eval_piecewise_3d!(Nₛ, x, y, z)
        N  = ξᵤ[:𝝭]
        𝑤  = ξₛ.𝑤

        n₁ = ξᵤ.n₁; n₂ = ξᵤ.n₂; n₃ = ξᵤ.n₃

        n₁₁ = ξᵤ.n₁₁; n₂₂ = ξᵤ.n₂₂; n₃₃ = ξᵤ.n₃₃
        n₁₂ = ξᵤ.n₁₂; n₂₃ = ξᵤ.n₂₃; n₁₃ = ξᵤ.n₁₃

        F₁₁ = 1.0; F₁₂ = 0.0; F₁₃ = 0.0
        F₂₁ = 0.0; F₂₂ = 1.0; F₂₃ = 0.0
        F₃₁ = 0.0; F₃₂ = 0.0; F₃₃ = 1.0

        for (i, xᵢ) in enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁; F₁₂ += B₂[i]*xᵢ.d₁; F₁₃ += B₃[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂; F₂₂ += B₂[i]*xᵢ.d₂; F₂₃ += B₃[i]*xᵢ.d₂
            F₃₁ += B₁[i]*xᵢ.d₃; F₃₂ += B₂[i]*xᵢ.d₃; F₃₃ += B₃[i]*xᵢ.d₃
        end

        P₁₁ = F₁₁*n₁₁ + F₂₁*n₁₂ + F₃₁*n₁₃
        P₁₂ = F₁₁*n₁₂ + F₂₁*n₂₂ + F₃₁*n₂₃
        P₁₃ = F₁₁*n₁₃ + F₂₁*n₂₃ + F₃₁*n₃₃

        P₂₁ = F₁₂*n₁₁ + F₂₂*n₁₂ + F₃₂*n₁₃
        P₂₂ = F₁₂*n₁₂ + F₂₂*n₂₂ + F₃₂*n₂₃
        P₂₃ = F₁₂*n₁₃ + F₂₂*n₂₃ + F₃₂*n₃₃

        P₃₁ = F₁₃*n₁₁ + F₂₃*n₁₂ + F₃₃*n₁₃
        P₃₂ = F₁₃*n₁₂ + F₂₃*n₂₂ + F₃₃*n₂₃
        P₃₃ = F₁₃*n₁₃ + F₂₃*n₂₃ + F₃₃*n₃₃

        for (i, xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            Nₛ_i = Nₛ[i]

            for (j, xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                N_j = N[j]

                w1_1 = P₁₁ * N_j
                w2_1 = P₂₁ * N_j
                w3_1 = P₃₁ * N_j

                k[6*I-5, 3*J-2] += Nₛ_i * n₁ * w1_1 * 𝑤
                k[6*I-4, 3*J-2] += Nₛ_i * n₂ * w2_1 * 𝑤
                k[6*I-3, 3*J-2] += Nₛ_i * n₃ * w3_1 * 𝑤
                k[6*I-2, 3*J-2] += Nₛ_i * (n₁ * w2_1 + n₂ * w1_1) * 𝑤
                k[6*I-1, 3*J-2] += Nₛ_i * (n₂ * w3_1 + n₃ * w2_1) * 𝑤
                k[6*I,   3*J-2] += Nₛ_i * (n₁ * w3_1 + n₃ * w1_1) * 𝑤

                w1_2 = P₁₂ * N_j
                w2_2 = P₂₂ * N_j
                w3_2 = P₃₂ * N_j

                k[6*I-5, 3*J-1] += Nₛ_i * n₁ * w1_2 * 𝑤
                k[6*I-4, 3*J-1] += Nₛ_i * n₂ * w2_2 * 𝑤
                k[6*I-3, 3*J-1] += Nₛ_i * n₃ * w3_2 * 𝑤
                k[6*I-2, 3*J-1] += Nₛ_i * (n₁ * w2_2 + n₂ * w1_2) * 𝑤
                k[6*I-1, 3*J-1] += Nₛ_i * (n₂ * w3_2 + n₃ * w2_2) * 𝑤
                k[6*I,   3*J-1] += Nₛ_i * (n₁ * w3_2 + n₃ * w1_2) * 𝑤

                w1_3 = P₁₃ * N_j
                w2_3 = P₂₃ * N_j
                w3_3 = P₃₃ * N_j

                k[6*I-5, 3*J] += Nₛ_i * n₁ * w1_3 * 𝑤
                k[6*I-4, 3*J] += Nₛ_i * n₂ * w2_3 * 𝑤
                k[6*I-3, 3*J] += Nₛ_i * n₃ * w3_3 * 𝑤
                k[6*I-2, 3*J] += Nₛ_i * (n₁ * w2_3 + n₂ * w1_3) * 𝑤
                k[6*I-1, 3*J] += Nₛ_i * (n₂ * w3_3 + n₃ * w2_3) * 𝑤
                k[6*I,   3*J] += Nₛ_i * (n₁ * w3_3 + n₃ * w1_3) * 𝑤
            end
        end
    end
end

function ∫∫ΔSδFngdxdy_HR(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖

    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        n₁₁ = ξᵤ.n₁₁
        n₁₂ = ξᵤ.n₁₂
        n₂₂ = ξᵤ.n₂₂
        g₁ = ξᵤ.g₁
        g₂ = ξᵤ.g₂

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
        end

        Δu₁ = g₁-u₁
        Δu₂ = g₂-u₂

      for (i,xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼

                k[2*I-1,2*J-1] += (N[i]*B₁[j]*(S₁₁*n₁+S₁₂*n₂)+N[i]*B₂[j]*(S₁₂*n₁+S₂₂*n₂))*n₁₁*𝑤
                k[2*I-1,2*J]   += (N[i]*B₁[j]*(S₁₁*n₁+S₁₂*n₂)+N[i]*B₂[j]*(S₁₂*n₁+S₂₂*n₂))*n₁₂*𝑤
                k[2*I,2*J-1]   += (N[i]*B₁[j]*(S₁₁*n₁+S₁₂*n₂)+N[i]*B₂[j]*(S₁₂*n₁+S₂₂*n₂))*n₁₂*𝑤
                k[2*I,2*J]     += (N[i]*B₁[j]*(S₁₁*n₁+S₁₂*n₂)+N[i]*B₂[j]*(S₁₂*n₁+S₂₂*n₂))*n₂₂*𝑤

            end

        end

                for (i,xᵢ) in enumerate(𝓒ᵤ)
                     I = xᵢ.𝐼

                        f[2*I-1] +=0.0
                        f[2*I]   +=0.0      

                 end
    end
end

function ∫∫∫ΔSδFngdxdydz_HR(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
   Nₛ = zeros(length(𝓒ₛ))
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        B₃ = ξᵤ[:∂𝝭∂z]
        𝑤 = ξₛ.𝑤

        x = ξᵤ.x
        y = ξᵤ.y
        z = ξᵤ.z
        eval_piecewise_3d!(Nₛ, x, y, z)
        N = ξᵤ[:𝝭]

        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        n₃ = ξᵤ.n₃

        n₁₁ = ξᵤ.n₁₁
        n₁₂ = ξᵤ.n₁₂
        n₁₃ = ξᵤ.n₁₃
        n₂₂ = ξᵤ.n₂₂
        n₂₃ = ξᵤ.n₂₃
        n₃₃ = ξᵤ.n₃₃

        g₁ = ξᵤ.g₁
        g₂ = ξᵤ.g₂
        g₃ = ξᵤ.g₃

        F₁₁ = 1.0; F₁₂ = 0.0; F₁₃ = 0.0
        F₂₁ = 0.0; F₂₂ = 1.0; F₂₃ = 0.0
        F₃₁ = 0.0; F₃₂ = 0.0; F₃₃ = 1.0

        u₁ = 0.0; u₂ = 0.0; u₃ = 0.0
        S₁₁ = 0.0; S₂₂ = 0.0; S₃₃ = 0.0
        S₁₂ = 0.0; S₂₃ = 0.0; S₁₃ = 0.0

        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₁₃ += B₃[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            F₂₃ += B₃[i]*xᵢ.d₂
            F₃₁ += B₁[i]*xᵢ.d₃
            F₃₂ += B₂[i]*xᵢ.d₃
            F₃₃ += B₃[i]*xᵢ.d₃

            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
            u₃  += N[i]*xᵢ.d₃
        end

        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₃₃ += Nₛ[i]*xᵢ.dₛ₃₃
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
           S₂₃ += Nₛ[i]*xᵢ.dₛ₂₃
           S₁₃ += Nₛ[i]*xᵢ.dₛ₁₃
        end

        Δu₁ = g₁ - u₁
        Δu₂ = g₂ - u₂
        Δu₃ = g₃ - u₃

      for (i,xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼

                k[3*I-2,3*J-2] += (N[i]*B₁[j]*(S₁₁*n₁ + S₁₂*n₂ + S₁₃*n₃) + N[i]*B₂[j]*(S₁₂*n₁ + S₂₂*n₂ + S₂₃*n₃) + N[i]*B₃[j]*(S₁₃*n₁ + S₂₃*n₂ + S₃₃*n₃))*n₁₁*𝑤
                k[3*I-2,3*J-1] += (N[i]*B₁[j]*(S₁₁*n₁ + S₁₂*n₂ + S₁₃*n₃) + N[i]*B₂[j]*(S₁₂*n₁ + S₂₂*n₂ + S₂₃*n₃) + N[i]*B₃[j]*(S₁₃*n₁ + S₂₃*n₂ + S₃₃*n₃))*n₁₂*𝑤
                k[3*I-2,3*J]   += (N[i]*B₁[j]*(S₁₁*n₁ + S₁₂*n₂ + S₁₃*n₃) + N[i]*B₂[j]*(S₁₂*n₁ + S₂₂*n₂ + S₂₃*n₃) + N[i]*B₃[j]*(S₁₃*n₁ + S₂₃*n₂ + S₃₃*n₃))*n₁₃*𝑤

                k[3*I-1,3*J-2] += (N[i]*B₁[j]*(S₁₁*n₁ + S₁₂*n₂ + S₁₃*n₃) + N[i]*B₂[j]*(S₁₂*n₁ + S₂₂*n₂ + S₂₃*n₃) + N[i]*B₃[j]*(S₁₃*n₁ + S₂₃*n₂ + S₃₃*n₃))*n₁₂*𝑤
                k[3*I-1,3*J-1] += (N[i]*B₁[j]*(S₁₁*n₁ + S₁₂*n₂ + S₁₃*n₃) + N[i]*B₂[j]*(S₁₂*n₁ + S₂₂*n₂ + S₂₃*n₃) + N[i]*B₃[j]*(S₁₃*n₁ + S₂₃*n₂ + S₃₃*n₃))*n₂₂*𝑤
                k[3*I-1,3*J]   += (N[i]*B₁[j]*(S₁₁*n₁ + S₁₂*n₂ + S₁₃*n₃) + N[i]*B₂[j]*(S₁₂*n₁ + S₂₂*n₂ + S₂₃*n₃) + N[i]*B₃[j]*(S₁₃*n₁ + S₂₃*n₂ + S₃₃*n₃))*n₂₃*𝑤

                k[3*I,3*J-2]   += (N[i]*B₁[j]*(S₁₁*n₁ + S₁₂*n₂ + S₁₃*n₃) + N[i]*B₂[j]*(S₁₂*n₁ + S₂₂*n₂ + S₂₃*n₃) + N[i]*B₃[j]*(S₁₃*n₁ + S₂₃*n₂ + S₃₃*n₃))*n₁₃*𝑤
                k[3*I,3*J-1]   += (N[i]*B₁[j]*(S₁₁*n₁ + S₁₂*n₂ + S₁₃*n₃) + N[i]*B₂[j]*(S₁₂*n₁ + S₂₂*n₂ + S₂₃*n₃) + N[i]*B₃[j]*(S₁₃*n₁ + S₂₃*n₂ + S₃₃*n₃))*n₂₃*𝑤
                k[3*I,3*J]     += (N[i]*B₁[j]*(S₁₁*n₁ + S₁₂*n₂ + S₁₃*n₃) + N[i]*B₂[j]*(S₁₂*n₁ + S₂₂*n₂ + S₂₃*n₃) + N[i]*B₃[j]*(S₁₃*n₁ + S₂₃*n₂ + S₃₃*n₃))*n₃₃*𝑤

            end

        end

        for (i,xᵢ) in enumerate(𝓒ᵤ)
             I = xᵢ.𝐼

                f[3*I-2] += 0.0
                f[3*I-1] += 0.0      
                f[3*I]   += 0.0      

         end
    end
end

function ∫∫ΔSδFngdxdy_HR_SaintVenantKirchhoff(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖

    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        n₁₁ = ξᵤ.n₁₁
        n₁₂ = ξᵤ.n₁₂
        n₂₂ = ξᵤ.n₂₂
        g₁ = ξᵤ.g₁
        g₂ = ξᵤ.g₂

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
        end

        Δu₁ = g₁-u₁
        Δu₂ = g₂-u₂

      for (i,xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼

                k[2*I-1,2*J-1] += (N[i]*B₁[j]*(S₁₁*n₁+S₁₂*n₂)+N[i]*B₂[j]*(S₁₂*n₁+S₂₂*n₂))*n₁₁*𝑤
                k[2*I-1,2*J]   += (N[i]*B₁[j]*(S₁₁*n₁+S₁₂*n₂)+N[i]*B₂[j]*(S₁₂*n₁+S₂₂*n₂))*n₁₂*𝑤
                k[2*I,2*J-1]   += (N[i]*B₁[j]*(S₁₁*n₁+S₁₂*n₂)+N[i]*B₂[j]*(S₁₂*n₁+S₂₂*n₂))*n₁₂*𝑤
                k[2*I,2*J]     += (N[i]*B₁[j]*(S₁₁*n₁+S₁₂*n₂)+N[i]*B₂[j]*(S₁₂*n₁+S₂₂*n₂))*n₂₂*𝑤

                k[2*I-1,2*J-1] += (N[j]*B₁[i]*(S₁₁*n₁+S₁₂*n₂)+N[j]*B₂[i]*(S₁₂*n₁+S₂₂*n₂))*n₁₁*𝑤
                k[2*I-1,2*J]   += (N[j]*B₁[i]*(S₁₁*n₁+S₁₂*n₂)+N[j]*B₂[i]*(S₁₂*n₁+S₂₂*n₂))*n₁₂*𝑤
                k[2*I,2*J-1]   += (N[j]*B₁[i]*(S₁₁*n₁+S₁₂*n₂)+N[j]*B₂[i]*(S₁₂*n₁+S₂₂*n₂))*n₁₂*𝑤
                k[2*I,2*J]     += (N[j]*B₁[i]*(S₁₁*n₁+S₁₂*n₂)+N[j]*B₂[i]*(S₁₂*n₁+S₂₂*n₂))*n₂₂*𝑤
            end

        end

                for (i,xᵢ) in enumerate(𝓒ᵤ)
                     I = xᵢ.𝐼

                        f[2*I-1] += ((Δu₁*B₁[i]*(S₁₁*n₁+S₁₂*n₂)+Δu₁*B₂[i]*(S₁₂*n₁+S₂₂*n₂))*n₁₁ + (Δu₂*B₁[i]*(S₁₁*n₁+S₁₂*n₂)+Δu₂*B₂[i]*(S₁₂*n₁+S₂₂*n₂))*n₁₂)*𝑤
                        f[2*I]   += ((Δu₁*B₁[i]*(S₁₁*n₁+S₁₂*n₂)+Δu₁*B₂[i]*(S₁₂*n₁+S₂₂*n₂))*n₁₂ + (Δu₂*B₁[i]*(S₁₁*n₁+S₁₂*n₂)+Δu₂*B₂[i]*(S₁₂*n₁+S₂₂*n₂))*n₂₂)*𝑤

                 end
    end
end

function ∫∫∫ΔSδFngdxdydz_HR_SaintVenantKirchhoff(aₛ::T, aᵤ::S, k::AbstractMatrix{Float64}, f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒; 𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒; 𝓖ᵤ = aᵤ.𝓖

    for (ξₛ, ξᵤ) in zip(𝓖ₛ, 𝓖ᵤ)

        B₁ = ξᵤ[:∂𝝭∂x]; B₂ = ξᵤ[:∂𝝭∂y]; B₃ = ξᵤ[:∂𝝭∂z]
        Nₛ = ξₛ[:𝝭]
        N  = ξᵤ[:𝝭]
        𝑤  = ξₛ.𝑤

        n₁ = ξᵤ.n₁; n₂ = ξᵤ.n₂; n₃ = ξᵤ.n₃

        n₁₁ = ξᵤ.n₁₁; n₂₂ = ξᵤ.n₂₂; n₃₃ = ξᵤ.n₃₃
        n₁₂ = ξᵤ.n₁₂; n₂₃ = ξᵤ.n₂₃; n₁₃ = ξᵤ.n₁₃

        S₁₁ = 0.0; S₂₂ = 0.0; S₃₃ = 0.0
        S₁₂ = 0.0; S₂₃ = 0.0; S₁₃ = 0.0

        for (i, xᵢ) in enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₃₃ += Nₛ[i]*xᵢ.dₛ₃₃
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
           S₂₃ += Nₛ[i]*xᵢ.dₛ₂₃
           S₁₃ += Nₛ[i]*xᵢ.dₛ₁₃
        end

        s₁ = S₁₁*n₁ + S₁₂*n₂ + S₁₃*n₃
        s₂ = S₁₂*n₁ + S₂₂*n₂ + S₂₃*n₃
        s₃ = S₁₃*n₁ + S₂₃*n₂ + S₃₃*n₃

        for (j, xⱼ) in enumerate(𝓒ᵤ)
            J = xⱼ.𝐼

            w_J = (B₁[j]*s₁ + B₂[j]*s₂ + B₃[j]*s₃) * 𝑤

            w_1 = w_J * n₁₁; w_2 = w_J * n₁₂; w_3 = w_J * n₁₃
            w_4 = w_J * n₁₂; w_5 = w_J * n₂₂; w_6 = w_J * n₂₃
            w_7 = w_J * n₁₃; w_8 = w_J * n₂₃; w_9 = w_J * n₃₃

            for (i, xᵢ) in enumerate(𝓒ᵤ)
                I = xᵢ.𝐼
                N_i = N[i]

                k[3*I-2, 3*J-2] += N_i * w_1
                k[3*I-2, 3*J-1] += N_i * w_2
                k[3*I-2, 3*J]   += N_i * w_3

                k[3*I-1, 3*J-2] += N_i * w_4
                k[3*I-1, 3*J-1] += N_i * w_5
                k[3*I-1, 3*J]   += N_i * w_6

                k[3*I,   3*J-2] += N_i * w_7
                k[3*I,   3*J-1] += N_i * w_8
                k[3*I,   3*J]   += N_i * w_9
            end
        end

    end
end

function ∫∫SFnδudxdy_HR(aₛ::T,aᵤ::S,f::AbstractVector{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖

    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        n₁₁ = ξᵤ.n₁₁
        n₁₂ = ξᵤ.n₁₂
        n₂₂ = ξᵤ.n₂₂
        g₁ = ξᵤ.g₁
        g₂ = ξᵤ.g₂

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
        end

                for (i,xᵢ) in enumerate(𝓒ᵤ)
                     I = xᵢ.𝐼

                        f[2*I-1] += ((N[i]*F₁₁*(S₁₁*n₁+S₁₂*n₂)+N[i]*F₁₂*(S₁₂*n₁+S₂₂*n₂))*n₁₁ + (N[i]*F₂₁*(S₁₁*n₁+S₁₂*n₂)+N[i]*F₂₂*(S₁₂*n₁+S₂₂*n₂))*n₁₂)*𝑤

                        f[2*I]   += ((N[i]*F₁₁*(S₁₁*n₁+S₁₂*n₂)+N[i]*F₁₂*(S₁₂*n₁+S₂₂*n₂))*n₁₂ + (N[i]*F₂₁*(S₁₁*n₁+S₁₂*n₂)+N[i]*F₂₂*(S₁₂*n₁+S₂₂*n₂))*n₂₂)*𝑤

                 end
    end
end

function ∫∫∫SFnδudxdydz_HR(aₛ::T,aᵤ::S,f::AbstractVector{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖

    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        B₃ = ξᵤ[:∂𝝭∂z]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]

        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        n₃ = ξᵤ.n₃

        n₁₁ = ξᵤ.n₁₁
        n₁₂ = ξᵤ.n₁₂
        n₁₃ = ξᵤ.n₁₃
        n₂₂ = ξᵤ.n₂₂
        n₂₃ = ξᵤ.n₂₃
        n₃₃ = ξᵤ.n₃₃

        g₁ = ξᵤ.g₁
        g₂ = ξᵤ.g₂
        g₃ = ξᵤ.g₃

        F₁₁ = 1.0; F₁₂ = 0.0; F₁₃ = 0.0
        F₂₁ = 0.0; F₂₂ = 1.0; F₂₃ = 0.0
        F₃₁ = 0.0; F₃₂ = 0.0; F₃₃ = 1.0

        u₁ = 0.0; u₂ = 0.0; u₃ = 0.0

        S₁₁ = 0.0; S₂₂ = 0.0; S₃₃ = 0.0
        S₁₂ = 0.0; S₂₃ = 0.0; S₁₃ = 0.0

        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₁₃ += B₃[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            F₂₃ += B₃[i]*xᵢ.d₂
            F₃₁ += B₁[i]*xᵢ.d₃
            F₃₂ += B₂[i]*xᵢ.d₃
            F₃₃ += B₃[i]*xᵢ.d₃

            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
            u₃  += N[i]*xᵢ.d₃
        end

        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₃₃ += Nₛ[i]*xᵢ.dₛ₃₃
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
           S₂₃ += Nₛ[i]*xᵢ.dₛ₂₃
           S₁₃ += Nₛ[i]*xᵢ.dₛ₁₃
        end

        for (i,xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼

            f[3*I-2] += ( ((N[i]*F₁₁*(S₁₁*n₁ + S₁₂*n₂ + S₁₃*n₃) + N[i]*F₁₂*(S₁₂*n₁ + S₂₂*n₂ + S₂₃*n₃) + N[i]*F₁₃*(S₁₃*n₁ + S₂₃*n₂ + S₃₃*n₃))*n₁₁) + 
                          ((N[i]*F₂₁*(S₁₁*n₁ + S₁₂*n₂ + S₁₃*n₃) + N[i]*F₂₂*(S₁₂*n₁ + S₂₂*n₂ + S₂₃*n₃) + N[i]*F₂₃*(S₁₃*n₁ + S₂₃*n₂ + S₃₃*n₃))*n₁₂) + 
                          ((N[i]*F₃₁*(S₁₁*n₁ + S₁₂*n₂ + S₁₃*n₃) + N[i]*F₃₂*(S₁₂*n₁ + S₂₂*n₂ + S₂₃*n₃) + N[i]*F₃₃*(S₁₃*n₁ + S₂₃*n₂ + S₃₃*n₃))*n₁₃) )*𝑤

            f[3*I-1] += ( ((N[i]*F₁₁*(S₁₁*n₁ + S₁₂*n₂ + S₁₃*n₃) + N[i]*F₁₂*(S₁₂*n₁ + S₂₂*n₂ + S₂₃*n₃) + N[i]*F₁₃*(S₁₃*n₁ + S₂₃*n₂ + S₃₃*n₃))*n₁₂) + 
                          ((N[i]*F₂₁*(S₁₁*n₁ + S₁₂*n₂ + S₁₃*n₃) + N[i]*F₂₂*(S₁₂*n₁ + S₂₂*n₂ + S₂₃*n₃) + N[i]*F₂₃*(S₁₃*n₁ + S₂₃*n₂ + S₃₃*n₃))*n₂₂) + 
                          ((N[i]*F₃₁*(S₁₁*n₁ + S₁₂*n₂ + S₁₃*n₃) + N[i]*F₃₂*(S₁₂*n₁ + S₂₂*n₂ + S₂₃*n₃) + N[i]*F₃₃*(S₁₃*n₁ + S₂₃*n₂ + S₃₃*n₃))*n₂₃) )*𝑤

            f[3*I]   += ( ((N[i]*F₁₁*(S₁₁*n₁ + S₁₂*n₂ + S₁₃*n₃) + N[i]*F₁₂*(S₁₂*n₁ + S₂₂*n₂ + S₂₃*n₃) + N[i]*F₁₃*(S₁₃*n₁ + S₂₃*n₂ + S₃₃*n₃))*n₁₃) + 
                          ((N[i]*F₂₁*(S₁₁*n₁ + S₁₂*n₂ + S₁₃*n₃) + N[i]*F₂₂*(S₁₂*n₁ + S₂₂*n₂ + S₂₃*n₃) + N[i]*F₂₃*(S₁₃*n₁ + S₂₃*n₂ + S₃₃*n₃))*n₂₃) + 
                          ((N[i]*F₃₁*(S₁₁*n₁ + S₁₂*n₂ + S₁₃*n₃) + N[i]*F₃₂*(S₁₂*n₁ + S₂₂*n₂ + S₂₃*n₃) + N[i]*F₃₃*(S₁₃*n₁ + S₂₃*n₂ + S₃₃*n₃))*n₃₃) )*𝑤
        end
    end
end
function ∫∫SFnδudxdy_HR_SaintVenantKirchhoff(aₛ::T,aᵤ::S,f::AbstractVector{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖

    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        n₁₁ = ξᵤ.n₁₁
        n₁₂ = ξᵤ.n₁₂
        n₂₂ = ξᵤ.n₂₂
        g₁ = ξᵤ.g₁
        g₂ = ξᵤ.g₂

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
        end

                for (i,xᵢ) in enumerate(𝓒ᵤ)
                     I = xᵢ.𝐼

                        f[2*I-1] += ((N[i]*F₁₁*(S₁₁*n₁+S₁₂*n₂)+N[i]*F₁₂*(S₁₂*n₁+S₂₂*n₂))*n₁₁ + (N[i]*F₂₁*(S₁₁*n₁+S₁₂*n₂)+N[i]*F₂₂*(S₁₂*n₁+S₂₂*n₂))*n₁₂)*𝑤

                        f[2*I]   += ((N[i]*F₁₁*(S₁₁*n₁+S₁₂*n₂)+N[i]*F₁₂*(S₁₂*n₁+S₂₂*n₂))*n₁₂ + (N[i]*F₂₁*(S₁₁*n₁+S₁₂*n₂)+N[i]*F₂₂*(S₁₂*n₁+S₂₂*n₂))*n₂₂)*𝑤

                 end
    end
end

function ∫∫∫SFnδudxdydz_HR_SaintVenantKirchhoff(aₛ::T, aᵤ::S, f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒; 𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒; 𝓖ᵤ = aᵤ.𝓖

    for (ξₛ, ξᵤ) in zip(𝓖ₛ, 𝓖ᵤ)

        B₁ = ξᵤ[:∂𝝭∂x]; B₂ = ξᵤ[:∂𝝭∂y]; B₃ = ξᵤ[:∂𝝭∂z]
        Nₛ = ξₛ[:𝝭]
        N  = ξᵤ[:𝝭]
        𝑤  = ξₛ.𝑤

        n₁ = ξᵤ.n₁; n₂ = ξᵤ.n₂; n₃ = ξᵤ.n₃

        n₁₁ = ξᵤ.n₁₁; n₂₂ = ξᵤ.n₂₂; n₃₃ = ξᵤ.n₃₃
        n₁₂ = ξᵤ.n₁₂; n₂₃ = ξᵤ.n₂₃; n₁₃ = ξᵤ.n₁₃

        g₁ = ξᵤ.g₁; g₂ = ξᵤ.g₂; g₃ = ξᵤ.g₃

        F₁₁ = 1.0; F₁₂ = 0.0; F₁₃ = 0.0
        F₂₁ = 0.0; F₂₂ = 1.0; F₂₃ = 0.0
        F₃₁ = 0.0; F₃₂ = 0.0; F₃₃ = 1.0

        u₁ = 0.0; u₂ = 0.0; u₃ = 0.0

        S₁₁ = 0.0; S₂₂ = 0.0; S₃₃ = 0.0
        S₁₂ = 0.0; S₂₃ = 0.0; S₁₃ = 0.0

        for (i, xᵢ) in enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁; F₁₂ += B₂[i]*xᵢ.d₁; F₁₃ += B₃[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂; F₂₂ += B₂[i]*xᵢ.d₂; F₂₃ += B₃[i]*xᵢ.d₂
            F₃₁ += B₁[i]*xᵢ.d₃; F₃₂ += B₂[i]*xᵢ.d₃; F₃₃ += B₃[i]*xᵢ.d₃
            u₁  += N[i]*xᵢ.d₁;  u₂  += N[i]*xᵢ.d₂;  u₃  += N[i]*xᵢ.d₃
        end

        for (i, xᵢ) in enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₃₃ += Nₛ[i]*xᵢ.dₛ₃₃
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
           S₂₃ += Nₛ[i]*xᵢ.dₛ₂₃
           S₁₃ += Nₛ[i]*xᵢ.dₛ₁₃
        end

        s₁ = S₁₁*n₁ + S₁₂*n₂ + S₁₃*n₃
        s₂ = S₁₂*n₁ + S₂₂*n₂ + S₂₃*n₃
        s₃ = S₁₃*n₁ + S₂₃*n₂ + S₃₃*n₃

        T₁ = F₁₁*s₁ + F₁₂*s₂ + F₁₃*s₃
        T₂ = F₂₁*s₁ + F₂₂*s₂ + F₂₃*s₃
        T₃ = F₃₁*s₁ + F₃₂*s₂ + F₃₃*s₃

        for (i, xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼

            f[3*I-2] += N[i] * (T₁*n₁₁ + T₂*n₁₂ + T₃*n₁₃) * 𝑤
            f[3*I-1] += N[i] * (T₁*n₁₂ + T₂*n₂₂ + T₃*n₂₃) * 𝑤
            f[3*I]   += N[i] * (T₁*n₁₃ + T₂*n₂₃ + T₃*n₃₃) * 𝑤
        end
    end
end

function Δ∫∫EᵢⱼSᵢⱼdxdy_SaintVenantKirchhoff(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    F = zeros(2,2)
    S = zeros(2,2)

    for ξ in 𝓖
        E_mod = ξ.Ē
        ν = ξ.ν̄ 

        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤

        λ = E_mod * ν / ((1 + ν) * (1 - 2 * ν))
        μ = E_mod / (2 * (1 + ν))

        fill!(F, 0.0)
        F[1,1] = 1.0; F[2,2] = 1.0

        for (i, xᵢ) in enumerate(𝓒)
            d = [xᵢ.d₁, xᵢ.d₂]

            F[1,1] += d[1] * B₁[i]; F[1,2] += d[1] * B₂[i]
            F[2,1] += d[2] * B₁[i]; F[2,2] += d[2] * B₂[i]
        end

        C = F' * F
        E_strain = 0.5 * (C - I)

        trE = tr(E_strain)
        S .= λ * trE * I + 2 * μ * E_strain

        for (i, xᵢ) in enumerate(𝓒)
            I_idx = xᵢ.𝐼
            ∇N_i = [B₁[i], B₂[i]] 

            for (j, xⱼ) in enumerate(𝓒)
                J_idx = xⱼ.𝐼
                ∇N_j = [B₁[j], B₂[j]] 

                s_geo = dot(∇N_i, S, ∇N_j) * 𝑤

                grad_dot = dot(∇N_i, ∇N_j)

                for a in 1:2 
                    F_a = view(F, a, :) 
                    Fa_dot_Ni = dot(F_a, ∇N_i)
                    Fa_dot_Nj = dot(F_a, ∇N_j)

                    for b in 1:2
                        F_b = view(F, b, :) 
                        Fb_dot_Nj = dot(F_b, ∇N_j)
                        Fb_dot_Ni = dot(F_b, ∇N_i)
                        Fa_dot_Fb = dot(F_a, F_b)

                        k_val = (λ * Fa_dot_Ni * Fb_dot_Nj + 
                                 μ * Fa_dot_Fb * grad_dot + 
                                 μ * Fb_dot_Ni * Fa_dot_Nj) * 𝑤

                        if a == b
                            k_val += s_geo
                        end

                        row = 2 * I_idx - (2 - a)
                        col = 2 * J_idx - (2 - b)
                        k[row, col] += k_val
                    end
                end
            end
        end
    end
end

function Δ∫∫∫EᵢⱼSᵢⱼdxdydz_SaintVenantKirchhoff(ap::T, k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    for ξ in 𝓖
        E = ξ.E
        ν = ξ.ν
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z]  
        𝑤 = ξ.𝑤

        λ = E*ν/((1+ν)*(1-2*ν))  
        μ = E/(2*(1+ν))          

       Cᵢᵢᵢᵢ = λ + 2*μ
       Cᵢᵢⱼⱼ = λ
       Cᵢⱼᵢⱼ = μ

        F₁₁ = 1.0; F₁₂ = 0.0; F₁₃ = 0.0
        F₂₁ = 0.0; F₂₂ = 1.0; F₂₃ = 0.0
        F₃₁ = 0.0; F₃₂ = 0.0; F₃₃ = 1.0

        for (i,xᵢ) in  enumerate(𝓒)

            F₁₁ += B₁[i]*xᵢ.d₁; F₁₂ += B₂[i]*xᵢ.d₁; F₁₃ += B₃[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂; F₂₂ += B₂[i]*xᵢ.d₂; F₂₃ += B₃[i]*xᵢ.d₂
            F₃₁ += B₁[i]*xᵢ.d₃; F₃₂ += B₂[i]*xᵢ.d₃; F₃₃ += B₃[i]*xᵢ.d₃
        end

        E₁₁ = 0.5*(F₁₁*F₁₁+F₂₁*F₂₁+F₃₁*F₃₁-1.0)
        E₂₂ = 0.5*(F₁₂*F₁₂+F₂₂*F₂₂+F₃₂*F₃₂-1.0)
        E₃₃ = 0.5*(F₁₃*F₁₃+F₂₃*F₂₃+F₃₃*F₃₃-1.0) 
        E₁₂ = 0.5*(F₁₁*F₁₂+F₂₁*F₂₂+F₃₁*F₃₂)
        E₁₃ = 0.5*(F₁₁*F₁₃+F₂₁*F₂₃+F₃₁*F₃₃)    
        E₂₃ = 0.5*(F₁₂*F₁₃+F₂₂*F₂₃+F₃₂*F₃₃)    

       S₁₁ = λ*(E₁₁ + E₂₂ + E₃₃) + 2*μ*E₁₁
       S₂₂ = λ*(E₁₁ + E₂₂ + E₃₃) + 2*μ*E₂₂
       S₃₃ = λ*(E₁₁ + E₂₂ + E₃₃) + 2*μ*E₃₃
       S₁₂ = 2*μ*E₁₂
       S₁₃ = 2*μ*E₁₃
       S₂₃ = 2*μ*E₂₃ 

       for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] +=((B₁[i]*F₁₁*B₁[j]*F₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂ + B₃[i]*F₁₃*B₃[j]*F₁₃)*Cᵢᵢᵢᵢ
                               +  (B₂[i]*F₁₂*B₁[j]*F₁₁ + B₃[i]*F₁₃*B₁[j]*F₁₁ + B₁[i]*F₁₁*B₂[j]*F₁₂ + B₃[i]*F₁₃*B₂[j]*F₁₂ + B₁[i]*F₁₁*B₃[j]*F₁₃ + B₂[i]*F₁₂*B₃[j]*F₁₃)*Cᵢᵢⱼⱼ
                               +  ((B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)+(B₃[i]*F₁₁+B₁[i]*F₁₃)*(B₃[j]*F₁₁+B₁[j]*F₁₃) + (B₃[i]*F₁₂+B₂[i]*F₁₃)*(B₃[j]*F₁₂+B₂[j]*F₁₃) )*Cᵢⱼᵢⱼ
                               +   B₁[i]*B₁[j]*S₁₁+B₂[i]*B₂[j]*S₂₂+B₃[i]*B₃[j]*S₃₃ + (B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂ + (B₁[i]*B₃[j]+B₃[i]*B₁[j])*S₁₃ + (B₃[i]*B₂[j]+B₂[i]*B₃[j])*S₂₃)*𝑤

                k[3*I-2,3*J-1] += ((B₁[i]*F₁₁*B₁[j]*F₂₁ + B₂[i]*F₁₂*B₂[j]*F₂₂ + B₃[i]*F₁₃*B₃[j]*F₂₃)*Cᵢᵢᵢᵢ
                               +   (B₂[i]*F₁₂*B₁[j]*F₂₁ + B₃[i]*F₁₃*B₁[j]*F₂₁ + B₁[i]*F₁₁*B₂[j]*F₂₂ + B₃[i]*F₁₃*B₂[j]*F₂₂ + B₁[i]*F₁₁*B₃[j]*F₂₃ + B₂[i]*F₁₂*B₃[j]*F₂₃)*Cᵢᵢⱼⱼ
                               +   ((B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁) + (B₃[i]*F₁₁+B₁[i]*F₁₃)*(B₃[j]*F₂₁+B₁[j]*F₂₃) + (B₃[i]*F₁₂+B₂[i]*F₁₃)*(B₃[j]*F₂₂+B₂[j]*F₂₃))*Cᵢⱼᵢⱼ)*𝑤

                k[3*I-2,3*J]   += ((B₁[i]*F₁₁*B₁[j]*F₃₁ + B₂[i]*F₁₂*B₂[j]*F₃₂ + B₃[i]*F₁₃*B₃[j]*F₃₃)*Cᵢᵢᵢᵢ
                               +   (B₂[i]*F₁₂*B₁[j]*F₃₁ + B₃[i]*F₁₃*B₁[j]*F₃₁ + B₁[i]*F₁₁*B₂[j]*F₃₂ + B₃[i]*F₁₃*B₂[j]*F₃₂ + B₁[i]*F₁₁*B₃[j]*F₃₃ + B₂[i]*F₁₂*B₃[j]*F₃₃)*Cᵢᵢⱼⱼ
                               +   ((B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₃₂+B₂[j]*F₃₁) + (B₃[i]*F₁₁+B₁[i]*F₁₃)*(B₃[j]*F₃₁+B₁[j]*F₃₃) + (B₃[i]*F₁₂+B₂[i]*F₁₃)*(B₃[j]*F₃₂+B₂[j]*F₃₃))*Cᵢⱼᵢⱼ)*𝑤

                k[3*I-1,3*J-2] +=((B₁[i]*F₂₁*B₁[j]*F₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂ + B₃[i]*F₂₃*B₃[j]*F₁₃)*Cᵢᵢᵢᵢ
                               +  (B₂[i]*F₂₂*B₁[j]*F₁₁ + B₃[i]*F₂₃*B₁[j]*F₁₁ + B₁[i]*F₂₁*B₂[j]*F₁₂ + B₃[i]*F₂₃*B₂[j]*F₁₂ + B₁[i]*F₂₁*B₃[j]*F₁₃ + B₂[i]*F₂₂*B₃[j]*F₁₃)*Cᵢᵢⱼⱼ
                               +  ((B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)+(B₃[j]*F₁₁+B₁[j]*F₁₃)*(B₃[i]*F₂₁+B₁[i]*F₂₃) + (B₃[j]*F₁₂+B₂[j]*F₁₃)*(B₃[i]*F₂₂+B₂[i]*F₂₃) )*Cᵢⱼᵢⱼ)*𝑤

                k[3*I-1,3*J-1] += ((B₁[i]*F₂₁*B₁[j]*F₂₁ + B₂[i]*F₂₂*B₂[j]*F₂₂ + B₃[i]*F₂₃*B₃[j]*F₂₃)*Cᵢᵢᵢᵢ
                               +  (B₂[i]*F₂₂*B₁[j]*F₂₁ + B₃[i]*F₂₃*B₁[j]*F₂₁ + B₁[i]*F₂₁*B₂[j]*F₂₂ + B₃[i]*F₂₃*B₂[j]*F₂₂ + B₁[i]*F₂₁*B₃[j]*F₂₃ + B₂[i]*F₂₂*B₃[j]*F₂₃)*Cᵢᵢⱼⱼ
                               +  ((B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁) + (B₃[i]*F₂₁+B₁[i]*F₂₃)*(B₃[j]*F₂₁+B₁[j]*F₂₃) + (B₃[i]*F₂₂+B₂[i]*F₂₃)*(B₃[j]*F₂₂+B₂[j]*F₂₃))*Cᵢⱼᵢⱼ
                               +  B₁[i]*B₁[j]*S₁₁+B₂[i]*B₂[j]*S₂₂+B₃[i]*B₃[j]*S₃₃ + (B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂ + (B₁[i]*B₃[j]+B₃[i]*B₁[j])*S₁₃ + (B₃[i]*B₂[j]+B₂[i]*B₃[j])*S₂₃)*𝑤

                k[3*I-1,3*J]   += ((B₁[i]*F₂₁*B₁[j]*F₃₁ + B₂[i]*F₂₂*B₂[j]*F₃₂ + B₃[i]*F₂₃*B₃[j]*F₃₃)*Cᵢᵢᵢᵢ
                               +   (B₂[i]*F₂₂*B₁[j]*F₃₁ + B₃[i]*F₂₃*B₁[j]*F₃₁ + B₁[i]*F₂₁*B₂[j]*F₃₂ + B₃[i]*F₂₃*B₂[j]*F₃₂ + B₁[i]*F₂₁*B₃[j]*F₃₃ + B₂[i]*F₂₂*B₃[j]*F₃₃)*Cᵢᵢⱼⱼ
                               +   ((B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₃₂+B₂[j]*F₃₁) + (B₃[i]*F₂₁+B₁[i]*F₂₃)*(B₃[j]*F₃₁+B₁[j]*F₃₃) + (B₃[i]*F₂₂+B₂[i]*F₂₃)*(B₃[j]*F₃₂+B₂[j]*F₃₃))*Cᵢⱼᵢⱼ)*𝑤

                k[3*I,3*J-2]   += ((B₁[i]*F₃₁*B₁[j]*F₁₁ + B₂[i]*F₃₂*B₂[j]*F₁₂ + B₃[i]*F₃₃*B₃[j]*F₁₃)*Cᵢᵢᵢᵢ
                               +   (B₂[i]*F₃₂*B₁[j]*F₁₁ + B₃[i]*F₃₃*B₁[j]*F₁₁ + B₁[i]*F₃₁*B₂[j]*F₁₂ + B₃[i]*F₃₃*B₂[j]*F₁₂ + B₁[i]*F₃₁*B₃[j]*F₁₃ + B₂[i]*F₃₂*B₃[j]*F₁₃)*Cᵢᵢⱼⱼ
                               +   ((B₁[i]*F₃₂+B₂[i]*F₃₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁) + (B₃[i]*F₃₁+B₁[i]*F₃₃)*(B₃[j]*F₁₁+B₁[j]*F₁₃) + (B₃[i]*F₃₂+B₂[i]*F₃₃)*(B₃[j]*F₁₂+B₂[j]*F₁₃))*Cᵢⱼᵢⱼ)*𝑤

                k[3*I,3*J-1]   += ((B₁[i]*F₃₁*B₁[j]*F₂₁ + B₂[i]*F₃₂*B₂[j]*F₂₂ + B₃[i]*F₃₃*B₃[j]*F₂₃)*Cᵢᵢᵢᵢ
                               +   (B₂[i]*F₃₂*B₁[j]*F₂₁ + B₃[i]*F₃₃*B₁[j]*F₂₁ + B₁[i]*F₃₁*B₂[j]*F₂₂ + B₃[i]*F₃₃*B₂[j]*F₂₂ + B₁[i]*F₃₁*B₃[j]*F₂₃ + B₂[i]*F₃₂*B₃[j]*F₂₃)*Cᵢᵢⱼⱼ
                               +   ((B₁[i]*F₃₂+B₂[i]*F₃₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁) + (B₃[i]*F₃₁+B₁[i]*F₃₃)*(B₃[j]*F₂₁+B₁[j]*F₂₃) + (B₃[i]*F₃₂+B₂[i]*F₃₃)*(B₃[j]*F₂₂+B₂[j]*F₂₃))*Cᵢⱼᵢⱼ)*𝑤

                k[3*I,3*J]   += ((B₁[i]*F₃₁*B₁[j]*F₃₁ + B₂[i]*F₃₂*B₂[j]*F₃₂ + B₃[i]*F₃₃*B₃[j]*F₃₃)*Cᵢᵢᵢᵢ
                             +   (B₂[i]*F₃₂*B₁[j]*F₃₁ + B₃[i]*F₃₃*B₁[j]*F₃₁ + B₁[i]*F₃₁*B₂[j]*F₃₂ + B₃[i]*F₃₃*B₂[j]*F₃₂ + B₁[i]*F₃₁*B₃[j]*F₃₃ + B₂[i]*F₃₂*B₃[j]*F₃₃)*Cᵢᵢⱼⱼ
                             +   ((B₁[i]*F₃₂+B₂[i]*F₃₁)*(B₁[j]*F₃₂+B₂[j]*F₃₁) + (B₃[i]*F₃₁+B₁[i]*F₃₃)*(B₃[j]*F₃₁+B₁[j]*F₃₃) + (B₃[i]*F₃₂+B₂[i]*F₃₃)*(B₃[j]*F₃₂+B₂[j]*F₃₃))*Cᵢⱼᵢⱼ
                             +   B₁[i]*B₁[j]*S₁₁+B₂[i]*B₂[j]*S₂₂+B₃[i]*B₃[j]*S₃₃ + (B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂ + (B₁[i]*B₃[j]+B₃[i]*B₁[j])*S₁₃ + (B₃[i]*B₂[j]+B₂[i]*B₃[j])*S₂₃)*𝑤

            end
        end
    end
end

function ∫∫EᵢⱼSᵢⱼdxdy_SaintVenantKirchhoff(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    for ξ in 𝓖
        E = ξ.Ē
        ν = ξ.ν̄ 
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

function ∫∫∫EᵢⱼSᵢⱼdxdydz_SaintVenantKirchhoff(ap::T, f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    for ξ in 𝓖
        E = ξ.E
        ν = ξ.ν
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z]  
        𝑤 = ξ.𝑤

        λ = E*ν/((1+ν)*(1-2*ν)) 
        μ = E/(2*(1+ν))           
        Cᵢᵢᵢᵢ = λ + 2*μ           

        F₁₁ = 1.0; F₁₂ = 0.0; F₁₃ = 0.0
        F₂₁ = 0.0; F₂₂ = 1.0; F₂₃ = 0.0
        F₃₁ = 0.0; F₃₂ = 0.0; F₃₃ = 1.0

        for (i,xᵢ) in enumerate(𝓒)
            F₁₁ += B₁[i]*xᵢ.d₁; F₁₂ += B₂[i]*xᵢ.d₁; F₁₃ += B₃[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂; F₂₂ += B₂[i]*xᵢ.d₂; F₂₃ += B₃[i]*xᵢ.d₂
            F₃₁ += B₁[i]*xᵢ.d₃; F₃₂ += B₂[i]*xᵢ.d₃; F₃₃ += B₃[i]*xᵢ.d₃
        end

        E₁₁ = 0.5*(F₁₁*F₁₁+F₂₁*F₂₁+F₃₁*F₃₁-1.0)
        E₂₂ = 0.5*(F₁₂*F₁₂+F₂₂*F₂₂+F₃₂*F₃₂-1.0)
        E₃₃ = 0.5*(F₁₃*F₁₃+F₂₃*F₂₃+F₃₃*F₃₃-1.0)
        E₁₂ = 0.5*(F₁₁*F₁₂+F₂₁*F₂₂+F₃₁*F₃₂)
        E₁₃ = 0.5*(F₁₁*F₁₃+F₂₁*F₂₃+F₃₁*F₃₃)
        E₂₃ = 0.5*(F₁₂*F₁₃+F₂₂*F₂₃+F₃₂*F₃₃)

        S₁₁ = Cᵢᵢᵢᵢ*E₁₁ + λ*E₂₂ + λ*E₃₃
        S₂₂ = λ*E₁₁ + Cᵢᵢᵢᵢ*E₂₂ + λ*E₃₃
        S₃₃ = λ*E₁₁ + λ*E₂₂ + Cᵢᵢᵢᵢ*E₃₃
        S₁₂ = 2.0*μ*E₁₂ 
        S₁₃ = 2.0*μ*E₁₃
        S₂₃ = 2.0*μ*E₂₃

        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼

            f[3*I-2] += (
                B₁[i] * (F₁₁*S₁₁ + F₁₂*S₁₂ + F₁₃*S₁₃) +
                B₂[i] * (F₁₁*S₁₂ + F₁₂*S₂₂ + F₁₃*S₂₃) +
               B₃[i] * (F₁₁*S₁₃ + F₁₂*S₂₃ + F₁₃*S₃₃)
            ) * 𝑤

            f[3*I-1] += (
                B₁[i] * (F₂₁*S₁₁ + F₂₂*S₁₂ + F₂₃*S₁₃) +
                B₂[i] * (F₂₁*S₁₂ + F₂₂*S₂₂ + F₂₃*S₂₃) +
                B₃[i] * (F₂₁*S₁₃ + F₂₂*S₂₃ + F₂₃*S₃₃)
            ) * 𝑤

            f[3*I] += (
               B₁[i] * (F₃₁*S₁₁ + F₃₂*S₁₂ + F₃₃*S₁₃) +
               B₂[i] * (F₃₁*S₁₂ + F₃₂*S₂₂ + F₃₃*S₂₃) +
               B₃[i] * (F₃₁*S₁₃ + F₃₂*S₂₃ + F₃₃*S₃₃)
            ) * 𝑤
        end
    end
end

function Δ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    for ξ in 𝓖
         E=ξ.Ē
         ν=ξ.ν̄ 
         λ=E*ν/(1+ν)/(1-2*ν)
         μ=E/(1+ν)/2
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

function Δ∫∫∫EᵢⱼSᵢⱼdxdydz_NeoHookean(ap::T, k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    F = zeros(3,3)
    C = zeros(3,3)
    Cinv = zeros(3,3)
    S = zeros(3,3)

    for ξ in 𝓖
        E = ξ.E
        ν = ξ.ν

        λ = E * ν / ((1 + ν) * (1 - 2 * ν))
        μ = E / ((1 + ν) * 2)

        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z] 
        𝑤 = ξ.𝑤

        fill!(F, 0.0)
        F[1,1] = 1.0; F[2,2] = 1.0; F[3,3] = 1.0

        for (i, xᵢ) in enumerate(𝓒)

            d₁ = xᵢ.d₁
            d₂ = xᵢ.d₂
            d₃ = xᵢ.d₃ 

            F[1,1] += B₁[i] * d₁;  F[1,2] += B₂[i] * d₁;  F[1,3] += B₃[i] * d₁
            F[2,1] += B₁[i] * d₂;  F[2,2] += B₂[i] * d₂;  F[2,3] += B₃[i] * d₂
            F[3,1] += B₁[i] * d₃;  F[3,2] += B₂[i] * d₃;  F[3,3] += B₃[i] * d₃
        end

        mul!(C, F', F) 

        detC = det(C)
        J = sqrt(detC) 

        Cinv .= inv(C)

        coeff_pres = λ * J * (J - 1.0)

        for i=1:3, j=1:3
            δ = (i == j) ? 1.0 : 0.0
            S[i,j] = coeff_pres * Cinv[i,j] + μ * (δ - Cinv[i,j])
        end

        c_vol = λ * J * (2.0 * J - 1.0)
        c_iso = 2.0 * (μ - λ * J * (J - 1.0)) 

        for (I, xI) in enumerate(𝓒)
            row_idx = xI.𝐼 

            bI = [B₁[I], B₂[I], B₃[I]] 

            for (J, xJ) in enumerate(𝓒)
                col_idx = xJ.𝐼
                bJ = [B₁[J], B₂[J], B₃[J]]

                s_geo = 0.0
                for p=1:3, q=1:3
                    s_geo += bI[p] * S[p,q] * bJ[q]
                end

                for m = 1:3
                    row = 3 * row_idx - 3 + m
                    for n = 1:3
                        col = 3 * col_idx - 3 + n

                        val = (m == n) ? s_geo : 0.0

                        term_mat = 0.0
                        for k=1:3, l=1:3, p=1:3, q=1:3

                            Isym = 0.5 * (Cinv[k,p]*Cinv[l,q] + Cinv[k,q]*Cinv[l,p])
                            C_ijkl = c_vol * Cinv[k,l] * Cinv[p,q] + c_iso * Isym 

                            term_mat += bI[k] * F[m,l] * C_ijkl * F[n,p] * bJ[q]
                        end

                        val += term_mat

                        k[row, col] += val * 𝑤
                    end
                end
            end
        end
    end
end

function Δ∫∫∫_NeoHookean_Dev(ap::T, k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    F = zeros(3,3); C = zeros(3,3); Cinv = zeros(3,3); S_dev = zeros(3,3)

    for ξ in 𝓖
        E, ν, 𝑤 = ξ.E, ξ.ν, ξ.𝑤
        μ = E / (2 * (1 + ν))
        B = [ξ[:∂𝝭∂x] ξ[:∂𝝭∂y] ξ[:∂𝝭∂z]] 

        fill!(F, 0.0); for i=1:3; F[i,i]=1.0; end
        for (i, xᵢ) in enumerate(𝓒)
            d = [xᵢ.d₁, xᵢ.d₂, xᵢ.d₃]
            F += d * B[i, :]' 
        end
        mul!(C, F', F)
        J = sqrt(det(C))
        Cinv .= inv(C)

        trC = tr(C)
        J23inv = J^(-2/3)
        for i=1:3, j=1:3
            δ = (i == j) ? 1.0 : 0.0
            S_dev[i,j] = μ * J23inv * (δ - (1/3) * trC * Cinv[i,j])
        end

        c_iso = 2.0 * μ * J23inv 

        assemble_stiffness!(k, 𝓒, B, F, Cinv, S_dev, c_iso, 0.0, 𝑤) 
    end
end

function Δ∫∫∫_NeoHookean_Vol(ap::T, k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖_low = ap.𝓖_reduced 
    F = zeros(3,3); C = zeros(3,3); Cinv = zeros(3,3); S_vol = zeros(3,3)

    for ξ in 𝓖_low
        E, ν, 𝑤 = ξ.E, ξ.ν, ξ.𝑤
        λ = E * ν / ((1 + ν) * (1 - 2 * ν))
        B = [ξ[:∂𝝭∂x] ξ[:∂𝝭∂y] ξ[:∂𝝭∂z]]

        fill!(F, 0.0); for i=1:3; F[i,i]=1.0; end
        for (i, xᵢ) in enumerate(𝓒)
            d = [xᵢ.d₁, xᵢ.d₂, xᵢ.d₃]
            F += d * B[i, :]'
        end
        detF = det(F)
        Cinv .= inv(F' * F)

        p = λ * (detF - 1.0) 
        S_vol .= (p * detF) .* Cinv

        c_vol = λ * detF * (2.0 * detF - 1.0)

        assemble_stiffness!(k, 𝓒, B, F, Cinv, S_vol, 0.0, c_vol, 𝑤)
    end
end

function Δ∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E=op.E
    ν=op.ν
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

function Δ∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E=op.E
    ν=op.ν
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

function ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    for ξ in 𝓖
        E=ξ.Ē
        ν=ξ.ν̄ 
        λ=E*ν/(1+ν)/(1-2*ν)
        μ=E/(1+ν)/2
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

function ∫∫∫EᵢⱼSᵢⱼdxdydz_NeoHookean(ap::T, f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    F = zeros(3,3)
    C = zeros(3,3)
    S = zeros(3,3)
    P = zeros(3,3) 

    for ξ in 𝓖
        E = ξ.E
        ν = ξ.ν

        λ = E * ν / ((1 + ν) * (1 - 2 * ν))
        μ = E / ((1 + ν) * 2)

        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z] 
        𝑤 = ξ.𝑤

        fill!(F, 0.0)
        F[1,1] = 1.0; F[2,2] = 1.0; F[3,3] = 1.0

        for (i, xᵢ) in enumerate(𝓒)
            d₁ = xᵢ.d₁
            d₂ = xᵢ.d₂
            d₃ = xᵢ.d₃ 

            F[1,1] += d₁ * B₁[i]; F[1,2] += d₁ * B₂[i]; F[1,3] += d₁ * B₃[i]
            F[2,1] += d₂ * B₁[i]; F[2,2] += d₂ * B₂[i]; F[2,3] += d₂ * B₃[i]
            F[3,1] += d₃ * B₁[i]; F[3,2] += d₃ * B₂[i]; F[3,3] += d₃ * B₃[i]
        end

        mul!(C, F', F) 

        J = det(F)       
        detC = det(C)    
        Cinv = inv(C)    

        coeff_1 = λ * J * (J - 1.0)

        for i in 1:3, j in 1:3
            δ = (i == j) ? 1.0 : 0.0
            S[i,j] = coeff_1 * Cinv[i,j] + μ * (δ - Cinv[i,j])
        end

        mul!(P, F, S)

        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 

            dN_dX = B₁[i]
            dN_dY = B₂[i]
            dN_dZ = B₃[i]

            val_x = P[1,1]*dN_dX + P[1,2]*dN_dY + P[1,3]*dN_dZ

            val_y = P[2,1]*dN_dX + P[2,2]*dN_dY + P[2,3]*dN_dZ

            val_z = P[3,1]*dN_dX + P[3,2]*dN_dY + P[3,3]*dN_dZ

            f[3*I-2] += val_x * 𝑤
            f[3*I-1] += val_y * 𝑤
            f[3*I]   += val_z * 𝑤
        end
    end
end

function ∫∫∫_NeoHookean_Force_Dev(ap::T, f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    F = zeros(3,3); C = zeros(3,3); S_dev = zeros(3,3); P_dev = zeros(3,3)

    for ξ in 𝓖
        E, ν, 𝑤 = ξ.E, ξ.ν, ξ.𝑤
        μ = E / (2 * (1 + ν))
        B₁ = ξ[:∂𝝭∂x]; B₂ = ξ[:∂𝝭∂y]; B₃ = ξ[:∂𝝭∂z]

        fill!(F, 0.0); for i=1:3; F[i,i]=1.0; end
        for (i, xᵢ) in enumerate(𝓒)
            F[1,1]+=xᵢ.d₁*B₁[i]; F[1,2]+=xᵢ.d₁*B₂[i]; F[1,3]+=xᵢ.d₁*B₃[i]
            F[2,1]+=xᵢ.d₂*B₁[i]; F[2,2]+=xᵢ.d₂*B₂[i]; F[2,3]+=xᵢ.d₂*B₃[i]
            F[3,1]+=xᵢ.d₃*B₁[i]; F[3,2]+=xᵢ.d₃*B₂[i]; F[3,3]+=xᵢ.d₃*B₃[i]
        end

        mul!(C, F', F)
        J = det(F)
        J23inv = J^(-2/3)
        Cinv = inv(C)
        trC = tr(C)

        for i=1:3, j=1:3
            δ = (i == j) ? 1.0 : 0.0
            S_dev[i,j] = μ * J23inv * (δ - (1/3) * trC * Cinv[i,j])
        end

        mul!(P_dev, F, S_dev)
        assemble_force!(f, 𝓒, P_dev, B₁, B₂, B₃, 𝑤)
    end
end

function ∫∫∫_NeoHookean_Force_Vol(ap::T, f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖_low = ap.𝓖_reduced 
    F = zeros(3,3); P_vol = zeros(3,3)

    for ξ in 𝓖_low
        E, ν, 𝑤 = ξ.E, ξ.ν, ξ.𝑤
        λ = E * ν / ((1 + ν) * (1 - 2 * ν))
        B₁ = ξ[:∂𝝭∂x]; B₂ = ξ[:∂𝝭∂y]; B₃ = ξ[:∂𝝭∂z]

        fill!(F, 0.0); for i=1:3; F[i,i]=1.0; end
        for (i, xᵢ) in enumerate(𝓒)
            F[1,1]+=xᵢ.d₁*B₁[i]; F[1,2]+=xᵢ.d₁*B₂[i]; F[1,3]+=xᵢ.d₁*B₃[i]
            F[2,1]+=xᵢ.d₂*B₁[i]; F[2,2]+=xᵢ.d₂*B₂[i]; F[2,3]+=xᵢ.d₂*B₃[i]
            F[3,1]+=xᵢ.d₃*B₁[i]; F[3,2]+=xᵢ.d₃*B₂[i]; F[3,3]+=xᵢ.d₃*B₃[i]
        end
        J = det(F)

        p = λ * (J - 1.0)
        FinvT = inv(F)' 
        P_vol .= (p * J) .* FinvT

        assemble_force!(f, 𝓒, P_vol, B₁, B₂, B₃, 𝑤)
    end
end

function assemble_force!(f, 𝓒, P, B₁, B₂, B₃, 𝑤)
    for (i, xᵢ) in enumerate(𝓒)
        I = xᵢ.𝐼

        val_x = P[1,1]*B₁[i] + P[1,2]*B₂[i] + P[1,3]*B₃[i]
        val_y = P[2,1]*B₁[i] + P[2,2]*B₂[i] + P[2,3]*B₃[i]
        val_z = P[3,1]*B₁[i] + P[3,2]*B₂[i] + P[3,3]*B₃[i]

        f[3*I-2] += val_x * 𝑤
        f[3*I-1] += val_y * 𝑤
        f[3*I]   += val_z * 𝑤
    end
end

function ∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E=op.E
    ν=op.ν
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

function ∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E=op.E
    ν=op.ν
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

function ∫vᵢuᵢds(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    for ξ in 𝓖
        α = ξ.α
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

function ∫vᵢuᵢdΓ(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    for ξ in 𝓖
        α   = ξ.α
        𝑤   = ξ.𝑤
        N   = ξ[:𝝭]

        n₁₁ = ξ.n₁₁
        n₂₂ = ξ.n₂₂
        n₃₃ = ξ.n₃₃
        n₁₂ = ξ.n₁₂
        n₁₃ = ξ.n₁₃
        n₂₃ = ξ.n₂₃

        g₁ = ξ.g₁
        g₂ = ξ.g₂
        g₃ = ξ.g₃

        u₁ = 0.0
        u₂ = 0.0
        u₃ = 0.0
        for (i, xᵢ) in enumerate(𝓒)
            u₁ += N[i] * xᵢ.d₁
            u₂ += N[i] * xᵢ.d₂
            u₃ += N[i] * xᵢ.d₃
        end

        Δu₁ = g₁ - u₁
        Δu₂ = g₂ - u₂
        Δu₃ = g₃ - u₃

        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼

            for (j, xⱼ) in enumerate(𝓒)
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

            f[3*I-2] += α * N[i] * (n₁₁*Δu₁ + n₁₂*Δu₂ + n₁₃*Δu₃) * 𝑤
            f[3*I-1] += α * N[i] * (n₁₂*Δu₁ + n₂₂*Δu₂ + n₂₃*Δu₃) * 𝑤
            f[3*I]   += α * N[i] * (n₁₃*Δu₁ + n₂₃*Δu₂ + n₃₃*Δu₃) * 𝑤
        end
    end
end

function ∫vᵢuᵢdΓ_center(ap::T, kα::AbstractMatrix{Float64}, fα::AbstractVector{Float64};
                          xc::Float64=0.5, yc::Float64=0.5, σ::Float64=0.02) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖

    for ξ in 𝓖
        N = ξ[:𝝭]
        wΓ = ξ.𝑤
        α   = ξ.α

        dx = ξ.x - xc
        dy = ξ.y - yc
        wloc = exp(-(dx*dx + dy*dy) / (σ*σ))

        ux = 0.0
        uy = 0.0
        for (i, xᵢ) in enumerate(𝓒)
            ux += N[i] * xᵢ.d₁
            uy += N[i] * xᵢ.d₂
        end

        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            Ni = N[i]
            coef = α * wloc * Ni * wΓ

            fα[3*I-2] += -coef * ux
            fα[3*I-1] += -coef * uy

            for (j, xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                Nj = N[j]
                kij = α * wloc * Ni * Nj * wΓ

                kα[3*I-2, 3*J-2] += kij

                kα[3*I-1, 3*J-1] += kij

            end
        end
    end
end

function Δ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean2(ap::T, k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        E=ξ.Ē
        ν=ξ.ν̄ 
        K=E/(1-2*ν)/3
        λ=E*ν/(1+ν)/(1-2*ν)
        μ=E/(1+ν)/2
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

        C₁₁₁₁=(-2/3*C⁻¹₁₁-2/3*C⁻¹₁₁)*J⁻²³*μ+(J^2*K+2/9*J⁻²³*I₁*μ)*C⁻¹₁₁*C⁻¹₁₁+(1/3*μ*J⁻²³*I₁-0.5*K*(J^2-1.0))*2*C⁻¹₁₁*C⁻¹₁₁

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

function Δ∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean2(ap::T, k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        E=ξ.Ē
        ν=ξ.ν̄ 
        K=E/(1-2*ν)/3
        λ=E*ν/(1+ν)/(1-2*ν)
        μ=E/(1+ν)/2
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
function Δ∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean2(ap::T, k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        E=ξ.Ē
        ν=ξ.ν̄ 
        K=E/(1-2*ν)/3
        λ=E*ν/(1+ν)/(1-2*ν)
        μ=E/(1+ν)/2
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

function ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean2(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        E=ξ.Ē
        ν=ξ.ν̄ 
        K=E/(1-2*ν)/3
        λ=E*ν/(1+ν)/(1-2*ν)
        μ=E/(1+ν)/2
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

function ∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean2(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        E=ξ.Ē
        ν=ξ.ν̄ 
        K=E/(1-2*ν)/3
        λ=E*ν/(1+ν)/(1-2*ν)
        μ=E/(1+ν)/2
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
function ∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean2(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        E=ξ.Ē
        ν=ξ.ν̄ 
        K=E/(1-2*ν)/3
        λ=E*ν/(1+ν)/(1-2*ν)
        μ=E/(1+ν)/2
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

function Δ∫∫∫EᵢⱼSᵢⱼdxdydz_NeoHookean2(ap::T, k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    F = zeros(3, 3)
    C = zeros(3, 3)
    S = zeros(3, 3)

    for ξ in 𝓖
        E = ξ.E
        ν = ξ.ν

        K_mod = E / (1 - 2*ν) / 3
        μ = E / (1 + ν) / 2

        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z] 
        𝑤 = ξ.𝑤

        fill!(F, 0.0)
        F[1,1] = 1.0; F[2,2] = 1.0; F[3,3] = 1.0

        for (i, xᵢ) in enumerate(𝓒)
            d₁ = xᵢ.d₁
            d₂ = xᵢ.d₂
            d₃ = xᵢ.d₃ 

            F[1,1] += d₁ * B₁[i]; F[1,2] += d₁ * B₂[i]; F[1,3] += d₁ * B₃[i]
            F[2,1] += d₂ * B₁[i]; F[2,2] += d₂ * B₂[i]; F[2,3] += d₂ * B₃[i]
            F[3,1] += d₃ * B₁[i]; F[3,2] += d₃ * B₂[i]; F[3,3] += d₃ * B₃[i]
        end

        mul!(C, F', F) 
        J = det(F)
        I₁ = tr(C)
        detC = J^2
        Cinv = inv(C) 
        J⁻²³ = cbrt(1.0 / detC)

        for i in 1:3, j in 1:3
            δ_ij = (i == j ? 1.0 : 0.0)
            S[i, j] = μ * J⁻²³ * (δ_ij - 1/3 * I₁ * Cinv[i, j]) + 0.5 * K_mod * (J^2 - 1.0) * Cinv[i, j]
        end

        t1 = -2/3 * μ * J⁻²³
        t2 = K_mod * J^2 + 2/9 * μ * J⁻²³ * I₁
        t3 = 1/3 * μ * J⁻²³ * I₁ - 0.5 * K_mod * (J^2 - 1.0)

        get_Cijkl = (i, j, m, n) -> begin
            δ_ij = (i == j ? 1.0 : 0.0)
            δ_mn = (m == n ? 1.0 : 0.0)
            return t1 * (Cinv[i,j]*δ_mn + δ_ij*Cinv[m,n]) + 
                   t2 * Cinv[i,j]*Cinv[m,n] + 
                   t3 * (Cinv[i,m]*Cinv[j,n] + Cinv[i,n]*Cinv[j,m])
        end

        for (i, xᵢ) in enumerate(𝓒)
            I_node = xᵢ.𝐼
            ∇Nᵢ = (B₁[i], B₂[i], B₃[i]) 

            for (j, xⱼ) in enumerate(𝓒)
                J_node = xⱼ.𝐼
                ∇Nⱼ = (B₁[j], B₂[j], B₃[j])

                k_geo = 0.0
                for m in 1:3, n in 1:3
                    k_geo += ∇Nᵢ[m] * S[m, n] * ∇Nⱼ[n]
                end

                for a in 1:3, b in 1:3
                    k_mat = 0.0
                    for k_idx in 1:3, l_idx in 1:3, p_idx in 1:3, q_idx in 1:3
                        C_val = get_Cijkl(l_idx, k_idx, p_idx, q_idx)
                        k_mat += ∇Nᵢ[k_idx] * F[a, l_idx] * C_val * F[b, p_idx] * ∇Nⱼ[q_idx]
                    end

                    total_k = k_mat + (a == b ? k_geo : 0.0)

                    k[(3*I_node-3) + a, (3*J_node-3) + b] += total_k * 𝑤
                end
            end
        end
    end
end

function ∫∫∫EᵢⱼSᵢⱼdxdydz_NeoHookean2(ap::T, f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    F = zeros(3, 3)
    C = zeros(3, 3)
    S = zeros(3, 3)
    P = zeros(3, 3) 

    for ξ in 𝓖
        E = ξ.E
        ν = ξ.ν

        K_mod = E / (1 - 2*ν) / 3
        μ = E / (1 + ν) / 2

        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z] 
        𝑤 = ξ.𝑤

        fill!(F, 0.0)
        F[1,1] = 1.0; F[2,2] = 1.0; F[3,3] = 1.0

        for (i, xᵢ) in enumerate(𝓒)
            d₁ = xᵢ.d₁
            d₂ = xᵢ.d₂
            d₃ = xᵢ.d₃ 

            F[1,1] += d₁ * B₁[i]; F[1,2] += d₁ * B₂[i]; F[1,3] += d₁ * B₃[i]
            F[2,1] += d₂ * B₁[i]; F[2,2] += d₂ * B₂[i]; F[2,3] += d₂ * B₃[i]
            F[3,1] += d₃ * B₁[i]; F[3,2] += d₃ * B₂[i]; F[3,3] += d₃ * B₃[i]
        end

        mul!(C, F', F) 

        J = det(F)       
        detC = det(C)    
        Cinv = inv(C)    

        I₁ = tr(C)
        J⁻²³ = cbrt(1.0 / detC)

        for i in 1:3, j in 1:3
            δ_ij = (i == j) ? 1.0 : 0.0
            S[i,j] = μ * J⁻²³ * (δ_ij - 1/3 * I₁ * Cinv[i,j]) + 0.5 * K_mod * (J^2 - 1.0) * Cinv[i,j]
        end

        mul!(P, F, S)

        for (i, xᵢ) in enumerate(𝓒)
            I_node = xᵢ.𝐼

            dN_dX = B₁[i]
            dN_dY = B₂[i]
            dN_dZ = B₃[i]

            val_x = P[1,1]*dN_dX + P[1,2]*dN_dY + P[1,3]*dN_dZ
            val_y = P[2,1]*dN_dX + P[2,2]*dN_dY + P[2,3]*dN_dZ
            val_z = P[3,1]*dN_dX + P[3,2]*dN_dY + P[3,3]*dN_dZ

            f[3*I_node-2] += val_x * 𝑤
            f[3*I_node-1] += val_y * 𝑤
            f[3*I_node]   += val_z * 𝑤
        end
    end
end

function get_neoHookean_S_and_H_components(E_vec::Vector{Float64},
                                           E::Float64, ν::Float64)
    E11, E22, E12 = E_vec

    λ = E*ν/((1+ν)*(1-2ν))
    μ = E/(2*(1+ν))

    C11 = 2*E11 + 1.0
    C22 = 2*E22 + 1.0
    C12 = 2*E12
    C33 = 1.0

    detC = C11*C22 - C12*C12
    if detC <= 0
        detC = 1e-12       
    end
    J = sqrt(detC)

    I1 = C11 + C22 + C33
    I2 = 0.5*(I1^2 - C11^2 - C22^2 - C33^2 - 2*C12^2)

    Cinv11 = 1.0/detC*(C11*C11 + C12*C12 - I1*C11 + I2)
    Cinv12 = 1.0/detC*(C11*C12 + C12*C22 - I1*C12)
    Cinv22 = 1.0/detC*(C12*C12 + C22*C22 - I1*C22 + I2)

    S11 = λ*J*(J-1.0)*Cinv11 + μ*(1.0 - Cinv11)
    S22 = λ*J*(J-1.0)*Cinv22 + μ*(1.0 - Cinv22)
    S12 = λ*J*(J-1.0)*Cinv12 - μ*Cinv12

    C1111 = λ*J*(2*J-1.0)*Cinv11*Cinv11 +
            (μ-λ*J*(J-1.0))*2*Cinv11*Cinv11

    C2222 = λ*J*(2*J-1.0)*Cinv22*Cinv22 +
            (μ-λ*J*(J-1.0))*2*Cinv22*Cinv22

    C1122 = λ*J*(2*J-1.0)*Cinv11*Cinv22 +
            (μ-λ*J*(J-1.0))*2*Cinv12*Cinv12

    C1112 = λ*J*(2*J-1.0)*Cinv11*Cinv12 +
            (μ-λ*J*(J-1.0))*2*Cinv11*Cinv12

    C2212 = λ*J*(2*J-1.0)*Cinv22*Cinv12 +
            (μ-λ*J*(J-1.0))*2*Cinv12*Cinv22

    C1212 = λ*J*(2*J-1.0)*Cinv12*Cinv12 +
            (μ-λ*J*(J-1.0))*(Cinv11*Cinv22 + Cinv12*Cinv12)

    S_vec = [S11, S22, S12]
    H_comp = (C1111, C2222, C1122, C1112, C2212, C1212)

    return S_vec, H_comp
end

function get_neoHookean_S_and_H(E_vec::Vector{Float64}, E::Float64, ν::Float64)

    E11, E22, E12 = E_vec

    λ = E * ν / ((1 + ν) * (1 - 2*ν))
    μ = E / (2 * (1 + ν))

    C11 = 2.0*E11 + 1.0
    C22 = 2.0*E22 + 1.0
    C12 = 2.0*E12
    C33 = 1.0

    C = @inbounds [
        C11  C12  0.0
        C12  C22  0.0
        0.0  0.0  C33
    ]

    detC = det(C)
    if detC <= 1e-14
        throw(DomainError(detC, "Non-physical state: det(C) <= 0"))
    end
    J = sqrt(detC)

    Cinv = inv(C)
    Cinv11 = Cinv[1,1]
    Cinv22 = Cinv[2,2]
    Cinv12 = Cinv[1,2]  

    vol = λ * J * (J - 1.0)

    S11 = vol*Cinv11 + μ*(1.0 - Cinv11)
    S22 = vol*Cinv22 + μ*(1.0 - Cinv22)
    S12 = vol*Cinv12 - μ*Cinv12

    S_vec = [S11, S22, S12]

    A = λ * J * (2.0*J - 1.0)
    B = μ - λ * J * (J - 1.0)

    C1111 = A*Cinv11*Cinv11 + 2.0*B*Cinv11*Cinv11
    C2222 = A*Cinv22*Cinv22 + 2.0*B*Cinv22*Cinv22
    C1122 = A*Cinv11*Cinv22 + 2.0*B*Cinv12*Cinv12
    C1112 = A*Cinv11*Cinv12 + 2.0*B*Cinv11*Cinv12
    C2212 = A*Cinv22*Cinv12 + 2.0*B*Cinv12*Cinv22
    C1212 = A*Cinv12*Cinv12 + B*(Cinv11*Cinv22 + Cinv12*Cinv12)

    H = zeros(3,3)
    H[1,1] = C1111
    H[1,2] = C1122
    H[1,3] = C1112

    H[2,1] = C1122
    H[2,2] = C2222
    H[2,3] = C2212

    H[3,1] = C1112
    H[3,2] = C2212
    H[3,3] = C1212

    return S_vec, H
end

function get_neoHookean2_S_and_H(E_vec::Vector{Float64}, E::Float64, ν::Float64)

    E11, E22, E12 = E_vec

    λ = E * ν / ((1 + ν) * (1 - 2*ν))
    μ = E / (2 * (1 + ν))

    C11 = 2.0*E11 + 1.0
    C22 = 2.0*E22 + 1.0
    C12 = 2.0*E12/2.0
    C33 = 1.0

    C = @inbounds [
        C11  C12  0.0
        C12  C22  0.0
        0.0  0.0  C33
    ]

    detC = det(C)
    if detC <= 1e-14
        throw(DomainError(detC, "Non-physical state: det(C) <= 0"))
    end
    J = sqrt(detC)

    Cinv = inv(C)
    Cinv11 = Cinv[1,1]
    Cinv22 = Cinv[2,2]
    Cinv12 = Cinv[1,2]  

    vol = λ * J * (J - 1.0)

    S11 = vol*Cinv11 + μ*(1.0 - Cinv11)
    S22 = vol*Cinv22 + μ*(1.0 - Cinv22)
    S12 = vol*Cinv12 - μ*Cinv12

    S_vec = [S11, S22, S12]

    A = λ * J * (2.0*J - 1.0)
    B = μ - λ * J * (J - 1.0)

    C1111 = A*Cinv11*Cinv11 + 2.0*B*Cinv11*Cinv11
    C2222 = A*Cinv22*Cinv22 + 2.0*B*Cinv22*Cinv22
    C1122 = A*Cinv11*Cinv22 + 2.0*B*Cinv12*Cinv12
    C1112 = A*Cinv11*Cinv12 + 2.0*B*Cinv11*Cinv12
    C2212 = A*Cinv22*Cinv12 + 2.0*B*Cinv12*Cinv22
    C1212 = A*Cinv12*Cinv12 + B*(Cinv11*Cinv22 + Cinv12*Cinv12)

    H = zeros(3,3)
    H[1,1] = C1111
    H[1,2] = C1122
    H[1,3] = C1112

    H[2,1] = C1122
    H[2,2] = C2222
    H[2,3] = C2212

    H[3,1] = C1112
    H[3,2] = C2212
    H[3,3] = C1212

    return S_vec, H
end

function get_neoHookean2_3D_S_and_H(E_vec::Vector{Float64}, E::Float64, ν::Float64)
    E11, E22, E33, E12, E23, E13 = E_vec

    λ = E * ν / ((1 + ν) * (1 - 2*ν))
    μ = E / (2 * (1 + ν))

    C11 = 2.0*E11 + 1.0
    C22 = 2.0*E22 + 1.0
    C33 = 2.0*E33 + 1.0
    C12 = E12
    C23 = E23
    C13 = E13

    detC = C11*(C22*C33 - C23^2) - C12*(C12*C33 - C13*C23) + C13*(C12*C23 - C13*C22)

    if detC <= 1e-14
        throw(DomainError(detC, "Non-physical state: det(C) <= 0"))
    end
    J = sqrt(detC)
    invDetC = 1.0 / detC

    c11 = (C22*C33 - C23*C23) * invDetC
    c22 = (C11*C33 - C13*C13) * invDetC
    c33 = (C11*C22 - C12*C12) * invDetC
    c12 = (C13*C23 - C12*C33) * invDetC
    c23 = (C12*C13 - C11*C23) * invDetC
    c13 = (C12*C23 - C13*C22) * invDetC

    vol = λ * J * (J - 1.0)

    S11 = vol*c11 + μ*(1.0 - c11)
    S22 = vol*c22 + μ*(1.0 - c22)
    S33 = vol*c33 + μ*(1.0 - c33)
    S12 = vol*c12 - μ*c12
    S23 = vol*c23 - μ*c23
    S13 = vol*c13 - μ*c13

    S_vec = [S11, S22, S33, S12, S23, S13]

    A = λ * J * (2.0*J - 1.0)
    B = μ - λ * J * (J - 1.0)

    H11 = A*c11*c11 + 2.0*B*c11*c11
    H22 = A*c22*c22 + 2.0*B*c22*c22
    H33 = A*c33*c33 + 2.0*B*c33*c33

    H12 = A*c11*c22 + 2.0*B*c12*c12
    H13 = A*c11*c33 + 2.0*B*c13*c13
    H23 = A*c22*c33 + 2.0*B*c23*c23

    H14 = A*c11*c12 + 2.0*B*c11*c12
    H15 = A*c11*c23 + 2.0*B*c12*c13
    H16 = A*c11*c13 + 2.0*B*c11*c13

    H24 = A*c22*c12 + 2.0*B*c22*c12
    H25 = A*c22*c23 + 2.0*B*c22*c23
    H26 = A*c22*c13 + 2.0*B*c12*c23

    H34 = A*c33*c12 + 2.0*B*c13*c23
    H35 = A*c33*c23 + 2.0*B*c33*c23
    H36 = A*c33*c13 + 2.0*B*c33*c13

    H44 = A*c12*c12 + B*(c11*c22 + c12*c12)
    H45 = A*c12*c23 + B*(c12*c23 + c13*c22)
    H46 = A*c12*c13 + B*(c11*c23 + c13*c12)

    H55 = A*c23*c23 + B*(c22*c33 + c23*c23)
    H56 = A*c23*c13 + B*(c12*c33 + c23*c13)

    H66 = A*c13*c13 + B*(c11*c33 + c13*c13)

    H = zeros(6,6)
    H[1,1]=H11; H[1,2]=H12; H[1,3]=H13; H[1,4]=H14; H[1,5]=H15; H[1,6]=H16
    H[2,1]=H12; H[2,2]=H22; H[2,3]=H23; H[2,4]=H24; H[2,5]=H25; H[2,6]=H26
    H[3,1]=H13; H[3,2]=H23; H[3,3]=H33; H[3,4]=H34; H[3,5]=H35; H[3,6]=H36
    H[4,1]=H14; H[4,2]=H24; H[4,3]=H34; H[4,4]=H44; H[4,5]=H45; H[4,6]=H46
    H[5,1]=H15; H[5,2]=H25; H[5,3]=H35; H[5,4]=H45; H[5,5]=H55; H[5,6]=H56
    H[6,1]=H16; H[6,2]=H26; H[6,3]=H36; H[6,4]=H46; H[6,5]=H56; H[6,6]=H66

    return S_vec, H
end

function get_neoHookean3_S_and_H(E_vec::Vector{Float64}, E::Float64, ν::Float64)
    E11, E22, E12 = E_vec

    λ = E * ν / ((1 + ν) * (1 - 2*ν))
    μ = E / (2 * (1 + ν))

    C11 = 2.0*E11 + 1.0
    C22 = 2.0*E22 + 1.0
    C12 = 2.0*E12
    C33 = 1.0

    C = @inbounds [
        C11  C12  0.0
        C12  C22  0.0
        0.0  0.0  C33
    ]

    detC = det(C)
    if detC <= 1e-14
        throw(DomainError(detC, "Non-physical state: det(C) <= 0"))
    end
    J = sqrt(detC)

    Cinv = inv(C)
    Cinv11 = Cinv[1,1]
    Cinv22 = Cinv[2,2]
    Cinv12 = Cinv[1,2]  

    factor = 0.5 * λ * (J^2 - 1.0)

    S11 = μ*(1.0 - Cinv11) + factor*Cinv11
    S22 = μ*(1.0 - Cinv22) + factor*Cinv22
    S12 =        - μ*Cinv12 + factor*Cinv12

    S_vec = [S11, S22, S12]

    A = λ * J^2
    B = μ - 0.5*λ*(J^2 - 1.0)

    C1111 = A*Cinv11*Cinv11 + 2.0*B*Cinv11*Cinv11
    C2222 = A*Cinv22*Cinv22 + 2.0*B*Cinv22*Cinv22
    C1122 = A*Cinv11*Cinv22 + 2.0*B*Cinv12*Cinv12
    C1112 = A*Cinv11*Cinv12 + 2.0*B*Cinv11*Cinv12
    C2212 = A*Cinv22*Cinv12 + 2.0*B*Cinv12*Cinv22
    C1212 = A*Cinv12*Cinv12 + B*(Cinv11*Cinv22 + Cinv12*Cinv12)

    H = zeros(3,3)
    H[1,1] = C1111
    H[1,2] = C1122
    H[1,3] = C1112*2.0

    H[2,1] = C1122
    H[2,2] = C2222
    H[2,3] = C2212*2.0

    H[3,1] = C1112
    H[3,2] = C2212
    H[3,3] = C1212*2.0

    return S_vec, H
end

function update_Econs!(ξ, S_vec, E_init::AbstractVector{<:Real}; tol_rel=1e-10, tol_abs=1e-12, maxiter=50, alpha_max=1.0)
    E_mod = ξ.Ē
    ν = ξ.ν̄

    λ = E_mod * ν / ((1.0 + ν) * (1.0 - 2.0*ν))
    μ = E_mod / (2.0 * (1.0 + ν))

    C_elastic = zeros(3, 3)
    C_elastic[1, 1] = λ + 2.0*μ
    C_elastic[1, 2] = λ
    C_elastic[2, 1] = λ
    C_elastic[2, 2] = λ + 2.0*μ
    C_elastic[3, 3] = μ

    D = inv(C_elastic)
    S_norm = norm(S_vec)
    tol_check = S_norm < 1.0e-12 ? tol_abs : max(tol_rel*S_norm, tol_abs)

    is_physical(E_vec) = begin
        C11 = 2.0*E_vec[1] + 1.0
        C22 = 2.0*E_vec[2] + 1.0
        C12 = E_vec[3]
        C11*C22 - C12*C12 > 1.0e-12
    end

    function residual_at(E_vec)
        dpsi_dE, H = get_neoHookean2_S_and_H(E_vec, E_mod, ν)
        r = S_vec .- dpsi_dE
        return r, H, norm(r)
    end

    candidates = Vector{Vector{Float64}}()
    if length(E_init) == 3 && all(isfinite, E_init)
        push!(candidates, Float64.(E_init))
    end
    push!(candidates, C_elastic \ S_vec)
    push!(candidates, zeros(3))

    E_vec = zeros(3)
    H_old = copy(C_elastic)
    r_old = copy(S_vec)
    resnorm_old = Inf

    for candidate in candidates
        if !is_physical(candidate)
            continue
        end
        try
            r_candidate, H_candidate, res_candidate = residual_at(candidate)
            if res_candidate < resnorm_old
                E_vec .= candidate
                r_old .= r_candidate
                H_old .= H_candidate
                resnorm_old = res_candidate
            end
        catch
            continue
        end
    end

    if !isfinite(resnorm_old)
        return E_vec, D
    end

    converged = resnorm_old < tol_check
    for it in 1:maxiter
        converged && break

        ΔE = H_old \ r_old
        alpha = alpha_max
        accepted = false

        for ls_it in 1:12
            E_trial = E_vec + alpha*ΔE
            if !is_physical(E_trial)
                alpha *= 0.5
                continue
            end

            try
                r_trial, H_trial, res_trial = residual_at(E_trial)
                if res_trial <= (1.0 - 1.0e-4*alpha)*resnorm_old || res_trial < tol_check
                    E_vec .= E_trial
                    r_old .= r_trial
                    H_old .= H_trial
                    resnorm_old = res_trial
                    accepted = true
                    break
                end
            catch
            end
            alpha *= 0.5
        end

        if !accepted
            break
        end
        converged = resnorm_old < tol_check
    end

    try
        D .= inv(H_old)
    catch
        D .= inv(C_elastic)
    end

    return E_vec, D
end

function update_Econs_3D!(ξ, S_vec, E_init::AbstractVector{<:Real}; tol_rel=1e-4, maxiter=30, alpha_max=1.0)
    E_mod = ξ.E
    ν     = ξ.ν

    λ = E_mod * ν / ((1 + ν) * (1 - 2 * ν))
    μ = E_mod / (2 * (1 + ν))

    C_elastic = zeros(6, 6)
    C_elastic[1,1] = λ + 2μ; C_elastic[1,2] = λ;      C_elastic[1,3] = λ
    C_elastic[2,1] = λ;      C_elastic[2,2] = λ + 2μ; C_elastic[2,3] = λ
    C_elastic[3,1] = λ;      C_elastic[3,2] = λ;      C_elastic[3,3] = λ + 2μ
    C_elastic[4,4] = μ;      C_elastic[5,5] = μ;      C_elastic[6,6] = μ

    E_vec = zeros(6)
    E_vec_old = copy(E_vec)
    D = copy(C_elastic)
    converged = false

    S_norm = norm(S_vec)
    tol_abs = 1e-10 

    local dpsi_dE_old = zeros(6)
    local H_old = zeros(6, 6)
    local r_old = zeros(6)
    local resnorm_old = 0.0

    try
        dpsi_dE_old, H_old = get_neoHookean2_3D_S_and_H(E_vec_old, E_mod, ν)
    catch e
        @warn "Initial guess failed in get_neoHookean3_S_and_H: $(e)"
        return E_vec, D 
    end

    r_old .= S_vec .- dpsi_dE_old
    resnorm_old = norm(r_old)

    tol_check = S_norm < 1e-12 ? tol_abs : max(tol_rel * S_norm, tol_abs)

    for it = 1:maxiter
        ΔE = H_old \ r_old

        alpha = alpha_max
        max_ls_iter = 10 

        local r_trial = zeros(6)
        local H_trial = zeros(6, 6)

        for ls_it = 1:max_ls_iter
            E_vec_trial = E_vec_old + alpha * ΔE

            C11_t = 2.0 * E_vec_trial[1] + 1.0
            C22_t = 2.0 * E_vec_trial[2] + 1.0
            C33_t = 2.0 * E_vec_trial[3] + 1.0
            C12_t = E_vec_trial[4] 
            C23_t = E_vec_trial[5]
            C13_t = E_vec_trial[6]

            detC_trial = C11_t*(C22_t*C33_t - C23_t^2) - C12_t*(C12_t*C33_t - C13_t*C23_t) + C13_t*(C12_t*C23_t - C13_t*C22_t)

            if detC_trial <= 1e-12
                alpha /= 2.0
                if alpha < 1e-5
                    @warn "Line search failed: Reached non-physical state (det(C)<=0)."
                    break
                end
                continue
            end

            dpsi_dE_trial, H_trial = get_neoHookean2_3D_S_and_H(E_vec_trial, E_mod, ν) 

            r_trial .= S_vec .- dpsi_dE_trial
            resnorm_trial = norm(r_trial)

            if resnorm_trial < resnorm_old
                E_vec .= E_vec_trial
                r_old .= r_trial
                H_old .= H_trial
                resnorm_old = resnorm_trial
                break
            else
                alpha /= 2.0
            end
        end 

        E_vec_old .= E_vec

        if resnorm_old < tol_check
            converged = true
            D .= inv(H_old)
            break
        end
    end 

    return E_vec, D
end

function Δ∫∫EᵢⱼSᵢⱼdxdy_HR_NeoHookean(e::Int,ap::T,k::AbstractMatrix{Float64},D_hist) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    for (g, ξ) in enumerate(𝓖)
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]

        S₁₁ = 0.0; S₂₂ = 0.0; S₁₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            S₁₁ += N[i]*xᵢ.dₛ₁₁
            S₂₂ += N[i]*xᵢ.dₛ₂₂
            S₁₂ += N[i]*xᵢ.dₛ₁₂
        end
        S_vec = [S₁₁, S₂₂, S₁₂]

        D = D_hist[e,g]
        D11 = D[1,1]; D12 = D[1,2]; D13 = D[1,3]
        D21 = D[2,1]; D22 = D[2,2]; D23 = D[2,3]
        D31 = D[3,1]; D32 = D[3,2]; D33 = D[3,3]

        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += N[i]*D11*N[j]*𝑤
                k[3*I-2,3*J-1] += N[i]*D12*N[j]*𝑤
                 k[3*I-2, 3*J] += N[i]*D13*N[j]*𝑤
                k[3*I-1,3*J-2] += N[i]*D21*N[j]*𝑤
                k[3*I-1,3*J-1] += N[i]*D22*N[j]*𝑤
                k[3*I-1,3*J] += N[i]*D23*N[j]*𝑤

                k[3*I  ,3*J-2 ] += N[i]*D31*N[j]*𝑤
                k[3*I  ,3*J-1 ] += N[i]*D32*N[j]*𝑤
                k[3*I  ,3*J  ] += N[i]*D33*N[j]*𝑤

            end
        end
    end
end

function Δ∫∫∫EᵢⱼSᵢⱼdxdydz_HR_NeoHookean(e::Int, ap::T, k::AbstractMatrix{Float64}, D_hist) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    for (g, ξ) in enumerate(𝓖)
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]

        D = D_hist[e,g]

        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            Ni_w = N[i] * 𝑤
            for (j, xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                fac = Ni_w * N[j]

                for r in 1:6
                    for c in 1:6
                        k[6*I-6+r, 6*J-6+c] += fac * D[r, c]
                    end
                end
            end
        end
    end
end

function ∫∫EᵢⱼSᵢⱼdxdy_HR_NeoHookean(e::Int, ap::T,f::AbstractVector{Float64},E_cons_hist, D_hist) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    for (g, ξ) in enumerate(𝓖)
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]

        S₁₁ = 0.0; S₂₂ = 0.0; S₁₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            S₁₁ += N[i]*xᵢ.dₛ₁₁
            S₂₂ += N[i]*xᵢ.dₛ₂₂
            S₁₂ += N[i]*xᵢ.dₛ₁₂
        end
        S_vec = [S₁₁, S₂₂, S₁₂]

        E_vec, D = update_Econs!(ξ, S_vec, E_cons_hist[e,g])
        E₁₁, E₂₂, E₁₂ = E_vec
        E_cons_hist[e,g] .= E_vec
        D_hist[e,g]      .= D

        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[3*I-2] += N[i]*E₁₁*𝑤
            f[3*I-1] += N[i]*E₂₂*𝑤

            f[3*I]   += N[i]*E₁₂*𝑤 

    end
end 
end

function ∫∫∫EᵢⱼSᵢⱼdxdydz_HR_NeoHookean(e::Int, ap::T, f::AbstractVector{Float64}, E_cons_hist, D_hist) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    for (g, ξ) in enumerate(𝓖)
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]

        S₁₁=0.0; S₂₂=0.0; S₃₃=0.0; S₁₂=0.0; S₂₃=0.0; S₁₃=0.0
        for (i, xᵢ) in enumerate(𝓒)
            S₁₁ += N[i]*xᵢ.dₛ₁₁
            S₂₂ += N[i]*xᵢ.dₛ₂₂
            S₃₃ += N[i]*xᵢ.dₛ₃₃
            S₁₂ += N[i]*xᵢ.dₛ₁₂
            S₂₃ += N[i]*xᵢ.dₛ₂₃
            S₁₃ += N[i]*xᵢ.dₛ₁₃
        end
        S_vec = [S₁₁, S₂₂, S₃₃, S₁₂, S₂₃, S₁₃]

        E_vec, D = update_Econs_3D!(ξ, S_vec, E_cons_hist[e,g])
        E₁₁, E₂₂, E₃₃, E₁₂, E₂₃, E₁₃ = E_vec

        E_cons_hist[e,g] .= E_vec
        D_hist[e,g]      .= D

        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[6*I-5] += N[i] * E₁₁ * 𝑤
            f[6*I-4] += N[i] * E₂₂ * 𝑤
            f[6*I-3] += N[i] * E₃₃ * 𝑤
            f[6*I-2] += N[i] * E₁₂ * 𝑤
            f[6*I-1] += N[i] * E₂₃ * 𝑤
            f[6*I]   += N[i] * E₁₃ * 𝑤
        end
    end
end

function ∫∫δSᵢⱼEᵢⱼdxdy_HR(aₛ::T,aᵤ::S, f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖

    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]

        𝑤 = ξₛ.𝑤

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

        E₁₁ = 0.5 * (F₁₁*F₁₁ + F₂₁*F₂₁ - 1.0)
        E₂₂ = 0.5 * (F₁₂*F₁₂ + F₂₂*F₂₂ - 1.0)
        E₁₂ = 0.5 * (F₁₁*F₁₂ + F₂₁*F₂₂) 

        for (i, xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼

            f[3*I-2] -= Nₛ[i] * E₁₁ * 𝑤

            f[3*I-1] -= Nₛ[i] * E₂₂ * 𝑤

            f[3*I]   -= 2.0*Nₛ[i]  * E₁₂ * 𝑤
        end
    end
end

function ∫∫∫δSᵢⱼEᵢⱼdxdydz_HR(aₛ::T, aᵤ::S, f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒; 𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒; 𝓖ᵤ = aᵤ.𝓖

    for (ξₛ, ξᵤ) in zip(𝓖ₛ, 𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        B₃ = ξᵤ[:∂𝝭∂z]
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        𝑤 = ξₛ.𝑤

        F₁₁ = 1.0; F₁₂ = 0.0; F₁₃ = 0.0
        F₂₁ = 0.0; F₂₂ = 1.0; F₂₃ = 0.0
        F₃₁ = 0.0; F₃₂ = 0.0; F₃₃ = 1.0

        for (i, xᵢ) in enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁; F₁₂ += B₂[i]*xᵢ.d₁; F₁₃ += B₃[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂; F₂₂ += B₂[i]*xᵢ.d₂; F₂₃ += B₃[i]*xᵢ.d₂
            F₃₁ += B₁[i]*xᵢ.d₃; F₃₂ += B₂[i]*xᵢ.d₃; F₃₃ += B₃[i]*xᵢ.d₃
        end

        E₁₁ = 0.5 * (F₁₁*F₁₁ + F₂₁*F₂₁ + F₃₁*F₃₁ - 1.0)
        E₂₂ = 0.5 * (F₁₂*F₁₂ + F₂₂*F₂₂ + F₃₂*F₃₂ - 1.0)
        E₃₃ = 0.5 * (F₁₃*F₁₃ + F₂₃*F₂₃ + F₃₃*F₃₃ - 1.0)
        E₁₂ = 0.5 * (F₁₁*F₁₂ + F₂₁*F₂₂ + F₃₁*F₃₂) 
        E₂₃ = 0.5 * (F₁₂*F₁₃ + F₂₂*F₂₃ + F₃₂*F₃₃) 
        E₁₃ = 0.5 * (F₁₁*F₁₃ + F₂₁*F₂₃ + F₃₁*F₃₃) 

        for (i, xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            f[6*I-5] -= Nₛ[i] * E₁₁ * 𝑤
            f[6*I-4] -= Nₛ[i] * E₂₂ * 𝑤
            f[6*I-3] -= Nₛ[i] * E₃₃ * 𝑤
            f[6*I-2] -= Nₛ[i] * 2.0 * E₁₂ * 𝑤
            f[6*I-1] -= Nₛ[i] * 2.0 * E₂₃ * 𝑤
            f[6*I]   -= Nₛ[i] * 2.0 * E₁₃ * 𝑤
        end
    end
end

function ∫∫δSᵢⱼΔEᵢⱼdxdy_HR(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
   𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]

        𝑤 = ξₛ.𝑤

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

        E₁₁ = 0.5 * (F₁₁*F₁₁ + F₂₁*F₂₁ - 1.0)
        E₂₂ = 0.5 * (F₁₂*F₁₂ + F₂₂*F₂₂ - 1.0)
        E₁₂ = 0.5 * (F₁₁*F₁₂ + F₂₁*F₂₂) 

        for (i, xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[3*I-2, 2*J-1] -= Nₛ[i] * F₁₁ * B₁[j] * 𝑤  
                k[3*I-2, 2*J]   -= Nₛ[i] * F₂₁ * B₁[j] * 𝑤  

                k[3*I-1, 2*J-1] -= Nₛ[i] * F₁₂ * B₂[j] * 𝑤
                k[3*I-1, 2*J]   -= Nₛ[i] * F₂₂ * B₂[j] * 𝑤

                k[3*I, 2*J-1] -= Nₛ[i] * (F₁₂* B₁[j]+F₁₁ * B₂[j])* 𝑤
                k[3*I, 2*J]   -= Nₛ[i] * (F₂₂* B₁[j]+F₂₁ * B₂[j]) * 𝑤

            end
        end
    end
end

function ∫∫∫δSᵢⱼΔEᵢⱼdxdydz_HR(aₛ::T, aᵤ::S, k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒; 𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒; 𝓖ᵤ = aᵤ.𝓖

    for (ξₛ, ξᵤ) in zip(𝓖ₛ, 𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]; B₂ = ξᵤ[:∂𝝭∂y]; B₃ = ξᵤ[:∂𝝭∂z]
        Nₛ = ξₛ[:𝝭]
        𝑤 = ξₛ.𝑤

        F₁₁ = 1.0; F₁₂ = 0.0; F₁₃ = 0.0
        F₂₁ = 0.0; F₂₂ = 1.0; F₂₃ = 0.0
        F₃₁ = 0.0; F₃₂ = 0.0; F₃₃ = 1.0

        for (i, xᵢ) in enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁; F₁₂ += B₂[i]*xᵢ.d₁; F₁₃ += B₃[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂; F₂₂ += B₂[i]*xᵢ.d₂; F₂₃ += B₃[i]*xᵢ.d₂
            F₃₁ += B₁[i]*xᵢ.d₃; F₃₂ += B₂[i]*xᵢ.d₃; F₃₃ += B₃[i]*xᵢ.d₃
        end

        for (i, xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j, xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼

                k[6*I-5, 3*J-2] -= Nₛ[i] * F₁₁ * B₁[j] * 𝑤
                k[6*I-5, 3*J-1] -= Nₛ[i] * F₂₁ * B₁[j] * 𝑤
                k[6*I-5, 3*J]   -= Nₛ[i] * F₃₁ * B₁[j] * 𝑤

                k[6*I-4, 3*J-2] -= Nₛ[i] * F₁₂ * B₂[j] * 𝑤
                k[6*I-4, 3*J-1] -= Nₛ[i] * F₂₂ * B₂[j] * 𝑤
                k[6*I-4, 3*J]   -= Nₛ[i] * F₃₂ * B₂[j] * 𝑤

                k[6*I-3, 3*J-2] -= Nₛ[i] * F₁₃ * B₃[j] * 𝑤
                k[6*I-3, 3*J-1] -= Nₛ[i] * F₂₃ * B₃[j] * 𝑤
                k[6*I-3, 3*J]   -= Nₛ[i] * F₃₃ * B₃[j] * 𝑤

                k[6*I-2, 3*J-2] -= Nₛ[i] * (F₁₁*B₂[j] + F₁₂*B₁[j]) * 𝑤
                k[6*I-2, 3*J-1] -= Nₛ[i] * (F₂₁*B₂[j] + F₂₂*B₁[j]) * 𝑤
                k[6*I-2, 3*J]   -= Nₛ[i] * (F₃₁*B₂[j] + F₃₂*B₁[j]) * 𝑤

                k[6*I-1, 3*J-2] -= Nₛ[i] * (F₁₂*B₃[j] + F₁₃*B₂[j]) * 𝑤
                k[6*I-1, 3*J-1] -= Nₛ[i] * (F₂₂*B₃[j] + F₂₃*B₂[j]) * 𝑤
                k[6*I-1, 3*J]   -= Nₛ[i] * (F₃₂*B₃[j] + F₃₃*B₂[j]) * 𝑤

                k[6*I,   3*J-2] -= Nₛ[i] * (F₁₁*B₃[j] + F₁₃*B₁[j]) * 𝑤
                k[6*I,   3*J-1] -= Nₛ[i] * (F₂₁*B₃[j] + F₂₃*B₁[j]) * 𝑤
                k[6*I,   3*J]   -= Nₛ[i] * (F₃₁*B₃[j] + F₃₃*B₁[j]) * 𝑤
            end
        end
    end
end

function ∫∫SᵢⱼδEᵢⱼdxdy_HR(aₛ::T,aᵤ::S, f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖

    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]

        𝑤 = ξₛ.𝑤

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

        E₁₁ = 0.5 * (F₁₁*F₁₁ + F₂₁*F₂₁ - 1.0)
        E₂₂ = 0.5 * (F₁₂*F₁₂ + F₂₂*F₂₂ - 1.0)
        E₁₂ = 0.5 * (F₁₁*F₁₂ + F₂₁*F₂₂) 

        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
        end

        P₁₁ = F₁₁ * S₁₁ + F₁₂ * S₁₂
        P₁₂ = F₁₁ * S₁₂ + F₁₂ * S₂₂
        P₂₁ = F₂₁ * S₁₁ + F₂₂ * S₁₂
        P₂₂ = F₂₁ * S₁₂ + F₂₂ * S₂₂

        for (i, xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
           f[2*I-1] -= (P₁₁ * B₁[i] + P₁₂ * B₂[i]) * 𝑤
           f[2*I]   -= (P₂₁ * B₁[i] + P₂₂ * B₂[i]) * 𝑤

        end
    end
end

function ∫∫∫SᵢⱼδEᵢⱼdxdydz_HR(aₛ::T, aᵤ::S, f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒; 𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒; 𝓖ᵤ = aᵤ.𝓖

    for (ξₛ, ξᵤ) in zip(𝓖ₛ, 𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]; B₂ = ξᵤ[:∂𝝭∂y]; B₃ = ξᵤ[:∂𝝭∂z]
        Nₛ = ξₛ[:𝝭]
        𝑤 = ξₛ.𝑤

        F₁₁ = 1.0; F₁₂ = 0.0; F₁₃ = 0.0
        F₂₁ = 0.0; F₂₂ = 1.0; F₂₃ = 0.0
        F₃₁ = 0.0; F₃₂ = 0.0; F₃₃ = 1.0

        for (i, xᵢ) in enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁; F₁₂ += B₂[i]*xᵢ.d₁; F₁₃ += B₃[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂; F₂₂ += B₂[i]*xᵢ.d₂; F₂₃ += B₃[i]*xᵢ.d₂
            F₃₁ += B₁[i]*xᵢ.d₃; F₃₂ += B₂[i]*xᵢ.d₃; F₃₃ += B₃[i]*xᵢ.d₃
        end

        S₁₁ = 0.0; S₂₂ = 0.0; S₃₃ = 0.0
        S₁₂ = 0.0; S₂₃ = 0.0; S₁₃ = 0.0
        for (i, xᵢ) in enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₃₃ += Nₛ[i]*xᵢ.dₛ₃₃
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
           S₂₃ += Nₛ[i]*xᵢ.dₛ₂₃
           S₁₃ += Nₛ[i]*xᵢ.dₛ₁₃
        end

        P₁₁ = F₁₁*S₁₁ + F₁₂*S₁₂ + F₁₃*S₁₃
        P₁₂ = F₁₁*S₁₂ + F₁₂*S₂₂ + F₁₃*S₂₃
        P₁₃ = F₁₁*S₁₃ + F₁₂*S₂₃ + F₁₃*S₃₃

        P₂₁ = F₂₁*S₁₁ + F₂₂*S₁₂ + F₂₃*S₁₃
        P₂₂ = F₂₁*S₁₂ + F₂₂*S₂₂ + F₂₃*S₂₃
        P₂₃ = F₂₁*S₁₃ + F₂₂*S₂₃ + F₂₃*S₃₃

        P₃₁ = F₃₁*S₁₁ + F₃₂*S₁₂ + F₃₃*S₁₃
        P₃₂ = F₃₁*S₁₂ + F₃₂*S₂₂ + F₃₃*S₂₃
        P₃₃ = F₃₁*S₁₃ + F₃₂*S₂₃ + F₃₃*S₃₃

        for (i, xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
            f[3*I-2] -= (P₁₁*B₁[i] + P₁₂*B₂[i] + P₁₃*B₃[i]) * 𝑤
            f[3*I-1] -= (P₂₁*B₁[i] + P₂₂*B₂[i] + P₂₃*B₃[i]) * 𝑤
            f[3*I]   -= (P₃₁*B₁[i] + P₃₂*B₂[i] + P₃₃*B₃[i]) * 𝑤
        end
    end
end

function ∫∫SᵢⱼδΔEᵢⱼdxdy_HR(aₛ::T,aᵤ::S, k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖

    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]

        𝑤 = ξₛ.𝑤

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

        E₁₁ = 0.5 * (F₁₁*F₁₁ + F₂₁*F₂₁ - 1.0)
        E₂₂ = 0.5 * (F₁₂*F₁₂ + F₂₂*F₂₂ - 1.0)
        E₁₂ = 0.5 * (F₁₁*F₁₂ + F₂₁*F₂₂) 

        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
        end

        P₁₁ = F₁₁ * S₁₁ + F₁₂ * S₁₂
        P₁₂ = F₁₁ * S₁₂ + F₁₂ * S₂₂
        P₂₁ = F₂₁ * S₁₁ + F₂₂ * S₁₂
        P₂₂ = F₂₁ * S₁₂ + F₂₂ * S₂₂

        for (i, xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
            b1_i = B₁[i]
            b2_i = B₂[i]
            for (j, xⱼ) in enumerate(𝓒ᵤ)
                 J = xⱼ.𝐼
                b1_j = B₁[j]
                b2_j = B₂[j]
                term1 = S₁₁ * b1_i * b1_j
                term2 = S₁₂ * (b1_i * b2_j + b2_i * b1_j) 
                term3 = S₂₂ * b2_i * b2_j
                g = (term1 + term2 + term3) * 𝑤
                k[2*I-1, 2*J-1] -= g
                k[2*I,   2*J]   -= g

            end

        end
    end
end

function ∫∫∫SᵢⱼδΔEᵢⱼdxdydz_HR(aₛ::T, aᵤ::S, k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒; 𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒; 𝓖ᵤ = aᵤ.𝓖

    for (ξₛ, ξᵤ) in zip(𝓖ₛ, 𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]; B₂ = ξᵤ[:∂𝝭∂y]; B₃ = ξᵤ[:∂𝝭∂z]
        Nₛ = ξₛ[:𝝭]
        𝑤 = ξₛ.𝑤

        S₁₁ = 0.0; S₂₂ = 0.0; S₃₃ = 0.0
        S₁₂ = 0.0; S₂₃ = 0.0; S₁₃ = 0.0
        for (i, xᵢ) in enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₃₃ += Nₛ[i]*xᵢ.dₛ₃₃
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
           S₂₃ += Nₛ[i]*xᵢ.dₛ₂₃
           S₁₃ += Nₛ[i]*xᵢ.dₛ₁₃
        end

        for (i, xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
            b1_i = B₁[i]; b2_i = B₂[i]; b3_i = B₃[i]

            for (j, xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                b1_j = B₁[j]; b2_j = B₂[j]; b3_j = B₃[j]

                term1 = S₁₁ * b1_i * b1_j
                term2 = S₂₂ * b2_i * b2_j
                term3 = S₃₃ * b3_i * b3_j
                term4 = S₁₂ * (b1_i * b2_j + b2_i * b1_j)
                term5 = S₂₃ * (b2_i * b3_j + b3_i * b2_j)
                term6 = S₁₃ * (b1_i * b3_j + b3_i * b1_j)

                g = (term1 + term2 + term3 + term4 + term5 + term6) * 𝑤

                k[3*I-2, 3*J-2] -= g
                k[3*I-1, 3*J-1] -= g
                k[3*I,   3*J]   -= g
            end
        end
    end
end

function ∫∫Stabilization_Operator_HR(
    aₛ::T, aᵤ::S,
    f_u::AbstractVector{Float64}, f_S::AbstractVector{Float64},
    k_uu::AbstractMatrix{Float64},
    k_Su::AbstractMatrix{Float64}, k_SS::AbstractMatrix{Float64}
) where {T<:AbstractElement, S<:AbstractElement}

    𝓒ₛ = aₛ.𝓒
    𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒
    𝓖ᵤ = aᵤ.𝓖

    for (ξₛ, ξᵤ) in zip(𝓖ₛ, 𝓖ᵤ)

        τ = ξₛ.τ
        𝑤 = ξₛ.𝑤

        B₁  = ξᵤ[:∂𝝭∂x]
        B₂  = ξᵤ[:∂𝝭∂y]
        B₁₁ = ξᵤ[:∂²𝝭∂x²]
        B₂₂ = ξᵤ[:∂²𝝭∂y²]
        B₁₂ = ξᵤ[:∂²𝝭∂x∂y]

        Nₛ  = ξₛ[:𝝭]
        Bₛ₁ = ξₛ[:∂𝝭∂x]
        Bₛ₂ = ξₛ[:∂𝝭∂y]

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0

        u1_11 = 0.0
        u1_12 = 0.0
        u1_22 = 0.0

        u2_11 = 0.0
        u2_12 = 0.0
        u2_22 = 0.0

        for (m, xₘ) in enumerate(𝓒ᵤ)
            u1 = xₘ.d₁
            u2 = xₘ.d₂

            F₁₁ += B₁[m] * u1
            F₁₂ += B₂[m] * u1
            F₂₁ += B₁[m] * u2
            F₂₂ += B₂[m] * u2

            u1_11 += B₁₁[m] * u1
            u1_12 += B₁₂[m] * u1
            u1_22 += B₂₂[m] * u1

            u2_11 += B₁₁[m] * u2
            u2_12 += B₁₂[m] * u2
            u2_22 += B₂₂[m] * u2
        end

        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0

        divS_1 = 0.0
        divS_2 = 0.0

        for (m, xₘ) in enumerate(𝓒ₛ)
            s11 = xₘ.dₛ₁₁
            s22 = xₘ.dₛ₂₂
            s12 = xₘ.dₛ₁₂

            S₁₁ += Nₛ[m] * s11
            S₂₂ += Nₛ[m] * s22
            S₁₂ += Nₛ[m] * s12

            divS_1 += Bₛ₁[m] * s11 + Bₛ₂[m] * s12
            divS_2 += Bₛ₁[m] * s12 + Bₛ₂[m] * s22
        end

        b₁ = ξₛ.b₁
        b₂ = ξₛ.b₂

        comp_1 = S₁₁*u1_11 + 2.0*S₁₂*u1_12 + S₂₂*u1_22
        comp_2 = S₁₁*u2_11 + 2.0*S₁₂*u2_12 + S₂₂*u2_22

        R1 = F₁₁*divS_1 + F₁₂*divS_2 + comp_1 + b₁
        R2 = F₂₁*divS_1 + F₂₂*divS_2 + comp_2 + b₂

        for (i, xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            N_i = Nₛ[i]

            vS11_1 = F₁₁*Bₛ₁[i] + u1_11*N_i
            vS11_2 = F₂₁*Bₛ₁[i] + u2_11*N_i

            vS22_1 = F₁₂*Bₛ₂[i] + u1_22*N_i
            vS22_2 = F₂₂*Bₛ₂[i] + u2_22*N_i

            vS12_1 = F₁₁*Bₛ₂[i] + F₁₂*Bₛ₁[i] + 2.0*u1_12*N_i
            vS12_2 = F₂₁*Bₛ₂[i] + F₂₂*Bₛ₁[i] + 2.0*u2_12*N_i

            f_S[3*I-2] -= τ * (vS11_1*R1 + vS11_2*R2) * 𝑤
            f_S[3*I-1] -= τ * (vS22_1*R1 + vS22_2*R2) * 𝑤
            f_S[3*I]   -= τ * (vS12_1*R1 + vS12_2*R2) * 𝑤

            for (j, xⱼ) in enumerate(𝓒ₛ)
                J = xⱼ.𝐼
                N_j = Nₛ[j]

                rS11_1 = F₁₁*Bₛ₁[j] + u1_11*N_j
                rS11_2 = F₂₁*Bₛ₁[j] + u2_11*N_j

                rS22_1 = F₁₂*Bₛ₂[j] + u1_22*N_j
                rS22_2 = F₂₂*Bₛ₂[j] + u2_22*N_j

                rS12_1 = F₁₁*Bₛ₂[j] + F₁₂*Bₛ₁[j] + 2.0*u1_12*N_j
                rS12_2 = F₂₁*Bₛ₂[j] + F₂₂*Bₛ₁[j] + 2.0*u2_12*N_j

                k_SS[3*I-2, 3*J-2] += τ * (vS11_1*rS11_1 + vS11_2*rS11_2) * 𝑤
                k_SS[3*I-2, 3*J-1] += τ * (vS11_1*rS22_1 + vS11_2*rS22_2) * 𝑤
                k_SS[3*I-2, 3*J]   += τ * (vS11_1*rS12_1 + vS11_2*rS12_2) * 𝑤

                k_SS[3*I-1, 3*J-2] += τ * (vS22_1*rS11_1 + vS22_2*rS11_2) * 𝑤
                k_SS[3*I-1, 3*J-1] += τ * (vS22_1*rS22_1 + vS22_2*rS22_2) * 𝑤
                k_SS[3*I-1, 3*J]   += τ * (vS22_1*rS12_1 + vS22_2*rS12_2) * 𝑤

                k_SS[3*I,   3*J-2] += τ * (vS12_1*rS11_1 + vS12_2*rS11_2) * 𝑤
                k_SS[3*I,   3*J-1] += τ * (vS12_1*rS22_1 + vS12_2*rS22_2) * 𝑤
                k_SS[3*I,   3*J]   += τ * (vS12_1*rS12_1 + vS12_2*rS12_2) * 𝑤
            end
        end

        for (i, xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼

            vu_i = B₁[i]*divS_1 + B₂[i]*divS_2 +
                   S₁₁*B₁₁[i] + 2.0*S₁₂*B₁₂[i] + S₂₂*B₂₂[i]

            f_u[2*I-1] -= τ * vu_i * R1 * 𝑤
            f_u[2*I]   -= τ * vu_i * R2 * 𝑤

            for (j, xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼

                ru_j = B₁[j]*divS_1 + B₂[j]*divS_2 +
                       S₁₁*B₁₁[j] + 2.0*S₁₂*B₁₂[j] + S₂₂*B₂₂[j]

                k_uu[2*I-1, 2*J-1] += τ * vu_i * ru_j * 𝑤
                k_uu[2*I,   2*J]   += τ * vu_i * ru_j * 𝑤
            end
        end

        for (i, xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            N_i = Nₛ[i]

            vS11_1 = F₁₁*Bₛ₁[i] + u1_11*N_i
            vS11_2 = F₂₁*Bₛ₁[i] + u2_11*N_i

            vS22_1 = F₁₂*Bₛ₂[i] + u1_22*N_i
            vS22_2 = F₂₂*Bₛ₂[i] + u2_22*N_i

            vS12_1 = F₁₁*Bₛ₂[i] + F₁₂*Bₛ₁[i] + 2.0*u1_12*N_i
            vS12_2 = F₂₁*Bₛ₂[i] + F₂₂*Bₛ₁[i] + 2.0*u2_12*N_i

            for (j, xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼

                ru_j = B₁[j]*divS_1 + B₂[j]*divS_2 +
                       S₁₁*B₁₁[j] + 2.0*S₁₂*B₁₂[j] + S₂₂*B₂₂[j]

                G11 = B₁[j]*Bₛ₁[i] + B₁₁[j]*N_i
                G22 = B₂[j]*Bₛ₂[i] + B₂₂[j]*N_i
                G12 = B₁[j]*Bₛ₂[i] + B₂[j]*Bₛ₁[i] + 2.0*B₁₂[j]*N_i

                k_Su[3*I-2, 2*J-1] += τ * (vS11_1*ru_j + G11*R1) * 𝑤
                k_Su[3*I-1, 2*J-1] += τ * (vS22_1*ru_j + G22*R1) * 𝑤
                k_Su[3*I,   2*J-1] += τ * (vS12_1*ru_j + G12*R1) * 𝑤

                k_Su[3*I-2, 2*J]   += τ * (vS11_2*ru_j + G11*R2) * 𝑤
                k_Su[3*I-1, 2*J]   += τ * (vS22_2*ru_j + G22*R2) * 𝑤
                k_Su[3*I,   2*J]   += τ * (vS12_2*ru_j + G12*R2) * 𝑤
            end
        end
    end
end

function ∫∫∫Stabilization_Operator_HR(
    aₛ::T, aᵤ::S,
    f_u::AbstractVector{Float64}, f_S::AbstractVector{Float64},
    k_uu::AbstractMatrix{Float64},
    k_Su::AbstractMatrix{Float64}, k_SS::AbstractMatrix{Float64}
) where {T<:AbstractElement, S<:AbstractElement}

    𝓒ₛ = aₛ.𝓒; 𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒; 𝓖ᵤ = aᵤ.𝓖

    for (ξₛ, ξᵤ) in zip(𝓖ₛ, 𝓖ᵤ)
        τ = ξₛ.τ
        B₁ = ξᵤ[:∂𝝭∂x]; B₂ = ξᵤ[:∂𝝭∂y]; B₃ = ξᵤ[:∂𝝭∂z]

        B₁₁ = ξᵤ[:∂²𝝭∂x²]
        B₂₂ = ξᵤ[:∂²𝝭∂y²]
        B₃₃ = ξᵤ[:∂²𝝭∂z²]
        B₁₂ = ξᵤ[:∂²𝝭∂x∂y]
        B₂₃ = ξᵤ[:∂²𝝭∂y∂z]
        B₁₃ = ξᵤ[:∂²𝝭∂x∂z]

        Nₛ  = ξₛ[:𝝭]
        Bₛ₁ = ξₛ[:∂𝝭∂x]
        Bₛ₂ = ξₛ[:∂𝝭∂y]
        Bₛ₃ = ξₛ[:∂𝝭∂z]

        𝑤 = ξₛ.𝑤

        F₁₁=1.0; F₁₂=0.0; F₁₃=0.0
        F₂₁=0.0; F₂₂=1.0; F₂₃=0.0
        F₃₁=0.0; F₃₂=0.0; F₃₃=1.0

        u1_11=0.0; u1_22=0.0; u1_33=0.0; u1_12=0.0; u1_23=0.0; u1_13=0.0
        u2_11=0.0; u2_22=0.0; u2_33=0.0; u2_12=0.0; u2_23=0.0; u2_13=0.0
        u3_11=0.0; u3_22=0.0; u3_33=0.0; u3_12=0.0; u3_23=0.0; u3_13=0.0

        for (m, xₘ) in enumerate(𝓒ᵤ)
            u1 = xₘ.d₁
            u2 = xₘ.d₂
            u3 = xₘ.d₃

            F₁₁ += B₁[m]*u1; F₁₂ += B₂[m]*u1; F₁₃ += B₃[m]*u1
            F₂₁ += B₁[m]*u2; F₂₂ += B₂[m]*u2; F₂₃ += B₃[m]*u2
            F₃₁ += B₁[m]*u3; F₃₂ += B₂[m]*u3; F₃₃ += B₃[m]*u3

            u1_11 += B₁₁[m]*u1; u1_22 += B₂₂[m]*u1; u1_33 += B₃₃[m]*u1
            u1_12 += B₁₂[m]*u1; u1_23 += B₂₃[m]*u1; u1_13 += B₁₃[m]*u1

            u2_11 += B₁₁[m]*u2; u2_22 += B₂₂[m]*u2; u2_33 += B₃₃[m]*u2
            u2_12 += B₁₂[m]*u2; u2_23 += B₂₃[m]*u2; u2_13 += B₁₃[m]*u2

            u3_11 += B₁₁[m]*u3; u3_22 += B₂₂[m]*u3; u3_33 += B₃₃[m]*u3
            u3_12 += B₁₂[m]*u3; u3_23 += B₂₃[m]*u3; u3_13 += B₁₃[m]*u3
        end

        S₁₁=0.0; S₂₂=0.0; S₃₃=0.0
        S₁₂=0.0; S₂₃=0.0; S₁₃=0.0

        divS_1=0.0; divS_2=0.0; divS_3=0.0

        for (m, xₘ) in enumerate(𝓒ₛ)
            s11 = xₘ.dₛ₁₁
            s22 = xₘ.dₛ₂₂
            s33 = xₘ.dₛ₃₃
            s12 = xₘ.dₛ₁₂
            s23 = xₘ.dₛ₂₃
            s13 = xₘ.dₛ₁₃

            S₁₁ += Nₛ[m]*s11
            S₂₂ += Nₛ[m]*s22
            S₃₃ += Nₛ[m]*s33
            S₁₂ += Nₛ[m]*s12
            S₂₃ += Nₛ[m]*s23
            S₁₃ += Nₛ[m]*s13

            divS_1 += Bₛ₁[m]*s11 + Bₛ₂[m]*s12 + Bₛ₃[m]*s13
            divS_2 += Bₛ₁[m]*s12 + Bₛ₂[m]*s22 + Bₛ₃[m]*s23
            divS_3 += Bₛ₁[m]*s13 + Bₛ₂[m]*s23 + Bₛ₃[m]*s33
        end

        b₁ = ξₛ.b₁
        b₂ = ξₛ.b₂
        b₃ = ξₛ.b₃

        comp_1 = S₁₁*u1_11 + S₂₂*u1_22 + S₃₃*u1_33 +
                 2.0*S₁₂*u1_12 + 2.0*S₂₃*u1_23 + 2.0*S₁₃*u1_13

        comp_2 = S₁₁*u2_11 + S₂₂*u2_22 + S₃₃*u2_33 +
                 2.0*S₁₂*u2_12 + 2.0*S₂₃*u2_23 + 2.0*S₁₃*u2_13

        comp_3 = S₁₁*u3_11 + S₂₂*u3_22 + S₃₃*u3_33 +
                 2.0*S₁₂*u3_12 + 2.0*S₂₃*u3_23 + 2.0*S₁₃*u3_13

        R1 = F₁₁*divS_1 + F₁₂*divS_2 + F₁₃*divS_3 + comp_1 + b₁
        R2 = F₂₁*divS_1 + F₂₂*divS_2 + F₂₃*divS_3 + comp_2 + b₂
        R3 = F₃₁*divS_1 + F₃₂*divS_2 + F₃₃*divS_3 + comp_3 + b₃

        for (i, xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼

            N_i = Nₛ[i]

            vS11_1 = F₁₁*Bₛ₁[i] + u1_11*N_i
            vS11_2 = F₂₁*Bₛ₁[i] + u2_11*N_i
            vS11_3 = F₃₁*Bₛ₁[i] + u3_11*N_i

            vS22_1 = F₁₂*Bₛ₂[i] + u1_22*N_i
            vS22_2 = F₂₂*Bₛ₂[i] + u2_22*N_i
            vS22_3 = F₃₂*Bₛ₂[i] + u3_22*N_i

            vS33_1 = F₁₃*Bₛ₃[i] + u1_33*N_i
            vS33_2 = F₂₃*Bₛ₃[i] + u2_33*N_i
            vS33_3 = F₃₃*Bₛ₃[i] + u3_33*N_i

            vS12_1 = F₁₁*Bₛ₂[i] + F₁₂*Bₛ₁[i] + 2.0*u1_12*N_i
            vS12_2 = F₂₁*Bₛ₂[i] + F₂₂*Bₛ₁[i] + 2.0*u2_12*N_i
            vS12_3 = F₃₁*Bₛ₂[i] + F₃₂*Bₛ₁[i] + 2.0*u3_12*N_i

            vS23_1 = F₁₂*Bₛ₃[i] + F₁₃*Bₛ₂[i] + 2.0*u1_23*N_i
            vS23_2 = F₂₂*Bₛ₃[i] + F₂₃*Bₛ₂[i] + 2.0*u2_23*N_i
            vS23_3 = F₃₂*Bₛ₃[i] + F₃₃*Bₛ₂[i] + 2.0*u3_23*N_i

            vS13_1 = F₁₁*Bₛ₃[i] + F₁₃*Bₛ₁[i] + 2.0*u1_13*N_i
            vS13_2 = F₂₁*Bₛ₃[i] + F₂₃*Bₛ₁[i] + 2.0*u2_13*N_i
            vS13_3 = F₃₁*Bₛ₃[i] + F₃₃*Bₛ₁[i] + 2.0*u3_13*N_i

            f_S[6*I-5] -= τ*(vS11_1*R1 + vS11_2*R2 + vS11_3*R3)*𝑤
            f_S[6*I-4] -= τ*(vS22_1*R1 + vS22_2*R2 + vS22_3*R3)*𝑤
            f_S[6*I-3] -= τ*(vS33_1*R1 + vS33_2*R2 + vS33_3*R3)*𝑤
            f_S[6*I-2] -= τ*(vS12_1*R1 + vS12_2*R2 + vS12_3*R3)*𝑤
            f_S[6*I-1] -= τ*(vS23_1*R1 + vS23_2*R2 + vS23_3*R3)*𝑤
            f_S[6*I]   -= τ*(vS13_1*R1 + vS13_2*R2 + vS13_3*R3)*𝑤

            for (j, xⱼ) in enumerate(𝓒ₛ)
                J = xⱼ.𝐼
                N_j = Nₛ[j]

                r11 = (
                    F₁₁*Bₛ₁[j] + u1_11*N_j,
                    F₂₁*Bₛ₁[j] + u2_11*N_j,
                    F₃₁*Bₛ₁[j] + u3_11*N_j
                )
                r22 = (
                    F₁₂*Bₛ₂[j] + u1_22*N_j,
                    F₂₂*Bₛ₂[j] + u2_22*N_j,
                    F₃₂*Bₛ₂[j] + u3_22*N_j
                )
                r33 = (
                    F₁₃*Bₛ₃[j] + u1_33*N_j,
                    F₂₃*Bₛ₃[j] + u2_33*N_j,
                    F₃₃*Bₛ₃[j] + u3_33*N_j
                )
                r12 = (
                    F₁₁*Bₛ₂[j] + F₁₂*Bₛ₁[j] + 2.0*u1_12*N_j,
                    F₂₁*Bₛ₂[j] + F₂₂*Bₛ₁[j] + 2.0*u2_12*N_j,
                    F₃₁*Bₛ₂[j] + F₃₂*Bₛ₁[j] + 2.0*u3_12*N_j
                )
                r23 = (
                    F₁₂*Bₛ₃[j] + F₁₃*Bₛ₂[j] + 2.0*u1_23*N_j,
                    F₂₂*Bₛ₃[j] + F₂₃*Bₛ₂[j] + 2.0*u2_23*N_j,
                    F₃₂*Bₛ₃[j] + F₃₃*Bₛ₂[j] + 2.0*u3_23*N_j
                )
                r13 = (
                    F₁₁*Bₛ₃[j] + F₁₃*Bₛ₁[j] + 2.0*u1_13*N_j,
                    F₂₁*Bₛ₃[j] + F₂₃*Bₛ₁[j] + 2.0*u2_13*N_j,
                    F₃₁*Bₛ₃[j] + F₃₃*Bₛ₁[j] + 2.0*u3_13*N_j
                )

                vlist = [
                    (vS11_1, vS11_2, vS11_3),
                    (vS22_1, vS22_2, vS22_3),
                    (vS33_1, vS33_2, vS33_3),
                    (vS12_1, vS12_2, vS12_3),
                    (vS23_1, vS23_2, vS23_3),
                    (vS13_1, vS13_2, vS13_3)
                ]

                rlist = [r11, r22, r33, r12, r23, r13]

                for aS in 1:6, bS in 1:6
                    k_SS[6*I-6+aS, 6*J-6+bS] +=
                        τ*(vlist[aS][1]*rlist[bS][1] +
                           vlist[aS][2]*rlist[bS][2] +
                           vlist[aS][3]*rlist[bS][3])*𝑤
                end
            end
        end

        for (i, xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼

            vu_i = B₁[i]*divS_1 + B₂[i]*divS_2 + B₃[i]*divS_3 +
                   S₁₁*B₁₁[i] + S₂₂*B₂₂[i] + S₃₃*B₃₃[i] +
                   2.0*S₁₂*B₁₂[i] + 2.0*S₂₃*B₂₃[i] + 2.0*S₁₃*B₁₃[i]

            f_u[3*I-2] -= τ*vu_i*R1*𝑤
            f_u[3*I-1] -= τ*vu_i*R2*𝑤
            f_u[3*I]   -= τ*vu_i*R3*𝑤

            for (j, xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼

                ru_j = B₁[j]*divS_1 + B₂[j]*divS_2 + B₃[j]*divS_3 +
                       S₁₁*B₁₁[j] + S₂₂*B₂₂[j] + S₃₃*B₃₃[j] +
                       2.0*S₁₂*B₁₂[j] + 2.0*S₂₃*B₂₃[j] + 2.0*S₁₃*B₁₃[j]

                k_uu[3*I-2, 3*J-2] += τ*vu_i*ru_j*𝑤
                k_uu[3*I-1, 3*J-1] += τ*vu_i*ru_j*𝑤
                k_uu[3*I,   3*J]   += τ*vu_i*ru_j*𝑤
            end

        end

        for (i, xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            N_i = Nₛ[i]

            vS11_1 = F₁₁*Bₛ₁[i] + u1_11*N_i
            vS11_2 = F₂₁*Bₛ₁[i] + u2_11*N_i
            vS11_3 = F₃₁*Bₛ₁[i] + u3_11*N_i

            vS22_1 = F₁₂*Bₛ₂[i] + u1_22*N_i
            vS22_2 = F₂₂*Bₛ₂[i] + u2_22*N_i
            vS22_3 = F₃₂*Bₛ₂[i] + u3_22*N_i

            vS33_1 = F₁₃*Bₛ₃[i] + u1_33*N_i
            vS33_2 = F₂₃*Bₛ₃[i] + u2_33*N_i
            vS33_3 = F₃₃*Bₛ₃[i] + u3_33*N_i

            vS12_1 = F₁₁*Bₛ₂[i] + F₁₂*Bₛ₁[i] + 2.0*u1_12*N_i
            vS12_2 = F₂₁*Bₛ₂[i] + F₂₂*Bₛ₁[i] + 2.0*u2_12*N_i
            vS12_3 = F₃₁*Bₛ₂[i] + F₃₂*Bₛ₁[i] + 2.0*u3_12*N_i

            vS23_1 = F₁₂*Bₛ₃[i] + F₁₃*Bₛ₂[i] + 2.0*u1_23*N_i
            vS23_2 = F₂₂*Bₛ₃[i] + F₂₃*Bₛ₂[i] + 2.0*u2_23*N_i
            vS23_3 = F₃₂*Bₛ₃[i] + F₃₃*Bₛ₂[i] + 2.0*u3_23*N_i

            vS13_1 = F₁₁*Bₛ₃[i] + F₁₃*Bₛ₁[i] + 2.0*u1_13*N_i
            vS13_2 = F₂₁*Bₛ₃[i] + F₂₃*Bₛ₁[i] + 2.0*u2_13*N_i
            vS13_3 = F₃₁*Bₛ₃[i] + F₃₃*Bₛ₁[i] + 2.0*u3_13*N_i

            for (j, xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼

                ru_j = B₁[j]*divS_1 + B₂[j]*divS_2 + B₃[j]*divS_3 +
                       S₁₁*B₁₁[j] + S₂₂*B₂₂[j] + S₃₃*B₃₃[j] +
                       2.0*S₁₂*B₁₂[j] + 2.0*S₂₃*B₂₃[j] + 2.0*S₁₃*B₁₃[j]

                G11 = B₁[j]*Bₛ₁[i] + B₁₁[j]*N_i
                G22 = B₂[j]*Bₛ₂[i] + B₂₂[j]*N_i
                G33 = B₃[j]*Bₛ₃[i] + B₃₃[j]*N_i
                G12 = B₁[j]*Bₛ₂[i] + B₂[j]*Bₛ₁[i] + 2.0*B₁₂[j]*N_i
                G23 = B₂[j]*Bₛ₃[i] + B₃[j]*Bₛ₂[i] + 2.0*B₂₃[j]*N_i
                G13 = B₁[j]*Bₛ₃[i] + B₃[j]*Bₛ₁[i] + 2.0*B₁₃[j]*N_i

                k_Su[6*I-5, 3*J-2] += τ*(vS11_1*ru_j + G11*R1)*𝑤
                k_Su[6*I-4, 3*J-2] += τ*(vS22_1*ru_j + G22*R1)*𝑤
                k_Su[6*I-3, 3*J-2] += τ*(vS33_1*ru_j + G33*R1)*𝑤
                k_Su[6*I-2, 3*J-2] += τ*(vS12_1*ru_j + G12*R1)*𝑤
                k_Su[6*I-1, 3*J-2] += τ*(vS23_1*ru_j + G23*R1)*𝑤
                k_Su[6*I,   3*J-2] += τ*(vS13_1*ru_j + G13*R1)*𝑤

                k_Su[6*I-5, 3*J-1] += τ*(vS11_2*ru_j + G11*R2)*𝑤
                k_Su[6*I-4, 3*J-1] += τ*(vS22_2*ru_j + G22*R2)*𝑤
                k_Su[6*I-3, 3*J-1] += τ*(vS33_2*ru_j + G33*R2)*𝑤
                k_Su[6*I-2, 3*J-1] += τ*(vS12_2*ru_j + G12*R2)*𝑤
                k_Su[6*I-1, 3*J-1] += τ*(vS23_2*ru_j + G23*R2)*𝑤
                k_Su[6*I,   3*J-1] += τ*(vS13_2*ru_j + G13*R2)*𝑤

                k_Su[6*I-5, 3*J] += τ*(vS11_3*ru_j + G11*R3)*𝑤
                k_Su[6*I-4, 3*J] += τ*(vS22_3*ru_j + G22*R3)*𝑤
                k_Su[6*I-3, 3*J] += τ*(vS33_3*ru_j + G33*R3)*𝑤
                k_Su[6*I-2, 3*J] += τ*(vS12_3*ru_j + G12*R3)*𝑤
                k_Su[6*I-1, 3*J] += τ*(vS23_3*ru_j + G23*R3)*𝑤
                k_Su[6*I,   3*J] += τ*(vS13_3*ru_j + G13*R3)*𝑤
            end
        end
    end
end

function ∫∫stabilization_Operator_HR(
    aₛ::T, aᵤ::S,
    f_u::AbstractVector{Float64}, f_S::AbstractVector{Float64},
    k_uu::AbstractMatrix{Float64}, 
    k_Su::AbstractMatrix{Float64}, 
    k_SS::AbstractMatrix{Float64}
) where {T<:AbstractElement, S<:AbstractElement}

    𝓒ₛ = aₛ.𝓒; 𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒; 𝓖ᵤ = aᵤ.𝓖

    for (ξₛ, ξᵤ) in zip(𝓖ₛ, 𝓖ᵤ)

        B₁ = ξᵤ[:∂𝝭∂x]; B₂ = ξᵤ[:∂𝝭∂y]
        Nₛ = ξₛ[:𝝭]
        𝑤 = ξₛ.𝑤

        F₁₁=1.0; F₁₂=0.0; F₂₁=0.0; F₂₂=1.0
        for (m, xₘ) in enumerate(𝓒ᵤ)
            u1 = xₘ.d₁; u2 = xₘ.d₂
            F₁₁ += B₁[m]*u1; F₁₂ += B₂[m]*u1
            F₂₁ += B₁[m]*u2; F₂₂ += B₂[m]*u2
        end

        S₁₁=0.0; S₂₂=0.0; S₁₂=0.0
        for (m, xₘ) in enumerate(𝓒ₛ)
            s11=xₘ.dₛ₁₁; s22=xₘ.dₛ₂₂; s12=xₘ.dₛ₁₂
            S₁₁ += Nₛ[m]*s11; S₂₂ += Nₛ[m]*s22; S₁₂ += Nₛ[m]*s12
        end

        R₁₁ = F₁₁*S₁₁ + F₁₂*S₁₂
        R₁₂ = F₁₁*S₁₂ + F₁₂*S₂₂
        R₂₁ = F₂₁*S₁₁ + F₂₂*S₁₂
        R₂₂ = F₂₁*S₁₂ + F₂₂*S₂₂

        for (i, xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            N_i = Nₛ[i]
             τ = xᵢ.β
              ℎ = xᵢ.ℎ
            a=1.0/ℎ^2

            f_S[3*I-2] -= a*τ * N_i * (F₁₁*R₁₁ + F₂₁*R₂₁) * 𝑤
            f_S[3*I-1] -= a*τ * N_i * (F₁₂*R₁₂ + F₂₂*R₂₂) * 𝑤
            f_S[3*I]   -= a*τ * N_i * (F₁₂*R₁₁ + F₁₁*R₁₂ + F₂₂*R₂₁ + F₂₁*R₂₂) * 𝑤

            for (j, xⱼ) in enumerate(𝓒ₛ)
                J = xⱼ.𝐼
                NN = N_i * Nₛ[j]

                term_11_11 = F₁₁^2 + F₂₁^2
                term_22_22 = F₁₂^2 + F₂₂^2
                term_cross = F₁₁*F₁₂ + F₂₁*F₂₂
                term_12_12 = F₁₂^2 + F₁₁^2 + F₂₂^2 + F₂₁^2

                k_SS[3*I-2, 3*J-2] += a*τ * NN * term_11_11 * 𝑤

                k_SS[3*I-2, 3*J]   += a*τ * NN * term_cross * 𝑤

                k_SS[3*I-1, 3*J-1] += a*τ * NN * term_22_22 * 𝑤
                k_SS[3*I-1, 3*J]   += a*τ * NN * term_cross * 𝑤

                k_SS[3*I,   3*J-2] += a*τ * NN * term_cross * 𝑤
                k_SS[3*I,   3*J-1] += a*τ * NN * term_cross * 𝑤
                k_SS[3*I,   3*J]   += a*τ * NN * term_12_12 * 𝑤
            end
        end

        for (i, xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
             τ = xᵢ.β
              ℎ = xᵢ.ℎ
              a=1.0/ℎ^2
            V1_i = B₁[i]*S₁₁ + B₂[i]*S₁₂
            V2_i = B₁[i]*S₁₂ + B₂[i]*S₂₂

            f_u[2*I-1] -= a*τ * (V1_i*R₁₁ + V2_i*R₁₂) * 𝑤
            f_u[2*I]   -= a*τ * (V1_i*R₂₁ + V2_i*R₂₂) * 𝑤

            for (j, xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                V1_j = B₁[j]*S₁₁ + B₂[j]*S₁₂
                V2_j = B₁[j]*S₁₂ + B₂[j]*S₂₂

                uu_val = a*τ * (V1_i*V1_j + V2_i*V2_j) * 𝑤

                k_uu[2*I-1, 2*J-1] += uu_val
                k_uu[2*I,   2*J]   += uu_val

            end

            for (j, xⱼ) in enumerate(𝓒ₛ)
                J = xⱼ.𝐼
                N_j = Nₛ[j]

                mat_x11 = V1_i * F₁₁ * N_j
                mat_x22 = V2_i * F₁₂ * N_j
                mat_x12 = V1_i * F₁₂ * N_j + V2_i * F₁₁ * N_j

                mat_y11 = V1_i * F₂₁ * N_j
                mat_y22 = V2_i * F₂₂ * N_j
                mat_y12 = V1_i * F₂₂ * N_j + V2_i * F₂₁ * N_j

                geo_x11 = B₁[i] * N_j * R₁₁
                geo_x22 = B₂[i] * N_j * R₁₂
                geo_x12 = B₂[i] * N_j * R₁₁ + B₁[i] * N_j * R₁₂

                geo_y11 = B₁[i] * N_j * R₂₁
                geo_y22 = B₂[i] * N_j * R₂₂
                geo_y12 = B₂[i] * N_j * R₂₁ + B₁[i] * N_j * R₂₂

                k_Su[3*J-2, 2*I-1] += a*τ * (mat_x11 + geo_x11) * 𝑤
                k_Su[3*J-1, 2*I-1] += a*τ * (mat_x22 + geo_x22) * 𝑤
                k_Su[3*J,   2*I-1] += a*τ * (mat_x12 + geo_x12) * 𝑤

                k_Su[3*J-2, 2*I]   += a*τ * (mat_y11 + geo_y11) * 𝑤
                k_Su[3*J-1, 2*I]   += a*τ * (mat_y22 + geo_y22) * 𝑤
                k_Su[3*J,   2*I]   += a*τ * (mat_y12 + geo_y12) * 𝑤
            end
        end
    end
end

function ∫∫∫stabilization_Operator_HR_New(
    aₛ::T, aᵤ::S,
    f_u::AbstractVector{Float64}, f_S::AbstractVector{Float64},
    k_uu::AbstractMatrix{Float64}, 
    k_Su::AbstractMatrix{Float64}, 
    k_SS::AbstractMatrix{Float64}
) where {T<:AbstractElement, S<:AbstractElement}

    𝓒ₛ = aₛ.𝓒; 𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒; 𝓖ᵤ = aᵤ.𝓖

    for (ξₛ, ξᵤ) in zip(𝓖ₛ, 𝓖ᵤ)
        τ = ξₛ.τ
        ℎ = ξₛ.ℎ
        a = 0.1/ ℎ^2

        B₁ = ξᵤ[:∂𝝭∂x]; B₂ = ξᵤ[:∂𝝭∂y]; B₃ = ξᵤ[:∂𝝭∂z]
        Nₛ = ξₛ[:𝝭]
        𝑤 = ξₛ.𝑤

        F₁₁=1.0; F₁₂=0.0; F₁₃=0.0
        F₂₁=0.0; F₂₂=1.0; F₂₃=0.0
        F₃₁=0.0; F₃₂=0.0; F₃₃=1.0

        for (m, xₘ) in enumerate(𝓒ᵤ)
            u1 = xₘ.d₁; u2 = xₘ.d₂; u3 = xₘ.d₃
            F₁₁ += B₁[m]*u1; F₁₂ += B₂[m]*u1; F₁₃ += B₃[m]*u1
            F₂₁ += B₁[m]*u2; F₂₂ += B₂[m]*u2; F₂₃ += B₃[m]*u2
            F₃₁ += B₁[m]*u3; F₃₂ += B₂[m]*u3; F₃₃ += B₃[m]*u3
        end

        S₁₁=0.0; S₂₂=0.0; S₃₃=0.0; S₁₂=0.0; S₂₃=0.0; S₁₃=0.0
        for (m, xₘ) in enumerate(𝓒ₛ)
            s11=xₘ.dₛ₁₁; s22=xₘ.dₛ₂₂; s33=xₘ.dₛ₃₃
            s12=xₘ.dₛ₁₂; s23=xₘ.dₛ₂₃; s13=xₘ.dₛ₁₃
            S₁₁ += Nₛ[m]*s11; S₂₂ += Nₛ[m]*s22; S₃₃ += Nₛ[m]*s33
            S₁₂ += Nₛ[m]*s12; S₂₃ += Nₛ[m]*s23; S₁₃ += Nₛ[m]*s13
        end

        R₁₁ = F₁₁*S₁₁ + F₁₂*S₁₂ + F₁₃*S₁₃
        R₁₂ = F₁₁*S₁₂ + F₁₂*S₂₂ + F₁₃*S₂₃
        R₁₃ = F₁₁*S₁₃ + F₁₂*S₂₃ + F₁₃*S₃₃

        R₂₁ = F₂₁*S₁₁ + F₂₂*S₁₂ + F₂₃*S₁₃
        R₂₂ = F₂₁*S₁₂ + F₂₂*S₂₂ + F₂₃*S₂₃
        R₂₃ = F₂₁*S₁₃ + F₂₂*S₂₃ + F₂₃*S₃₃

        R₃₁ = F₃₁*S₁₁ + F₃₂*S₁₂ + F₃₃*S₁₃
        R₃₂ = F₃₁*S₁₂ + F₃₂*S₂₂ + F₃₃*S₂₃
        R₃₃ = F₃₁*S₁₃ + F₃₂*S₂₃ + F₃₃*S₃₃

        T₁₁ = F₁₁*R₁₁ + F₂₁*R₂₁ + F₃₁*R₃₁
        T₂₂ = F₁₂*R₁₂ + F₂₂*R₂₂ + F₃₂*R₃₂
        T₃₃ = F₁₃*R₁₃ + F₂₃*R₂₃ + F₃₃*R₃₃
        T₁₂ = F₁₁*R₁₂ + F₂₁*R₂₂ + F₃₁*R₃₂
        T₂₁ = F₁₂*R₁₁ + F₂₂*R₂₁ + F₃₂*R₃₁
        T₂₃ = F₁₂*R₁₃ + F₂₂*R₂₃ + F₃₂*R₃₃
        T₃₂ = F₁₃*R₁₂ + F₂₃*R₂₂ + F₃₃*R₃₂
        T₁₃ = F₁₁*R₁₃ + F₂₁*R₂₃ + F₃₁*R₃₃
        T₃₁ = F₁₃*R₁₁ + F₂₃*R₂₁ + F₃₃*R₃₁

        C₁₁ = F₁₁^2 + F₂₁^2 + F₃₁^2
        C₂₂ = F₁₂^2 + F₂₂^2 + F₃₂^2
        C₃₃ = F₁₃^2 + F₂₃^2 + F₃₃^2
        C₁₂ = F₁₁*F₁₂ + F₂₁*F₂₂ + F₃₁*F₃₂
        C₂₃ = F₁₂*F₁₃ + F₂₂*F₂₃ + F₃₂*F₃₃
        C₁₃ = F₁₁*F₁₃ + F₂₁*F₂₃ + F₃₁*F₃₃

        for (i, xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            N_i = Nₛ[i]

            f_S[6*I-5] -= a*τ * N_i * T₁₁ * 𝑤
            f_S[6*I-4] -= a*τ * N_i * T₂₂ * 𝑤
            f_S[6*I-3] -= a*τ * N_i * T₃₃ * 𝑤
            f_S[6*I-2] -= a*τ * N_i * (T₁₂ + T₂₁) * 𝑤
            f_S[6*I-1] -= a*τ * N_i * (T₂₃ + T₃₂) * 𝑤
            f_S[6*I]   -= a*τ * N_i * (T₁₃ + T₃₁) * 𝑤

            for (j, xⱼ) in enumerate(𝓒ₛ)
                J = xⱼ.𝐼
                NN = N_i * Nₛ[j]

                k_SS[6*I-5, 6*J-5] += a*τ * NN * C₁₁ * 𝑤
                k_SS[6*I-5, 6*J-2] += a*τ * NN * C₁₂ * 𝑤
                k_SS[6*I-5, 6*J]   += a*τ * NN * C₁₃ * 𝑤

                k_SS[6*I-4, 6*J-4] += a*τ * NN * C₂₂ * 𝑤
                k_SS[6*I-4, 6*J-2] += a*τ * NN * C₁₂ * 𝑤
                k_SS[6*I-4, 6*J-1] += a*τ * NN * C₂₃ * 𝑤

                k_SS[6*I-3, 6*J-3] += a*τ * NN * C₃₃ * 𝑤
                k_SS[6*I-3, 6*J-1] += a*τ * NN * C₂₃ * 𝑤
                k_SS[6*I-3, 6*J]   += a*τ * NN * C₁₃ * 𝑤

                k_SS[6*I-2, 6*J-5] += a*τ * NN * C₁₂ * 𝑤
                k_SS[6*I-2, 6*J-4] += a*τ * NN * C₁₂ * 𝑤
                k_SS[6*I-2, 6*J-2] += a*τ * NN * (C₁₁ + C₂₂) * 𝑤
                k_SS[6*I-2, 6*J-1] += a*τ * NN * C₁₃ * 𝑤
                k_SS[6*I-2, 6*J]   += a*τ * NN * C₂₃ * 𝑤

                k_SS[6*I-1, 6*J-4] += a*τ * NN * C₂₃ * 𝑤
                k_SS[6*I-1, 6*J-3] += a*τ * NN * C₂₃ * 𝑤
                k_SS[6*I-1, 6*J-2] += a*τ * NN * C₁₃ * 𝑤
                k_SS[6*I-1, 6*J-1] += a*τ * NN * (C₂₂ + C₃₃) * 𝑤
                k_SS[6*I-1, 6*J]   += a*τ * NN * C₁₂ * 𝑤

                k_SS[6*I,   6*J-5] += a*τ * NN * C₁₃ * 𝑤
                k_SS[6*I,   6*J-3] += a*τ * NN * C₁₃ * 𝑤
                k_SS[6*I,   6*J-2] += a*τ * NN * C₂₃ * 𝑤
                k_SS[6*I,   6*J-1] += a*τ * NN * C₁₂ * 𝑤
                k_SS[6*I,   6*J]   += a*τ * NN * (C₁₁ + C₃₃) * 𝑤
            end
        end

        for (i, xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
            V1_i = B₁[i]*S₁₁ + B₂[i]*S₁₂ + B₃[i]*S₁₃
            V2_i = B₁[i]*S₁₂ + B₂[i]*S₂₂ + B₃[i]*S₂₃
            V3_i = B₁[i]*S₁₃ + B₂[i]*S₂₃ + B₃[i]*S₃₃

            f_u[3*I-2] -= a*τ * (V1_i*R₁₁ + V2_i*R₁₂ + V3_i*R₁₃) * 𝑤
            f_u[3*I-1] -= a*τ * (V1_i*R₂₁ + V2_i*R₂₂ + V3_i*R₂₃) * 𝑤
            f_u[3*I]   -= a*τ * (V1_i*R₃₁ + V2_i*R₃₂ + V3_i*R₃₃) * 𝑤

            for (j, xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                V1_j = B₁[j]*S₁₁ + B₂[j]*S₁₂ + B₃[j]*S₁₃
                V2_j = B₁[j]*S₁₂ + B₂[j]*S₂₂ + B₃[j]*S₂₃
                V3_j = B₁[j]*S₁₃ + B₂[j]*S₂₃ + B₃[j]*S₃₃

                uu_val = a*τ * (V1_i*V1_j + V2_i*V2_j + V3_i*V3_j) * 𝑤

                k_uu[3*I-2, 3*J-2] += uu_val
                k_uu[3*I-1, 3*J-1] += uu_val
                k_uu[3*I,   3*J]   += uu_val

            end

            for (j, xⱼ) in enumerate(𝓒ₛ)
                J = xⱼ.𝐼
                N_j = Nₛ[j]

                mat_x11 = V1_i * F₁₁ * N_j; geo_x11 = B₁[i] * N_j * R₁₁
                mat_x22 = V2_i * F₁₂ * N_j; geo_x22 = B₂[i] * N_j * R₁₂
                mat_x33 = V3_i * F₁₃ * N_j; geo_x33 = B₃[i] * N_j * R₁₃
                mat_x12 = (V2_i * F₁₁ + V1_i * F₁₂) * N_j; geo_x12 = (B₁[i] * R₁₂ + B₂[i] * R₁₁) * N_j
                mat_x23 = (V3_i * F₁₂ + V2_i * F₁₃) * N_j; geo_x23 = (B₂[i] * R₁₃ + B₃[i] * R₁₂) * N_j
                mat_x13 = (V3_i * F₁₁ + V1_i * F₁₃) * N_j; geo_x13 = (B₁[i] * R₁₃ + B₃[i] * R₁₁) * N_j

                k_Su[6*J-5, 3*I-2] += a*τ * (mat_x11 + geo_x11) * 𝑤
                k_Su[6*J-4, 3*I-2] += a*τ * (mat_x22 + geo_x22) * 𝑤
                k_Su[6*J-3, 3*I-2] += a*τ * (mat_x33 + geo_x33) * 𝑤
                k_Su[6*J-2, 3*I-2] += a*τ * (mat_x12 + geo_x12) * 𝑤
                k_Su[6*J-1, 3*I-2] += a*τ * (mat_x23 + geo_x23) * 𝑤
                k_Su[6*J,   3*I-2] += a*τ * (mat_x13 + geo_x13) * 𝑤

                mat_y11 = V1_i * F₂₁ * N_j; geo_y11 = B₁[i] * N_j * R₂₁
                mat_y22 = V2_i * F₂₂ * N_j; geo_y22 = B₂[i] * N_j * R₂₂
                mat_y33 = V3_i * F₂₃ * N_j; geo_y33 = B₃[i] * N_j * R₂₃
                mat_y12 = (V2_i * F₂₁ + V1_i * F₂₂) * N_j; geo_y12 = (B₁[i] * R₂₂ + B₂[i] * R₂₁) * N_j
                mat_y23 = (V3_i * F₂₂ + V2_i * F₂₃) * N_j; geo_y23 = (B₂[i] * R₂₃ + B₃[i] * R₂₂) * N_j
                mat_y13 = (V3_i * F₂₁ + V1_i * F₂₃) * N_j; geo_y13 = (B₁[i] * R₂₃ + B₃[i] * R₂₁) * N_j

                k_Su[6*J-5, 3*I-1] += a*τ * (mat_y11 + geo_y11) * 𝑤
                k_Su[6*J-4, 3*I-1] += a*τ * (mat_y22 + geo_y22) * 𝑤
                k_Su[6*J-3, 3*I-1] += a*τ * (mat_y33 + geo_y33) * 𝑤
                k_Su[6*J-2, 3*I-1] += a*τ * (mat_y12 + geo_y12) * 𝑤
                k_Su[6*J-1, 3*I-1] += a*τ * (mat_y23 + geo_y23) * 𝑤
                k_Su[6*J,   3*I-1] += a*τ * (mat_y13 + geo_y13) * 𝑤

                mat_z11 = V1_i * F₃₁ * N_j; geo_z11 = B₁[i] * N_j * R₃₁
                mat_z22 = V2_i * F₃₂ * N_j; geo_z22 = B₂[i] * N_j * R₃₂
                mat_z33 = V3_i * F₃₃ * N_j; geo_z33 = B₃[i] * N_j * R₃₃
                mat_z12 = (V2_i * F₃₁ + V1_i * F₃₂) * N_j; geo_z12 = (B₁[i] * R₃₂ + B₂[i] * R₃₁) * N_j
                mat_z23 = (V3_i * F₃₂ + V2_i * F₃₃) * N_j; geo_z23 = (B₂[i] * R₃₃ + B₃[i] * R₃₂) * N_j
                mat_z13 = (V3_i * F₃₁ + V1_i * F₃₃) * N_j; geo_z13 = (B₁[i] * R₃₃ + B₃[i] * R₃₁) * N_j

                k_Su[6*J-5, 3*I]   += a*τ * (mat_z11 + geo_z11) * 𝑤
                k_Su[6*J-4, 3*I]   += a*τ * (mat_z22 + geo_z22) * 𝑤
                k_Su[6*J-3, 3*I]   += a*τ * (mat_z33 + geo_z33) * 𝑤
                k_Su[6*J-2, 3*I]   += a*τ * (mat_z12 + geo_z12) * 𝑤
                k_Su[6*J-1, 3*I]   += a*τ * (mat_z23 + geo_z23) * 𝑤
                k_Su[6*J,   3*I]   += a*τ * (mat_z13 + geo_z13) * 𝑤
            end
        end
    end
end

function gp_key(ξ; ndigits=12)
    return (
        round(ξ.x, digits=ndigits),
        round(ξ.y, digits=ndigits),
        round(ξ.z, digits=ndigits),
    )
end

function build_gp_map(𝓖; ndigits=12)
    m = Dict{Tuple{Float64,Float64,Float64}, Int}()
    for (q, ξ) in enumerate(𝓖)
        key = gp_key(ξ, ndigits=ndigits)
        m[key] = q
    end
    return m
end
function eval_piecewise_3d!(N::AbstractVector, x, y, z)
    np = length(N)

    if np == 4

        N[1] = 1.0
        N[2] = x
        N[3] = y
        N[4] = z

    elseif np == 10

        N[1]  = 1.0
        N[2]  = x
        N[3]  = y
        N[4]  = z
        N[5]  = x^2
        N[6]  = x*y
        N[7]  = x*z
        N[8]  = y^2
        N[9]  = y*z
        N[10] = z^2

    else
        error("eval_piecewise_3d!: unsupported number of basis functions np = $np")
    end
end

end
