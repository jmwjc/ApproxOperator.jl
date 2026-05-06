module Timoshenko

using ..ApproxOperator: AbstractElement

function ∫κκdΩ(ap::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = ap.E
    ν = ap.ν
    h = ap.h
    Dᵢᵢᵢᵢ = E*h^3/12/(1-ν^2)
    Dᵢᵢⱼⱼ = E*ν*h^3/12/(1-ν^2)
    Dᵢⱼᵢⱼ = E*h^3/24/(1+ν)
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += (Dᵢᵢᵢᵢ*B₁[i]*B₁[j] + Dᵢⱼᵢⱼ*B₂[i]*B₂[j])*𝑤
                k[2*I-1,2*J]   += (Dᵢᵢⱼⱼ*B₁[i]*B₂[j] + Dᵢⱼᵢⱼ*B₂[i]*B₁[j])*𝑤
                k[2*I,2*J-1]   += (Dᵢᵢⱼⱼ*B₂[i]*B₁[j] + Dᵢⱼᵢⱼ*B₁[i]*B₂[j])*𝑤
                k[2*I,2*J]     += (Dᵢᵢᵢᵢ*B₂[i]*B₂[j] + Dᵢⱼᵢⱼ*B₁[i]*B₁[j])*𝑤
            end
        end
    end
end

function ∫wwdΩ(ap::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    𝑤 = ap.𝑤
    E = ap.E
    ν = ap.ν
    h = ap.h
    D = 5/6*h*E/2/(1+ν)
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += D*(B₁[i]*B₁[j] + B₂[i]*B₂[j])*𝑤
            end
        end
    end
end

function ∫φφdΩ(ap::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    𝑤 = ap.𝑤
    E = ap.E
    ν = ap.ν
    h = ap.h
    D = 5/6*h*E/2/(1+ν)
    for ξ in 𝓖
        N = ξ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += D*N[i]*N[j]*𝑤
                k[2*I,2*J]     += D*N[i]*N[j]*𝑤
            end
        end
    end
end

function ∫φwdΩ(ap::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    𝑤 = ap.𝑤
    E = ap.E
    ν = ap.ν
    h = ap.h
    D = 5/6*h*E/2/(1+ν)
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        N = ξ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,J] -= D*N[i]*B₁[j]*𝑤
                k[2*I,J]   -= D*N[i]*B₂[j]*𝑤
            end
        end
    end
end

function ∫φwdΩ(a₁::T,a₂::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒₁ = a₁.𝓒; 𝓖₁ = a₁.𝓖
    𝓒₂ = a₂.𝓒; 𝓖₂ = a₂.𝓖
    𝑤 = a₁.𝑤
    E = a₁.E
    ν = a₁.ν
    h = a₁.h
    D = 5/6*h*E/2/(1+ν)
    for (ξ₁,ξ₂) in zip(𝓖₁,𝓖₂)
        B₁ = ξ₁[:∂𝝭∂x]
        B₂ = ξ₁[:∂𝝭∂y]
        N = ξ₂[:𝝭]
        for (i,xᵢ) in enumerate(𝓒₁)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒₂)
                J = xⱼ.𝐼
                k[2*I-1,J] -= D*N[i]*B₁[j]*𝑤
                k[2*I,J]   -= D*N[i]*B₂[j]*𝑤
            end
        end
    end
end

function ∫wwGdΩ(ap::T, k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        B = ξ[:∂𝝭∂x]
        𝑤 = ξ.𝑤
        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j, xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I, J] += B[i] * B[j] * 𝑤
            end
        end
    end
end

function ∫wwGdΩ2D(ap::T, k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        σ₁₁ = ξ.σ₁₁
        σ₂₂ = ξ.σ₂₂
        σ₁₂ = ξ.σ₁₂
        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j, xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3 * I - 2, 3 * J - 2] += (
                    σ₁₁ * B₁[i] * B₁[j] +
                    σ₂₂ * B₂[i] * B₂[j] +
                    σ₁₂ * (B₁[i] * B₂[j] + B₂[i] * B₁[j])
                ) * 𝑤
            end
        end
    end
end

# SS exact displacement: w = w_b + w_s (報告 2.2.4(b), SS)
function w_exact_ss(x::Float64, E::Float64, I::Float64, κ::Float64, G::Float64, A::Float64, L::Float64, q::Float64)
    w_b = q * x * (L^3 - 2.0 * L * x^2 + x^3) / (24.0 * E * I)
    w_s = q * x * (L - x) / (2.0 * κ * G * A)
    return w_b + w_s
end

# CF exact displacement: w = w_b + w_s (報告 2.2.4(b), CF)
function w_exact_cf(x::Float64, E::Float64, I::Float64, κ::Float64, G::Float64, A::Float64, L::Float64, q::Float64)
    w_b = q * x^2 * (6.0 * L^2 - 4.0 * L * x + x^2) / (24.0 * E * I)
    w_s = q * x * (2.0 * L - x) / (2.0 * κ * G * A)
    return w_b + w_s
end

# SS exact rotation φ(x): 對應 SS 位移解的轉角閉式
function φ_exact_ss(x::Float64, E::Float64, I::Float64, L::Float64, q::Float64)
    return q * (L^3 - 6.0 * L * x^2 + 4.0 * x^3) / (24.0 * E * I)
end

# CF exact rotation φ(x): 對應 CF 位移解的轉角閉式
function φ_exact_cf(x::Float64, E::Float64, I::Float64, L::Float64, q::Float64)
    return q * x * (3.0 * L^2 - 3.0 * L * x + x^2) / (6.0 * E * I)
end

end
