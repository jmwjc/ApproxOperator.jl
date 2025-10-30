module Hamilton
    
using ..ApproxOperator: AbstractElement

function ∫∫q̇mṗqkpdxdt(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        m = ξ.m
        kᶜ = ξ.k
        B = ξ[:∂𝝭∂x]
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (B[i]*m*B[j] - N[i]*kᶜ*N[j])*𝑤
            end
        end
    end
end

function ∫∫∇q∇pdxdt(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        ρA = ξ.ρA
        EA = ξ.EA
        Bₓ = ξ[:∂𝝭∂x]
        Bₜ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (Bₜ[i]*ρA*Bₜ[j] - Bₓ[i]*EA*Bₓ[j])*𝑤
            end
        end
    end
end

function ∫∫∇q∇pdxdt(a₁::T,a₂::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒₁ = a₁.𝓒; 𝓖₁ = a₁.𝓖
    𝓒₂ = a₂.𝓒; 𝓖₂ = a₂.𝓖
    for (ξ₁,ξ₂) in zip(𝓖₁,𝓖₂)
        B̄ₓ = ξ₁[:∂𝝭∂x]
        B̄ₜ = ξ₁[:∂𝝭∂y]
        ρA = ξ₂.ρA
        EA = ξ₂.EA
        Bₓ = ξ₂[:∂𝝭∂x]
        Bₜ = ξ₂[:∂𝝭∂y]
        𝑤 = ξ₂.𝑤
        for (i,xᵢ) in enumerate(𝓒₁)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒₂)
                J = xⱼ.𝐼
                k[I,J] += (-B̄ₜ[i]*ρA*Bₜ[j] + B̄ₓ[i]*EA*Bₓ[j])*𝑤
            end
        end
    end
end

function ∫q∇𝑛pds(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        Bₓ = ξ[:∂𝝭∂x]
        Bₜ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        nₓ = ξ.n₁
        nₜ = ξ.n₂
        ρA = ξ.ρA
        EA = ξ.EA
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] -= (-N[i]*ρA*Bₜ[j]*nₜ + N[i]*EA*Bₓ[j]*nₓ)*𝑤
            end
            # f[I] -= (-ρA*Bₜ[i]*nₜ + EA*Bₓ[i]*nₓ)*g*𝑤
        end
    end
end

function ∫∫αqṗdxdt(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    α = ap.α
    for ξ in 𝓖
        B = ξ[:∂𝝭∂y]
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += α*(B[i]*N[j] + N[i]*B[j])*𝑤
            end
        end
    end
end


function stabilization_bar_LSG(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        ρA = ξ.ρA
        EA = ξ.EA
        α = ξ.α
        Bₓₓ = ξ[:∂²𝝭∂x²]
        Bₜₜ = ξ[:∂²𝝭∂y²]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += α*(ρA*Bₜₜ[i] - EA*Bₓₓ[i])*(ρA*Bₜₜ[j] - EA*Bₓₓ[j])*𝑤
            end
        end
    end
end

function stabilization_bar_LSG_Γ(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        ρA = ξ.ρA
        EA = ξ.EA
        α = ξ.α
        Bₓ = ξ[:∂𝝭∂x]
        Bₜ = ξ[:∂𝝭∂y]
        nₓ = ξ.n₁
        nₜ = ξ.n₂
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += α*(ρA*Bₜ[i]*nₜ - EA*Bₓ[i]*nₓ)*(ρA*Bₜ[j]*nₜ - EA*Bₓ[j]*nₓ)*𝑤
            end
        end
    end
end

function truncation_error(ap::T,fₓ::AbstractVector{Float64},fₜ::AbstractVector{Float64},fₓₓ::AbstractVector{Float64},fₜₜ::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        c = ξ.c
        Bₓ = ξ[:∂𝝭∂x]
        Bₜ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        x = ξ.x
        t = ξ.y
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                fₓ[I] += Bₓ[i]*Bₓ[j]*(xⱼ.x-xᵢ.x - c*(xⱼ.y-xᵢ.y))*𝑤
                fₜ[I] += Bₜ[i]*Bₜ[j]*(xⱼ.x-xᵢ.x - c*(xⱼ.y-xᵢ.y))*𝑤
                fₓₓ[I] += Bₓ[i]*Bₓ[j]*(xⱼ.x-xᵢ.x - c*(xⱼ.y-xᵢ.y))^2*𝑤
                fₜₜ[I] += Bₜ[i]*Bₜ[j]*(xⱼ.x-xᵢ.x - c*(xⱼ.y-xᵢ.y))^2*𝑤
            end
        end
    end
end

function truncation_error(aps::Vector{T},nₚ::Int) where T<:AbstractElement
    fₓ = zeros(nₚ)
    fₜ = zeros(nₚ)
    fₓₓ = zeros(nₚ)
    fₜₜ = zeros(nₚ)
    for ap in aps
        truncation_error(ap,fₓ,fₜ,fₓₓ,fₜₜ)
    end
    return fₓ,fₜ,fₓₓ,fₜₜ
end

function ∫pudΩ(a₁::T,a₂::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒₁ = a₁.𝓒; 𝓖₁ = a₁.𝓖
    𝓒₂ = a₂.𝓒; 𝓖₂ = a₂.𝓖
    for (ξ₁,ξ₂) in zip(𝓖₁,𝓖₂)
        Bₜ = ξ₂[:∂𝝭∂y]
        N = ξ₁[:𝝭]
        𝑤 = ξ₂.𝑤
        for (i,xᵢ) in enumerate(𝓒₁)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒₂)
                J = xⱼ.𝐼
                k[I,J] += N[i]*Bₜ[j]*𝑤
            end
        end
    end
end

function ∫ppdΩ(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        ρA = ξ.ρA
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += -(1/ρA)N[i]*N[j]*𝑤
            end
        end
    end
end

function ∫uudΩ(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        EA = ξ.EA
        Bₓ = ξ[:∂𝝭∂x]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += -Bₓ[i]*EA*Bₓ[j]*𝑤
            end
        end
    end
end


end