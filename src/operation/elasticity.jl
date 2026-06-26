module Elasticity
    
using ..ApproxOperator: AbstractElement
function ∫εᵢⱼσᵢⱼdΩ(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z]
        𝑤 = ξ.𝑤
        E = ξ.E
        ν = ξ.ν
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
function ∫∫ρvᵢuᵢdxdy(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        ρ = ξ.ρ
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

function ∫∫εᵢⱼσᵢⱼdxdy(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        E = ξ.E
        ν = ξ.ν
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

function ∫∫εᵢⱼσᵢⱼdxdy_PlaneStrian(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        E = ξ.E
        ν = ξ.ν
        Cᵢᵢᵢᵢ = E*(1-ν)/(1-2*ν)/(1+ν)
        Cᵢᵢⱼⱼ = E*ν/(1-2*ν)/(1+ν)
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

function ∫∫εᵢⱼσᵢⱼdxdy(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        Bᵤ₁ = ξᵤ[:∂𝝭∂x]
        Bᵤ₂ = ξᵤ[:∂𝝭∂y]
        Bₛ₁ = ξₛ[:∂𝝭∂x]
        Bₛ₂ = ξₛ[:∂𝝭∂y]
        𝑤 = ξᵤ.𝑤
        E = ξᵤ.E
        ν = ξᵤ.ν
        Cᵢᵢᵢᵢ = E/(1-ν^2)
        Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += ( Cᵢᵢᵢᵢ*Bₛ₁[i]*Bᵤ₁[j]+Cᵢⱼᵢⱼ*Bₛ₂[i]*Bᵤ₂[j])*𝑤
                k[2*I-1,2*J]   += (Cᵢᵢⱼⱼ*Bₛ₁[i]*Bᵤ₂[j]+Cᵢⱼᵢⱼ*Bₛ₂[i]*Bᵤ₁[j])*𝑤
                k[2*I,2*J-1]   += (Cᵢᵢⱼⱼ*Bₛ₂[i]*Bᵤ₁[j]+Cᵢⱼᵢⱼ*Bₛ₁[i]*Bᵤ₂[j])*𝑤
                k[2*I,2*J]     += (Cᵢᵢᵢᵢ*Bₛ₂[i]*Bᵤ₂[j]+Cᵢⱼᵢⱼ*Bₛ₁[i]*Bᵤ₁[j])*𝑤
            end
        end
    end
end

function ∫∫εᵢⱼσᵢⱼdxdy_HR(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        Bᵤ₁ = ξᵤ[:∂𝝭∂x]
        Bᵤ₂ = ξᵤ[:∂𝝭∂y]
        N = ξₛ[:𝝭]
        𝑤 = ξᵤ.𝑤
       
     
        for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
            

                k[3*I-2,2*J-1] -= Bᵤ₁[j]*N[i]*𝑤
                k[3*I-1,2*J]   -= Bᵤ₂[j]*N[i]*𝑤
                k[3*I,2*J-1]   -= Bᵤ₂[j]*N[i]*𝑤
                k[3*I,2*J]     -= Bᵤ₁[j]*N[i]*𝑤
            end
        end
    end
end
function ∫εᵢⱼσᵢⱼdΩ(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        B₃ = ξᵤ[:∂𝝭∂z]
        N = ξₛ[:𝝭]
       
        𝑤 = ξᵤ.𝑤
        E = ξᵤ.E
        ν = ξᵤ.ν
        Cᵢᵢᵢᵢ = E*(1-ν)/(1+ν)/(1-2ν)
        Cᵢᵢⱼⱼ = E*ν/(1+ν)/(1-2ν)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
       
        for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
              

                k[6*I-5,3*J-2] -=  N[i]*B₁[j]*𝑤
                k[6*I-4,3*J-1] -=  N[i]*B₂[j]*𝑤
                k[6*I-3,3*J]   -=  N[i]*B₃[j]*𝑤
                k[6*I-2,3*J-2] -=  N[i]*B₂[j]*𝑤
                k[6*I-2,3*J-1] -=  N[i]*B₁[j]*𝑤
                k[6*I-1,3*J-1] -=  N[i]*B₃[j]*𝑤
                k[6*I-1,3*J]   -=  N[i]*B₂[j]*𝑤
                k[6*I,3*J-2]   -=  N[i]*B₃[j]*𝑤
                k[6*I,3*J]     -=  N[i]*B₁[j]*𝑤


              
            end
        end
    end
end
function ∫∫εᵢⱼσᵢⱼdxdy_PlaneStrian(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    for (ξᵤ,ξₛ) in zip(𝓖ᵤ,𝓖ₛ)
        Bᵤ₁ = ξᵤ[:∂𝝭∂x]
        Bᵤ₂ = ξᵤ[:∂𝝭∂y]
        Bₛ₁ = ξₛ[:∂𝝭∂x]
        Bₛ₂ = ξₛ[:∂𝝭∂y]
        𝑤 = ξᵤ.𝑤
        E = ξₛ.E
        ν = ξₛ.ν
        Cᵢᵢᵢᵢ = E*(1-ν)/(1-2*ν)/(1+ν)
        Cᵢᵢⱼⱼ = E*ν/(1-2*ν)/(1+ν)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += ( Cᵢᵢᵢᵢ*Bₛ₁[i]*Bᵤ₁[j]+Cᵢⱼᵢⱼ*Bₛ₂[i]*Bᵤ₂[j])*𝑤
                k[2*I-1,2*J]   += (Cᵢᵢⱼⱼ*Bₛ₁[i]*Bᵤ₂[j]+Cᵢⱼᵢⱼ*Bₛ₂[i]*Bᵤ₁[j])*𝑤
                k[2*I,2*J-1]   += (Cᵢᵢⱼⱼ*Bₛ₂[i]*Bᵤ₁[j]+Cᵢⱼᵢⱼ*Bₛ₁[i]*Bᵤ₂[j])*𝑤
                k[2*I,2*J]     += (Cᵢᵢᵢᵢ*Bₛ₂[i]*Bᵤ₂[j]+Cᵢⱼᵢⱼ*Bₛ₁[i]*Bᵤ₁[j])*𝑤
            end
        end
    end
end
function ∫∫εᵛᵢⱼσᵛᵢⱼdxdy(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        E = ξ.E
        ν = ξ.ν
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

function ∫∫εᵛᵢⱼσᵛᵢⱼdxdy_bbar(ap::T, k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    n = length(𝓒)

    B̄₁ = zeros(Float64, n)
    B̄₂ = zeros(Float64, n)
    𝑤sum = 0.0

    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤

        for i in 1:n
            B̄₁[i] += B₁[i] * 𝑤
            B̄₂[i] += B₂[i] * 𝑤
        end
        𝑤sum += 𝑤
    end

    for i in 1:n
        B̄₁[i] /= 𝑤sum
        B̄₂[i] /= 𝑤sum
    end

    for ξ in 𝓖
        𝑤 = ξ.𝑤
        E = ξ.E
        ν = ξ.ν
        Cᵛ = E / (1 - 2*ν)

        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j, xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += Cᵛ/3 * B̄₁[i] * B̄₁[j] * 𝑤
                k[2*I-1,2*J]   += Cᵛ/3 * B̄₁[i] * B̄₂[j] * 𝑤
                k[2*I,2*J-1]   += Cᵛ/3 * B̄₂[i] * B̄₁[j] * 𝑤
                k[2*I,2*J]     += Cᵛ/3 * B̄₂[i] * B̄₂[j] * 𝑤
            end
        end
    end
end

function ∫∫εᵈᵢⱼσᵈᵢⱼdxdy(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        E = ξ.E
        ν = ξ.ν
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
function ∫εᵈᵢⱼσᵈᵢⱼdΩ(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z]
        𝑤 = ξ.𝑤
        E = ξ.E
        ν = ξ.ν
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


function ∫εᵈᵢⱼσᵈᵢⱼdΩ(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    for (ξᵤ,ξₛ) in zip(𝓖ᵤ,𝓖ₛ)
        Bᵤ₁ = ξᵤ[:∂𝝭∂x]
        Bᵤ₂ = ξᵤ[:∂𝝭∂y]
        Bₛ₁ = ξₛ[:∂𝝭∂x]
        Bₛ₂ = ξₛ[:∂𝝭∂y]
        𝑤 = ξᵤ.𝑤
        E = ξₛ.E
        ν = ξₛ.ν
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



function ∫∫εᵈᵢⱼσᵈᵢⱼdxdy(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    for (ξᵤ,ξₛ) in zip(𝓖ᵤ,𝓖ₛ)
        Bᵤ₁ = ξᵤ[:∂𝝭∂x]
        Bᵤ₂ = ξᵤ[:∂𝝭∂y]
        Bₛ₁ = ξₛ[:∂𝝭∂x]
        Bₛ₂ = ξₛ[:∂𝝭∂y]
        𝑤 = ξᵤ.𝑤
        E = ξₛ.E
        ν = ξₛ.ν
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

function ∫∫σᵢⱼσₖₗdxdy(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        E = ξ.E
        ν = ξ.ν
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

function ∫∫uᵢⱼuₖₗdxdy(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
   
       
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += N[i]*N[j]*𝑤
                # k[2*I-1,2*J] += N[i]*N[j]*𝑤
                # k[2*I,2*J-1] += N[i]*N[j]*𝑤
                k[2*I,2*J] += N[i]*N[j]*𝑤
              
            end
        end
    end
end

function ∫∫σᵢⱼσₖₗdΩ(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        E = ξ.E
        ν = ξ.ν
        C⁻¹ᵢᵢᵢᵢ = 1/E
        C⁻¹ᵢᵢⱼⱼ = -ν/E
        C⁻¹ᵢⱼᵢⱼ = 2*(1+ν)/E
        
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                

                k[6*I-5,6*J-5] += N[i]*C⁻¹ᵢᵢᵢᵢ*N[j]*𝑤
                k[6*I-5,6*J-4] += N[i]*C⁻¹ᵢᵢⱼⱼ*N[j]*𝑤
                k[6*I-5,6*J-3] += N[i]*C⁻¹ᵢᵢⱼⱼ*N[j]*𝑤

                k[6*I-4,6*J-5] += N[i]*C⁻¹ᵢᵢⱼⱼ*N[j]*𝑤
                k[6*I-4,6*J-4] += N[i]*C⁻¹ᵢᵢᵢᵢ*N[j]*𝑤
                k[6*I-4,6*J-3] += N[i]*C⁻¹ᵢᵢⱼⱼ*N[j]*𝑤

                k[6*I-3,6*J-5] += N[i]*C⁻¹ᵢᵢⱼⱼ*N[j]*𝑤
                k[6*I-3,6*J-4] += N[i]*C⁻¹ᵢᵢⱼⱼ*N[j]*𝑤
                k[6*I-3,6*J-3] += N[i]*C⁻¹ᵢᵢᵢᵢ*N[j]*𝑤

                k[6*I-2,6*J-2] += N[i]*C⁻¹ᵢⱼᵢⱼ*N[j]*𝑤
                k[6*I-1,6*J-1] += N[i]*C⁻¹ᵢⱼᵢⱼ*N[j]*𝑤
                k[6*I,6*J]     += N[i]*C⁻¹ᵢⱼᵢⱼ*N[j]*𝑤
            end
        end
    end
end

# function ∫∫σᵢⱼσₖₗdxdy(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
#     𝓒 = ap.𝓒;𝓖 = ap.𝓖
#     for ξ in 𝓖
#         N = ξ[:𝝭]
#         𝑤 = ξ.𝑤
#         E = ξ.E
#         ν = ξ.ν
#         C⁻¹ᵢᵢᵢᵢ = 1/E
#         C⁻¹ᵢᵢⱼⱼ = -ν/E
#         C⁻¹ᵢⱼᵢⱼ = 2*(1+ν)/E
#         # C⁻¹ᵢⱼᵢⱼ = (1+ν)/E
#         for (i,xᵢ) in enumerate(𝓒)
#             I = xᵢ.𝐼
#             for (j,xⱼ) in enumerate(𝓒)
#                 J = xⱼ.𝐼
#                 k[3*I-2,3*J-2] += N[i]*C⁻¹ᵢᵢᵢᵢ*N[j]*𝑤
#                 k[3*I-2,3*J-1] += N[i]*C⁻¹ᵢᵢⱼⱼ*N[j]*𝑤
#                 k[3*I-1,3*J-2] += N[i]*C⁻¹ᵢᵢⱼⱼ*N[j]*𝑤
#                 k[3*I-1,3*J-1] += N[i]*C⁻¹ᵢᵢᵢᵢ*N[j]*𝑤
#                 k[3*I,3*J]     += N[i]*C⁻¹ᵢⱼᵢⱼ*N[j]*𝑤
#             end
#         end
#     end
# end


function ∫∫Mσσdxdy(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        E = ξ.E
        ν = ξ.ν
        C⁻¹ᵢᵢᵢᵢ = 1/E
        C⁻¹ᵢᵢⱼⱼ = -ν/E
        C⁻¹ᵢⱼᵢⱼ = 2*(1+ν)/E
        # C⁻¹ᵢⱼᵢⱼ = (1+ν)/E
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
function ∫∫εᵢⱼCᵢⱼₖₗEₖₗdxdy(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        N = ξₛ[:𝝭]
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        E = ξₛ.E
        ν = ξₛ.ν
        Cᵢᵢᵢᵢ = E/(1-ν^2)
        Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                
                k[3*I-2,2*J-1] -= N[i]*Cᵢᵢᵢᵢ*B₁[j]*𝑤
                k[3*I-2,2*J]   -= N[i]*Cᵢᵢⱼⱼ*B₂[j]*𝑤
                k[3*I-1,2*J-1] -= N[i]*Cᵢᵢⱼⱼ*B₁[j]*𝑤
                k[3*I-1,2*J]   -= N[i]*Cᵢᵢᵢᵢ*B₂[j]*𝑤
                k[3*I,2*J-1]   -= N[i]*Cᵢⱼᵢⱼ*B₂[j]*𝑤
                k[3*I,2*J]     -= N[i]*Cᵢⱼᵢⱼ*B₁[j]*𝑤
            end
        end
    end
end

function ∫∫εᵢⱼCᵢⱼₖₗεₖₗdxdy(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        E = ξ.E
        ν = ξ.ν
        Cᵢᵢᵢᵢ = E/(1-ν^2)
        Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += N[i]*Cᵢᵢᵢᵢ*N[j]*𝑤
                k[3*I-2,3*J-1] += N[i]*Cᵢᵢⱼⱼ*N[j]*𝑤
                k[3*I-1,3*J-2] += N[i]*Cᵢᵢⱼⱼ*N[j]*𝑤
                k[3*I-1,3*J-1] += N[i]*Cᵢᵢᵢᵢ*N[j]*𝑤
                k[3*I,3*J]     += N[i]*Cᵢⱼᵢⱼ*N[j]*𝑤
            end
        end
    end
end

function ∫∫εᵈᵢⱼCᵢⱼₖₗεᵈₖₗdxdy(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        Ē = ξ.Ē
        ν̄  = ξ.ν̄ 
       
        Cᵈ = Ē/(1+ν̄ )
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += 2/3*N[i]*Cᵈ*N[j]*𝑤
                k[3*I-2,3*J-1] += -1/3*N[i]*Cᵈ*N[j]*𝑤
                
                k[3*I-1,3*J-2] += -1/3*N[i]*Cᵈ*N[j]*𝑤
                k[3*I-1,3*J-1] += 2/3*N[i]*Cᵈ*N[j]*𝑤
                k[3*I,3*J]     += 1/2*N[i]*Cᵈ*N[j]*𝑤

            end
        end
    end
end

function ∫∫σᵈᵢⱼσᵈₖₗdxdy(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        Ē = ξ.Ē
        ν̄  = ξ.ν̄ 
        Cᵈ = Ē/(1+ν̄ )
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += 2/Cᵈ*N[i]*N[j]*𝑤
                k[3*I-2,3*J-1] += 1/Cᵈ*N[i]*Cᵈ*N[j]*𝑤
                
                k[3*I-1,3*J-2] += 1/Cᵈ*N[i]*Cᵈ*N[j]*𝑤
                k[3*I-1,3*J-1] += 2/Cᵈ*N[i]*Cᵈ*N[j]*𝑤
                k[3*I,3*J]     += 2/Cᵈ*N[i]*Cᵈ*N[j]*𝑤

            end
        end
    end
end
function ∫∫σᵛᵢⱼσᵛₖₗdxdy(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        Ē = ξ.Ē
        ν̄  = ξ.ν̄ 
        Cᵛ = Ē/(1-2*ν̄ )
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += 3/4/Cᵛ*N[i]*N[j]*𝑤
                k[3*I-2,3*J-1] += 3/4/Cᵛ*N[i]*N[j]*𝑤
                k[3*I-1,3*J-2] += 3/4/Cᵛ*N[i]*N[j]*𝑤
                k[3*I-1,3*J-1] += 3/4/Cᵛ*N[i]*N[j]*𝑤
                
            end
        end
    end
end

function ∫∫εᵛᵢⱼCᵢⱼₖₗεᵛₖₗdxdy(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        Ē = ξ.Ē
        ν̄  = ξ.ν̄ 
        
        Cᵛ = Ē/(1-2*ν̄ )
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += 1/3*N[i]*Cᵛ*N[j]*𝑤
                k[3*I-2,3*J-1] += 1/3*N[i]*Cᵛ*N[j]*𝑤
                k[3*I-1,3*J-2] += 1/3*N[i]*Cᵛ*N[j]*𝑤
                k[3*I-1,3*J-1] += 1/3*N[i]*Cᵛ*N[j]*𝑤
                
            end
        end
    end
end


function ∫∫σᵢⱼσₖₗdxdy_PlaneStrian(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        E = ξ.E
        ν = ξ.ν
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

# function ∫σᵢⱼnⱼgᵢdΓ(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
#     𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
#     𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
#     for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
#         # 𝑤 = ξₛ.𝑤
#         𝑤 = ξᵤ.𝑤
    

#         N = ξₛ[:𝝭]
#         N̄ = ξᵤ[:𝝭]
#         n₁ = ξᵤ.n₁
#         n₂ = ξᵤ.n₂
#         n₃ = ξᵤ.n₃
#         n₁₁ = ξᵤ.n₁₁
#         n₁₂ = ξᵤ.n₁₂
#         n₁₃ = ξᵤ.n₁₃
#         n₂₂ = ξᵤ.n₂₂
#         n₂₃ = ξᵤ.n₂₃
#         n₃₃ = ξᵤ.n₃₃
#         g₁ = ξᵤ.g₁
#         g₂ = ξᵤ.g₂
#         g₃ = ξᵤ.g₃
#     #     @assert abs(ξₛ.x - ξᵤ.x) < 1e-12 "boundary Gauss x mismatch : ξₛ.x=$(ξₛ.x), ξᵤ.x=$(ξᵤ.x)"
#     # @assert abs(ξₛ.y - ξᵤ.y) < 1e-12 "boundary Gauss y mismatch : ξₛ.y=$(ξₛ.y), ξᵤ.y=$(ξᵤ.y)"
#     # @assert abs(ξₛ.z - ξᵤ.z) < 1e-12 "boundary Gauss z mismatch at : ξₛ.z=$(ξₛ.z), ξᵤ.z=$(ξᵤ.z)"
#     # @assert abs(ξₛ.𝑤 - ξᵤ.𝑤) < 1e-12 "boundary Gauss weight mismatch at : ξₛ.w=$(ξₛ.𝑤), ξᵤ.w=$(ξᵤ.𝑤)"


#         for (i,xᵢ) in enumerate(𝓒ₛ)
#             I = xᵢ.𝐼
#             for (j,xⱼ) in enumerate(𝓒ᵤ)
#                 J = xⱼ.𝐼
            
#                 k[6*I-5,3*J-2] += N[i]*n₁*n₁₁*N̄[j]*𝑤
#                 k[6*I-5,3*J-1] += N[i]*n₁*n₁₂*N̄[j]*𝑤
#                 k[6*I-5,3*J]   += N[i]*n₁*n₁₃*N̄[j]*𝑤

#                 k[6*I-4,3*J-2] += N[i]*n₂*n₁₂*N̄[j]*𝑤
#                 k[6*I-4,3*J-1] += N[i]*n₂*n₂₂*N̄[j]*𝑤
#                 k[6*I-4,3*J]   += N[i]*n₂*n₂₃*N̄[j]*𝑤

#                 k[6*I-3,3*J-2] += N[i]*n₃*n₁₃*N̄[j]*𝑤
#                 k[6*I-3,3*J-1] += N[i]*n₃*n₂₃*N̄[j]*𝑤
#                 k[6*I-3,3*J]   += N[i]*n₃*n₃₃*N̄[j]*𝑤

#                 k[6*I-2,3*J-2] += N[i]*(n₁*n₁₂ + n₂*n₁₁)*N̄[j]*𝑤
#                 k[6*I-2,3*J-1] += N[i]*(n₁*n₂₂ + n₂*n₁₂)*N̄[j]*𝑤
#                 k[6*I-2,3*J]   += N[i]*(n₁*n₂₃ + n₂*n₁₃)*N̄[j]*𝑤

#                 k[6*I-1,3*J-2] += N[i]*(n₂*n₁₃ + n₃*n₁₂)*N̄[j]*𝑤
#                 k[6*I-1,3*J-1] += N[i]*(n₂*n₂₃ + n₃*n₂₂)*N̄[j]*𝑤
#                 k[6*I-1,3*J]   += N[i]*(n₂*n₃₃ + n₃*n₂₃)*N̄[j]*𝑤

#                 k[6*I,3*J-2] += N[i]*(n₁*n₁₃ + n₃*n₁₁)*N̄[j]*𝑤
#                 k[6*I,3*J-1] += N[i]*(n₁*n₂₃ + n₃*n₁₂)*N̄[j]*𝑤
#                 k[6*I,3*J]   += N[i]*(n₁*n₃₃ + n₃*n₁₃)*N̄[j]*𝑤

#             end
#             f[6*I-5] += N[i]*(n₁*n₁₁*g₁ + n₁*n₁₂*g₂ + n₁*n₁₃*g₃)*𝑤
#             f[6*I-4] += N[i]*(n₂*n₁₂*g₁ + n₂*n₂₂*g₂ + n₂*n₂₃*g₃)*𝑤
#             f[6*I-3] += N[i]*(n₃*n₁₃*g₁ + n₃*n₂₃*g₂ + n₃*n₃₃*g₃)*𝑤
#             f[6*I-2] += N[i]*((n₁*n₁₂+n₂*n₁₁)*g₁ + (n₁*n₂₂+n₂*n₁₂)*g₂ + (n₁*n₂₃+n₂*n₁₃)*g₃)*𝑤 
#             f[6*I-1] += N[i]*((n₃*n₁₂+n₂*n₁₃)*g₁ + (n₂*n₂₃+n₃*n₂₂)*g₂ + (n₂*n₃₃+n₃*n₂₃)*g₃)*𝑤 
#             f[6*I]   += N[i]*((n₁*n₁₃+n₃*n₁₁)*g₁ + (n₁*n₂₃+n₃*n₁₂)*g₂ + (n₁*n₃₃+n₃*n₁₃)*g₃)*𝑤

          
#         end
#     end
# end

# function ∫σᵢⱼnⱼgᵢdΓ(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
#     𝓒ₛ = aₛ.𝓒
#     𝓖ₛ = aₛ.𝓖
#     𝓒ᵤ = aᵤ.𝓒
#     𝓖ᵤ = aᵤ.𝓖

#     gpmap = build_gp_map(𝓖ᵤ, ndigits=12)

#     for ξₛ in 𝓖ₛ
#         key = gp_key(ξₛ, ndigits=12)
#         @assert haskey(gpmap, key) "No matching Gauss point found for key=$key"

#         ξᵤ = 𝓖ᵤ[gpmap[key]]

#         𝑤 = ξᵤ.𝑤

#         N  = ξₛ[:𝝭]
#         N̄ = ξᵤ[:𝝭]

#         n₁ = ξᵤ.n₁
#         n₂ = ξᵤ.n₂
#         n₃ = ξᵤ.n₃

#         n₁₁ = ξᵤ.n₁₁
#         n₁₂ = ξᵤ.n₁₂
#         n₁₃ = ξᵤ.n₁₃
#         n₂₂ = ξᵤ.n₂₂
#         n₂₃ = ξᵤ.n₂₃
#         n₃₃ = ξᵤ.n₃₃

#         g₁ = ξᵤ.g₁
#         g₂ = ξᵤ.g₂
#         g₃ = ξᵤ.g₃

#         for (i,xᵢ) in enumerate(𝓒ₛ)
#             I = xᵢ.𝐼
#             for (j,xⱼ) in enumerate(𝓒ᵤ)
#                 J = xⱼ.𝐼

#                 k[6*I-5,3*J-2] += N[i]*n₁*n₁₁*N̄[j]*𝑤
#                 k[6*I-5,3*J-1] += N[i]*n₁*n₁₂*N̄[j]*𝑤
#                 k[6*I-5,3*J]   += N[i]*n₁*n₁₃*N̄[j]*𝑤

#                 k[6*I-4,3*J-2] += N[i]*n₂*n₁₂*N̄[j]*𝑤
#                 k[6*I-4,3*J-1] += N[i]*n₂*n₂₂*N̄[j]*𝑤
#                 k[6*I-4,3*J]   += N[i]*n₂*n₂₃*N̄[j]*𝑤

#                 k[6*I-3,3*J-2] += N[i]*n₃*n₁₃*N̄[j]*𝑤
#                 k[6*I-3,3*J-1] += N[i]*n₃*n₂₃*N̄[j]*𝑤
#                 k[6*I-3,3*J]   += N[i]*n₃*n₃₃*N̄[j]*𝑤

#                 k[6*I-2,3*J-2] += N[i]*(n₁*n₁₂ + n₂*n₁₁)*N̄[j]*𝑤
#                 k[6*I-2,3*J-1] += N[i]*(n₁*n₂₂ + n₂*n₁₂)*N̄[j]*𝑤
#                 k[6*I-2,3*J]   += N[i]*(n₁*n₂₃ + n₂*n₁₃)*N̄[j]*𝑤

#                 k[6*I-1,3*J-2] += N[i]*(n₂*n₁₃ + n₃*n₁₂)*N̄[j]*𝑤
#                 k[6*I-1,3*J-1] += N[i]*(n₂*n₂₃ + n₃*n₂₂)*N̄[j]*𝑤
#                 k[6*I-1,3*J]   += N[i]*(n₂*n₃₃ + n₃*n₂₃)*N̄[j]*𝑤

#                 k[6*I,3*J-2] += N[i]*(n₁*n₁₃ + n₃*n₁₁)*N̄[j]*𝑤
#                 k[6*I,3*J-1] += N[i]*(n₁*n₂₃ + n₃*n₁₂)*N̄[j]*𝑤
#                 k[6*I,3*J]   += N[i]*(n₁*n₃₃ + n₃*n₁₃)*N̄[j]*𝑤
#             end

#             f[6*I-5] += N[i]*(n₁*n₁₁*g₁ + n₁*n₁₂*g₂ + n₁*n₁₃*g₃)*𝑤
#             f[6*I-4] += N[i]*(n₂*n₁₂*g₁ + n₂*n₂₂*g₂ + n₂*n₂₃*g₃)*𝑤
#             f[6*I-3] += N[i]*(n₃*n₁₃*g₁ + n₃*n₂₃*g₂ + n₃*n₃₃*g₃)*𝑤
#             f[6*I-2] += N[i]*((n₁*n₁₂+n₂*n₁₁)*g₁ + (n₁*n₂₂+n₂*n₁₂)*g₂ + (n₁*n₂₃+n₂*n₁₃)*g₃)*𝑤
#             f[6*I-1] += N[i]*((n₃*n₁₂+n₂*n₁₃)*g₁ + (n₂*n₂₃+n₃*n₂₂)*g₂ + (n₂*n₃₃+n₃*n₂₃)*g₃)*𝑤
#             f[6*I]   += N[i]*((n₁*n₁₃+n₃*n₁₁)*g₁ + (n₁*n₂₃+n₃*n₁₂)*g₂ + (n₁*n₃₃+n₃*n₁₃)*g₃)*𝑤
#         end
#     end
# end



function ∫σᵢⱼnⱼgᵢdΓ(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒
    𝓒ᵤ = aᵤ.𝓒
    𝓖ᵤ = aᵤ.𝓖

    N = zeros(length(𝓒ₛ))

    for ξᵤ in 𝓖ᵤ
        𝑤 = ξᵤ.𝑤
        x = ξᵤ.x
        y = ξᵤ.y
        z = ξᵤ.z

         eval_piecewise_3d!(N, x, y, z)

        N̄ = ξᵤ[:𝝭]

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

        for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼

                k[6*I-5,3*J-2] += N[i]*n₁*n₁₁*N̄[j]*𝑤
                k[6*I-5,3*J-1] += N[i]*n₁*n₁₂*N̄[j]*𝑤
                k[6*I-5,3*J]   += N[i]*n₁*n₁₃*N̄[j]*𝑤

                k[6*I-4,3*J-2] += N[i]*n₂*n₁₂*N̄[j]*𝑤
                k[6*I-4,3*J-1] += N[i]*n₂*n₂₂*N̄[j]*𝑤
                k[6*I-4,3*J]   += N[i]*n₂*n₂₃*N̄[j]*𝑤

                k[6*I-3,3*J-2] += N[i]*n₃*n₁₃*N̄[j]*𝑤
                k[6*I-3,3*J-1] += N[i]*n₃*n₂₃*N̄[j]*𝑤
                k[6*I-3,3*J]   += N[i]*n₃*n₃₃*N̄[j]*𝑤

                k[6*I-2,3*J-2] += N[i]*(n₁*n₁₂ + n₂*n₁₁)*N̄[j]*𝑤
                k[6*I-2,3*J-1] += N[i]*(n₁*n₂₂ + n₂*n₁₂)*N̄[j]*𝑤
                k[6*I-2,3*J]   += N[i]*(n₁*n₂₃ + n₂*n₁₃)*N̄[j]*𝑤

                k[6*I-1,3*J-2] += N[i]*(n₂*n₁₃ + n₃*n₁₂)*N̄[j]*𝑤
                k[6*I-1,3*J-1] += N[i]*(n₂*n₂₃ + n₃*n₂₂)*N̄[j]*𝑤
                k[6*I-1,3*J]   += N[i]*(n₂*n₃₃ + n₃*n₂₃)*N̄[j]*𝑤

                k[6*I,3*J-2] += N[i]*(n₁*n₁₃ + n₃*n₁₁)*N̄[j]*𝑤
                k[6*I,3*J-1] += N[i]*(n₁*n₂₃ + n₃*n₁₂)*N̄[j]*𝑤
                k[6*I,3*J]   += N[i]*(n₁*n₃₃ + n₃*n₁₃)*N̄[j]*𝑤
            end

            f[6*I-5] += N[i]*(n₁*n₁₁*g₁ + n₁*n₁₂*g₂ + n₁*n₁₃*g₃)*𝑤
            f[6*I-4] += N[i]*(n₂*n₁₂*g₁ + n₂*n₂₂*g₂ + n₂*n₂₃*g₃)*𝑤
            f[6*I-3] += N[i]*(n₃*n₁₃*g₁ + n₃*n₂₃*g₂ + n₃*n₃₃*g₃)*𝑤
            f[6*I-2] += N[i]*((n₁*n₁₂+n₂*n₁₁)*g₁ + (n₁*n₂₂+n₂*n₁₂)*g₂ + (n₁*n₂₃+n₂*n₁₃)*g₃)*𝑤
            f[6*I-1] += N[i]*((n₃*n₁₂+n₂*n₁₃)*g₁ + (n₂*n₂₃+n₃*n₂₂)*g₂ + (n₂*n₃₃+n₃*n₂₃)*g₃)*𝑤
            f[6*I]   += N[i]*((n₁*n₁₃+n₃*n₁₁)*g₁ + (n₁*n₂₃+n₃*n₁₂)*g₂ + (n₁*n₃₃+n₃*n₁₃)*g₃)*𝑤
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
        # Linear3D
        N[1] = 1.0
        N[2] = x
        N[3] = y
        N[4] = z

    elseif np == 10
        # Quadratic3D
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

function ∫σᵢⱼnⱼuᵢds(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        # 𝑤 = ξₛ.𝑤
        𝑤 = ξᵤ.𝑤

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



function ∫σᵢⱼnⱼuᵢdΓ(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
  
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
      
        # 𝑤 = ξₛ.𝑤
        𝑤 = ξᵤ.𝑤
    
        N = ξₛ[:𝝭]
        N̄ = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        n₃ = ξᵤ.n₃
     @assert abs(ξₛ.x - ξᵤ.x) < 1e-12 "boundary Gauss x mismatch : ξₛ.x=$(ξₛ.x), ξᵤ.x=$(ξᵤ.x)"
    @assert abs(ξₛ.y - ξᵤ.y) < 1e-12 "boundary Gauss y mismatch : ξₛ.y=$(ξₛ.y), ξᵤ.y=$(ξᵤ.y)"
    @assert abs(ξₛ.z - ξᵤ.z) < 1e-12 "boundary Gauss z mismatch at : ξₛ.z=$(ξₛ.z), ξᵤ.z=$(ξᵤ.z)"
    @assert abs(ξₛ.𝑤 - ξᵤ.𝑤) < 1e-12 "boundary Gauss weight mismatch at : ξₛ.w=$(ξₛ.𝑤), ξᵤ.w=$(ξᵤ.𝑤)"

        
        for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[6*I-5,3*J-2] -= N[i]*n₁*N̄[j]*𝑤
                k[6*I-4,3*J-1] -= N[i]*n₂*N̄[j]*𝑤
                k[6*I-3,3*J]   -= N[i]*n₃*N̄[j]*𝑤

                k[6*I-2,3*J-2] -= N[i]*n₂*N̄[j]*𝑤
                k[6*I-2,3*J-1] -= N[i]*n₁*N̄[j]*𝑤

                k[6*I-1,3*J-1] -= N[i]*n₃*N̄[j]*𝑤
                k[6*I-1,3*J]   -= N[i]*n₂*N̄[j]*𝑤

                k[6*I,3*J-2]   -= N[i]*n₃*N̄[j]*𝑤
                k[6*I,3*J]     -= N[i]*n₁*N̄[j]*𝑤
            
            end

        end
   


    end
end



function ∫∫∇σᵢⱼuᵢdxdy(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        # 𝑤 = ξₛ.𝑤
        𝑤 = ξᵤ.𝑤
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

function ∫∫∇σᵢⱼuᵢdΩ(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        # 𝑤 = ξₛ.𝑤
        𝑤 = ξᵤ.𝑤
        
        B₁ = ξₛ[:∂𝝭∂x]
        B₂ = ξₛ[:∂𝝭∂y]
        B₃ = ξₛ[:∂𝝭∂z]
        N = ξᵤ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
               
                k[6*I-5,3*J-2] += B₁[i]*N[j]*𝑤
                k[6*I-4,3*J-1] += B₂[i]*N[j]*𝑤
                k[6*I-3,3*J]   += B₃[i]*N[j]*𝑤
                k[6*I-2,3*J-2] += B₂[i]*N[j]*𝑤
                k[6*I-2,3*J-1] += B₁[i]*N[j]*𝑤
                k[6*I-1,3*J-1] += B₃[i]*N[j]*𝑤
                k[6*I-1,3*J]   += B₂[i]*N[j]*𝑤
                k[6*I,3*J-2]   += B₃[i]*N[j]*𝑤
                k[6*I,3*J]     += B₁[i]*N[j]*𝑤

               

            end
        end
    end
end
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

function ∫qpdΩ(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        E = ξ.E
        ν = ξ.ν
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

function ∫∫sᵢⱼsᵢⱼdxdy(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        E = ξ.E
        ν = ξ.ν
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

function ∫σᵢⱼnⱼgᵢds(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
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
        E = ξ.E
        ν = ξ.ν
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

function ∫σᵢⱼnⱼgᵢvᵢgᵢds(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
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
        E = ξ.E
        ν = ξ.ν
        α = ξ.α
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

function ∫vᵢgᵢds(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        n₁₁ = ξ.n₁₁
        n₂₂ = ξ.n₂₂
        n₁₂ = ξ.n₁₂
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        α = ξ.α
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

function ∫vᵢgᵢdΓ(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
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
        α = ξ.α
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



function ∫∫τ∇σᵢⱼ∇σᵢₖdxdy(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        τ = ξ.τ
        b₁ = ξ.b₁
        b₂ = ξ.b₂
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        E = ξ.E
        ν = ξ.ν
        C⁻¹ᵢᵢᵢᵢ = 1/E
        C⁻¹ᵢᵢⱼⱼ = -ν/E
        C⁻¹ᵢⱼᵢⱼ = 2*(1+ν)/E
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            # τ = xᵢ.β 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                
                k[3*I-2,3*J-2] += τ*(B₁[i]*B₁[j])*𝑤
                k[3*I-2,3*J]   += τ*B₁[i]*B₂[j]*𝑤
                k[3*I-1,3*J-1] += τ*(B₂[i]*B₂[j])*𝑤
                k[3*I-1,3*J]   += τ*B₂[i]*B₁[j]*𝑤
                k[3*I,3*J-2]   += τ*B₂[i]*B₁[j]*𝑤
                k[3*I,3*J-1]   += τ*B₁[i]*B₂[j]*𝑤
                k[3*I,3*J]     += τ*(B₁[i]*B₁[j] + B₂[i]*B₂[j])*𝑤
              
            end
            f[3*I-2] += τ*(B₁[i]*b₁)*𝑤
            f[3*I-1] += τ*(B₂[i]*b₂ )*𝑤
            f[3*I]   += τ*(B₁[i]*b₂ + B₂[i]*b₁)*𝑤

        end
    end
end



function ∫∫τ∇σᵢⱼ∇σᵢₖbdxdy(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        τ = ξ.τ
        ℎ = ξ.ℎ
        b₁ = ξ.b₁
        b₂ = ξ.b₂
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        E = ξ.E
        ν = ξ.ν
        C⁻¹ᵢᵢᵢᵢ = 1/E
        C⁻¹ᵢᵢⱼⱼ = -ν/E
        C⁻¹ᵢⱼᵢⱼ = 2*(1+ν)/E
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            # τ = xᵢ.β 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                gij = (0.5/ℎ^2*N[i]*N[j]) * 𝑤
                k[3*I-2,3*J-2] += τ*(B₁[i]*B₁[j]+gij)*𝑤
                k[3*I-2,3*J]   += τ*B₁[i]*B₂[j]*𝑤
                k[3*I-1,3*J-1] += τ*(B₂[i]*B₂[j]+gij)*𝑤
                k[3*I-1,3*J]   += τ*B₂[i]*B₁[j]*𝑤
                k[3*I,3*J-2]   += τ*B₂[i]*B₁[j]*𝑤
                k[3*I,3*J-1]   += τ*B₁[i]*B₂[j]*𝑤
                k[3*I,3*J]     += τ*(B₁[i]*B₁[j] + B₂[i]*B₂[j]+gij)*𝑤
              
            end
            f[3*I-2] += τ*(B₁[i]*b₁)*𝑤
            f[3*I-1] += τ*(B₂[i]*b₂ )*𝑤
            f[3*I]   += τ*(B₁[i]*b₂ + B₂[i]*b₁)*𝑤

        end
    end
end


function ∫∫τσᵢⱼσᵢₖdxdy(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖

    
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        τ = ξ.τ
        ℎ = ξ.ℎ 
        b₁ = ξ.b₁     
        b₂ = ξ.b₂     
        N  = ξ[:𝝭]    
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        E = ξ.E
        ν = ξ.ν
        G = E/(2*(1+ν))            
       

        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            #  τ = xᵢ.β 
            #  ℎ= xᵢ.ℎ 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                gij = (0.008*ℎ^2*N[i]*N[j]) * 𝑤
                k[3*I-2,3*J-2] += τ * gij   
                k[3*I-1,3*J-1] += τ * gij  
                k[3*I,3*J] += τ * gij   
            end
        end
    end
end



function ∫∫τ∇σᵢⱼ∇σᵢₖdΩ(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        τ = ξ.τ
        b₁ = ξ.b₁
        b₂ = ξ.b₂
        b₃ = ξ.b₃
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                # σ11
                k[6*I-5,6*J-5] += τ*B₁[i]*B₁[j]*𝑤
                k[6*I-5,6*J-2] += τ*B₁[i]*B₂[j]*𝑤
                k[6*I-5,6*J]   += τ*B₁[i]*B₃[j]*𝑤
                # σ22
                k[6*I-4,6*J-4] += τ*B₂[i]*B₂[j]*𝑤
                k[6*I-4,6*J-2] += τ*B₂[i]*B₁[j]*𝑤
                k[6*I-4,6*J-1] += τ*B₂[i]*B₃[j]*𝑤
                # σ33
                k[6*I-3,6*J-3] += τ*B₃[i]*B₃[j]*𝑤
                k[6*I-3,6*J-1] += τ*B₃[i]*B₂[j]*𝑤
                k[6*I-3,6*J]   += τ*B₃[i]*B₁[j]*𝑤
                # σ12
                k[6*I-2,6*J-5] += τ*B₂[i]*B₁[j]*𝑤
                k[6*I-2,6*J-4] += τ*B₁[i]*B₂[j]*𝑤
                k[6*I-2,6*J-2] += τ*(B₁[i]*B₁[j] + B₂[i]*B₂[j])*𝑤
                k[6*I-2,6*J-1] += τ*B₁[i]*B₃[j]*𝑤
                k[6*I-2,6*J]   += τ*B₂[i]*B₃[j]*𝑤
                # σ23
                k[6*I-1,6*J-4] += τ*B₃[i]*B₂[j]*𝑤
                k[6*I-1,6*J-3] += τ*B₂[i]*B₃[j]*𝑤
                k[6*I-1,6*J-2] += τ*B₃[i]*B₁[j]*𝑤
                k[6*I-1,6*J-1] += τ*(B₃[i]*B₃[j] + B₂[i]*B₂[j])*𝑤
                k[6*I-1,6*J]   += τ*B₂[i]*B₁[j]*𝑤
                # σ13
                k[6*I,6*J-5] += τ*B₃[i]*B₁[j]*𝑤
                k[6*I,6*J-3] += τ*B₁[i]*B₃[j]*𝑤
                k[6*I,6*J-2] += τ*B₃[i]*B₂[j]*𝑤
                k[6*I,6*J-1] += τ*B₁[i]*B₂[j]*𝑤
                k[6*I,6*J]   += τ*(B₃[i]*B₃[j] + B₁[i]*B₁[j])*𝑤
            end
            f[6*I-5] += τ*(B₁[i]*b₁)*𝑤
            f[6*I-4] += τ*(B₂[i]*b₂)*𝑤
            f[6*I-3] += τ*(B₃[i]*b₃)*𝑤
            f[6*I-2] += τ*(B₁[i]*b₂ + B₂[i]*b₁)*𝑤
            f[6*I-1] += τ*(B₃[i]*b₂ + B₂[i]*b₃)*𝑤
            f[6*I]   += τ*(B₃[i]*b₁ + B₁[i]*b₃)*𝑤
        end
    end
end


function ∫∫τ∇σᵢⱼ∇σᵢₖbdΩ(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        τ = ξ.τ
        b₁ = ξ.b₁
        b₂ = ξ.b₂
        b₃ = ξ.b₃
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z]
        E = ξ.E
        ν = ξ.ν
        C⁻¹ᵢᵢᵢᵢ = 1/E
        C⁻¹ᵢᵢⱼⱼ = -ν/E
        C⁻¹ᵢⱼᵢⱼ = 2*(1+ν)/E
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                gij = (1.0/ℎ^2*N[i]*N[j]) * 𝑤
                k[6*I-5,6*J-5] += τ*(B₁[i]*B₁[j]+ gij)*𝑤
                k[6*I-5,6*J-2] += τ*B₁[i]*B₂[j]*𝑤
                k[6*I-5,6*J]   += τ*B₁[i]*B₃[j]*𝑤

                k[6*I-4,6*J-4] += τ*(B₂[i]*B₂[j]+ gij)*𝑤
                k[6*I-4,6*J-2] += τ*B₂[i]*B₁[j]*𝑤
                k[6*I-4,6*J-1] += τ*B₂[i]*B₃[j]*𝑤
               
                k[6*I-3,6*J-3] += τ*(B₃[i]*B₃[j]+ gij)*𝑤
                k[6*I-3,6*J-1] += τ*B₃[i]*B₂[j]*𝑤
                k[6*I-3,6*J]   += τ*B₃[i]*B₁[j]*𝑤

                k[6*I-2,6*J-5] += τ*B₂[i]*B₁[j]*𝑤
                k[6*I-2,6*J-4] += τ*B₁[i]*B₂[j]*𝑤
                k[6*I-2,6*J-2] += τ*(B₁[i]*B₁[j]+B₂[i]*B₂[j]+ gij)*𝑤
                k[6*I-2,6*J-1] += τ*B₁[i]*B₃[j]*𝑤
                k[6*I-2,6*J]   += τ*B₂[i]*B₃[j]*𝑤

                k[6*I-1,6*J-4] += τ*B₃[i]*B₂[j]*𝑤
                k[6*I-1,6*J-3] += τ*B₂[i]*B₃[j]*𝑤
                k[6*I-1,6*J-2] += τ*B₃[i]*B₁[j]*𝑤
                k[6*I-1,6*J-1] += τ*(B₃[i]*B₃[j]+B₂[i]*B₂[j]+ gij)*𝑤
                k[6*I-1,6*J]   += τ*B₂[i]*B₁[j]*𝑤

                k[6*I-1,6*J-5] += τ*B₃[i]*B₁[j]*𝑤
                k[6*I-1,6*J-3] += τ*B₁[i]*B₃[j]*𝑤
                k[6*I-1,6*J-2] += τ*B₃[i]*B₂[j]*𝑤
                k[6*I-1,6*J-1] += τ*B₁[i]*B₂[j]*𝑤
                k[6*I-1,6*J]   += τ*(B₃[i]*B₃[j]+B₁[i]*B₁[j]+ gij)*𝑤

            end
            f[6*I-5] += τ*(B₁[i]*b₁)*𝑤
            f[6*I-4] += τ*(B₂[i]*b₂)*𝑤
            f[6*I-3] += τ*(B₃[i]*b₃)*𝑤
            f[6*I-2] += τ*(B₁[i]*b₂ + B₂[i]*b₁)*𝑤
            f[6*I-1] += τ*(B₃[i]*b₂ + B₂[i]*b₃)*𝑤
            f[6*I]   += τ*(B₃[i]*b₁ + B₁[i]*b₃)*𝑤
        end
    end
end


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

function L₂(ap::T) where T<:AbstractElement
    Δu²= 0
    ū² = 0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        ū₁ = ξ.u
        ū₂ = ξ.v
        u₁ = 0.
        u₂ = 0.
        for (i,xᵢ) in enumerate(ap.𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
        end
        Δu² += ((u₁ - ū₁)^2 + (u₂ - ū₂)^2)*𝑤
        ū² += (ū₁^2 + ū₂^2)*𝑤
    end
    return Δu², ū²
end

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


function Hₑ_PlaneStress(ap::T) where T<:AbstractElement
    ΔW²= 0
    W̄² = 0
    Δu²= 0
    ū² = 0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        E = ξ.E
        ν = ξ.ν
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

function Hₑ_PlaneStrain_Deviatoric(ap::T) where T<:AbstractElement
    ΔW²= 0
    W̄² = 0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        # E = ξ.E
        # ν = ξ.ν

        E = ξ.Ē
        ν = ξ.ν̄ 
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

function Hₑ_PlaneStrain_Dil(ap::T) where T<:AbstractElement
    ΔW²= 0
    W̄² = 0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        E = ξ.Ē
        ν = ξ.ν̄ 
        Cᵢᵢᵢᵢ = E/(1+ν)/(1-2*ν)*(1-ν)
        Cᵢᵢⱼⱼ = E/(1+ν)/(1-2*ν)*ν
        Cᵢⱼᵢⱼ = E/(1+ν)/2
        
        # Cᵢᵢᵢᵢ = E/(1-ν^2)
        # Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
        # Cᵢⱼᵢⱼ = E/2/(1+ν)

        K=E/3/(1-2ν)
        G=E/2/(1+ν)
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
        p̄ = (σ̄₁₁ + σ̄₂₂)/2
        ε₁₁ = 0.
        ε₂₂ = 0.
        ε₁₂ = 0.
        for (i,xᵢ) in enumerate(ap.𝓒)
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += B₂[i]*xᵢ.d₁ + B₁[i]*xᵢ.d₂
        end
        σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂
        σ₂₂ = Cᵢᵢⱼⱼ*ε₁₁ + Cᵢᵢᵢᵢ*ε₂₂ 
        p = (σ₁₁ + σ₂₂)/2
        
        ΔW² += (3*(p-p̄)^2/2/K)*𝑤
        W̄² += (3*p̄^2/2/K)*𝑤
    end
    return ΔW², W̄²
end

function Hₑ_PlaneStrain_Dil(aps::Vector{T}) where T<:AbstractElement
    HₑNorm_ΔW²= 0.0
    HₑNorm_W̄² = 0.0
    for ap in aps
        ΔW², W̄² = Hₑ_PlaneStrain_Dil(ap)
        HₑNorm_ΔW² += ΔW²
        HₑNorm_W̄²  += W̄²
    end
    return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5
end


function 𝐿₂_PlaneStrain_Pressure(ap::T) where T<:AbstractElement
    ΔW²= 0
    W̄² = 0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        E = ξ.Ē
        ν = ξ.ν̄ 
        K=E/3/(1-2ν)
        G=E/2/(1+ν)
        Cᵢᵢᵢᵢ = E/(1+ν)/(1-2*ν)*(1-ν)
        Cᵢᵢⱼⱼ = E/(1+ν)/(1-2*ν)*ν
        Cᵢⱼᵢⱼ = E/(1+ν)/2
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
        p̄ = (σ̄₁₁ + σ̄₂₂)/2
        ε₁₁ = 0.
        ε₂₂ = 0.
        ε₁₂ = 0.
        for (i,xᵢ) in enumerate(ap.𝓒)
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += B₂[i]*xᵢ.d₁ + B₁[i]*xᵢ.d₂
        end
        σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂
        σ₂₂ = Cᵢᵢⱼⱼ*ε₁₁ + Cᵢᵢᵢᵢ*ε₂₂ 
        p = (σ₁₁ + σ₂₂)/2
        
        ΔW² += ((p-p̄)^2)*𝑤
        W̄² += (p̄^2)*𝑤
    end
    return ΔW², W̄²
end

function 𝐿₂_PlaneStrain_Pressure(aps::Vector{T}) where T<:AbstractElement
    HₑNorm_ΔW²= 0.0
    HₑNorm_W̄² = 0.0
    for ap in aps
        ΔW², W̄² = 𝐿₂_PlaneStrain_Pressure(ap)
        HₑNorm_ΔW² += ΔW²
        HₑNorm_W̄²  += W̄²
    end
    return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5
end


function 𝐿₂_PlaneStrain_Pressure_HR(ap::T) where T<:AbstractElement
    ΔW²= 0
    W̄² = 0
    for ξ in ap.𝓖
        xc = ξ.x
        yc = ξ.y
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        E = ξ.Ē
        ν = ξ.ν̄ 
        K=E/3/(1-2ν)
        G=E/2/(1+ν)
        Cᵢᵢᵢᵢ = E/(1+ν)/(1-2*ν)*(1-ν)
        Cᵢᵢⱼⱼ = E/(1+ν)/(1-2*ν)*ν
        Cᵢⱼᵢⱼ = E/(1+ν)/2
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
        p̄ = (σ̄₁₁ + σ̄₂₂)/2
        ε₁₁ = 0.
        ε₂₂ = 0.
        ε₁₂ = 0.
        
        𝓒 = ap.𝓒
        σ₁₁ = 𝓒[1].dₛ₁₁+𝓒[2].dₛ₁₁*xc+𝓒[3].dₛ₁₁*yc
        σ₂₂ = 𝓒[1].dₛ₂₂+𝓒[2].dₛ₂₂*xc+𝓒[3].dₛ₂₂*yc
        σ₁₂ = 𝓒[1].dₛ₁₂+𝓒[2].dₛ₁₂*xc+𝓒[3].dₛ₁₂*yc
        
        # σ₁₁ = 0.
        # σ₂₂ = 0.
        # σ₁₂ = 0.
        #  for (i,xᵢ) in  enumerate(𝓒)
        #    σ₁₁ += N[i]*xᵢ.dₛ₁₁
        #    σ₂₂ += N[i]*xᵢ.dₛ₂₂
        #    σ₁₂ += N[i]*xᵢ.dₛ₁₂
        # end
        # σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂
        # σ₂₂ = Cᵢᵢⱼⱼ*ε₁₁ + Cᵢᵢᵢᵢ*ε₂₂ 
        p = (σ₁₁ + σ₂₂)/2
        
        ΔW² += ((p-p̄)^2)*𝑤
        W̄² += (p̄^2)*𝑤
    end
    return ΔW², W̄²
end

function 𝐿₂_PlaneStrain_Pressure_HR(aps::Vector{T}) where T<:AbstractElement
    HₑNorm_ΔW²= 0.0
    HₑNorm_W̄² = 0.0
    for ap in aps
        ΔW², W̄² = 𝐿₂_PlaneStrain_Pressure_HR(ap)
        HₑNorm_ΔW² += ΔW²
        HₑNorm_W̄²  += W̄²
    end
    return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5
end

function 𝐿₂_3D_Pressure(ap::T) where T<:AbstractElement
    ΔW² = 0.0
    W̄² = 0.0

    for ξ in ap.𝓖
        𝑤 = ξ.𝑤

        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z]

        E = ξ.E
        ν = ξ.ν

        K = E / 3.0 / (1.0 - 2.0*ν)

        # ---------------------------------------------------------
        # exact pressure from exact displacement gradient
        # small-strain 3D pressure:
        # p̄ = K * div(ū)
        # ---------------------------------------------------------
        ∂ū₁∂x = ξ.∂u₁∂x
        ∂ū₂∂y = ξ.∂u₂∂y
        ∂ū₃∂z = ξ.∂u₃∂z

        divū = ∂ū₁∂x + ∂ū₂∂y + ∂ū₃∂z
        p̄ = K * divū

        # ---------------------------------------------------------
        # numerical pressure from numerical displacement gradient
        # ---------------------------------------------------------
        ∂u₁∂x = 0.0
        ∂u₂∂y = 0.0
        ∂u₃∂z = 0.0

        for (i, xᵢ) in enumerate(ap.𝓒)
            ∂u₁∂x += B₁[i] * xᵢ.d₁
            ∂u₂∂y += B₂[i] * xᵢ.d₂
            ∂u₃∂z += B₃[i] * xᵢ.d₃
        end

        divu = ∂u₁∂x + ∂u₂∂y + ∂u₃∂z
        p = K * divu

        ΔW² += (p - p̄)^2 * 𝑤
        W̄² += p̄^2 * 𝑤
    end

    return ΔW², W̄²
end


function 𝐿₂_3D_Pressure(aps::Vector{T}) where T<:AbstractElement
    ΔW²_total = 0.0
    W̄²_total = 0.0

    for ap in aps
        ΔW², W̄² = 𝐿₂_3D_Pressure(ap)
        ΔW²_total += ΔW²
        W̄²_total += W̄²
    end

    return sqrt(ΔW²_total / W̄²_total)
end

function 𝐿₂_3D_Pressure_SVK(ap::T) where T<:AbstractElement
    ΔW² = 0.0
    W̄² = 0.0

    for ξ in ap.𝓖
        𝑤 = ξ.𝑤

        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z]

        E = ξ.E
        ν = ξ.ν

        K = E / 3.0 / (1.0 - 2.0*ν)

        # ---------------------------------------------------------
        # exact F
        # ---------------------------------------------------------
        F̄11 = 1.0 + ξ.∂u₁∂x
        F̄12 =       ξ.∂u₁∂y
        F̄13 =       ξ.∂u₁∂z

        F̄21 =       ξ.∂u₂∂x
        F̄22 = 1.0 + ξ.∂u₂∂y
        F̄23 =       ξ.∂u₂∂z

        F̄31 =       ξ.∂u₃∂x
        F̄32 =       ξ.∂u₃∂y
        F̄33 = 1.0 + ξ.∂u₃∂z

        Ē11 = 0.5 * (F̄11^2 + F̄21^2 + F̄31^2 - 1.0)
        Ē22 = 0.5 * (F̄12^2 + F̄22^2 + F̄32^2 - 1.0)
        Ē33 = 0.5 * (F̄13^2 + F̄23^2 + F̄33^2 - 1.0)

        p̄ = K * (Ē11 + Ē22 + Ē33)

        # ---------------------------------------------------------
        # numerical F
        # ---------------------------------------------------------
        u1_x = 0.0
        u1_y = 0.0
        u1_z = 0.0

        u2_x = 0.0
        u2_y = 0.0
        u2_z = 0.0

        u3_x = 0.0
        u3_y = 0.0
        u3_z = 0.0

        for (i, xᵢ) in enumerate(ap.𝓒)
            u1_x += B₁[i] * xᵢ.d₁
            u1_y += B₂[i] * xᵢ.d₁
            u1_z += B₃[i] * xᵢ.d₁

            u2_x += B₁[i] * xᵢ.d₂
            u2_y += B₂[i] * xᵢ.d₂
            u2_z += B₃[i] * xᵢ.d₂

            u3_x += B₁[i] * xᵢ.d₃
            u3_y += B₂[i] * xᵢ.d₃
            u3_z += B₃[i] * xᵢ.d₃
        end

        F11 = 1.0 + u1_x
        F12 =       u1_y
        F13 =       u1_z

        F21 =       u2_x
        F22 = 1.0 + u2_y
        F23 =       u2_z

        F31 =       u3_x
        F32 =       u3_y
        F33 = 1.0 + u3_z

        E11 = 0.5 * (F11^2 + F21^2 + F31^2 - 1.0)
        E22 = 0.5 * (F12^2 + F22^2 + F32^2 - 1.0)
        E33 = 0.5 * (F13^2 + F23^2 + F33^2 - 1.0)

        p = K * (E11 + E22 + E33)

        ΔW² += (p - p̄)^2 * 𝑤
        W̄² += p̄^2 * 𝑤
    end

    return ΔW², W̄²
end


function 𝐿₂_3D_Pressure_SVK(aps::Vector{T}) where T<:AbstractElement
    ΔW²_total = 0.0
    W̄²_total = 0.0

    for ap in aps
        ΔW², W̄² = 𝐿₂_3D_Pressure_SVK(ap)
        ΔW²_total += ΔW²
        W̄²_total += W̄²
    end

    return sqrt(ΔW²_total / W̄²_total)
end
# function Hₑ(ap::T) where T<:AbstractElement
#     ΔW²= 0.0
#     W̄² = 0.0
#     Δu²= 0.0
#     ū² = 0.0
#     for ξ in ap.𝓖
#         𝑤 = ξ.𝑤
#         #  
#         N = ξ[:𝝭]
#         B₁ = ξ[:∂𝝭∂x]
#         B₂ = ξ[:∂𝝭∂y]
#         B₃ = ξ[:∂𝝭∂z]
#         E = ξ.E
#         ν = ξ.ν
#         Cᵢᵢᵢᵢ = E*(1-ν)/(1-2*ν)/(1+ν)
#         Cᵢᵢⱼⱼ = E*ν/(1-2*ν)/(1+ν)
#         Cᵢⱼᵢⱼ = E/2/(1+ν)
#         ū₁ = ξ.u₁
#         ū₂ = ξ.u₂
#         ū₃ = ξ.u₃
#         ∂ū₁∂x = ξ.∂u₁∂x
#         ∂ū₁∂y = ξ.∂u₁∂y
#         ∂ū₁∂z = ξ.∂u₁∂z
#         ∂ū₂∂x = ξ.∂u₂∂x
#         ∂ū₂∂y = ξ.∂u₂∂y
#         ∂ū₂∂z = ξ.∂u₂∂z
#         ∂ū₃∂x = ξ.∂u₃∂x
#         ∂ū₃∂y = ξ.∂u₃∂y
#         ∂ū₃∂z = ξ.∂u₃∂z
#         ε̄₁₁ = ∂ū₁∂x
#         ε̄₂₂ = ∂ū₂∂y
#         ε̄₃₃ = ∂ū₃∂z
#         ε̄₁₂ = ∂ū₁∂y + ∂ū₂∂x
#         ε̄₁₃ = ∂ū₁∂z + ∂ū₃∂x
#         ε̄₂₃ = ∂ū₂∂z + ∂ū₃∂y
#         σ̄₁₁ = Cᵢᵢᵢᵢ*ε̄₁₁ + Cᵢᵢⱼⱼ*ε̄₂₂ + Cᵢᵢⱼⱼ*ε̄₃₃
#         σ̄₂₂ = Cᵢᵢⱼⱼ*ε̄₁₁ + Cᵢᵢᵢᵢ*ε̄₂₂ + Cᵢᵢⱼⱼ*ε̄₃₃ 
#         σ̄₃₃ = Cᵢᵢⱼⱼ*ε̄₁₁ + Cᵢᵢⱼⱼ*ε̄₂₂ + Cᵢᵢᵢᵢ*ε̄₃₃ 
#         σ̄₁₂ = Cᵢⱼᵢⱼ*ε̄₁₂
#         σ̄₁₃ = Cᵢⱼᵢⱼ*ε̄₁₃
#         σ̄₂₃ = Cᵢⱼᵢⱼ*ε̄₂₃
#         u₁ = 0.
#         u₂ = 0.
#         u₃ = 0.
#         ε₁₁ = 0.
#         ε₂₂ = 0.
#         ε₃₃ = 0.
#         ε₁₂ = 0.
#         ε₁₃ = 0.
#         ε₂₃ = 0.
#         for (i,xᵢ) in enumerate(ap.𝓒)
#             u₁ += N[i]*xᵢ.d₁
#             u₂ += N[i]*xᵢ.d₂
#             u₃ += N[i]*xᵢ.d₃   
#             ε₁₁ += B₁[i]*xᵢ.d₁
#             ε₂₂ += B₂[i]*xᵢ.d₂
#             ε₃₃ += B₃[i]*xᵢ.d₃
#             ε₁₂ += B₂[i]*xᵢ.d₁ + B₁[i]*xᵢ.d₂
#             ε₁₃ += B₃[i]*xᵢ.d₁ + B₁[i]*xᵢ.d₃
#             ε₂₃ += B₃[i]*xᵢ.d₂ + B₂[i]*xᵢ.d₃
#         end
#         σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂ + Cᵢᵢⱼⱼ*ε₃₃
#         σ₂₂ = Cᵢᵢⱼⱼ*ε₁₁ + Cᵢᵢᵢᵢ*ε₂₂ + Cᵢᵢⱼⱼ*ε₃₃ 
#         σ₃₃ = Cᵢᵢⱼⱼ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂ + Cᵢᵢᵢᵢ*ε₃₃ 
#         σ₁₂ = Cᵢⱼᵢⱼ*ε₁₂
#         σ₁₃ = Cᵢⱼᵢⱼ*ε₁₃
#         σ₂₃ = Cᵢⱼᵢⱼ*ε₂₃
#         ΔW² += 0.5*((σ₁₁-σ̄₁₁)*(ε₁₁-ε̄₁₁) + (σ₂₂-σ̄₂₂)*(ε₂₂-ε̄₂₂) + (σ₃₃-σ̄₃₃)*(ε₃₃-ε̄₃₃) + (σ₁₂-σ̄₁₂)*(ε₁₂-ε̄₁₂) + (σ₁₃-σ̄₁₃)*(ε₁₃-ε̄₁₃) + (σ₂₃-σ̄₂₃)*(ε₂₃-ε̄₂₃))*𝑤
#         W̄² += 0.5*(σ̄₁₁*ε̄₁₁ + σ̄₂₂*ε̄₂₂ + σ̄₃₃*ε̄₃₃ + σ̄₁₂*ε̄₁₂ + σ̄₁₃*ε̄₁₃ + σ̄₂₃*ε̄₂₃)*𝑤
#         Δu² += ((u₁ - ū₁)^2 + (u₂ - ū₂)^2 + (u₃ - ū₃)^2)*𝑤
#         ū² += (ū₁^2 + ū₂^2 + ū₃^2)*𝑤
#     end
#     return ΔW², W̄², Δu², ū²
# end

# function Hₑ(aps::Vector{T}) where T<:AbstractElement
#     HₑNorm_ΔW²= 0.0
#     HₑNorm_W̄² = 0.0
#     L₂Norm_Δu²= 0.0
#     L₂Norm_ū² = 0.0
#     for ap in aps
#         ΔW², W̄², Δu², ū² = Hₑ(ap)
#         HₑNorm_ΔW² += ΔW²
#         HₑNorm_W̄²  += W̄²
#         L₂Norm_Δu² += Δu²
#         L₂Norm_ū²  += ū²
#     end
#     return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
# end


function Hₑ(ap::T) where T<:AbstractElement
    ΔW² = 0.0
    W̄² = 0.0
    Δu² = 0.0
    ū² = 0.0

    for ξ in ap.𝓖
        w = ξ.𝑤

        N  = ξ[:𝝭]
        Bx = ξ[:∂𝝭∂x]
        By = ξ[:∂𝝭∂y]
        Bz = ξ[:∂𝝭∂z]

        E = ξ.E
        ν = ξ.ν

        λ = E * ν / ((1 + ν) * (1 - 2ν))
        G = E / (2 * (1 + ν))

        # exact displacement
        ū1 = ξ.u₁
        ū2 = ξ.u₂
        ū3 = ξ.u₃

        # exact strain, engineering shear strain convention
        ε̄11 = ξ.∂u₁∂x
        ε̄22 = ξ.∂u₂∂y
        ε̄33 = ξ.∂u₃∂z
        ε̄12 = ξ.∂u₁∂y + ξ.∂u₂∂x
        ε̄13 = ξ.∂u₁∂z + ξ.∂u₃∂x
        ε̄23 = ξ.∂u₂∂z + ξ.∂u₃∂y

        # FE displacement and strain
        u1 = 0.0
        u2 = 0.0
        u3 = 0.0
        ε11 = 0.0
        ε22 = 0.0
        ε33 = 0.0
        ε12 = 0.0
        ε13 = 0.0
        ε23 = 0.0

        for (i, xᵢ) in enumerate(ap.𝓒)
            d1 = xᵢ.d₁
            d2 = xᵢ.d₂
            d3 = xᵢ.d₃

            Ni  = N[i]
            Bxi = Bx[i]
            Byi = By[i]
            Bzi = Bz[i]

            u1 = muladd(Ni, d1, u1)
            u2 = muladd(Ni, d2, u2)
            u3 = muladd(Ni, d3, u3)

            ε11 = muladd(Bxi, d1, ε11)
            ε22 = muladd(Byi, d2, ε22)
            ε33 = muladd(Bzi, d3, ε33)

            ε12 = muladd(Byi, d1, ε12)
            ε12 = muladd(Bxi, d2, ε12)

            ε13 = muladd(Bzi, d1, ε13)
            ε13 = muladd(Bxi, d3, ε13)

            ε23 = muladd(Bzi, d2, ε23)
            ε23 = muladd(Byi, d3, ε23)
        end

        # strain error
        Δε11 = ε11 - ε̄11
        Δε22 = ε22 - ε̄22
        Δε33 = ε33 - ε̄33
        Δε12 = ε12 - ε̄12
        Δε13 = ε13 - ε̄13
        Δε23 = ε23 - ε̄23

        # exact energy density: 1/2 * ε̄^T D ε̄
        trε̄ = ε̄11 + ε̄22 + ε̄33
        σ̄11 = λ * trε̄ + 2G * ε̄11
        σ̄22 = λ * trε̄ + 2G * ε̄22
        σ̄33 = λ * trε̄ + 2G * ε̄33
        σ̄12 = G * ε̄12
        σ̄13 = G * ε̄13
        σ̄23 = G * ε̄23

        # error energy density: 1/2 * Δε^T D Δε
        trΔε = Δε11 + Δε22 + Δε33
        Δσ11 = λ * trΔε + 2G * Δε11
        Δσ22 = λ * trΔε + 2G * Δε22
        Δσ33 = λ * trΔε + 2G * Δε33
        Δσ12 = G * Δε12
        Δσ13 = G * Δε13
        Δσ23 = G * Δε23

        ΔW² = muladd(0.5 * (
            Δσ11 * Δε11 + Δσ22 * Δε22 + Δσ33 * Δε33 +
            Δσ12 * Δε12 + Δσ13 * Δε13 + Δσ23 * Δε23
        ), w, ΔW²)

        W̄² = muladd(0.5 * (
            σ̄11 * ε̄11 + σ̄22 * ε̄22 + σ̄33 * ε̄33 +
            σ̄12 * ε̄12 + σ̄13 * ε̄13 + σ̄23 * ε̄23
        ), w, W̄²)

        du1 = u1 - ū1
        du2 = u2 - ū2
        du3 = u3 - ū3

        Δu² = muladd(du1, du1, Δu²)
        Δu² = muladd(du2, du2, Δu²)
        Δu² = muladd(du3, du3, Δu²)
        Δu² *= 1.0
        Δu² = Δu²  # keep style explicit
        # multiply by weight without losing previous sum
        Δu² = Δu²

        # safer weighted accumulation
        # undo temporary style above by using separate scalar:
        # (done below for clarity)
        ū_loc = ū1*ū1 + ū2*ū2 + ū3*ū3
        Δu_loc = du1*du1 + du2*du2 + du3*du3

        # correct weighted sums
        Δu² -= Δu_loc  # cancel local temp contribution
        Δu² = muladd(Δu_loc, w, Δu²)
        ū² = muladd(ū_loc, w, ū²)
    end

    return ΔW², W̄², Δu², ū²
end

function Hₑ(aps::Vector{T}) where T<:AbstractElement
    HₑNorm_ΔW² = 0.0
    HₑNorm_W̄² = 0.0
    L₂Norm_Δu² = 0.0
    L₂Norm_ū² = 0.0

    for ap in aps
        ΔW², W̄², Δu², ū² = Hₑ(ap)
        HₑNorm_ΔW² += ΔW²
        HₑNorm_W̄² += W̄²
        L₂Norm_Δu² += Δu²
        L₂Norm_ū² += ū²
    end

    eH = HₑNorm_W̄² > eps(Float64) ? sqrt(max(HₑNorm_ΔW² / HₑNorm_W̄², 0.0)) :
                                    sqrt(max(HₑNorm_ΔW², 0.0))

    eL2 = L₂Norm_ū² > eps(Float64) ? sqrt(max(L₂Norm_Δu² / L₂Norm_ū², 0.0)) :
                                      sqrt(max(L₂Norm_Δu², 0.0))

    return eH, eL2
end

function max_error(ap::T) where T<:AbstractElement
    max_u = 0.0
    max_ε = 0.0

    for ξ in ap.𝓖
        N  = ξ[:𝝭]
        Bx = ξ[:∂𝝭∂x]
        By = ξ[:∂𝝭∂y]
        Bz = ξ[:∂𝝭∂z]

        u1 = 0.0; u2 = 0.0; u3 = 0.0
        ε11 = 0.0; ε22 = 0.0; ε33 = 0.0
        ε12 = 0.0; ε13 = 0.0; ε23 = 0.0

        for (i, xᵢ) in enumerate(ap.𝓒)
            u1 += N[i]*xᵢ.d₁
            u2 += N[i]*xᵢ.d₂
            u3 += N[i]*xᵢ.d₃

            ε11 += Bx[i]*xᵢ.d₁
            ε22 += By[i]*xᵢ.d₂
            ε33 += Bz[i]*xᵢ.d₃
            ε12 += By[i]*xᵢ.d₁ + Bx[i]*xᵢ.d₂
            ε13 += Bz[i]*xᵢ.d₁ + Bx[i]*xᵢ.d₃
            ε23 += Bz[i]*xᵢ.d₂ + By[i]*xᵢ.d₃
        end

        max_u = max(max_u,
            abs(u1 - ξ.u₁),
            abs(u2 - ξ.u₂),
            abs(u3 - ξ.u₃)
        )

        max_ε = max(max_ε,
            abs(ε11 - ξ.∂u₁∂x),
            abs(ε22 - ξ.∂u₂∂y),
            abs(ε33 - ξ.∂u₃∂z),
            abs(ε12 - (ξ.∂u₁∂y + ξ.∂u₂∂x)),
            abs(ε13 - (ξ.∂u₁∂z + ξ.∂u₃∂x)),
            abs(ε23 - (ξ.∂u₂∂z + ξ.∂u₃∂y))
        )
    end

    return max_u, max_ε
end


end