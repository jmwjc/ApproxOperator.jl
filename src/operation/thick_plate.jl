module MindlinPlate 
    
using ..ApproxOperator: AbstractElement

function ∫κMγQdΩ(ap::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = ap.E
    ν = ap.ν
    h = ap.h
    for ξ in 𝓖
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        Dᵇᵢᵢᵢᵢ = E*h^3/12/(1-ν^2)
        Dᵇᵢᵢⱼⱼ = E*ν*h^3/12/(1-ν^2)
        Dᵇᵢⱼᵢⱼ = E*h^3/24/(1+ν)
        Dˢ =  5/6*h*E/2/(1+ν)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += ( Dˢ*B₁[i]*B₁[j] + Dˢ*B₂[i]*B₂[j])*𝑤
                k[3*I-2,3*J-1] += (-Dˢ*B₁[i]*N[j])*𝑤
                k[3*I-2,3*J]   += (-Dˢ*B₂[i]*N[j])*𝑤
                k[3*I-1,3*J-2] += (-Dˢ*N[i]*B₁[j])*𝑤
                k[3*I-1,3*J-1] += (-Dᵇᵢᵢᵢᵢ*B₁[i]*B₁[j] - Dᵇᵢⱼᵢⱼ*B₂[i]*B₂[j] + Dˢ*N[i]*N[j])*𝑤
                k[3*I-1,3*J]   += (-Dᵇᵢᵢⱼⱼ*B₁[i]*B₂[j] - Dᵇᵢⱼᵢⱼ*B₂[i]*B₁[j])*𝑤
                k[3*I,3*J-2]   += (-Dˢ*N[i]*B₂[j])*𝑤
                k[3*I,3*J-1]   += (-Dᵇᵢᵢⱼⱼ*B₂[i]*B₁[j] - Dᵇᵢⱼᵢⱼ*B₁[i]*B₂[j])*𝑤
                k[3*I,3*J]     += (-Dᵇᵢᵢᵢᵢ*B₂[i]*B₂[j] - Dᵇᵢⱼᵢⱼ*B₁[i]*B₁[j] + Dˢ*N[i]*N[j])*𝑤
            end
        end
    end
end

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

function ∫wqdΩ(ap::T,f::AbstractVector) where T<:AbstractElement
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

function ∫φmdΩ(ap::T,f::AbstractVector) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        m₁ = ξ.m₁
        m₂ = ξ.m₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += N[i]*m₁*𝑤
            f[2*I]   += N[i]*m₂*𝑤
        end
    end
end

function ∫wVdΓ(ap::T,f::AbstractVector) where T<:AbstractElement
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

function ∫φMdΓ(ap::T,f::AbstractVector) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        M₁ = ξ.M₁
        M₂ = ξ.M₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] -= N[i]*M₁*𝑤
            f[2*I]   -= N[i]*M₂*𝑤
        end
    end
end

function ∫Q∇wdΩ(ap::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,J] += N[i]*B₁[j]*𝑤
                k[2*I,J]   += N[i]*B₂[j]*𝑤
            end
        end
    end
end

function ∫Q∇wdΩ(a::T,b::S,k::AbstractMatrix) where {T<:AbstractElement,S<:AbstractElement}
    𝓒₁ = a.𝓒; 𝓖₁ = a.𝓖
    𝓒₂ = b.𝓒; 𝓖₂ = b.𝓖
    for (ξ₁,ξ₂) in zip(𝓖₁,𝓖₂)
        𝑤 = ξ₁.𝑤
        N = ξ₁[:𝝭]
        B̄₁ = ξ₂[:∂𝝭∂x]
        B̄₂ = ξ₂[:∂𝝭∂y]
        for (i,xᵢ) in enumerate(𝓒₁)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒₂)
                J = xⱼ.𝐼
                k[2*I-1,J] += N[i]*B̄₁[j]*𝑤
                k[2*I,J]   += N[i]*B̄₂[j]*𝑤
            end
        end
    end
end

function ∫∇QwdΩ(ap::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,J] -= B₁[i]*N[j]*𝑤
                k[2*I,J]   -= B₂[i]*N[j]*𝑤
            end
        end
    end
end

function ∫∇QwdΩ(a::T,b::S,k::AbstractMatrix) where {T<:AbstractElement,S<:AbstractElement}
    𝓒₁ = a.𝓒; 𝓖₁ = a.𝓖
    𝓒₂ = b.𝓒; 𝓖₂ = b.𝓖
    for (ξ₁,ξ₂) in zip(𝓖₁,𝓖₂)
        𝑤 = ξ₁.𝑤
        B₁ = ξ₁[:∂𝝭∂x]
        B₂ = ξ₁[:∂𝝭∂y]
        N̄ = ξ₂[:𝝭]
        for (i,xᵢ) in enumerate(𝓒₁)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒₂)
                J = xⱼ.𝐼
                k[2*I-1,J] -= B₁[i]*N̄[j]*𝑤
                k[2*I,J]   -= B₂[i]*N̄[j]*𝑤
            end
        end
    end
end

function ∫QφdΩ(ap::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] -= N[i]*N[j]*𝑤
                k[2*I,2*J]     -= N[i]*N[j]*𝑤
            end
        end
    end
end

function ∫QφdΩ(a::T,b::S,k::AbstractMatrix) where {T<:AbstractElement,S<:AbstractElement}
    𝓒₁ = a.𝓒; 𝓖₁ = a.𝓖
    𝓒₂ = b.𝓒; 𝓖₂ = b.𝓖
    for (ξ₁,ξ₂) in zip(𝓖₁,𝓖₂)
        𝑤 = ξ₁.𝑤
        N = ξ₁[:𝝭]
        N̄ = ξ₂[:𝝭]
        for (i,xᵢ) in enumerate(𝓒₁)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒₂)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] -= N[i]*N̄[j]*𝑤
                k[2*I,2*J]     -= N[i]*N̄[j]*𝑤
            end
        end
    end
end

function ∫QQdΩ(ap::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = ap.E
    ν = ap.ν
    h = ap.h
    D = 5/6*h*E/2/(1+ν)
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] -= N[i]*N[j]/D*𝑤
                k[2*I,2*J]     -= N[i]*N[j]/D*𝑤
            end
        end
    end
end

function ∫QwdΓ(ap::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,J] += n₁*N[i]*N[j]*𝑤
                k[2*I,J]   += n₂*N[i]*N[j]*𝑤
            end
        end
    end
end

function ∫QwdΓ(a::T,b::S,k::AbstractMatrix,f::AbstractVector) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵛ = a.𝓒; 𝓖ᵛ= a.𝓖
    𝓒ʷ = b.𝓒; 𝓖ʷ= b.𝓖
    for (ξᵛ,ξʷ) in zip(𝓖ᵛ,𝓖ʷ)
        𝑤 = ξᵛ.𝑤
        Nᵛ = ξᵛ[:𝝭]
        Nʷ = ξʷ[:𝝭]
        n₁ = ξᵛ.n₁
        n₂ = ξᵛ.n₂
        g = ξʷ.g
        for (i,xᵢ) in enumerate(𝓒ᵛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ʷ)
                J = xⱼ.𝐼
                k[2*I-1,J] -= n₁*Nᵛ[i]*Nʷ[j]*𝑤
                k[2*I,J]   -= n₂*Nᵛ[i]*Nʷ[j]*𝑤
            end
            f[2*I-1] -= n₁*Nᵛ[i]*g*𝑤
            f[2*I]   -= n₂*Nᵛ[i]*g*𝑤
        end
    end
end

function ∫QwdΓ(a::T,b::S,k::AbstractMatrix) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵛ = a.𝓒; 𝓖ᵛ= a.𝓖
    𝓒ʷ = b.𝓒; 𝓖ʷ= b.𝓖
    for (ξᵛ,ξʷ) in zip(𝓖ᵛ,𝓖ʷ)
        𝑤 = ξᵛ.𝑤
        Nᵛ = ξᵛ[:𝝭]
        Nʷ = ξʷ[:𝝭]
        n₁ = ξᵛ.n₁
        n₂ = ξᵛ.n₂
        for (i,xᵢ) in enumerate(𝓒ᵛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ʷ)
                J = xⱼ.𝐼
                k[2*I-1,J] += n₁*Nᵛ[i]*Nʷ[j]*𝑤
                k[2*I,J]   += n₂*Nᵛ[i]*Nʷ[j]*𝑤
            end
        end
    end
end

function ∫QwdΓ(ap::T,k::AbstractMatrix,f::AbstractVector) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,J] -= n₁*N[i]*N[j]*𝑤
                k[2*I,J]   -= n₂*N[i]*N[j]*𝑤
            end
            f[2*I-1] -= n₁*N[i]*g*𝑤
            f[2*I]   -= n₂*N[i]*g*𝑤
        end
    end
end

function ∫MMdΩ(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        E = ξ.E
        ν = ξ.ν
        h = ξ.h
        C⁻¹ᵢᵢᵢᵢ = 1/E*12/h^3
        C⁻¹ᵢᵢⱼⱼ = -ν/E*12/h^3
        C⁻¹ᵢⱼᵢⱼ = 2*(1+ν)/E*12/h^3
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] -= N[i]*C⁻¹ᵢᵢᵢᵢ*N[j]*𝑤
                k[3*I-2,3*J-1] -= N[i]*C⁻¹ᵢᵢⱼⱼ*N[j]*𝑤
                k[3*I-1,3*J-2] -= N[i]*C⁻¹ᵢᵢⱼⱼ*N[j]*𝑤
                k[3*I-1,3*J-1] -= N[i]*C⁻¹ᵢᵢᵢᵢ*N[j]*𝑤
                k[3*I,3*J]     -= N[i]*C⁻¹ᵢⱼᵢⱼ*N[j]*𝑤
            end
        end
    end
end

function ∫MφdΓ(aₘ::T,aᵩ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₘ = aₘ.𝓒;𝓖ₘ = aₘ.𝓖
    𝓒ᵩ = aᵩ.𝓒;𝓖ᵩ = aᵩ.𝓖
    for (ξₘ,ξᵩ) in zip(𝓖ₘ,𝓖ᵩ)
        𝑤 = ξₘ.𝑤
        N = ξₘ[:𝝭]
        N̄ = ξᵩ[:𝝭]
        n₁ = ξᵩ.n₁
        n₂ = ξᵩ.n₂
        for (i,xᵢ) in enumerate(𝓒ₘ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵩ)
                J = xⱼ.𝐼
                k[3*I-2,2*J-1] -= N[i]*n₁*N̄[j]*𝑤
                k[3*I-1,2*J]   -= N[i]*n₂*N̄[j]*𝑤
                k[3*I,2*J-1]   -= N[i]*n₂*N̄[j]*𝑤
                k[3*I,2*J]     -= N[i]*n₁*N̄[j]*𝑤
            end
        end
    end
end

function ∫MφdΓ(aₘ::T,aᵩ::S,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₘ = aₘ.𝓒;𝓖ₘ = aₘ.𝓖
    𝓒ᵩ = aᵩ.𝓒;𝓖ᵩ = aᵩ.𝓖
    for (ξₘ,ξᵩ) in zip(𝓖ₘ,𝓖ᵩ)
        𝑤 = ξₘ.𝑤
        N = ξₘ[:𝝭]
        N̄ = ξᵩ[:𝝭]
        n₁ = ξᵩ.n₁
        n₂ = ξᵩ.n₂
        n₁₁ = ξᵩ.n₁₁
        n₁₂ = ξᵩ.n₁₂
        n₂₂ = ξᵩ.n₂₂
        g₁ = ξᵩ.g₁
        g₂ = ξᵩ.g₂
        for (i,xᵢ) in enumerate(𝓒ₘ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵩ)
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

function ∫∇MφdΩ(aₘ::T,aᵩ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₘ = aₘ.𝓒;𝓖ₘ = aₘ.𝓖
    𝓒ᵩ = aᵩ.𝓒;𝓖ᵩ = aᵩ.𝓖
    for (ξₘ,ξᵩ) in zip(𝓖ₘ,𝓖ᵩ)
        𝑤 = ξₘ.𝑤
        B₁ = ξₘ[:∂𝝭∂x]
        B₂ = ξₘ[:∂𝝭∂y]
        N = ξᵩ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒ₘ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵩ)
                J = xⱼ.𝐼
                k[3*I-2,2*J-1] += B₁[i]*N[j]*𝑤
                k[3*I-1,2*J]   += B₂[i]*N[j]*𝑤
                k[3*I,2*J-1]   += B₂[i]*N[j]*𝑤
                k[3*I,2*J]     += B₁[i]*N[j]*𝑤
            end
        end
    end
end

function ∫αwwdΓ(ap::T,k::AbstractMatrix,f::AbstractVector) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    α = ap.α
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += α*N[i]*N[j]*𝑤
            end
            f[I] += α*N[i]*g*𝑤
        end
    end
end

function ∫αφφdΓ(ap::T,k::AbstractMatrix,f::AbstractVector) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    α = ap.α
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        n₁₁ = ξ.n₁₁
        n₁₂ = ξ.n₁₂
        n₂₂ = ξ.n₂₂
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

function L₂Q(ap::T) where T<:AbstractElement
    ΔQ²= 0
    Q̄² = 0
    E = ap.E
    ν = ap.ν
    h = ap.h
    Dˢ = 5/6*E*h/(2*(1+ν))
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        Q̄₁ = ξ.Q₁
        Q̄₂ = ξ.Q₂
        u = 0
        Q₁ = 0
        Q₂ = 0
        for (i,xᵢ) in enumerate(ap.𝓒)
            Q₁ += N[i]*xᵢ.q₁
            Q₂ += N[i]*xᵢ.q₂
        end
        ΔQ² +=((Q₁ - Q̄₁)^2 + (Q₂ - Q̄₂)^2)/Dˢ*𝑤
        Q̄²  += (Q̄₁^2 + Q̄₂^2)/Dˢ*𝑤
    end
    return ΔQ², Q̄²
end

function L₂Q(aps::Vector{T}) where T<:AbstractElement
    L₂Norm_ΔQ²= 0
    L₂Norm_Q̄² = 0
    for ap in aps
        ΔQ², Q̄² = L₂Q(ap)
        L₂Norm_ΔQ² += ΔQ²
        L₂Norm_Q̄²  += Q̄²
    end
    return (L₂Norm_ΔQ²/L₂Norm_Q̄²)^0.5
end


function L₂φ(ap::T) where T<:AbstractElement
    Δφ²= 0
    φ̄² = 0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        φ̄₁ = ξ.φ₁
        φ̄₂ = ξ.φ₂
        u = 0
        φ₁ = 0
        φ₂ = 0
        for (i,xᵢ) in enumerate(ap.𝓒)
            φ₁ += N[i]*xᵢ.d₁
            φ₂ += N[i]*xᵢ.d₂
        end
        Δφ² +=((φ₁ - φ̄₁)^2 + (φ₂ - φ̄₂)^2)*𝑤
        φ̄²  += (φ̄₁^2 + φ̄₂^2)*𝑤
    end
    return Δφ², φ̄²
end

function L₂φ(aps::Vector{T}) where T<:AbstractElement
    L₂Norm_Δφ²= 0
    L₂Norm_φ̄² = 0
    for ap in aps
        Δφ², φ̄² = L₂φ(ap)
        L₂Norm_Δφ² += Δφ²
        L₂Norm_φ̄²  += φ̄²
    end
    return (L₂Norm_Δφ²/L₂Norm_φ̄²)^0.5
end

function L₂(ap::T) where T<:AbstractElement
    Δu²= 0
    ū² = 0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        ū = ξ.u
        u = 0
        for (i,xᵢ) in enumerate(ap.𝓒)
            u += N[i]*xᵢ.d
        end
        Δu² +=(u - ū)^2*𝑤
        ū²  += ū^2*𝑤
    end
    return Δu², ū²
end

function L₂(aps::Vector{T}) where T<:AbstractElement
    L₂Norm_Δu²= 0
    L₂Norm_ū² = 0
    for ap in aps
        Δu², ū² = L₂(ap)
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end

end