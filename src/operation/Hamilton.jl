function (op::Operator{:∫q̇mpqkpdx})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    m = op.m
    kᶜ = op.kᶜ
    for ξ in 𝓖
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

function (op::Operator{:∫∫q̇mpqkpdx})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    ρA = op.ρA
    EA = op.EA
    for ξ in 𝓖
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

function (op::Operator{:∫∫q̇mΨqkΨdx})(a::T;b::S;k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒₁ = a.𝓒; 𝓖₁ = a.𝓖
    𝓒₂ = b.𝓒; 𝓖₂ = b.𝓖
    ρA = op.ρA
    EA = op.EA
    for (ξ₁,ξ₂) in zip(𝓖₁,𝓖₂)
        Bₓ = ξ₁[:∂𝝭∂x]
        Bₜ = ξ₁[:∂𝝭∂y]
        Ψₜ = ξ₂[:∂𝝭∂y]
        Ψₓ = ξ₂[:∂𝝭∂x]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒₁)
            I = xᵢ.𝐼
            for (k,xₖ) in enumerate(𝓒₂)
                K = xₖ.𝐼
                k[I,K] += (Bₜ[i]*ρA*Ψₜ[k] - Bₓ[i]*EA*Ψₓ[k])*𝑤
            end
        end
    end
end

# f[I] += -(ρ*A*Bₜ[i])*N[j] + ((N[i]*b) + N[j]*P)*𝑤
#
function (op::Operator{:∫𝑃δudx})(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        𝑃 = ξ.𝑃
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] -= N[i]*𝑃*𝑤
        end
    end
end

function (op::Operator{:∫qmpdΩ})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    ρA = op.ρA
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += N[i]*ρA*N[j]*𝑤
            end    
        end
    end
end

function (op::Operator{:∫qkpdΩ})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    EA = op.EA
    for ξ in 𝓖
        Bₓ = ξ[:∂𝝭∂x]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] +=  Bₓ[i]*EA*Bₓ[j]*𝑤
            end    
        end
    end
end