module MindlinPlate 
    
using ..ApproxOperator: AbstractElement

function (op::Operator{:∫κMγQdΩ})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = op.E
    ν = op.ν
    h = op.h
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

function (op::Operator{:∫κMdΩ})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = op.E
    ν = op.ν
    h = op.h
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        Dᵇᵢᵢᵢᵢ = E*h^3/12/(1-ν^2)
        Dᵇᵢᵢⱼⱼ = E*ν*h^3/12/(1-ν^2)
        Dᵇᵢⱼᵢⱼ = E*h^3/24/(1+ν)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-1,3*J-1] += ( Dᵇᵢᵢᵢᵢ*B₁[i]*B₁[j] + Dᵇᵢⱼᵢⱼ*B₂[i]*B₂[j])*𝑤
                k[3*I-1,3*J]   += ( Dᵇᵢᵢⱼⱼ*B₁[i]*B₂[j] + Dᵇᵢⱼᵢⱼ*B₂[i]*B₁[j])*𝑤
                k[3*I,3*J-1]   += ( Dᵇᵢᵢⱼⱼ*B₂[i]*B₁[j] + Dᵇᵢⱼᵢⱼ*B₁[i]*B₂[j])*𝑤
                k[3*I,3*J]     += ( Dᵇᵢᵢᵢᵢ*B₂[i]*B₂[j] + Dᵇᵢⱼᵢⱼ*B₁[i]*B₁[j])*𝑤
            end
        end
    end
end

function (op::Operator{:∫γQdΩ})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = op.E
    ν = op.ν
    h = op.h
    for ξ in 𝓖
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        Dˢ = 5/6*h*E/2/(1+ν)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += ( Dˢ*B₁[i]*B₁[j] + Dˢ*B₂[i]*B₂[j])*𝑤
                k[3*I-2,3*J-1] += (-Dˢ*B₁[i]*N[j])*𝑤
                k[3*I-2,3*J]   += (-Dˢ*B₂[i]*N[j])*𝑤
                k[3*I-1,3*J-2] += (-Dˢ*N[i]*B₁[j])*𝑤
                k[3*I-1,3*J-1] += ( Dˢ*N[i]*N[j])*𝑤
                k[3*I,3*J-2]   += (-Dˢ*N[i]*B₂[j])*𝑤
                k[3*I,3*J]     += ( Dˢ*N[i]*N[j])*𝑤
            end
        end
    end
end

function (op::Operator{:∫wqdΩ})(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        q = ξ.q
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[3*I-2] += N[i]*q*𝑤
        end
    end
end

function (op::Operator{:∫∇MQdΩ})(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        Q₁ = ξ.Q₁
        Q₂ = ξ.Q₂
        M₁ᵢᵢ = ξ.M₁ᵢᵢ
        M₂ᵢᵢ = ξ.M₂ᵢᵢ
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[3*I-1] -= N[i]*(M₁ᵢᵢ+Q₁)*𝑤
            f[3*I] -= N[i]*(M₂ᵢᵢ+Q₂)*𝑤
        end
    end
end

function (op::Operator{:∫wVdΓ})(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        Q₁ = ξ.Q₁
        Q₂ = ξ.Q₂
        V = Q₁*n₁+Q₂*n₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[3*I-2] += N[i]*V*𝑤
        end
    end
end

function (op::Operator{:∫θMdΓ})(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        M₁₁ = ξ.M₁₁
        M₁₂ = ξ.M₁₂
        M₂₂ = ξ.M₂₂
        M₁ = M₁₁*n₁+M₁₂*n₂
        M₂ = M₁₂*n₁+M₂₂*n₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[3*I-1] -= N[i]*M₁*𝑤
            f[3*I]   -= N[i]*M₂*𝑤
        end
    end
end

function (op::Operator{:∫θM₁dΓ})(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        M₁₁ = ξ.M₁₁
        M₁₂ = ξ.M₁₂
        M₁ = M₁₁*n₁+M₁₂*n₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[3*I-1] -= N[i]*M₁*𝑤
        end
    end
end

function (op::Operator{:∫θM₂dΓ})(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        M₁₂ = ξ.M₁₂
        M₂₂ = ξ.M₂₂
        M₂ = M₁₂*n₁+M₂₂*n₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[3*I]   -= N[i]*M₂*𝑤
        end
    end
end

function (op::Operator{:∫vwdΓ})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    α = op.α
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += α*N[i]*N[j]*𝑤
            end
            f[3*I-2] += α*N[i]*g*𝑤
        end
    end
end
function (op::Operator{:∫vθdΓ})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    α = op.α
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        θ₁ = ξ.θ₁
        θ₂ = ξ.θ₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-1,3*J-1] += α*N[i]*N[j]*𝑤
                k[3*I,3*J]     += α*N[i]*N[j]*𝑤
            end
            f[3*I-1] += α*N[i]*θ₁*𝑤
            f[3*I]   += α*N[i]*θ₂*𝑤
        end
    end
end

function (op::Operator{:∫vθ₁dΓ})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    α = op.α
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        θ₁ = ξ.θ₁
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-1,3*J-1] += α*N[i]*N[j]*𝑤
            end
            f[3*I-1] += α*N[i]*θ₁*𝑤
        end
    end
end

function (op::Operator{:∫vθ₂dΓ})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    α = op.α
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        θ₂ = ξ.θ₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I,3*J]     += α*N[i]*N[j]*𝑤
            end
            f[3*I]   += α*N[i]*θ₂*𝑤
        end
    end
end

function (op::Operator{:∫wQdΩ})(a::T,b::S;k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒₁ = a.𝓒; 𝓖₁ = a.𝓖
    𝓒₂ = b.𝓒; 𝓖₂ = b.𝓖
    for (ξ₁,ξ₂) in zip(𝓖₁,𝓖₂)
        𝑤 = ξ₁.𝑤
        N = ξ₁[:𝝭]
        B₁ = ξ₁[:∂𝝭∂x]
        B₂ = ξ₁[:∂𝝭∂y]
        Ñ = ξ₂[:𝝭]
        for (i,xᵢ) in enumerate(𝓒₁)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒₂)
                J = xⱼ.𝐼
                k[3*I-2,2*J-1] += B₁[i]*Ñ[j]*𝑤
                k[3*I-2,2*J]   += B₂[i]*Ñ[j]*𝑤
                k[3*I-1,2*J-1] -=  N[i]*Ñ[j]*𝑤
                k[3*I,2*J]     -=  N[i]*Ñ[j]*𝑤
            end
        end
    end
end

function (op::Operator{:∫QQdΩ})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = op.E
    ν = op.ν
    h = op.h
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        Dˢ = 5/6*h*E/2/(1+ν)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] -= N[i]*N[j]/Dˢ*𝑤
                k[2*I,2*J]     -= N[i]*N[j]/Dˢ*𝑤
            end
        end
    end
end


function (op::Operator{:L₂_ThickPlate})(ap::T) where T<:AbstractElement
    Δu²= 0
    ū² = 0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        ū₁ = ξ.u
        ū₂ = ξ.θ₁
        ū₃ = ξ.θ₂
        # ū₄ = ξ.Q₁
        # ū₅ = ξ.Q₂
        u₁ = 0
        u₂ = 0
        u₃ = 0
        # u₄ = 0
        # u₅ = 0
        for (i,xᵢ) in enumerate(ap.𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            u₃ += N[i]*xᵢ.d₃
            # u₄ += N[i]*xᵢ.q₁
            # u₅ += N[i]*xᵢ.q₂
        end
        Δu² +=((u₁ - ū₁)^2 + (u₂ - ū₂)^2 + (u₃ - ū₃)^2)*𝑤
        ū²  += (ū₁^2 + ū₂^2 + ū₃^2)*𝑤
        # Δu² +=((u₄ - ū₄)^2 + (u₅ - ū₅)^2)*𝑤
        # ū²  += (ū₄^2 + ū₅^2 )*𝑤
    end
    return Δu², ū²
end

function (op::Operator{:L₂_ThickPlate})(aps::Vector{T}) where T<:AbstractElement
    L₂Norm_Δu²= 0
    L₂Norm_ū² = 0
    for ap in aps
        Δu², ū² = op(ap)
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end

end