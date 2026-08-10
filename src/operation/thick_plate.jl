module MindlinPlate 
    
using ..ApproxOperator: AbstractElement, Element

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
    Dᵢᵢᵢᵢ = E/12/(1-ν^2)
    Dᵢᵢⱼⱼ = E*ν/12/(1-ν^2)
    Dᵢⱼᵢⱼ = E/24/(1+ν)
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
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += N[i]*N[j]*𝑤
            end
        end
    end
end

function ∫∇w∇wdΩ(ap::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = ap.E
    ν = ap.ν
    h = ap.h
    D = 5/6/h^2*E/2/(1+ν)
    for ξ in 𝓖
        𝑤 = ξ.𝑤
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

# function ∫wwdΩ(ap::T,k::AbstractMatrix) where T<:AbstractElement
#     𝓒 = ap.𝓒; 𝓖 = ap.𝓖
#     for ξ in 𝓖
#         𝑤 = ξ.𝑤
#         N = ξ[:𝝭]
#         for (i,xᵢ) in enumerate(𝓒)
#             I = xᵢ.𝐼 
#             for (j,xⱼ) in enumerate(𝓒)
#                 J = xⱼ.𝐼
#                 k[I,J] += N[i]*N[j]*𝑤
#             end
#         end
#     end
# end

function ∫φφdΩ(ap::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = ap.E
    ν = ap.ν
    h = ap.h
    D = 5/6/h^2*E/2/(1+ν)
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        𝑤 = ξ.𝑤
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

function ∫φwdΩ(a₁::T,a₂::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒₁ = a₁.𝓒; 𝓖₁ = a₁.𝓖
    𝓒₂ = a₂.𝓒; 𝓖₂ = a₂.𝓖
    E = a₁.E
    ν = a₁.ν
    h = a₁.h
    D = 5/6*h*E/2/(1+ν)
    for (ξ₁,ξ₂) in zip(𝓖₁,𝓖₂)
        𝑤 = ξ₁.𝑤
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

function ∫φwdΩ(ap::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = ap.E
    ν = ap.ν
    h = ap.h
    D = 5/6/h^2*E/2/(1+ν)
    for ξ in 𝓖
        𝑤 = ξ.𝑤
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
    E = a₁.E
    ν = a₁.ν
    h = a₁.h
    D = 5/6/h^2*E/2/(1+ν)
    for (ξ₁,ξ₂) in zip(𝓖₁,𝓖₂)
        𝑤 = ξ₁.𝑤
        B₁ = ξ₁[:∂𝝭∂x]
        B₂ = ξ₁[:∂𝝭∂y]
        N = ξ₂[:𝝭]
        for (i,xᵢ) in enumerate(𝓒₁)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,J] -= D*N[i]*B₁[j]*𝑤
                k[2*I,J]   -= D*N[i]*B₂[j]*𝑤
            end
        end
    end
end

function ∫κκdΩBui(ap::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = ap.E
    ν = ap.ν
    h = ap.h
    Dᵢᵢᵢᵢ = E*h^3/12/(1-ν^2)
    Dᵢᵢⱼⱼ = E*ν*h^3/12/(1-ν^2)
    Dᵢⱼᵢⱼ = E*h^3/24/(1+ν)
    for ξ in 𝓖
        B₁₁ = ξ[:∂²𝝭∂x²]
        B₁₂ = ξ[:∂²𝝭∂x∂y]
        B₂₂ = ξ[:∂²𝝭∂y²]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += (Dᵢᵢᵢᵢ*B₁₁[i]*B₁₁[j] + Dᵢⱼᵢⱼ*B₁₂[i]*B₁₂[j])*𝑤
                k[2*I-1,2*J]   += (Dᵢᵢⱼⱼ*B₁₁[i]*B₂₂[j] + Dᵢⱼᵢⱼ*B₁₂[i]*B₁₂[j])*𝑤
                k[2*I,2*J-1]   += (Dᵢᵢⱼⱼ*B₂₂[i]*B₁₁[j] + Dᵢⱼᵢⱼ*B₁₂[i]*B₁₂[j])*𝑤
                k[2*I,2*J]     += (Dᵢᵢᵢᵢ*B₂₂[i]*B₂₂[j] + Dᵢⱼᵢⱼ*B₁₂[i]*B₁₂[j])*𝑤
            end
        end
    end
end

function ∫ψψdΩBui(ap::T,k::AbstractMatrix) where T<:AbstractElement
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
                k[2*I-1,2*J-1] += D*B₁[i]*B₁[j]*𝑤
                k[2*I,2*J]     += D*B₂[i]*B₂[j]*𝑤
            end
        end
    end
end

function ∫ψwdΩBui(ap::T,k::AbstractMatrix) where T<:AbstractElement
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
                k[2*I-1,J] += D*B₁[i]*B₁[j]*𝑤
                k[2*I,J]   += D*B₂[i]*B₂[j]*𝑤
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
    D = 5/6/h^2*E/2/(1+ν)
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
        C⁻¹ᵢᵢᵢᵢ = 1/E*12
        C⁻¹ᵢᵢⱼⱼ = -ν/E*12
        C⁻¹ᵢⱼᵢⱼ = 2*(1+ν)/E*12
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

function ∫M∇φdΩ(aₘ::T,aᵩ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₘ = aₘ.𝓒;𝓖ₘ = aₘ.𝓖
    𝓒ᵩ = aᵩ.𝓒;𝓖ᵩ = aᵩ.𝓖
    for (ξₘ,ξᵩ) in zip(𝓖ₘ,𝓖ᵩ)
        𝑤 = ξₘ.𝑤
        N = ξₘ[:𝝭]
        B₁ = ξᵩ[:∂𝝭∂x]
        B₂ = ξᵩ[:∂𝝭∂y]
        for (i,xᵢ) in enumerate(𝓒ₘ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵩ)
                J = xⱼ.𝐼
                k[3*I-2,2*J-1] -= N[i]*B₁[j]*𝑤
                k[3*I-1,2*J]   -= N[i]*B₂[j]*𝑤
                k[3*I,2*J-1]   -= N[i]*B₂[j]*𝑤
                k[3*I,2*J]     -= N[i]*B₁[j]*𝑤
            end
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

function ∫αwwdΓ(ap::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    α = ap.α
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += α*N[i]*N[j]*𝑤
            end
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

function ∫αφφdΓ(ap::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    α = ap.α
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
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
        end
    end
end

function ∫αwwdΓ(ap::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    α = ap.α
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += α*N[i]*N[j]*𝑤
            end
        end
    end
end

function ∫αwwdΓ(ap::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    α = ap.α
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += α*N[i]*N[j]*𝑤
            end
        end
    end
end

function ∫αφφdΓ(ap::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    α = ap.α
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
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
        end
    end
end

function ∫∇w∇wdΩ_MITC(ap::Element{:Quad4},k::AbstractMatrix)
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = ap.E
    ν = ap.ν
    h = ap.h
    D = 5/6/h^2*E/2/(1+ν)

    x₁,x₂,x₃,x₄ = (xᵢ.x for xᵢ in 𝓒)
    y₁,y₂,y₃,y₄ = (xᵢ.y for xᵢ in 𝓒)
    ∂Nᵃ∂ξ = (-0.5, 0.5, 0.0, 0.0)
    ∂Nᵇ∂η = ( 0.0,-0.5, 0.5, 0.0)
    ∂Nᶜ∂ξ = ( 0.0, 0.0, 0.5,-0.5)
    ∂Nᵈ∂η = (-0.5, 0.0, 0.0, 0.5)
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        ξ̄ = ξ.ξ
        η̄ = ξ.η
        ∂N₁∂ξ = - 0.25*(1-η̄)
        ∂N₂∂ξ =   0.25*(1-η̄)
        ∂N₃∂ξ =   0.25*(1+η̄)
        ∂N₄∂ξ = - 0.25*(1+η̄)
        ∂N₁∂η = - 0.25*(1-ξ̄)
        ∂N₂∂η = - 0.25*(1+ξ̄)
        ∂N₃∂η =   0.25*(1+ξ̄)
        ∂N₄∂η =   0.25*(1-ξ̄)
        ∂x∂ξ = ∂N₁∂ξ*x₁ + ∂N₂∂ξ*x₂ + ∂N₃∂ξ*x₃ + ∂N₄∂ξ*x₄
        ∂x∂η = ∂N₁∂η*x₁ + ∂N₂∂η*x₂ + ∂N₃∂η*x₃ + ∂N₄∂η*x₄
        ∂y∂ξ = ∂N₁∂ξ*y₁ + ∂N₂∂ξ*y₂ + ∂N₃∂ξ*y₃ + ∂N₄∂ξ*y₄
        ∂y∂η = ∂N₁∂η*y₁ + ∂N₂∂η*y₂ + ∂N₃∂η*y₃ + ∂N₄∂η*y₄
        detJ = ∂x∂ξ*∂y∂η - ∂x∂η*∂y∂ξ
        ∂ξ∂x =   ∂y∂η/detJ
        ∂η∂x = - ∂y∂ξ/detJ
        ∂ξ∂y = - ∂x∂η/detJ
        ∂η∂y =   ∂x∂ξ/detJ
        ∂N̄∂ξ = 0.5*(1-η̄).*∂Nᵃ∂ξ .+ 0.5*(1+η̄).*∂Nᶜ∂ξ
        ∂N̄∂η = 0.5*(1-ξ̄).*∂Nᵈ∂η .+ 0.5*(1+ξ̄).*∂Nᵇ∂η
        B̄₁ = ∂ξ∂x.*∂N̄∂ξ .+ ∂η∂x.*∂N̄∂η
        B̄₂ = ∂ξ∂y.*∂N̄∂ξ .+ ∂η∂y.*∂N̄∂η
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += D*(B̄₁[i]*B̄₁[j] + B̄₂[i]*B̄₂[j])*𝑤
            end
        end
    end
end


function ∫φφdΩ_MITC(ap::Element{:Quad4},k::AbstractMatrix)
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = ap.E
    ν = ap.ν
    h = ap.h
    D = 5/6/h^2*E/2/(1+ν)
    Nᵃ = (0.5,0.5,0.0,0.0)
    Nᵇ = (0.0,0.5,0.5,0.0)
    Nᶜ = (0.0,0.0,0.5,0.5)
    Nᵈ = (0.5,0.0,0.0,0.5)

    for ξ in 𝓖
        𝑤 = ξ.𝑤
        ξ̄ = ξ.ξ
        η̄ = ξ.η
        N̄₁ = 0.5*(1-η̄).*Nᵃ .+ 0.5*(1+η̄).*Nᶜ
        N̄₂ = 0.5*(1-ξ̄).*Nᵈ .+ 0.5*(1+ξ̄).*Nᵇ
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += D*N̄₁[i]*N̄₁[j]*𝑤
                k[2*I,2*J]     += D*N̄₂[i]*N̄₂[j]*𝑤
            end
        end
    end
end

function ∫φwdΩ_MITC(ap::Element{:Quad4},k::AbstractMatrix)
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = ap.E
    ν = ap.ν
    h = ap.h
    D = 5/6/h^2*E/2/(1+ν)

    x₁,x₂,x₃,x₄ = (xᵢ.x for xᵢ in 𝓒)
    y₁,y₂,y₃,y₄ = (xᵢ.y for xᵢ in 𝓒)
    Nᵃ = (0.5,0.5,0.0,0.0)
    Nᵇ = (0.0,0.5,0.5,0.0)
    Nᶜ = (0.0,0.0,0.5,0.5)
    Nᵈ = (0.5,0.0,0.0,0.5)
    ∂Nᵃ∂ξ = (-0.5, 0.5, 0.0, 0.0)
    ∂Nᵇ∂η = ( 0.0,-0.5, 0.5, 0.0)
    ∂Nᶜ∂ξ = ( 0.0, 0.0, 0.5,-0.5)
    ∂Nᵈ∂η = (-0.5, 0.0, 0.0, 0.5)
 
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        ξ̄ = ξ.ξ
        η̄ = ξ.η
        N̄₁ = 0.5*(1-η̄).*Nᵃ .+ 0.5*(1+η̄).*Nᶜ
        N̄₂ = 0.5*(1-ξ̄).*Nᵈ .+ 0.5*(1+ξ̄).*Nᵇ
        ∂N₁∂ξ = - 0.25*(1-η̄)
        ∂N₂∂ξ =   0.25*(1-η̄)
        ∂N₃∂ξ =   0.25*(1+η̄)
        ∂N₄∂ξ = - 0.25*(1+η̄)
        ∂N₁∂η = - 0.25*(1-ξ̄)
        ∂N₂∂η = - 0.25*(1+ξ̄)
        ∂N₃∂η =   0.25*(1+ξ̄)
        ∂N₄∂η =   0.25*(1-ξ̄)
        ∂x∂ξ = ∂N₁∂ξ*x₁ + ∂N₂∂ξ*x₂ + ∂N₃∂ξ*x₃ + ∂N₄∂ξ*x₄
        ∂x∂η = ∂N₁∂η*x₁ + ∂N₂∂η*x₂ + ∂N₃∂η*x₃ + ∂N₄∂η*x₄
        ∂y∂ξ = ∂N₁∂ξ*y₁ + ∂N₂∂ξ*y₂ + ∂N₃∂ξ*y₃ + ∂N₄∂ξ*y₄
        ∂y∂η = ∂N₁∂η*y₁ + ∂N₂∂η*y₂ + ∂N₃∂η*y₃ + ∂N₄∂η*y₄
        detJ = ∂x∂ξ*∂y∂η - ∂x∂η*∂y∂ξ
        ∂ξ∂x =   ∂y∂η/detJ
        ∂η∂x = - ∂y∂ξ/detJ
        ∂ξ∂y = - ∂x∂η/detJ
        ∂η∂y =   ∂x∂ξ/detJ
        ∂N̄∂ξ = 0.5*(1-η̄).*∂Nᵃ∂ξ .+ 0.5*(1+η̄).*∂Nᶜ∂ξ
        ∂N̄∂η = 0.5*(1-ξ̄).*∂Nᵈ∂η .+ 0.5*(1+ξ̄).*∂Nᵇ∂η
        B̄₁ = ∂ξ∂x.*∂N̄∂ξ .+ ∂η∂x.*∂N̄∂η
        B̄₂ = ∂ξ∂y.*∂N̄∂ξ .+ ∂η∂y.*∂N̄∂η
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,J] -= D*N̄₁[i]*B̄₁[j]*𝑤
                k[2*I,J]   -= D*N̄₂[i]*B̄₂[j]*𝑤
            end
        end
    end
end

function ∫φφdΩ_DSG(ap::Element{:Tri3},k::AbstractMatrix)
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = ap.E
    ν = ap.ν
    h = ap.h
    D = 5/6/h^2*E/2/(1+ν)
    x₁,x₂,x₃ = (xᵢ.x for xᵢ in 𝓒)
    y₁,y₂,y₃ = (xᵢ.y for xᵢ in 𝓒)
    𝐽 = ap.𝐽
    𝑤 = 𝐽/2
    a = x₂-x₁
    b = y₂-y₁
    c = y₃-y₁
    d = x₃-x₁
    B₁₁ = 0.5/𝐽.*(  𝐽, a*c,-b*d)
    B₁₂ = 0.5/𝐽.*(0.0, b*c,-b*c)
    B₂₁ = 0.5/𝐽.*(0.0,-a*d, a*d)
    B₂₂ = 0.5/𝐽.*(  𝐽,-b*d, a*c)

    for (i,xᵢ) in enumerate(𝓒)
        I = xᵢ.𝐼
        for (j,xⱼ) in enumerate(𝓒)
            J = xⱼ.𝐼
            k[2*I-1,2*J-1] += D*(B₁₁[i]*B₁₁[j]+B₂₁[i]*B₂₁[j])*𝑤
            k[2*I-1,2*J]   += D*(B₁₁[i]*B₁₂[j]+B₂₁[i]*B₂₂[j])*𝑤
            k[2*I,2*J-1]   += D*(B₁₂[i]*B₁₁[j]+B₂₂[i]*B₂₁[j])*𝑤
            k[2*I,2*J]     += D*(B₁₂[i]*B₁₂[j]+B₂₂[i]*B₂₂[j])*𝑤
        end
    end
end

function ∫φwdΩ_DSG(ap::Element{:Tri3},k::AbstractMatrix)
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = ap.E
    ν = ap.ν
    h = ap.h
    D = 5/6/h^2*E/2/(1+ν)
    x₁,x₂,x₃ = (xᵢ.x for xᵢ in 𝓒)
    y₁,y₂,y₃ = (xᵢ.y for xᵢ in 𝓒)
    𝐽 = ap.𝐽
    a = x₂-x₁
    b = y₂-y₁
    c = y₃-y₁
    d = x₃-x₁
    Bᵠ₁₁ = 0.5/𝐽.*(  𝐽, a*c,-b*d)
    Bᵠ₁₂ = 0.5/𝐽.*(0.0, b*c,-b*c)
    Bᵠ₂₁ = 0.5/𝐽.*(0.0,-a*d, a*d)
    Bᵠ₂₂ = 0.5/𝐽.*(  𝐽,-b*d, a*c)
 
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,J] -= D*(Bᵠ₁₁[i]*B₁[j]+Bᵠ₂₁[i]*B₂[j])*𝑤
                k[2*I,J]   -= D*(Bᵠ₁₂[i]*B₁[j]+Bᵠ₂₂[i]*B₂[j])*𝑤
            end
        end
    end
end

function L₂Q(ap::T) where T<:AbstractElement
    ΔQ²= BigFloat(0.0)
    Q̄² = BigFloat(0.0)
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        Q̄₁ = ξ.Q₁
        Q̄₂ = ξ.Q₂
        Q₁ = 0.0
        Q₂ = 0.0
        for (i,xᵢ) in enumerate(ap.𝓒)
            Q₁ += N[i]*xᵢ.q₁
            Q₂ += N[i]*xᵢ.q₂
        end
        ΔQ² +=((Q₁ - Q̄₁)^2 + (Q₂ - Q̄₂)^2)*𝑤
        Q̄²  += (Q̄₁^2 + Q̄₂^2)*𝑤
    end
    return ΔQ², Q̄²
end

function L₂Q(aps::Vector{T}) where T<:AbstractElement
    L₂Norm_ΔQ²= BigFloat(0.0)
    L₂Norm_Q̄² = BigFloat(0.0)
    for ap in aps
        ΔQ², Q̄² = L₂Q(ap)
        L₂Norm_ΔQ² += ΔQ²
        L₂Norm_Q̄²  += Q̄²
    end
    return (L₂Norm_ΔQ²/L₂Norm_Q̄²)^0.5
end

function L₂γ(ap::T) where T<:AbstractElement
    Δγ²= BigFloat(0.0)
    γ̄² = BigFloat(0.0)
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        γ̄₁ = ξ.γ₁
        γ̄₂ = ξ.γ₂
        γ₁ = 0.0
        γ₂ = 0.0
        for (i,xᵢ) in enumerate(ap.𝓒)
            γ₁ += B₁[i]*xᵢ.d - N[i]*xᵢ.d₁
            γ₂ += B₂[i]*xᵢ.d - N[i]*xᵢ.d₂
        end
        Δγ² += ((γ₁ - γ̄₁)^2 + (γ₂ - γ̄₂)^2)*𝑤
        γ̄²  += (γ̄₁^2 + γ̄₂^2)*𝑤
    end
    return Δγ², γ̄²
end

function L₂γ(aps::Vector{T}) where T<:AbstractElement
    L₂Norm_Δγ²= BigFloat(0.0)
    L₂Norm_γ̄² = BigFloat(0.0)
    for ap in aps
        Δγ², γ̄² = L₂γ(ap)
        L₂Norm_Δγ² += Δγ²
        L₂Norm_γ̄²  += γ̄²
    end
    return (L₂Norm_Δγ²/L₂Norm_γ̄²)^0.5
end

function L₂γ(a::T,b::S) where {T,S<:AbstractElement}
    Δγ²= BigFloat(0.0)
    γ̄² = BigFloat(0.0)
    for (ξʷ,ξᵠ) in zip(a.𝓖,b.𝓖)
        𝑤 = ξʷ.𝑤
        N = ξᵠ[:𝝭]
        B₁ = ξʷ[:∂𝝭∂x]
        B₂ = ξʷ[:∂𝝭∂y]
        γ̄₁ = ξʷ.γ₁
        γ̄₂ = ξʷ.γ₂
        γ₁ = 0.0
        γ₂ = 0.0
        for (i,xᵢ) in enumerate(a.𝓒)
            γ₁ += B₁[i]*xᵢ.d
            γ₂ += B₂[i]*xᵢ.d
        end
        for (i,xᵢ) in enumerate(b.𝓒)
            γ₁ -= N[i]*xᵢ.d₁
            γ₂ -= N[i]*xᵢ.d₂
        end
        Δγ² += ((γ₁ - γ̄₁)^2 + (γ₂ - γ̄₂)^2)*𝑤
        γ̄²  += (γ̄₁^2 + γ̄₂^2)*𝑤
    end
    return Δγ², γ̄²
end

function L₂γ(as::Vector{T},bs::Vector{S}) where {T,S<:AbstractElement}
    L₂Norm_Δγ²= BigFloat(0.0)
    L₂Norm_γ̄² = BigFloat(0.0)
    for (a,b) in zip(as,bs)
        Δγ², γ̄² = L₂γ(a,b)
        L₂Norm_Δγ² += Δγ²
        L₂Norm_γ̄²  += γ̄²
    end
    return (L₂Norm_Δγ²/L₂Norm_γ̄²)^0.5
end


function L₂φ(ap::T) where T<:AbstractElement
    Δφ²= BigFloat(0.0)
    φ̄² = BigFloat(0.0)
    Δφ²= BigFloat(0.0)
    φ̄² = BigFloat(0.0)
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        φ̄₁ = ξ.φ₁
        φ̄₂ = ξ.φ₂
        φ₁ = 0.0
        φ₂ = 0.0
        φ₁ = 0.0
        φ₂ = 0.0
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
    L₂Norm_Δφ²= BigFloat(0.0)
    L₂Norm_φ̄² = BigFloat(0.0)
    for ap in aps
        Δφ², φ̄² = L₂φ(ap)
        L₂Norm_Δφ² += Δφ²
        L₂Norm_φ̄²  += φ̄²
    end
    return (L₂Norm_Δφ²/L₂Norm_φ̄²)^0.5
end

function H₁φ(ap::T) where T<:AbstractElement
    Δ∇φ²= 0.0
    ∇φ̄² = 0.0
    Δφ²= 0.0
    φ̄² = 0.0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        φ̄₁ = ξ.φ₁
        ∂φ̄₁∂x = ξ.∂φ₁∂x
        ∂φ̄₁∂y = ξ.∂φ₁∂y
        φ̄₂ = ξ.φ₂
        ∂φ̄₂∂x = ξ.∂φ₂∂x
        ∂φ̄₂∂y = ξ.∂φ₂∂y
        φ₁ = 0.0
        ∂φ₁∂x = 0.0
        ∂φ₁∂y = 0.0
        φ₂ = 0.0
        ∂φ₂∂x = 0.0
        ∂φ₂∂y = 0.0
        for (i,xᵢ) in enumerate(ap.𝓒)
            φ₁ += N[i]*xᵢ.d₁
            ∂φ₁∂x += B₁[i]*xᵢ.d₁
            ∂φ₁∂y += B₂[i]*xᵢ.d₁
            φ₂ += N[i]*xᵢ.d₂
            ∂φ₂∂x += B₁[i]*xᵢ.d₂
            ∂φ₂∂y += B₂[i]*xᵢ.d₂
        end
        Δ∇φ² += ((∂φ₁∂x - ∂φ̄₁∂x)^2 + (∂φ₁∂y - ∂φ̄₁∂y)^2 + (∂φ₂∂x - ∂φ̄₂∂x)^2 + (∂φ₂∂y - ∂φ̄₂∂y)^2)*𝑤
        ∇φ̄²  += (∂φ̄₁∂x^2 + ∂φ̄₁∂y^2 + ∂φ̄₂∂x^2 + ∂φ̄₂∂y^2)*𝑤
        Δφ² +=((φ₁ - φ̄₁)^2 + (φ₂ - φ̄₂)^2)*𝑤
        φ̄²  += (φ̄₁^2 + φ̄₂^2)*𝑤
    end
    return Δ∇φ², ∇φ̄², Δφ², φ̄²
end

function L₂w(ap::T) where T<:AbstractElement
    Δw²= BigFloat(0.0)
    w̄² = BigFloat(0.0)
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        w̄ = ξ.w
        w = 0.0
        for (i,xᵢ) in enumerate(ap.𝓒)
            w += N[i]*xᵢ.d
        end
        Δw² +=(w - w̄)^2*𝑤
        w̄²  += w̄^2*𝑤
    end
    return Δw², w̄²
end

function H₁w(ap::T) where T<:AbstractElement
    Δ∇w²= BigFloat(0.0)
    ∇w̄² = BigFloat(0.0)
    Δw²= BigFloat(0.0)
    w̄² = BigFloat(0.0)
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        w̄ = ξ.w
        ∂w̄∂x = ξ.∂w∂x
        ∂w̄∂y = ξ.∂w∂y
        w = 0.0
        ∂w∂x = 0.0
        ∂w∂y = 0.0
        for (i,xᵢ) in enumerate(ap.𝓒)
            w += N[i]*xᵢ.d
            ∂w∂x += B₁[i]*xᵢ.d
            ∂w∂y += B₂[i]*xᵢ.d
        end
        Δw² += (w - w̄)^2*𝑤
        w̄²  += w̄^2*𝑤
        Δ∇w² += ((∂w∂x - ∂w̄∂x)^2 + (∂w∂y - ∂w̄∂y)^2)*𝑤
        ∇w̄² += (∂w̄∂x^2 + ∂w̄∂y^2)*𝑤
    end
    return Δ∇w², ∇w̄², Δw², w̄²
end

function L₂w(aps::Vector{T}) where T<:AbstractElement
    L₂Norm_Δw²= BigFloat(0.0)
    L₂Norm_w̄² = BigFloat(0.0)
    for ap in aps
        Δw², w̄² = L₂w(ap)
        L₂Norm_Δw² += Δw²
        L₂Norm_w̄²  += w̄²
    end
    return (L₂Norm_Δw²/L₂Norm_w̄²)^0.5
end

function H₁(aps_w::Vector{Tʷ},aps_φ::Vector{Tᵠ}) where {Tʷ,Tᵠ<:AbstractElement}
    H₁Norm_Δu²= 0.0
    H₁Norm_ū² = 0.0
    L₂Norm_Δu²= 0.0
    L₂Norm_ū² = 0.0
    for ap in aps_w
        Δ∇w², ∇w̄², Δw², w̄² = H₁w(ap)
        H₁Norm_Δu² += Δ∇w²
        H₁Norm_ū²  += ∇w̄²
        L₂Norm_Δu² += Δw²
        L₂Norm_ū²  += w̄²
    end
    for ap in aps_φ
        Δ∇φ², ∇φ̄², Δφ², φ̄² = H₁φ(ap)
        H₁Norm_Δu² += Δ∇φ²
        H₁Norm_ū²  += ∇φ̄²
        L₂Norm_Δu² += Δφ²
        L₂Norm_ū²  += φ̄²
    end
    return (H₁Norm_Δu²/H₁Norm_ū²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end

function Hₑ(aps::Vector{T}) where T<:AbstractElement
    HₑNorm_Δu² = 0.0
    HₑNorm_ū² = 0.0
    for ap in aps
        Δu², ū² = Hₑ(ap)
        HₑNorm_Δu² += Δu²
        HₑNorm_ū² += ū²
    end
    return (HₑNorm_Δu²/HₑNorm_ū²)^0.5
end

function Hₑ(ap::T) where T<:AbstractElement
    Δu² = BigFloat(0.0)
    ū² = BigFloat(0.0)
    E = ap.E
    ν = ap.ν
    h = ap.h
    Dᵇᵢᵢᵢᵢ = E*h^3/12/(1-ν^2)
    Dᵇᵢᵢⱼⱼ = E*ν*h^3/12/(1-ν^2)
    Dᵇᵢⱼᵢⱼ = E*h^3/24/(1+ν)
    Dˢ =  5/6*h*E/2/(1+ν)
    for ξ in ap.𝓖
        κ̄₁₁ = -ξ.∂φ₁∂x
        κ̄₂₂ = -ξ.∂φ₂∂y
        κ̄₁₂ = -ξ.∂φ₁∂y-ξ.∂φ₂∂x
        γ̄₁ = ξ.∂w∂x - ξ.φ₁
        γ̄₂ = ξ.∂w∂y - ξ.φ₂
        𝑤 = ξ.𝑤

        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        κ₁₁ = 0.0
        κ₂₂ = 0.0
        κ₁₂ = 0.0
        γ₁ = 0.0
        γ₂ = 0.0
        for (i,xᵢ) in enumerate(ap.𝓒)
            γ₁ += B₁[i]*xᵢ.d - N[i]*xᵢ.d₁
            γ₂ += B₂[i]*xᵢ.d - N[i]*xᵢ.d₂
            κ₁₁ -= B₁[i]*xᵢ.d₁
            κ₂₂ -= B₂[i]*xᵢ.d₂
            κ₁₂ -= B₁[i]*xᵢ.d₂+B₂[i]*xᵢ.d₁
        end
        Δu² += 0.5*(Dᵇᵢᵢᵢᵢ*(κ₁₁-κ̄₁₁)^2+Dᵇᵢᵢᵢᵢ*(κ₂₂-κ̄₂₂)^2+2*Dᵇᵢᵢⱼⱼ*(κ₁₁-κ̄₁₁)*(κ₂₂-κ̄₂₂)+2*Dᵇᵢⱼᵢⱼ*(κ₁₂-κ̄₁₂)^2 + ((γ₁-γ̄₁)^2 + (γ₂-γ̄₂)^2)*Dˢ)*𝑤
        ū² += 0.5*(Dᵇᵢᵢᵢᵢ*κ̄₁₁^2+Dᵇᵢᵢᵢᵢ*κ̄₂₂^2+2*Dᵇᵢᵢⱼⱼ*κ̄₁₁*κ̄₂₂+2*Dᵇᵢⱼᵢⱼ*κ̄₁₂^2 + (γ̄₁^2 + γ̄₂^2)*Dˢ)*𝑤
    end
    return Δu², ū²
end

function Hₑ(aps_w::Vector{Tʷ},aps_φ::Vector{Tᵠ},aps_Q::Vector{Tˢ}) where {Tʷ,Tᵠ,Tˢ<:AbstractElement}
    HₑNorm_Δu² = 0.0
    HₑNorm_ū² = 0.0
    for (ap_w,ap_φ,ap_Q) in zip(aps_w,aps_φ,aps_Q)
        Δu², ū² = Hₑ(ap_w,ap_φ,ap_Q)
        HₑNorm_Δu² += Δu²
        HₑNorm_ū² += ū²
    end
    return (HₑNorm_Δu²/HₑNorm_ū²)^0.5
end

function Hₑ(ap_w::Tʷ,ap_φ::Tᵠ,ap_Q::Tˢ) where {Tʷ,Tᵠ,Tˢ<:AbstractElement}
    Δu² = BigFloat(0.0)
    ū² = BigFloat(0.0)
    E = ap_φ.E
    ν = ap_φ.ν
    h = ap_φ.h
    Dᵇᵢᵢᵢᵢ = E*h^3/12/(1-ν^2)
    Dᵇᵢᵢⱼⱼ = E*ν*h^3/12/(1-ν^2)
    Dᵇᵢⱼᵢⱼ = E*h^3/24/(1+ν)
    Dˢ =  5/6*h*E/2/(1+ν)
    for (ξʷ,ξᵠ,ξˢ) in zip(ap_w.𝓖,ap_φ.𝓖,ap_Q.𝓖)
        κ̄₁₁ = -ξᵠ.∂φ₁∂x
        κ̄₂₂ = -ξᵠ.∂φ₂∂y
        κ̄₁₂ = -ξᵠ.∂φ₁∂y-ξᵠ.∂φ₂∂x
        γ̄₁ = ξʷ.∂w∂x - ξᵠ.φ₁
        γ̄₂ = ξʷ.∂w∂y - ξᵠ.φ₂
        q̄₁ = Dˢ*γ̄₁
        q̄₂ = Dˢ*γ̄₂
        q̃₁ = ξˢ.Q₁
        q̃₂ = ξˢ.Q₂
        𝑤 = ξᵠ.𝑤
        if abs(q̄₁-q̃₁)>1e3*eps() || abs(q̄₂-q̃₂)>1e3*eps()
            error("Some problem in exact solution! q̄₁: $q̄₁, q̃₁: $q̃₁, q̄₂: $q̄₂, q̃₂: $q̃₂")
        end

        Nʷ = ξʷ[:𝝭]
        Bʷ₁ = ξʷ[:∂𝝭∂x]
        Bʷ₂ = ξʷ[:∂𝝭∂y]
        Nᵠ = ξᵠ[:𝝭]
        Bᵠ₁ = ξᵠ[:∂𝝭∂x]
        Bᵠ₂ = ξᵠ[:∂𝝭∂y]
        Nˢ = ξˢ[:𝝭]
        κ₁₁ = 0.0
        κ₂₂ = 0.0
        κ₁₂ = 0.0
        γ₁ = 0.0
        γ₂ = 0.0
        q₁ = 0.0
        q₂ = 0.0
        t₁₀₀ = 0.0
        t₁₁₀ = 0.0
        t₁₀₁ = 0.0
        for (i,xᵢ) in enumerate(ap_w.𝓒)
            γ₁ += Bʷ₁[i]*xᵢ.d
            γ₂ += Bʷ₂[i]*xᵢ.d
        end
        for (i,xᵢ) in enumerate(ap_φ.𝓒)
            κ₁₁ -= Bᵠ₁[i]*xᵢ.d₁
            κ₂₂ -= Bᵠ₂[i]*xᵢ.d₂
            κ₁₂ -= Bᵠ₁[i]*xᵢ.d₂+Bᵠ₂[i]*xᵢ.d₁
            γ₁ -= Nᵠ[i]*xᵢ.d₁
            γ₂ -= Nᵠ[i]*xᵢ.d₂
            t₁₀₀ += Bᵠ₁[i]
            t₁₁₀ += Bᵠ₁[i]*xᵢ.t
            t₁₀₁ += Bᵠ₂[i]*xᵢ.t
        end
        # if abs(t₁₀₀) > 1e2*eps() || abs(t₁₁₀-1.0) > 1e2*eps() || abs(t₁₀₁) > 1e2*eps()
        #     error("Some problem in consistency condition! t₁₀₀: $t₁₀₀, t₁₁₀: $t₁₁₀, t₁₀₁: $t₁₀₁")
        # end
        for (i,xᵢ) in enumerate(ap_Q.𝓒)
            q₁ += Nˢ[i]*xᵢ.q₁
            q₂ += Nˢ[i]*xᵢ.q₂
        end

        Δu² += 0.5*(Dᵇᵢᵢᵢᵢ*(κ₁₁-κ̄₁₁)^2+Dᵇᵢᵢᵢᵢ*(κ₂₂-κ̄₂₂)^2+2*Dᵇᵢᵢⱼⱼ*(κ₁₁-κ̄₁₁)*(κ₂₂-κ̄₂₂)+2*Dᵇᵢⱼᵢⱼ*(κ₁₂-κ̄₁₂)^2 + ((γ₁-γ̄₁)^2 + (γ₂-γ̄₂)^2)*Dˢ)*𝑤
        ū² += 0.5*(Dᵇᵢᵢᵢᵢ*κ̄₁₁^2+Dᵇᵢᵢᵢᵢ*κ̄₂₂^2+2*Dᵇᵢᵢⱼⱼ*κ̄₁₁*κ̄₂₂+2*Dᵇᵢⱼᵢⱼ*κ̄₁₂^2 + (γ̄₁^2 + γ̄₂^2)*Dˢ)*𝑤
    end
    return Δu², ū²
end

function ∫∇wσ∇wdΩ(ap::T, k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    h = ap.h
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
                k[I, J] += h*(
                    σ₁₁ * B₁[i] * B₁[j] +
                    σ₂₂ * B₂[i] * B₂[j] +
                    σ₁₂ * (B₁[i] * B₂[j] + B₂[i] * B₁[j])
                ) * 𝑤
            end
        end
    end
end

function ∫∇φσ∇φdΩ(ap::T, k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    h = ap.h
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
                k[2*I-1, 2*J-1] += h^3 / 12 * (
                    σ₁₁ * B₁[i] * B₁[j] +
                    σ₂₂ * B₂[i] * B₂[j] +
                    σ₁₂ * (B₁[i] * B₂[j] + B₂[i] * B₁[j])
                ) * 𝑤
                k[2*I-1, 2*J]   += h^3 / 12 * (
                    σ₁₁ * B₁[i] * B₁[j] +
                    σ₂₂ * B₂[i] * B₂[j] +
                    σ₁₂ * (B₁[i] * B₂[j] + B₂[i] * B₁[j])
                ) * 𝑤
                k[2*I, 2*J-1]   += h^3 / 12 * (
                    σ₁₁ * B₁[i] * B₁[j] +
                    σ₂₂ * B₂[i] * B₂[j] +
                    σ₁₂ * (B₁[i] * B₂[j] + B₂[i] * B₁[j])
                ) * 𝑤
                k[2*I, 2*J]     += h^3 / 12 * (
                    σ₁₁ * B₁[i] * B₁[j] +
                    σ₂₂ * B₂[i] * B₂[j] +
                    σ₁₂ * (B₁[i] * B₂[j] + B₂[i] * B₁[j])
                ) * 𝑤
            end
        end
    end
end

end
