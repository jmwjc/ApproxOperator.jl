function (op::Operator{:∫vᵢθᵢds})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    α = op.α
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        n₁₁ = ξ.n₁₁
        n₂₂ = ξ.n₂₂
        n₁₂ = ξ.n₁₂
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        g₃ = ξ.g₃
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += α*N[i]*n₁₁*N[j]*𝑤
                k[3*I-2,3*J-1] += α*N[i]*n₁₂*N[j]*𝑤
                k[3*I-1,3*J-2] += α*N[i]*n₁₂*N[j]*𝑤
                k[3*I-1,3*J-1] += α*N[i]*n₂₂*N[j]*𝑤
                k[3*I,3*J]     += α*N[i]*N[j]*𝑤
            end
            f[3*I-2] += α*N[i]*(n₁₁*g₁+n₁₂*g₂)*𝑤
            f[3*I-1] += α*N[i]*(n₁₂*g₁+n₂₂*g₂)*𝑤
            f[3*I]   += α*N[i]*g₃*𝑤
        end
    end
end

function (op::Operator{:∫κεγds})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    EI = op.EI
    EA = op.EA
    kGA = op.kGA
    R = op.R
    for ξ in 𝓖
        N = ξ[:𝝭]
        B = ξ[:∂𝝭∂x]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += (N[i]*kGA/R^2*N[j]+B[i]*EA*B[j])*𝑤
                k[3*I-2,3*J-1] += (N[i]*kGA/R*B[j]-B[i]*EA/R*N[j])*𝑤
                k[3*I-2,3*J]   += (-N[i]*kGA/R*N[j])*𝑤
                k[3*I-1,3*J-2] += (B[i]*kGA/R*N[j]-N[i]*EA/R*B[j])*𝑤
                k[3*I-1,3*J-1] += (B[i]*kGA*B[j]+N[i]*EA/R^2*N[j])*𝑤
                k[3*I-1,3*J]   += (-B[i]*kGA*N[j])*𝑤
                k[3*I,3*J-2]   += (-N[i]*kGA/R*N[j])*𝑤
                k[3*I,3*J-1]   += (-N[i]*kGA*B[j])*𝑤
                k[3*I,3*J]     += (B[i]*EI*B[j]+N[i]*kGA*N[j])*𝑤
            end
        end
    end
end