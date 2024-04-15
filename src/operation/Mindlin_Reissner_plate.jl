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
        Dˢ = h*E/2/(1+ν)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += ( Dˢ*B₁[i]*B₁[j] + Dˢ*B₂[i]*B₂[j])*𝑤
                k[3*I-2,3*J-1] += (-Dˢ*N[i]*B₁[j])*𝑤
                k[3*I-2,3*J]   += (-Dˢ*B₂[i]*N[j])*𝑤
                k[3*I-1,3*J-2] += (-Dˢ*N[i]*B₁[j])*𝑤
                k[3*I-1,3*J-1] += ( Dᵇᵢᵢᵢᵢ*B₁[i]*B₁[j] + Dᵇᵢⱼᵢⱼ*B₂[i]*B₂[j] + Dˢ*N[i]*N[j])*𝑤
                k[3*I-1,3*J]   += ( Dᵇᵢᵢⱼⱼ*B₁[i]*B₂[j] + Dᵇᵢⱼᵢⱼ*B₂[i]*B₁[j])*𝑤
                k[3*I,3*J-2]   += (-Dˢ*N[i]*B₂[j])*𝑤
                k[3*I,3*J-1]   += ( Dᵇᵢᵢⱼⱼ*B₂[i]*B₁[j] + Dᵇᵢⱼᵢⱼ*B₁[i]*B₂[j])*𝑤
                k[3*I,3*J]     += ( Dᵇᵢᵢᵢᵢ*B₂[i]*B₂[j] + Dᵇᵢⱼᵢⱼ*B₁[i]*B₁[j] + Dˢ*N[i]*N[j])*𝑤
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
                k[3*I-2,3*J-1] += (-Dˢ*N[i]*B₁[j])*𝑤
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

function (op::Operator{:∫∇M-QdΩ})(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
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
            f[3*I-1] += N[i]*(M₁ᵢᵢ-Q₁)*𝑤
            f[3*I-2] += N[i]*(M₂ᵢᵢ-Q₂)*𝑤
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
            f[3*I-2] += 0
            f[3*I-1] += α*N[i]*θ₁*𝑤
            f[3*I]   += α*N[i]*θ₂*𝑤
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