function (op::Operator{:∫∫ρvᵢuᵢdxdy})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    ρ = op.ρ
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
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
function (op::Operator{:∫∫εᵢⱼσᵢⱼvᵢbᵢdxdy})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = op.E
    ν = op.ν
    for ξ in 𝓖
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        Cᵢᵢᵢᵢ = E/(1-ν^2)
        Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        b₁ = ξ.b₁
        b₂ = ξ.b₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += (Cᵢᵢᵢᵢ*B₁[i]*B₁[j] + Cᵢⱼᵢⱼ*B₂[i]*B₂[j])*𝑤
                k[2*I-1,2*J]   += (Cᵢᵢⱼⱼ*B₁[i]*B₂[j] + Cᵢⱼᵢⱼ*B₂[i]*B₁[j])*𝑤
                k[2*I,2*J-1]   += (Cᵢᵢⱼⱼ*B₂[i]*B₁[j] + Cᵢⱼᵢⱼ*B₁[i]*B₂[j])*𝑤
                k[2*I,2*J]     += (Cᵢᵢᵢᵢ*B₂[i]*B₂[j] + Cᵢⱼᵢⱼ*B₁[i]*B₁[j])*𝑤
            end
            f[2*I-1] += N[i]*b₁*𝑤
            f[2*I]   += N[i]*b₂*𝑤
        end
    end
end

function (op::Operator{:∫∫εᵢⱼσᵢⱼdxdy})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = op.E
    ν = op.ν
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
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


function (op::Operator{:∫∫εᵛᵢⱼσᵛᵢⱼdxdy})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = op.E
    ν = op.ν
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
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

function (op::Operator{:∫∫εᵈᵢⱼσᵈᵢⱼdxdy})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = op.E
    ν = op.ν
    Cᵈ = E/(1+ν)
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
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

function (op::Operator{:∫∫p∇vdxdy})(apu::T,app::S;k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒u = apu.𝓒
    𝓒p = app.𝓒
    𝓖u = apu.𝓖
    𝓖p = app.𝓖
    for (ξu,ξp) in zip(𝓖u,𝓖p)
        N = ξp[:𝝭]
        B₁ = ξu[:∂𝝭∂x]
        B₂ = ξu[:∂𝝭∂y]
        𝑤 = ξu.𝑤
        for (i,xᵢ) in enumerate(𝓒u)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒p)
                J = xⱼ.𝐼
                k[2*I-1,J] += B₁[i]*N[j]*𝑤
                k[2*I,J]   += B₂[i]*N[j]*𝑤
            end
        end
    end
end

function (op::Operator{:∫∫qpdxdy})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = op.E
    ν = op.ν
    K = E/3/(1-2*ν)
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] -= N[i]*N[j]/K*𝑤
            end
        end
    end
end

function (op::Operator{:∫∫vᵢbᵢdxdy})(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
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

function (op::Operator{:∫vᵢtᵢds})(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
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

function (op::Operator{:g₂})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64},dof::Symbol) where T<:AbstractElement{:Poi1}
    x, = ap.𝓒
    if dof == :d₁
        j = 2*x.𝐼-1
    else
        j = 2*x.𝐼
    end
    g = getproperty(x,dof)
    for i in 1:length(f)
        f[i] -= k[i,j]*g
    end
    k[j,:] .= 0.
    k[:,j] .= 0.
    k[j,j] = 1.
    f[j] = g
end

function (op::Operator{:∫λᵢgᵢds})(ap1::T,ap2::S;g::AbstractMatrix{Float64},q::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    for j in 1:length(ap1.𝓖)
        ξ₁ = ap1.𝓖[j]
        ξ₂ = ap2.𝓖[j]
        𝑤 = ξ₁.𝑤
        N = ξ₁[:𝝭]
        N̄ = ξ₂[:𝝭]
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        n₁₁ = ξ.n₁₁
        n₂₂ = ξ.n₂₂
        n₁₂ = ξ.n₁₂
        for (k,xₖ) in enumerate(ap2.𝓒)
            K = xₖ.𝐼
            for (i,xᵢ) in enumerate(ap1.𝓒)
                I = xᵢ.𝐼
                N̄ₖNᵢ = N̄[k]*N[i]
                g[2*K-1,2*I-1] -= n₁₁*N̄ₖNᵢ*𝑤
                g[2*K-1,2*I]   -= n₁₂*N̄ₖNᵢ*𝑤
                g[2*K,2*I-1]   -= n₁₂*N̄ₖNᵢ*𝑤
                g[2*K,2*I]     -= n₂₂*N̄ₖNᵢ*𝑤
            end
            q[2*K-1] -= N̄[k]*(n₁₁*g₁+n₁₂*g₂)*𝑤
            q[2*K]   -= N̄[k]*(n₁₂*g₁+n₂₂*g₂)*𝑤
        end
    end
end

function (op::Operator{:∫σᵢⱼnⱼgᵢds})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = op.E
    ν = op.ν
    Cᵢᵢᵢᵢ = E/(1-ν^2)
    Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
    Cᵢⱼᵢⱼ = E/2/(1+ν)
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

function (op::Operator{:∫σᵈᵢⱼnⱼgᵢds})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = op.E
    ν = op.ν
    Cᵈ = E/(1+ν)
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
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            Tᵢ₁₁ =  2/3*n₁*B₁[i]+0.5*n₂*B₂[i]
            Tᵢ₁₂ = -1/3*n₁*B₂[i]+0.5*n₂*B₁[i]
            Tᵢ₂₁ = -1/3*n₂*B₁[i]+0.5*n₁*B₂[i]
            Tᵢ₂₂ =  2/3*n₂*B₂[i]+0.5*n₁*B₁[i]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                Tⱼ₁₁ =  2/3*n₁*B₁[j]+0.5*n₂*B₂[j]
                Tⱼ₁₂ = -1/3*n₁*B₂[j]+0.5*n₂*B₁[j]
                Tⱼ₂₁ = -1/3*n₂*B₁[j]+0.5*n₁*B₂[j]
                Tⱼ₂₂ =  2/3*n₂*B₂[j]+0.5*n₁*B₁[j]
                k[2*I-1,2*J-1] -= Cᵈ*(N[i]*(n₁₁*Tⱼ₁₁+n₁₂*Tⱼ₂₁)+(n₁₁*Tᵢ₁₁+n₁₂*Tᵢ₂₁)*N[j])*𝑤
                k[2*I-1,2*J]   -= Cᵈ*(N[i]*(n₁₁*Tⱼ₁₂+n₁₂*Tⱼ₂₂)+(n₁₁*Tᵢ₁₂+n₁₂*Tᵢ₂₂)*N[j])*𝑤
                k[2*I,2*J-1]   -= Cᵈ*(N[i]*(n₁₂*Tⱼ₁₁+n₂₂*Tⱼ₂₁)+(n₁₂*Tᵢ₁₁+n₂₂*Tᵢ₂₁)*N[j])*𝑤
                k[2*I,2*J]     -= Cᵈ*(N[i]*(n₁₂*Tⱼ₁₂+n₂₂*Tⱼ₂₂)+(n₁₂*Tᵢ₁₂+n₂₂*Tᵢ₂₂)*N[j])*𝑤
            end
            f[2*I-1] -= Cᵈ*((n₁₁*Tᵢ₁₁+n₁₂*Tᵢ₂₁)*g₁+(n₁₁*Tᵢ₁₂+n₁₂*Tᵢ₂₂)*g₂)*𝑤
            f[2*I]   -= Cᵈ*((n₁₂*Tᵢ₁₁+n₂₂*Tᵢ₂₁)*g₁+(n₁₂*Tᵢ₁₂+n₂₂*Tᵢ₂₂)*g₂)*𝑤
        end
    end
end

function (op::Operator{:∫σᵢⱼnⱼgᵢvᵢgᵢds})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = op.E
    ν = op.ν
    α = op.α
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        Cᵢᵢᵢᵢ = E/(1-ν^2)
        Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        n₁₁ = ξ.n₁₁
        n₁₂ = ξ.n₁₂
        n₂₂ = ξ.n₂₂
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        n₁ = ξ.n₁
        n₂ = ξ.n₂
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

function (op::Operator{:∫σ̄ᵢⱼnⱼgᵢds})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = op.E
    ν = op.ν
    Cᵢᵢᵢᵢ = E/(1-ν^2)
    Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
    Cᵢⱼᵢⱼ = E/2/(1+ν)
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x_]
        B₂ = ξ[:∂𝝭∂y_]
        n₁₁ = ξ.n₁₁
        n₁₂ = ξ.n₁₂
        n₂₂ = ξ.n₂₂
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        n₁ = ξ.n₁
        n₂ = ξ.n₂
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
                k[2*I-1,2*J-1] += (C₁₁₁*B₁[i]*N[j] + C₁₂₁*B₂[i]*N[j])*𝑤
                k[2*I-1,2*J]   += (C₁₁₂*B₁[i]*N[j] + C₁₂₂*B₂[i]*N[j])*𝑤
                k[2*I,2*J-1]   += (C₁₂₁*B₁[i]*N[j] + C₂₂₁*B₂[i]*N[j])*𝑤
                k[2*I,2*J]     += (C₁₂₂*B₁[i]*N[j] + C₂₂₂*B₂[i]*N[j])*𝑤
            end
            f[2*I-1] += ((C₁₁₁*B₁[i] + C₁₂₁*B₂[i])*g₁ + (C₁₁₂*B₁[i] + C₁₂₂*B₂[i])*g₂)*𝑤
            f[2*I]   += ((C₁₂₁*B₁[i] + C₂₂₁*B₂[i])*g₁ + (C₁₂₂*B₁[i] + C₂₂₂*B₂[i])*g₂)*𝑤
        end
    end
end

function (op::Operator{:∫σ̃̄ᵢⱼnⱼgᵢds})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = op.E
    ν = op.ν
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B̄₁ = ξ[:∂𝝭∂x_]
        B̄₂ = ξ[:∂𝝭∂y_]
        Cᵢᵢᵢᵢ = E/(1-ν^2)
        Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        n₁₁ = ξ.n₁₁
        n₁₂ = ξ.n₁₂
        n₂₂ = ξ.n₂₂
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        C₁₁₁ = Cᵢᵢᵢᵢ*n₁*n₁₁+Cᵢᵢⱼⱼ*n₂*n₁₂
        C₁₁₂ = Cᵢᵢᵢᵢ*n₁*n₁₂+Cᵢᵢⱼⱼ*n₂*n₂₂
        C₂₂₁ = Cᵢᵢⱼⱼ*n₁*n₁₁+Cᵢᵢᵢᵢ*n₂*n₁₂
        C₂₂₂ = Cᵢᵢⱼⱼ*n₁*n₁₂+Cᵢᵢᵢᵢ*n₂*n₂₂
        C₁₂₁ = Cᵢⱼᵢⱼ*(n₁*n₁₂+n₂*n₁₁)
        C₁₂₂ = Cᵢⱼᵢⱼ*(n₂*n₁₂+n₁*n₂₂)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            ΔB₁ = B₁[i]-B̄₁[i]
            ΔB₂ = B₂[i]-B̄₂[i]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] -= (C₁₁₁*(N[i]*B₁[j]+ΔB₁*N[j]) + C₁₂₁*(N[i]*B₂[j]+ΔB₂*N[j]))*𝑤
                k[2*I-1,2*J]   -= (C₁₂₁*N[i]*B₁[j] + C₁₁₂*ΔB₁*N[j] + C₂₂₁*N[i]*B₂[j] + C₁₂₂*ΔB₂*N[j])*𝑤
                k[2*I,2*J-1]   -= (C₁₁₂*N[i]*B₁[j] + C₁₂₁*ΔB₁*N[j] + C₁₂₂*N[i]*B₂[j] + C₂₂₁*ΔB₂*N[j])*𝑤
                k[2*I,2*J]     -= (C₁₂₂*(N[i]*B₁[j]+ΔB₁*N[j]) + C₂₂₂*(N[i]*B₂[j]+ΔB₂*N[j]))*𝑤
            end
            f[2*I-1] -= ((C₁₁₁*ΔB₁+C₁₂₁*ΔB₂)*g₁ + (C₁₁₂*ΔB₁+C₁₂₂*ΔB₂)*g₂)*𝑤
            f[2*I]   -= ((C₁₂₁*ΔB₁+C₂₂₁*ΔB₂)*g₁ + (C₁₂₂*ΔB₁+C₂₂₂*ΔB₂)*g₂)*𝑤
        end
    end
end

function (op::Operator{:∫vᵢgᵢds})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
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
