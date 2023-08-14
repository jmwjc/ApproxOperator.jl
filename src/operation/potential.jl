
function (op::Operator{:∫vudΩ})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += N[i]*N[j]*𝑤
            end
        end
    end
end

function (op::Operator{:∫∇v∇uvbdΩ})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    kᶜ = op.k
    for ξ in 𝓖
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z]
        𝑤 = ξ.𝑤
        b = ξ.b
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += kᶜ*(B₁[i]*B₁[j] + B₂[i]*B₂[j] + B₃[i]*B₃[j])*𝑤
            end
            f[I] += N[i]*b*𝑤
        end
    end
end

function (op::Operator{:∫vₓuₓdx})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    EA = op.EA
    for ξ in 𝓖
        B = ξ[:∂𝝭∂x]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += B[i]*EA*B[j]*𝑤
            end
        end
    end
end

function (op::Operator{:∫∫∇v∇udxdy})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    kᶜ = op.k
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += kᶜ*(B₁[i]*B₁[j] + B₂[i]*B₂[j])*𝑤
            end
        end
    end
end

function (op::Operator{:∫∇v∇udΩ})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    kᶜ = op.k
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += kᶜ*(B₁[i]*B₁[j] + B₂[i]*B₂[j] + B₃[i]*B₃[j])*𝑤
            end
        end
    end
end

function (op::Operator{:∫vbdΩ})(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        b = ξ.b
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] += N[i]*b*𝑤
        end
    end
end

function (op::Operator{:∫vtdΓ})(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        t = ξ.t
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] += N[i]*t*𝑤
        end
    end
end

function (op::Operator{:∫vgdΓ})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
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
                k[I,J] += α*N[i]*N[j]*𝑤
            end
            f[I] += α*N[i]*g*𝑤
        end
    end
end

function (op::Operator{:∫λgdΓ})(ap1::T,ap2::S,g::AbstractMatrix{Float64},q::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    for j in 1:length(ap1.𝓖)
        ξ₁ = ap1.𝓖[j]
        ξ₂ = ap2.𝓖[j]
        𝑤 = ξ₁.𝑤
        N = ξ₁[:𝝭]
        N̄ = ξ₂[:𝝭]
        ḡ = ξ₁.g
        for (k,xₖ) in enumerate(ap2.𝓒)
            K = xₖ.𝐼
            for (i,xᵢ) in enumerate(ap1.𝓒)
                I = xᵢ.𝐼
                g[I,K] -= N[i]*N̄[k]*𝑤
            end
            q[K] -= N̄[k]*ḡ*𝑤
        end
    end
end

function (op::Operator{:∫λₙgdΓ})(ap1::T,ap2::S,g::AbstractMatrix{Float64},q::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    for j in 1:length(ap1.𝓖)
        ξ₁ = ap1.𝓖[j]
        ξ₂ = ap2.𝓖[j]
        𝑤 = ξ₁.𝑤
        N = ξ₁[:𝝭]
        N̄ = ξ₂[:𝝭]
        ḡ = ξ₁.g
        sn = sign(ξ₁.n₁ + ξ₂.n₂)
        for (k,xₖ) in enumerate(ap2.𝓒)
            K = xₖ.𝐼
            for (i,xᵢ) in enumerate(ap1.𝓒)
                I = xᵢ.𝐼
                g[I,K] -= sn*N[i]*N̄[k]*𝑤
            end
            q[K] -= sn*N̄[k]*ḡ*𝑤
        end
    end
end

function (op::Operator{:∫∇𝑛vgdΓ})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    kᶜ = op.k
    α = op.α
    for ξ in 𝓖
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z]
        𝑤 = ξ.𝑤
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        n₃ = ξ.n₃
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (-kᶜ*((B₁[i]*n₁+B₂[i]*n₂+B₃[i]*n₃)*N[j]+N[i]*(B₁[j]*n₁+B₂[j]*n₂+B₃[j]*n₃)) + α*N[i]*N[j])*𝑤
            end
            f[I] += (-kᶜ*(B₁[i]*n₁+B₂[i]*n₂+B₃[i]*n₃) + α*N[i])*g*𝑤
        end
    end
end

function (op::Operator{:∫∇𝑛vgds})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    kᶜ = op.k
    α = op.α
    for ξ in 𝓖
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (-kᶜ*((B₁[i]*n₁+B₂[i]*n₂)*N[j]+N[i]*(B₁[j]*n₁+B₂[j]*n₂)) + α*N[i]*N[j])*𝑤
            end
            f[I] += (-kᶜ*(B₁[i]*n₁+B₂[i]*n₂) + α*N[i])*g*𝑤
        end
    end
end

function (op::Operator{:∫∇̃𝑛vgdx})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    kᶜ = op.k
    for ξ in 𝓖
        N = ξ[:𝝭]
        B = ξ[:∂𝝭∂x]
        𝑤 = ξ.𝑤
        n₁ = ξ.n₁
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] -= kᶜ*(B[i]*n₁*N[j]+N[i]*B[j]*n₁)*𝑤
            end
            f[I] -= kᶜ*B[i]*n₁*g*𝑤
        end
    end
end

function (op::Operator{:∫∇̄𝑛vgdx})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    kᶜ = op.k
    for ξ in 𝓖
        N = ξ[:𝝭]
        B = ξ[:∂𝝭∂x]
        𝑤 = ξ.𝑤
        n₁ = ξ.n₁
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += kᶜ*B[i]*n₁*N[j]*𝑤
            end
            f[I] += kᶜ*B[i]*n₁*g*𝑤
        end
    end
end

function (op::Operator{:∫∇̄𝑛vgdΓ})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    kᶜ = op.k
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z]
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        n₃ = ξ.n₃
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += kᶜ*(B₁[i]*n₁+B₂[i]*n₂+B₃[i]*n₃)*N[j]*𝑤
            end
            f[I] += kᶜ*(B₁[i]*n₁+B₂[i]*n₂+B₃[i]*n₃)*g*𝑤
        end
    end
end

function (op::Operator{:g})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64},dof::Symbol=:d) where T<:AbstractElement{:Poi1}
    x = ap.𝓒[1]
    j = x.𝐼
    g = getproperty(x,dof)
    for i in 1:length(f)
        f[i] -= k[i,j]*g
    end
    k[j,:] .= 0.
    k[:,j] .= 0.
    k[j,j] = 1.
    f[j] = g
end

function (op::Operator{:∫vᵢnᵢuds})(a₁::T,a₂::S;k::AbstractMatrix{Float64}) where {T,S<:AbstractElement}
    𝓖 = zip(a₁.𝓖,a₂.𝓖)
    kᶜ = op.k
    for (ξ₁,ξ₂) in 𝓖
        N = ξ₂[:𝝭]
        B₁ = ξ₁[:∂𝝭∂x]
        B₂ = ξ₁[:∂𝝭∂y]
        𝑤 = ξ₁.𝑤
        n₁ = ξ₁.n₁
        n₂ = ξ₁.n₂
        if ξ₁.𝐺 == 9 || ξ₁.𝐺 == 10
            𝒙₁,𝒙₂ = a₂.𝓒
            x₁ = 𝒙₁.x
            x₂ = 𝒙₂.x
            y₁ = 𝒙₁.y
            y₂ = 𝒙₂.y
            g₁ = 1.0+2x₁+3y₁
            g₂ = 1.0+2x₂+3y₂
            # println(N[1]*g₁ + N[2]*g₂)
            # println(𝑤)
            # println(n₁)
            # println(n₂)
            # println(B₁[1])
            # println(B₁[2])
            # println(B₁[3])
            # println(B₂[1])
            # println(B₂[2])
            # println(B₂[3])
            # println(g₁)
            # println(g₂)
            # println((B₁[3]*n₁ + B₂[3]*n₂)*N[1]*𝑤)
            # println((B₁[3]*n₁ + B₂[3]*n₂)*N[2]*𝑤)
            # println((B₁[3]*n₁ + B₂[3]*n₂)*(N[1]*g₁ + N[2]*g₂)*𝑤)
        end
        for (i,xᵢ) in enumerate(a₁.𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(a₂.𝓒)
                J = xⱼ.𝐼
                k[I,J] += kᶜ*(B₁[i]*n₁+B₂[i]*n₂)*N[j]*𝑤
            end
        end
    end
end

function (op::Operator{:∫vᵢnᵢgds})(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    kᶜ = op.k
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        g = ξ.g
        if ξ.𝐺 == 9 || ξ.𝐺 == 10
            # println(B₁[1])
            # println(B₁[2])
            # println(B₁[3])
            # println(B₂[1])
            # println(B₂[2])
            # println(B₂[3])
            # println((B₁[3]*n₁ + B₂[3]*n₂)*g*𝑤)
        end
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            if I == 1
                # println(f[I])
                # println((B₁[i]*n₁+B₂[i]*n₂)*g*𝑤)
            end
            f[I] += kᶜ*(B₁[i]*n₁+B₂[i]*n₂)*g*𝑤
        end
    end
end

function (op::Operator{:∫uds})(aps::Vector{T}) where T<:AbstractElement
    u = zeros(length(aps))
    for (c,ap) in enumerate(aps)
        𝓖 = ap.𝓖
        for ξ in 𝓖
            𝑤 = ξ.𝑤
            u[c] += ξ.u*𝑤
        end
        u[c] /= ap.𝐿
    end
    return u
end
