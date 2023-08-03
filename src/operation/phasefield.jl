
function (op::Operator{:∫v²uₓuₓdx})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    η = op.η
    for ξ in 𝓖
        N = ξ[:𝝭]
        B = ξ[:∂𝝭∂x]
        v = sum(N[i]*xᵢ.v for (i,xᵢ) in enumerate(𝓒))
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (v^2+η)*B[i]*B[j]*𝑤
            end
        end
    end
end

function (op::Operator{:∫vₓvₓvvdx_hard_device})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    kc = op.k
    l = op.l
    for ξ in 𝓖
        N = ξ[:𝝭]
        B = ξ[:∂𝝭∂x]
        ℋ = ξ.ℋ
        ε = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            ε += B[i]*xᵢ.u
        end
        ℋₜ = max(ℋ,(ε-1)^2)
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (kc*(2*l*B[i]*B[j] + N[i]*N[j]/2/l) + ℋₜ*N[i]*N[j])*𝑤
            end
            f[I] += N[i]*kc/2/l*𝑤
        end
    end
end

function (op::Operator{:∫vₓvₓvvdx})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    kc = op.k
    l = op.l
    η = op.η
    for ξ in 𝓖
        N = ξ[:𝝭]
        B = ξ[:∂𝝭∂x]
        ℋ = ξ.ℋ
        ε = 0.0
        v = 0.0
        ∂v∂x = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            ε += B[i]*xᵢ.u
            v += N[i]*xᵢ.v
            ∂v∂x += B[i]*xᵢ.v
        end
        ℋₜ = max(ℋ,(v^2+η)*ε^2)
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (kc*(2*l*B[i]*B[j] + N[i]*N[j]/2/l) + ℋₜ*N[i]*N[j])*𝑤
            end
            f[I] += N[i]*(kc/2/l - η*ℋₜ)*𝑤
        end
    end
end

function (op::Operator{:UPDATE_PFM_1D})(ap::T) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    η = op.η
    for ξ in 𝓖
        N = ξ[:𝝭]
        B = ξ[:∂𝝭∂x]
        ℋ = ξ.ℋ
        for (i,xᵢ) in enumerate(𝓒)
            ε += B[i]*xᵢ.u
            v += N[i]*xᵢ.v
        end
        ξ.ℋ = max(ℋ,(v^2+η)*ε^2)
    end
end


function (op::Operator{:∫∫∇v∇vvvdxdy})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    kc = op.k
    l = op.l
    η = op.η
    for ξ in 𝓖
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        ℋ = ξ.ℋ
        ε₁₁ = 0.0
        ε₂₂ = 0.0
        ε₁₂ = 0.0
        σ₁₁ = ξ.σ₁₁
        σ₂₂ = ξ.σ₂₂
        σ₁₂ = ξ.σ₁₂
        for (i,xᵢ) in enumerate(𝓒)
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += B₁[i]*xᵢ.d₂ + B₂[i]*xᵢ.d₁
        end
        ℋₜ = max(ℋ,0.5*(ε₁₁*σ₁₁ + ε₂₂*σ₂₂ + ε₁₂*σ₁₂))
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (kc*(2*l*(B₁[i]*B₁[j] + B₂[i]*B₂[j]) + N[i]*N[j]/2/l) + ℋₜ*N[i]*N[j])*𝑤
            end
            f[I] += N[i]*(kc/2/l - η*ℋₜ)*𝑤
        end
    end
end

function (op::Operator{:UPDATE_PFM_2D})(ap::T) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    η = op.η
    for ξ in 𝓖
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂x]
        σ₁₁ = ξ.σ₁₁
        σ₂₂ = ξ.σ₂₂
        σ₁₂ = ξ.σ₁₂
        ε₁₁ = 0.0
        ε₂₂ = 0.0
        ε₁₂ = 0.0
        ℋ = ξ.ℋ
        for (i,xᵢ) in enumerate(𝓒)
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += B₁[i]*xᵢ.d₂ + B₂[i]*xᵢ.d₁
            v += N[i]*xᵢ.v
        end
        ξ.ℋ = max(ℋ,0.5*(ε₁₁*σ₁₁ + ε₂₂*σ₂₂ + ε₁₂*σ₁₂))
    end
end