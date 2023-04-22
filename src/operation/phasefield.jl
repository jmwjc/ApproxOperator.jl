
function (op::Operator{:∫v²uₓuₓdx})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    EA = op.EA
    η = op.η
    for ξ in 𝓖
        N = ξ[:𝝭]
        B = ξ[:∂𝝭∂x]
        v = sum(N[i]*xᵢ.v for (i,xᵢ) in enumerate(𝓒))
        # println(v)
        b = ξ.b
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (v^2+η)*EA*B[i]*B[j]*𝑤
            end
            f[I] += N[i]*b*𝑤
        end
    end
end

function (op::Operator{:∫vₓvₓvvdx})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    kc = op.k
    l = op.l
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
        ℋₜ = max(ℋ,ε^2)
        # println(ℋₜ)
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += kc*(2*l*B[i]*B[j] + N[i]*N[j]/2/l)*𝑤
            end
            f[I] += N[i]*(kc/2/l - ℋₜ)*𝑤
            # println(f[I])
            # f[I] = N[i]*(kc/2/l - ℋₜ - kc*(2*l*∂v∂x + v/2/l))*𝑤
        end
    end
end

function (op::Operator{:UPDATE_PFM_1D})(ap::T) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B = ξ[:∂𝝭∂x]
        ℋ = ξ.ℋ
        ε = sum(B[i]*xᵢ.u for (i,xᵢ) in enumerate(𝓒))
        ξ.ℋ = max(ℋ,ε^2)
    end
end