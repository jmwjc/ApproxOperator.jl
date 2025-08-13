module Stokes

using ..ApproxOperator: AbstractElement

function ∫∫μ∇u∇vdxdy(a::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = a.𝓒; 𝓖 = a.𝓖
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        μ = ξ.μ
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2I-1,2J-1] += μ * (B₁[i]*B₁[j] + B₂[i]*B₂[j]) * 𝑤
                # k[2I,2J-1]   += μ * (B₁[i]*B₂[j] + B₂[i]*B₁[j]) * 𝑤 / 2
                # k[2I-1,2J]   += μ * (B₂[i]*B₁[j] + B₁[i]*B₂[j]) * 𝑤 / 2
                k[2I,2J]     += μ * (B₁[i]*B₁[j] + B₂[i]*B₂[j]) * 𝑤
            end
        end
    end
end

end