
module Test

using ..ApproxOperator: AbstractElement

function cc𝝭(ap::T) where T<:AbstractElement
    Δu²= 0
    ū² = 0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        ūᵢ = ξ.u
        uᵢ = 0
        for (i,xᵢ) in enumerate(ap.𝓒)
            uᵢ += N[i]*xᵢ.d
            # print("N$i = ")
            # println(N[i])
            # print("d$i = ")
            # println(xᵢ.d)
        end
        # print("ūᵢ = ")
        # println(ūᵢ)
        # print("uᵢ = ")
        # println(uᵢ)
        Δu² += (uᵢ - ūᵢ)^2*𝑤
        ū²  += ūᵢ^2*𝑤
    end
    return Δu², ū²
end

function cc𝝭(aps::Vector{T}) where T<:AbstractElement
    L₂Norm_Δu²= 0.
    L₂Norm_ū² = 0.
    for ap in aps
        Δu², ū² = cc𝝭(ap)
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end

function cc∇𝝭(ap::T) where T<:AbstractElement
    Δ∇u²= 0
    ∇ū² = 0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        ∂ūᵢ∂x = ξ.∂u∂x
        ∂ūᵢ∂y = ξ.∂u∂y
        ∂uᵢ∂x = 0
        ∂uᵢ∂y = 0
        for (i,xᵢ) in enumerate(ap.𝓒)
            ∂uᵢ∂x += B₁[i]*xᵢ.d
            ∂uᵢ∂y += B₂[i]*xᵢ.d
        end
        Δ∇u² += ((∂uᵢ∂x - ∂ūᵢ∂x)^2 + (∂uᵢ∂y - ∂ūᵢ∂y)^2)*𝑤
        ∇ū²  += (∂ūᵢ∂x^2 + ∂ūᵢ∂y^2)*𝑤
    end
    return Δ∇u², ∇ū²
end

function cc∇𝝭(aps::Vector{T}) where T<:AbstractElement
    H₁Norm_Δu²= 0.
    H₁Norm_ū² = 0.
    for ap in aps
        Δu², ū² = cc𝝭(ap)
        H₁Norm_Δu² += Δu²
        H₁Norm_ū²  += ū²
    end
    return (H₁Norm_Δu²/H₁Norm_ū²)^0.5
end

end