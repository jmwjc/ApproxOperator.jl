
struct GRKGradientSmoothing{𝑝,𝑠,𝜙,T}<:AbstractReproducingKernel{𝑠,𝜙,T}
    𝓒::Tuple{Int,Int,Vector{Node{(:𝐼,),1}}}
    𝓒ᵖ::Tuple{Int,Int,Vector{Node{(:𝐼,),1}}}
    𝓖::Tuple{Int,Int,Vector{Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}}}
    𝓖ˢ::Tuple{Int,Int,Vector{Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}}}
    𝓖ᵖ::Tuple{Int,Int,Vector{Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}}}
    𝗚::Matrix{Float64}
    𝗴₁::Matrix{Float64}
    𝗴₂::Matrix{Float64}
end

function Base.getproperty(a::GRKGradientSmoothing,s::Symbol)
    if s∈(:𝓒,:𝓒ᵖ,:𝓖,:𝓖ˢ,:𝓖ᵖ)
        𝓐 =  getfield(a,s)
        return (𝓐[3][𝓐[1]+i] for i in 1:𝓐[2])
    elseif s∈(:𝗚,:𝗴₁,:𝗴₂)
        return getfield(a,s)
    else
        𝓖 = getfield(a,:𝓖)
        ξ = 𝓖[3][𝓖[1]+1]
        return getproperty(ξ,s)
    end
end

function cal𝗠!(aps::Vector{T}) where T<:GRKGradientSmoothing
    cal𝗠!.(aps)
    𝗚 = aps[1].𝗚
    𝗴₁ = aps[1].𝗴₁
    𝗴₂ = aps[1].𝗴₂
    𝗴₁ .= 𝗚\𝗴₁
    𝗴₂ .= 𝗚\𝗴₂
end

function cal𝗠!(ap::GRKGradientSmoothing)
    𝗚 = ap.𝗚
    𝗴₁ = ap.𝗴₁
    𝗴₂ = ap.𝗴₂
    𝓒 = ap.𝓒
    𝓒ᵖ = ap.𝓒ᵖ
    𝓖 = ap.𝓖
    𝓖ˢ = ap.𝓖ˢ
    𝓖ᵖ = ap.𝓖ᵖ
    𝓖ᵗ = zip(𝓖ˢ,𝓖ᵖ)
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒ᵖ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵖ)
                J = xⱼ.𝐼
                𝗚[I,J] += N[i]*N[j]*𝑤
            end
        end
    end
    for (ξˢ,ξᵖ) in 𝓖ᵗ
        D₁ = ξˢ.D₁
        D₂ = ξˢ.D₂
        wᵇ = ξˢ.wᵇ
        𝑤 = ξˢ.𝑤
        Nᵖ = ξᵖ[:𝝭]
        B₁ᵖ = ξᵖ[:∂𝝭∂x]
        B₂ᵖ = ξᵖ[:∂𝝭∂y]
        Nˢ = ξˢ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒ᵖ)
            I = xᵢ.𝐼
            for (k,xₖ) in enumerate(𝓒)
                K = xₖ.𝐼
                𝗴₁[I,K] += Nᵖ[i]*Nˢ[k]*D₁*wᵇ - B₁ᵖ[i]*Nˢ[k]*𝑤
                𝗴₂[I,K] += Nᵖ[i]*Nˢ[k]*D₂*wᵇ - B₂ᵖ[i]*Nˢ[k]*𝑤
            end
        end
    end
end

function set∇𝝭!(ap::GRKGradientSmoothing{𝒑,𝑠,𝜙,:Tri3}) where {𝒑,𝑠,𝜙}
    𝓒 = ap.𝓒
    𝓒ᵖ = ap.𝓒ᵖ
    𝓖 = ap.𝓖
    𝗴₁ = ap.𝗴₁
    𝗴₂ = ap.𝗴₂
    for ξ in 𝓖
        N = ξ[:𝝭]
        ∂𝝭∂x = ξ[:∂𝝭∂x]
        ∂𝝭∂y = ξ[:∂𝝭∂y]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵖ)
                J = xⱼ.𝐼
                ∂𝝭∂x[i] += N[j]*𝗴₁[J,I]
                ∂𝝭∂y[i] += N[j]*𝗴₂[J,I]
            end
        end
    end
end