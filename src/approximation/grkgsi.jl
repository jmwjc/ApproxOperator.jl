
struct GRKGradientSmoothing{𝑝,𝑠,𝜙,T}<:AbstractReproducingKernel{𝑠,𝜙}
    𝓒::Tuple{Int,Int,Vector{Node{(:𝐼,),1}}}
    𝓒ᵘ::Tuple{Int,Int,Vector{Node{(:𝐼,),1}}}
    𝓒ᵖ::Tuple{Int,Int,Vector{Node{(:𝐼,),1}}}
    𝓖::Tuple{Int,Int,Vector{Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}}}
    𝓖ᵖ::Tuple{Int,Int,Vector{Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}}}
    𝓖ˢ::Tuple{Int,Int,Vector{Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}}}
    𝓖ˢᵖ::Tuple{Int,Int,Vector{Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}}}
    𝗚::Matrix{Float64}
    𝗴₁::Matrix{Float64}
    𝗴₂::Matrix{Float64}
end

function Base.getproperty(a::GRKGradientSmoothing,s::Symbol)
    if s∈(:𝓒,:𝓒ᵘ,:𝓒ᵖ,:𝓖,:𝓖ᵖ,:𝓖ˢ,:𝓖ˢᵖ)
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
    𝗚 = aps[1].𝗚
    𝗴₁ = aps[1].𝗴₁
    𝗴₂ = aps[1].𝗴₂
    fill!(𝗚,0.0)
    fill!(𝗴₁,0.0)
    fill!(𝗴₂,0.0)
    cal𝗠!.(aps)
    𝗴₁ .= 𝗚\𝗴₁
    𝗴₂ .= 𝗚\𝗴₂
end

function cal𝗠!(ap::GRKGradientSmoothing)
    𝗚 = ap.𝗚
    𝗴₁ = ap.𝗴₁
    𝗴₂ = ap.𝗴₂
    𝓒ᵘ= ap.𝓒ᵘ
    𝓒ᵖ = ap.𝓒ᵖ
    𝓖ˢ = ap.𝓖ˢ
    𝓖ᵖ = ap.𝓖ᵖ
    𝓖ˢᵖ = ap.𝓖ˢᵖ
    𝓖ᵗ = zip(𝓖ˢ,𝓖ˢᵖ)
    for ξ in 𝓖ᵖ
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
    for (ξ,ξᵖ) in 𝓖ᵗ
        D₁ = ξ.D₁
        D₂ = ξ.D₂
        wᵇ = ξ.wᵇ
        𝑤 = ξ.𝑤
        Nᵖ = ξᵖ[:𝝭]
        B₁ᵖ = ξᵖ[:∂𝝭∂x]
        B₂ᵖ = ξᵖ[:∂𝝭∂y]
        N = ξ[:𝝭]

        for (i,xᵢ) in enumerate(𝓒ᵖ)
            I = xᵢ.𝐼
            for (k,xₖ) in enumerate(𝓒ᵘ)
                K = xₖ.𝐼
                𝗴₁[I,K] += Nᵖ[i]*N[k]*D₁*wᵇ - B₁ᵖ[i]*N[k]*𝑤
                𝗴₂[I,K] += Nᵖ[i]*N[k]*D₂*wᵇ - B₂ᵖ[i]*N[k]*𝑤
            end
        end
    end
end

function set∇𝝭!(ap::GRKGradientSmoothing{𝒑,𝑠,𝜙,:Tri3}) where {𝒑,𝑠,𝜙}
    𝓒 = ap.𝓒
    𝓒ᵖ = ap.𝓒ᵖ
    𝓖 = ap.𝓖
    𝓖ᵖ = ap.𝓖ᵖ
    𝗴₁ = ap.𝗴₁
    𝗴₂ = ap.𝗴₂
    𝓖ᵗ = zip(𝓖,𝓖ᵖ)
    for (ξ,ξᵖ) in 𝓖ᵗ
        N = ξᵖ[:𝝭]
        ∂𝝭∂x = ξ[:∂𝝭∂x]
        ∂𝝭∂y = ξ[:∂𝝭∂y]
        for i in 1:length(𝓒)
            ∂𝝭∂x[i] = 0.0
            ∂𝝭∂y[i] = 0.0
        end

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

function set∇𝝭!(aps::Vector{T}) where T<:GRKGradientSmoothing
    cal𝗠!(aps)
    set∇𝝭!.(aps)
end

struct FRKGradientSmoothing{𝑝,𝑠,𝜙,T}<:AbstractReproducingKernel{𝑠,𝜙}
    𝓒::Tuple{Int,Int,Vector{Node{(:𝐼,),1}}}
    𝓒ᵐ::Tuple{Int,Int,Vector{Node{(:𝐼,),1}}}
    𝓒ᶠ::Tuple{Int,Int,Vector{Node{(:𝐼,),1}}}
    𝓖::Tuple{Int,Int,Vector{Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}}}
    𝓖ˢ::Tuple{Int,Int,Vector{Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}}}
    𝓖ᶠ::Tuple{Int,Int,Vector{Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}}}
    𝗚::Matrix{Float64}
    𝗴₁::Matrix{Float64}
    𝗴₂::Matrix{Float64}
end

function Base.getproperty(a::FRKGradientSmoothing,s::Symbol)
    if s∈(:𝓒,:𝓒ᵐ,:𝓒ᶠ,:𝓖,:𝓖ˢ,:𝓖ᶠ)
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

function cal𝗠!(aps::Vector{T}) where T<:FRKGradientSmoothing
    𝗚 = aps[1].𝗚
    𝗴₁ = aps[1].𝗴₁
    𝗴₂ = aps[1].𝗴₂
    fill!(𝗚,0.0)
    fill!(𝗴₁,0.0)
    fill!(𝗴₂,0.0)
    cal𝗠!.(aps)
    𝗴₁ .= 𝗚\𝗴₁
    𝗴₂ .= 𝗚\𝗴₂
end

function cal𝗠!(ap::FRKGradientSmoothing)
    𝗚 = ap.𝗚
    𝗴₁ = ap.𝗴₁
    𝗴₂ = ap.𝗴₂
    𝓒ᵐ = ap.𝓒ᵐ
    𝓒ᶠ = ap.𝓒ᶠ
    (v₁,v₂,v₃) = 𝓒ᶠ
    D₁₁ = v₃.y-v₂.y
    D₁₂ = v₂.x-v₃.x
    D₂₁ = v₁.y-v₃.y
    D₂₂ = v₃.x-v₁.x
    D₃₁ = v₂.y-v₁.y
    D₃₂ = v₁.x-v₂.x
    𝐴 = ap.𝐴
    𝓖 = ap.𝓖
    𝓖ˢ = ap.𝓖ˢ
    𝓖ᶠ = ap.𝓖ᶠ
    𝓖ᵗ = zip(𝓖ˢ,𝓖ᶠ)
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        ξ₁ = ξ.ξ
        ξ₂ = ξ.η
        ξ₃ = 1.0-ξ₁-ξ₂
        N = (ξ₁,ξ₂,ξ₃)
        for (i,xᵢ) in enumerate(𝓒ᶠ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᶠ)
                J = xⱼ.𝐼
                𝗚[I,J] += N[i]*N[j]*𝑤
            end
        end
    end
    for (ξ,ξᶠ) in 𝓖ᵗ
        D₁ = ξ.D₁
        D₂ = ξ.D₂
        wᵇ = ξ.wᵇ
        𝑤 = ξ.𝑤
        Nᶠ = ξᶠ[:𝝭]
        B₁ᶠ = ξᶠ[:∂𝝭∂x]
        B₂ᶠ = ξᶠ[:∂𝝭∂y]
        # Nᶠ = (ξ.ξ,ξ.η,1-ξ.ξ-ξ.η)
        # B₁ᶠ = (-D₁₁/2/𝐴,-D₂₁/2/𝐴,-D₃₁/2/𝐴)
        # B₂ᶠ = (-D₁₂/2/𝐴,-D₂₂/2/𝐴,-D₃₂/2/𝐴)
        N = ξ[:𝝭]

        for (i,xᵢ) in enumerate(𝓒ᶠ)
            I = xᵢ.𝐼
            for (k,xₖ) in enumerate(𝓒ᵐ)
                K = xₖ.𝐼
                𝗴₁[I,K] += Nᶠ[i]*N[k]*D₁*wᵇ - B₁ᶠ[i]*N[k]*𝑤
                𝗴₂[I,K] += Nᶠ[i]*N[k]*D₂*wᵇ - B₂ᶠ[i]*N[k]*𝑤
            end
        end
    end
end

function set∇𝝭!(ap::FRKGradientSmoothing{𝒑,𝑠,𝜙,:Tri3}) where {𝒑,𝑠,𝜙}
    𝓒 = ap.𝓒
    𝓒ᶠ = ap.𝓒ᶠ
    𝓖 = ap.𝓖
    𝗴₁ = ap.𝗴₁
    𝗴₂ = ap.𝗴₂
    for ξ in 𝓖
        ξ₁ = ξ.ξ
        ξ₂ = ξ.η
        ξ₃ = 1.0-ξ₁-ξ₂
        N = (ξ₁,ξ₂,ξ₃)
        ∂𝝭∂x = ξ[:∂𝝭∂x]
        ∂𝝭∂y = ξ[:∂𝝭∂y]

        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᶠ)
                J = xⱼ.𝐼
                ∂𝝭∂x[i] += N[j]*𝗴₁[J,I]
                ∂𝝭∂y[i] += N[j]*𝗴₂[J,I]
            end
        end
    end
end

function set∇𝝭!(aps::Vector{T}) where T<:FRKGradientSmoothing
    cal𝗠!(aps)
    set∇𝝭!.(aps)
end