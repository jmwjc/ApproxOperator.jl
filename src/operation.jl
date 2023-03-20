"""
Operator
"""
struct Operator{T}
    data::Dict{Symbol,Float64}
end
Operator{T}(d::Pair{Symbol,D}...) where T = Operator{T}(Dict(d))


getproperty(op::Operator,f::Symbol) = getfield(op,:data)[f]

(op::Operator)(aps::Vector{T},gps::Vector{S},k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement} = op.(aps,gps,k=k,f=f)
(op::Operator)(aps::Vector{T},k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement = op.(aps,k=k,f=f)
(op::Operator)(aps::Vector{T},k::AbstractMatrix{Float64}) where T<:AbstractElement = op.(aps,k=k)
(op::Operator)(aps::Vector{T},f::AbstractVector{Float64}) where T<:AbstractElement = op.(aps,f=f)
(op::Operator)(aps::Vector{T}) where T<:AbstractElement = op.(ap)

function prescribe!(ap::T,sf::Pair{Symbol,F}) where {T<:AbstractElement,F<:Function}
    𝓖 = ap.𝓖
    s,f = sf
    for ξ in 𝓖
        𝒙 = (ξ.x,ξ.y,ξ.z)
        if applicable(f,𝒙...)
            v = f(𝒙...)
        elseif applicable(f,𝒙...,ξ.n₁)
            v = f(𝒙...,ξ.n₁)
        elseif applicable(f,𝒙...,ξ.n₁,ξ.n₂)
            v = f(𝒙...,ξ.n₁,ξ.n₂)
        elseif applicable(f,𝒙...,ξ.n₁,ξ.n₂,ξ.n₃)
            v = f(𝒙...,ξ.n₁,ξ.n₂,ξ.n₃)
        end
        setproperty!(ξ,s,v)
    end
end

function prescribe!(aps::Vector{T},sf::Pair{Symbol,F}) where {T<:AbstractElement,F<:Function}
    s,f = sf
    n = length(getfield(aps[1].𝓖[1],:data)[:x][2])
    haskey(getfield(aps[1].𝓖[1],:data),s) ? nothing : push!(getfield(aps[1].𝓖[1],:data),s=>(2,zeros(n)))
    for ap in aps
        prescribe!(ap,sf)
    end
end

"""
Phase field modeling fracture
"""
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

"""
Error Estimates
"""
function (op::Operator{:L₂})(ap::T) where T<:AbstractElement
    Δu²= 0
    ū² = 0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        ūᵢ = ξ.u
        uᵢ = 0
        for (i,xᵢ) in enumerate(ap.𝓒)
            uᵢ += N[i]*xᵢ.d
        end
        Δu² += (uᵢ - ūᵢ)^2*𝑤
        ū²  += ūᵢ^2*𝑤
    end
    return Δu², ū²
end

function (op::Operator{:L₂})(aps::Vector{T}) where T<:AbstractElement
    L₂Norm_Δu²= 0
    L₂Norm_ū² = 0
    for ap in aps
        Δu², ū² = op(ap)
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end

function (op::Operator{:H₁})(ap::T) where T<:AbstractElement
    Δ∇u²= 0
    ∇ū² = 0
    Δu²= 0
    ū² = 0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z]
        ūᵢ = ξ.u
        ∂ūᵢ∂x = ξ.∂u∂x
        ∂ūᵢ∂y = ξ.∂u∂y
        ∂ūᵢ∂z = ξ.∂u∂z
        uᵢ = 0.
        ∂uᵢ∂x = 0.
        ∂uᵢ∂y = 0.
        ∂uᵢ∂z = 0.
        for (i,xᵢ) in enumerate(ap.𝓒)
            uᵢ += N[i]*xᵢ.d
            ∂uᵢ∂x += B₁[i]*xᵢ.d
            ∂uᵢ∂y += B₂[i]*xᵢ.d
            ∂uᵢ∂z += B₃[i]*xᵢ.d
        end
        Δ∇u² += ((∂uᵢ∂x - ∂ūᵢ∂x)^2 + (∂uᵢ∂y - ∂ūᵢ∂y)^2 + (∂uᵢ∂z - ∂ūᵢ∂z)^2)*𝑤
        ∇ū² += (∂ūᵢ∂x^2 + ∂ūᵢ∂y^2 + ∂ūᵢ∂z^2)*𝑤
        Δu² += (uᵢ - ūᵢ)^2*𝑤
        ū² += ūᵢ^2*𝑤
    end
    return Δ∇u², ∇ū², Δu², ū²
end

function (op::Operator{:H₁})(aps::Vector{T}) where T<:AbstractElement
    H₁Norm_Δu²= 0
    H₁Norm_ū² = 0
    L₂Norm_Δu²= 0
    L₂Norm_ū² = 0
    for ap in aps
        Δ∇u², ∇ū², Δu², ū² = op(ap)
        H₁Norm_Δu² += Δu² + Δ∇u²
        H₁Norm_ū²  += ū² + ∇ū²
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (H₁Norm_Δu²/H₁Norm_ū²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end

function (op::Operator{:Hₑ_PlaneStress})(ap::T) where T<:AbstractElement
    ΔW²= 0
    W̄² = 0
    Δu²= 0
    ū² = 0
    E = op.E
    ν = op.ν
    Cᵢᵢᵢᵢ = E/(1-ν^2)
    Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
    Cᵢⱼᵢⱼ = E/2/(1+ν)
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        ū₁ = ξ.u
        ū₂ = ξ.v
        ∂ū₁∂x = ξ.∂u∂x
        ∂ū₁∂y = ξ.∂u∂y
        ∂ū₂∂x = ξ.∂v∂x
        ∂ū₂∂y = ξ.∂v∂y
        ε̄₁₁ = ∂ū₁∂x
        ε̄₂₂ = ∂ū₂∂y
        ε̄₁₂ = ∂ū₁∂y + ∂ū₂∂x
        σ̄₁₁ = Cᵢᵢᵢᵢ*ε̄₁₁ + Cᵢᵢⱼⱼ*ε̄₂₂
        σ̄₂₂ = Cᵢᵢᵢᵢ*ε̄₂₂ + Cᵢᵢⱼⱼ*ε̄₁₁
        σ̄₁₂ = Cᵢⱼᵢⱼ*ε̄₁₂
        u₁ = 0.
        u₂ = 0.
        ε₁₁ = 0.
        ε₂₂ = 0.
        ε₁₂ = 0.
        for (i,xᵢ) in enumerate(ap.𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += B₂[i]*xᵢ.d₁ + B₁[i]*xᵢ.d₂
        end
        σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂
        σ₂₂ = Cᵢᵢᵢᵢ*ε₂₂ + Cᵢᵢⱼⱼ*ε₁₁
        σ₁₂ = Cᵢⱼᵢⱼ*ε₁₂
        ΔW² += 0.5*((σ₁₁-σ̄₁₁)*(ε₁₁-ε̄₁₁) + (σ₂₂-σ̄₂₂)*(ε₂₂-ε̄₂₂) + (σ₁₂-σ̄₁₂)*(ε₁₂-ε̄₁₂))*𝑤
        W̄² += 0.5*(σ₁₁*ε₁₁ + σ₂₂*ε₂₂ + σ₁₂*ε₁₂)*𝑤
        Δu² += ((u₁ - ū₁)^2 + (u₂ - ū₂)^2)*𝑤
        ū² += (ū₁^2 + ū₂^2)*𝑤
    end
    return ΔW², W̄², Δu², ū²
end

function (op::Operator{:Hₑ_PlaneStress})(aps::Vector{T}) where T<:AbstractElement
    HₑNorm_ΔW²= 0.0
    HₑNorm_W̄² = 0.0
    L₂Norm_Δu²= 0.0
    L₂Norm_ū² = 0.0
    for ap in aps
        ΔW², W̄², Δu², ū² = op(ap)
        HₑNorm_ΔW² += ΔW²
        HₑNorm_W̄²  += W̄²
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end

function set∇𝑢!(ap::T) where T<:AbstractElement
    for ξ in ap.𝓖
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x_]
        B₂ = ξ[:∂𝝭∂y_]
        B₃ = ξ[:∂𝝭∂z_]
        𝒙 = (ξ.x,ξ.y,ξ.z)
        u = 0.
        ∂u∂x = 0.
        ∂u∂y = 0.
        ∂u∂z = 0.
        for (i,xᵢ) in enumerate(ap.𝓒)
            u += N[i]*x.d
            ∂u∂x += B₁[i]*x.d
            ∂u∂y += B₂[i]*x.d
            ∂u∂z += B₃[i]*x.d
        end
        ξ.x = 𝒙[1]
        ξ.y = 𝒙[2]
        ξ.z = 𝒙[3]
        ξ.u = u
        ξ.∂u∂x = ∂u∂x
        ξ.∂u∂y = ∂u∂y
        ξ.∂u∂z = ∂u∂z
    end
end

function (op::Operator{:H₃})(aps::Vector{T}) where T<:AbstractElement
    H₃Norm_Δu²= 0
    H₃Norm_ū² = 0
    H₂Norm_Δu²= 0
    H₂Norm_ū² = 0
    H₁Norm_Δu²= 0
    H₁Norm_ū² = 0
    L₂Norm_Δu²= 0
    L₂Norm_ū² = 0
    for ap in aps
        Δ∇³u², ∇³ū²,Δ∇²u², ∇²ū²,Δ∇u², ∇ū², Δu², ū² = op(ap)
        H₃Norm_Δu² += Δu² + Δ∇u² + Δ∇²u² + Δ∇³u²
        H₃Norm_ū²  += ū² + ∇ū² + ∇²ū² + ∇³ū²
        H₂Norm_Δu² += Δu² + Δ∇u² + Δ∇²u²
        H₂Norm_ū²  += ū² + ∇ū² + ∇²ū²
        H₁Norm_Δu² += Δu² + Δ∇u²
        H₁Norm_ū²  += ū² + ∇ū²
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (H₃Norm_Δu²/H₃Norm_ū²)^0.5, (H₂Norm_Δu²/H₂Norm_ū²)^0.5, (H₁Norm_Δu²/H₁Norm_ū²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end
function (op::Operator{:H₃})(ap::T) where T<:AbstractElement
    Δ∇³u²= 0
    ∇³ū² = 0
    Δ∇²u²= 0
    ∇²ū² = 0
    Δ∇u²= 0
    ∇ū² = 0
    Δu²= 0
    ū² = 0
    for ξ in ap.𝓖
        𝑤 = get𝑤(ap,ξ)
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₁₁ = ξ[:∂²𝝭∂x²]
        B₁₂ = ξ[:∂²𝝭∂x∂y]
        B₂₂ = ξ[:∂²𝝭∂y²]
        B₁₁₁ = ξ[:∂³𝝭∂x³]
        B₁₁₂ = ξ[:∂³𝝭∂x²∂y]
        B₁₂₂ = ξ[:∂³𝝭∂x∂y²]
        B₂₂₂ = ξ[:∂³𝝭∂y³]
        ūᵢ = ξ.u
        ∂ūᵢ∂x = ξ.∂u∂x
        ∂ūᵢ∂y = ξ.∂u∂y
        ∂²ūᵢ∂x² = ξ.∂²u∂x²
        ∂²ūᵢ∂x∂y = ξ.∂²u∂x∂y
        ∂²ūᵢ∂y² = ξ.∂²u∂y²
        ∂³ūᵢ∂x³ = ξ.∂³u∂x³
        ∂³ūᵢ∂x²∂y = ξ.∂³u∂x²∂y
        ∂³ūᵢ∂x∂y² = ξ.∂³u∂x∂y²
        ∂³ūᵢ∂y³ = ξ.∂³u∂y³
        uᵢ = 0.
        ∂uᵢ∂x = 0.
        ∂uᵢ∂y = 0.
        ∂²uᵢ∂x² = 0.
        ∂²uᵢ∂x∂y = 0.
        ∂²uᵢ∂y² = 0.
        ∂³uᵢ∂x³ = 0.
        ∂³uᵢ∂x²∂y = 0.
        ∂³uᵢ∂x∂y² = 0.
        ∂³uᵢ∂y³ = 0.
        for i in 1:length(ap.𝓒)
            xᵢ = ap.𝓒[i]
            I = xᵢ.id
            uᵢ += N[i]*xᵢ.d
            ∂uᵢ∂x += B₁[i]*xᵢ.d
            ∂uᵢ∂y += B₂[i]*xᵢ.d
            ∂²uᵢ∂x² += B₁₁[i]*xᵢ.d
            ∂²uᵢ∂x∂y += B₁₂[i]*xᵢ.d
            ∂²uᵢ∂y² += B₂₂[i]*xᵢ.d
            ∂³uᵢ∂x³ += B₁₁₁[i]*xᵢ.d
            ∂³uᵢ∂x²∂y += B₁₁₂[i]*xᵢ.d
            ∂³uᵢ∂x∂y² += B₁₂₂[i]*xᵢ.d
            ∂³uᵢ∂y³ += B₂₂₂[i]*xᵢ.d
        end
        Δ∇³u² += ((∂³uᵢ∂x³ - ∂³ūᵢ∂x³)^2 + (∂³uᵢ∂x²∂y - ∂³ūᵢ∂x²∂y)^2 + (∂³uᵢ∂x∂y² - ∂³ūᵢ∂x∂y²)^2 + (∂³uᵢ∂y³ - ∂³ūᵢ∂y³)^2)*𝑤
        ∇³ū² += (∂³ūᵢ∂x³^2 + ∂³ūᵢ∂x²∂y^2  + ∂³ūᵢ∂x∂y²^2+ ∂³ūᵢ∂y³^2)*𝑤
        Δ∇²u² += ((∂²uᵢ∂x² - ∂²ūᵢ∂x²)^2 + (∂²uᵢ∂x∂y - ∂²ūᵢ∂x∂y)^2 + (∂²uᵢ∂y² - ∂²ūᵢ∂y²)^2)*𝑤
        ∇²ū² += (∂²ūᵢ∂x²^2 + ∂²ūᵢ∂x∂y^2 + ∂²ūᵢ∂y²^2)*𝑤
        Δ∇u² += ((∂uᵢ∂x - ∂ūᵢ∂x)^2 + (∂uᵢ∂y - ∂ūᵢ∂y)^2)*𝑤
        ∇ū² += (∂ūᵢ∂x^2 + ∂ūᵢ∂y^2)*𝑤
        Δu² += (uᵢ - ūᵢ)^2*𝑤
        ū² += ūᵢ^2*𝑤
    end
    return Δ∇³u², ∇³ū², Δ∇²u², ∇²ū², Δ∇u², ∇ū², Δu², ū²
end

function set∇𝑢!(aps::Vector{T}) where T<:AbstractElement
    for ap in aps
        set∇𝑢!(ap)
    end
end

function (op::Operator{:∫udΓ})(aps::Vector{T}) where T<:AbstractElement
    d = zeros(length(aps))
    for (i,ap) in enumerate(aps)
        d[i] = op(ap)
    end
    return d
end

function (op::Operator{:∫udΓ})(ap::T) where T<:AbstractElement{:Seg2}
    𝓖 = ap.𝓖
    d = sum(ξ.u*ξ.w for ξ in 𝓖)/2
    return d
end

function (op::Operator{:∫∇udΓ})(aps::Vector{T}) where T<:AbstractElement
    ∂u∂x = zeros(length(aps))
    ∂u∂y = zeros(length(aps))
    for (i,ap) in enumerate(aps)
        ∂u∂x_,∂u∂y_ = op(ap)
        ∂u∂x[i] = ∂u∂x_
        ∂u∂y[i] = ∂u∂y_
    end
    return ∂u∂x,∂u∂y
end

function (op::Operator{:∫∇udΓ})(ap::T) where T<:AbstractElement
    𝓖 = ap.𝓖
    ∂u∂x = sum(ξ.∂u∂x*ξ.w for ξ in 𝓖)/2
    ∂u∂y = sum(ξ.∂u∂y*ξ.w for ξ in 𝓖)/2
    return ∂u∂x,∂u∂y
end

# function (op::Operator{:∫udΓ})(aps::Vector{T},k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
#     for ap in aps
#         op(ap,k,f)
#     end
# end

function (op::Operator{:∫udΓ})(ap::AbstractElement,k::AbstractMatrix{Float64},f::AbstractVector{Float64})
    ξ = ap.𝓖[1]
    j = ξ.𝐶
    g = op(ap)
    for i in 1:length(f)
        f[i] -= k[i,j]*g
    end
    k[j,:] .= 0.
    k[:,j] .= 0.
    k[j,j] = 1.
    f[j] = g
end

# function (op::Operator{:∫udΓ})(ap::DBelement{:Seg2},f::AbstractVector{Float64})
#     x = ap.𝓒[3]
#     j = x.𝐼
#     g = op(ap)
#     for i in 1:length(f)
#         f[i] -= k[i,j]*g
#     end
#     k[j,:] .= 0.
#     k[:,j] .= 0.
#     k[j,j] = 1.
#     f[j] = g
# end

function (op::Operator{:Δ∫vtdΓ})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                d = xⱼ.d
                f[I] += op.k*(B₁[i]*B₁[j] + B₂[i]*B₂[j] + B₃[i]*B₃[j])*d*𝑤
            end
        end
    end
end

function (op::Operator{:∫vudΓ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    α = op.α
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += α*N[i]*N[j]*sign(n₁+n₂)*𝑤
                # k[I,J] += α*N[i]*N[j]*𝑤
            end
            f[I] += α*N[i]*g*sign(n₁+n₂)*𝑤
            # f[I] += α*N[i]*g*𝑤
        end
    end
end

function (op::Operator{:∫uλdΓ})(ap::T,g::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        N̄ = ξ[:𝝭̄]
        for (k,xₖ) in enumerate(𝓒)
            K = xₖ.𝑖
            for (i,xᵢ) in enumerate(𝓒)
                I = xᵢ.𝐼
                g[I,K] -= N[i]*N̄[k]*𝑤
            end
        end
    end
end

function (op::Operator{:∫uλ̄dΓ})(aps::Vector{T},g::AbstractMatrix{Float64}) where T<:AbstractElement
    for (K,ap) in enumerate(aps)
        𝓒 = ap.𝓒
        𝓖 = ap.𝓖
        for ξ in 𝓖
            𝑤 = ξ.𝑤
            N = ξ[:𝝭]
            n₁ = ξ.n₁
            n₂ = ξ.n₂
            for (i,xᵢ) in enumerate(𝓒)
                I = xᵢ.𝐼
                g[I,K] -= sign(n₁+n₂)*N[i]*𝑤
            end
        end
    end
end

function (op::Operator{:∫∇𝑛uvdΓ})(ap::T,g::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭̄]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (k,xₖ) in enumerate(𝓒)
                K = xₖ.𝑖
                g[I,K] -= (B₁[i]*n₁+B₂[i]*n₂)*N[k]*𝑤
            end
        end
    end
end

function (op::Operator{:∫∇𝑛ugdΓ})(ap::T,g::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        ḡ = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (k,xₖ) in enumerate(𝓒)
                K = xₖ.𝑖
                g[I,K] += (B₁[i]*n₁+B₂[i]*n₂)*N[k]*𝑤
            end
            f[I] += (B₁[i]*n₁+B₂[i]*n₂)*ḡ*𝑤
        end
    end
end

function (op::Operator{:∫∇v̄∇udΩ})(ap::T,g::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B̄₁ = ξ[:∂𝝭̄∂x]
        B̄₂ = ξ[:∂𝝭̄∂y]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (k,xₖ) in enumerate(𝓒)
                K = xₖ.𝑖
                g[I,K] += (B₁[i]*B̄₁[k]+B₂[i]*B̄₂[k])*𝑤
            end
        end
    end
end

function (op::Operator{:∫vu𝑛dΓ})(ap::T,g::AbstractMatrix{Float64},q::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        N̄ = ξ[:𝝭̄]
        sn = sign(ξ.n₁+ξ.n₂)
        u = ξ.u
        for (k,xₖ) in enumerate(𝓒)
            K = xₖ.𝑖
            for (i,xᵢ) in enumerate(𝓒)
                I = xᵢ.𝐼
                g[I,K] -= sn*N[i]*N̄[k]*𝑤
            end
            q[K] -= sn*N̄[k]*u*𝑤
        end
    end
end

# function (op::Operator{:∫sᵢnᵢudΓ})(ap::T,g::AbstractMatrix{Float64}) where T<:AbstractElement
#     𝓒 = ap.𝓒
#     𝓖 = ap.𝓖
#     for ξ in 𝓖
#         𝑤 = ξ.𝑤
#         N = ξ[:𝝭]
#         N̄ = ξ[:𝝭̄]
#         n₁ = ξ.n₁
#         n₂ = ξ.n₂
#         for (k,xₖ) in enumerate(𝓒)
#             K = xₖ.𝑖
#             for (i,xᵢ) in enumerate(𝓒)
#                 I = xᵢ.𝐼
#                 g[I,2*K-1] -= N̄[k]*N[i]*n₁*𝑤
#                 g[I,2*K]   -= N̄[k]*N[i]*n₂*𝑤
#             end
#         end
#     end
# end

function (op::Operator{:∫sᵢnᵢudΓ})(ap::T,g::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        C = ξ.𝐶
        for (k,xₖ) in enumerate(𝓒)
            K = xₖ.𝐼
            g[2*C-1,K] -= N[k]*n₁*𝑤
            g[2*C,K]   -= N[k]*n₂*𝑤
        end
    end
end

# function (op::Operator{:∫sᵢnᵢgdΓ})(ap::T,q::AbstractVector{Float64}) where T<:AbstractElement
#     𝓒 = ap.𝓒
#     𝓖 = ap.𝓖
#     for ξ in 𝓖
#         𝑤 = ξ.𝑤
#         N̄ = ξ[:𝝭̄]
#         n₁ = ξ.n₁
#         n₂ = ξ.n₂
#         ḡ = ξ.g
#         for (k,xₖ) in enumerate(𝓒)
#             K = xₖ.𝑖
#             q[2*K-1] -= N̄[k]*ḡ*n₁*𝑤
#             q[2*K]   -= N̄[k]*ḡ*n₂*𝑤
#         end
#     end
# end

function (op::Operator{:∫sᵢnᵢgdΓ})(ap::T,g::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        C = ξ.𝐶
        ḡ = ξ.g
        for (k,xₖ) in enumerate(𝓒)
            K = xₖ.𝐼
            g[2*C-1,K] += N[k]*n₁*𝑤
            g[2*C,K]   += N[k]*n₂*𝑤
        end
        f[2*C-1] += n₁*ḡ*𝑤 
        f[2*C]   += n₂*ḡ*𝑤 
    end
end

function (op::Operator{:∫usᵢnᵢdΓ})(ap::T,g::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        C = ξ.𝐶
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            g[I,2*C-1] -= N[i]*n₁*𝑤
            g[I,2*C]   -= N[i]*n₂*𝑤
        end
    end
end

function (op::Operator{:∫gsᵢnᵢdΓ})(ap::T,q::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        C = ξ.𝐶
        g = ξ.g
        q[2*C-1] -= n₁*g*𝑤 
        q[2*C]   -= n₂*g*𝑤 
    end
end

function (op::Operator{:∫tdΓ})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        b𝑤 = ξ.b*ξ.𝑤
        for x in 𝓒
            I = x.𝐼
            f[I] += b𝑤 
        end
    end
end

function (op::Operator{:∫bdΩ})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        b𝑤 = ξ.b*ξ.𝑤
        for x in 𝓒
            I = x.𝐼
            f[I] += b𝑤 
        end
    end
end

"""
pure meshfree Potential
"""
function (op::Operator{:∫∇̃v∇uvbdΩ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    kᶜ = op.k
    for (K,ξ) in enumerate(𝓖)
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B̃₁ = ξ[:∂𝝭∂x_]
        # B₂ = ξ[:∂𝝭∂y]
        # B₃ = ξ[:∂𝝭∂z]
        𝑤 = ξ.𝑤
        b = ξ.b
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += kᶜ*B̃₁[i]*B₁[j]*𝑤
                # k[I,J] += kᶜ*(B₁[i]*B₁[j] + B₂[i]*B₂[j] + B₃[i]*B₃[j])*𝑤
            end
        end
        I = 𝓒[K].𝐼
        f[I] += b*𝑤
    end
end

function (op::Operator{:∫ṽtdΓ})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    x = 𝓒[1]
    I = x.𝐼
    ξ = ap.𝓖[1]
    t = ξ.t
    𝑤 = ξ.𝑤
    f[I] += t*𝑤
    # for ξ in 𝓖
    #     𝑤 = ξ.𝑤
    #     t = ξ.t
    #     f[I] += t*𝑤
    # end
end

function (op::Operator{:∫ṽgdΓ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    α = op.α
    for (i,ξ) in enumerate(𝓖)
        I = 𝓒[i].𝐼
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        g = ξ.g
        for (j,xⱼ) in enumerate(𝓒)
            J = xⱼ.𝐼
            k[I,J] += α*N[j]*𝑤
        end
        f[I] += α*g*𝑤
    end
end

