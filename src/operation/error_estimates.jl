
function (op::Operator{:T_error})(ap::T) where T<:AbstractElement
    ΔT²= 0
    T² = 0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        Tᵢ = ξ.T
        Tʰᵢ = 0
        for (i,xᵢ) in enumerate(ap.𝓒)
            Tʰᵢ += N[i]*xᵢ.T
        end
        ΔT² += (Tʰᵢ -Tᵢ)^2*𝑤
        T²  += Tᵢ^2*𝑤
    end
    return ΔT², T²
end

function (op::Operator{:T_error})(aps::Vector{T}) where T<:AbstractElement
    L₂Norm_ΔT²= 0
    L₂Norm_T² = 0
    for ap in aps
        ΔT², T² = op(ap)
        L₂Norm_ΔT² += ΔT²
        L₂Norm_T²  += T²
    end
    return (L₂Norm_ΔT²/L₂Norm_T²)^0.5
end

function (op::Operator{:L₂_heat_flux})(ap::T) where T<:AbstractElement
    Δu²= 0
    ū² = 0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        ū₁ = ξ.p₁
        ū₂ = ξ.p₂
        u₁ = 0.
        u₂ = 0.
        for (i,xᵢ) in enumerate(ap.𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
        end
        # println(ū₁)
        # println(u₁)
        Δu² += ((u₁ - ū₁)^2 + (u₂ - ū₂)^2)*𝑤
        ū² += (ū₁^2 + ū₂^2)*𝑤
    end
    return Δu², ū²
end
function (op::Operator{:L₂_heat_flux_u})(ap::T) where T<:AbstractElement
    Δu²= 0
    ū² = 0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        ū₁ = ξ.p₁
        ū₂ = ξ.p₂
        u₁ = 0.
        u₂ = 0.
        for (i,xᵢ) in enumerate(ap.𝓒)
            u₁ -= B₁[i]*xᵢ.d
            u₂ -= B₂[i]*xᵢ.d
        end
        # println(ū₁)
        # println(u₁)
        Δu² += ((u₁ - ū₁)^2 + (u₂ - ū₂)^2)*𝑤
        ū² += (ū₁^2 + ū₂^2)*𝑤
    end
    return Δu², ū²
end

function (op::Operator{:L₂_heat_flux_u})(aps::Vector{T}) where T<:AbstractElement
    L₂Norm_Δu²= 0.0
    L₂Norm_ū² = 0.0
    for ap in aps
      Δu², ū² = op(ap)
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end
function (op::Operator{:L₂_heat_flux})(aps::Vector{T}) where T<:AbstractElement
    L₂Norm_Δu²= 0.0
    L₂Norm_ū² = 0.0
    for ap in aps
      Δu², ū² = op(ap)
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end
function (op::Operator{:Hₑ_ThinShell})(ap::T) where T<:AbstractElement
    ΔW²= 0
    W̄² = 0
    Δu²= 0
    ū² = 0
    E = op.E
    ν = op.ν
    h = op.h
    Dᵐ = E*h/(1-ν^2)
    Dᵇ = E*h^3/12/(1-ν^2)
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₁₁ = ξ[:∂²𝝭∂x²]
        B₁₂ = ξ[:∂²𝝭∂x∂y]
        B₂₂ = ξ[:∂²𝝭∂y²]
        a¹¹ = ξ.a¹¹
        a¹² = ξ.a¹²
        a²² = ξ.a²²
        𝒂₁₍₁₎ = ξ.𝒂₁₍₁₎
        𝒂₁₍₂₎ = ξ.𝒂₁₍₂₎
        𝒂₁₍₃₎ = ξ.𝒂₁₍₃₎
        𝒂₂₍₁₎ = ξ.𝒂₂₍₁₎
        𝒂₂₍₂₎ = ξ.𝒂₂₍₂₎
        𝒂₂₍₃₎ = ξ.𝒂₂₍₃₎
        𝒂₃₍₁₎ = ξ.𝒂₃₍₁₎
        𝒂₃₍₂₎ = ξ.𝒂₃₍₂₎
        𝒂₃₍₃₎ = ξ.𝒂₃₍₃₎
        Γ¹₁₁ = ξ.Γ¹₁₁
        Γ¹₂₂ = ξ.Γ¹₂₂
        Γ¹₁₂ = ξ.Γ¹₁₂
        Γ²₁₁ = ξ.Γ²₁₁
        Γ²₂₂ = ξ.Γ²₂₂
        Γ²₁₂ = ξ.Γ²₁₂
        ū₁ = ξ.u₁
        ū₂ = ξ.u₂
        ū₃ = ξ.u₃
        ε̄₁₁ = ξ.ε₁₁
        ε̄₂₂ = ξ.ε₂₂
        ε̄₁₂ = ξ.ε₁₂
        κ̄₁₁ = ξ.κ₁₁
        κ̄₂₂ = ξ.κ₂₂
        κ̄₁₂ = ξ.κ₁₂
        𝑪 = @SArray [                  a¹¹*a¹¹ ν*a¹¹*a²² + (1-ν)*a¹²*a¹²                            a¹¹*a¹²;
                     ν*a¹¹*a²² + (1-ν)*a¹²*a¹²                   a²²*a²²                            a²²*a¹²;
                                       a¹¹*a¹²                   a²²*a¹² 0.5*((1-ν)*a¹¹*a²² + (1+ν)*a¹²*a¹²)]
        𝜺̄ = @SArray [ε̄₁₁,ε̄₂₂,2*ε̄₁₂]
        𝜿̄ = @SArray [κ̄₁₁,κ̄₂₂,2*κ̄₁₂]
        u₁ = 0.
        u₂ = 0.
        u₃ = 0.
        ε₁₁ = 0.
        ε₂₂ = 0.
        ε₁₂ = 0.
        κ₁₁ = 0.
        κ₂₂ = 0.
        κ₁₂ = 0.
        for (i,xᵢ) in enumerate(ap.𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            u₃ += N[i]*xᵢ.d₃
            ε₁₁ += 𝒂₁₍₁₎*B₁[i]*xᵢ.d₁ + 𝒂₁₍₂₎*B₁[i]*xᵢ.d₂ + 𝒂₁₍₃₎*B₁[i]*xᵢ.d₃
            ε₂₂ += 𝒂₂₍₁₎*B₂[i]*xᵢ.d₁ + 𝒂₂₍₂₎*B₂[i]*xᵢ.d₂ + 𝒂₂₍₃₎*B₂[i]*xᵢ.d₃
            ε₁₂ += (𝒂₁₍₁₎*B₂[i]+𝒂₂₍₁₎*B₁[i])*xᵢ.d₁ + (𝒂₁₍₂₎*B₂[i]+𝒂₂₍₂₎*B₁[i])*xᵢ.d₂ + (𝒂₁₍₃₎*B₂[i]+𝒂₂₍₃₎*B₁[i])*xᵢ.d₃
            Bᵇ₁₁ = Γ¹₁₁*B₁[i]+Γ²₁₁*B₂[i]-B₁₁[i]
            Bᵇ₂₂ = Γ¹₂₂*B₁[i]+Γ²₂₂*B₂[i]-B₂₂[i]
            Bᵇ₁₂ = Γ¹₁₂*B₁[i]+Γ²₁₂*B₂[i]-B₁₂[i]
            κ₁₁ += Bᵇ₁₁*𝒂₃₍₁₎*xᵢ.d₁ + Bᵇ₁₁*𝒂₃₍₂₎*xᵢ.d₂ + Bᵇ₁₁*𝒂₃₍₃₎*xᵢ.d₃
            κ₂₂ += Bᵇ₂₂*𝒂₃₍₁₎*xᵢ.d₁ + Bᵇ₂₂*𝒂₃₍₂₎*xᵢ.d₂ + Bᵇ₂₂*𝒂₃₍₃₎*xᵢ.d₃
            κ₁₂ += 2*(Bᵇ₁₂*𝒂₃₍₁₎*xᵢ.d₁ + Bᵇ₁₂*𝒂₃₍₂₎*xᵢ.d₂ + Bᵇ₁₂*𝒂₃₍₃₎*xᵢ.d₃)
        end
        # println(u₃/ū₃)
        # println(u₂)
        # println(ū₂)
        # println(u₃)
        # println(ū₃)
        𝜺 = @SArray [ε₁₁,ε₂₂,ε₁₂]
        𝜿 = @SArray [κ₁₁,κ₂₂,κ₁₂]
        ΔW² += 0.5*(Dᵐ*((𝜺-𝜺̄)'*𝑪*(𝜺-𝜺̄)) + Dᵇ*((𝜿-𝜿̄)'*𝑪*(𝜿-𝜿̄)))*𝑤
        W̄² += 0.5*(Dᵐ*(𝜺̄'*𝑪*𝜺̄) + Dᵇ*(𝜿̄'*𝑪*𝜿̄))*𝑤
        Δu² += ((u₁ - ū₁)^2 + (u₂ - ū₂)^2 + (u₃ - ū₃)^2)*𝑤
        ū² += (ū₁^2 + ū₂^2 + ū₃^2)*𝑤
    end
    return ΔW², W̄², Δu², ū²
end

function (op::Operator{:Hₑ_ThinShell})(aps::Vector{T}) where T<:AbstractElement
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

function (op::Operator{:Hₑ_up_fem})(apu::T) where {T<:AbstractElement,S<:AbstractElement}
    ΔW²= 0
    W̄² = 0
    ΔW²_dil= 0.0
    ΔW²_dev= 0.0
    W̄²_dil = 0.0
    W̄²_dev = 0.0
    Δu²= 0
    Δp²= 0
    ū² = 0
    p̄² =0
    p² =0
    E = op.E
    ν = op.ν
    K=E/3/(1-2ν)
    G=E/2/(1+ν)
    𝓖u = apu.𝓖
    for ξu in 𝓖u
    𝑤 = ξu.𝑤
        Nᵘ= ξu[:𝝭]
        B₁ = ξu[:∂𝝭∂x]
        B₂ = ξu[:∂𝝭∂y] 
        ū₁ = ξu.u
        ū₂ = ξu.v
        ∂ū₁∂x = ξu.∂u∂x
        ∂ū₁∂y = ξu.∂u∂y
        ∂ū₂∂x = ξu.∂v∂x
        ∂ū₂∂y = ξu.∂v∂y
        
        ε̄₁₁ = ∂ū₁∂x
        ε̄₂₂ = ∂ū₂∂y
        ε̄₁₂ = ∂ū₁∂y + ∂ū₂∂x
        σ̄₁₁ = K*(ε̄₁₁+ε̄₂₂)+2G*(2/3*ε̄₁₁-1/3*ε̄₂₂)
        σ̄₂₂ = K*(ε̄₁₁+ε̄₂₂)+2G*(-1/3*ε̄₁₁+2/3*ε̄₂₂)
        σ̄₁₂ = 2G*ε̄₁₂
        σ̄₃₃ = K*(ε̄₁₁+ε̄₂₂)+2G*(-1/3*ε̄₁₁-1/3*ε̄₂₂)

        p̄   = K*(ε̄₁₁+ε̄₂₂)
        s̄₁₁ = 2G*(2/3*ε̄₁₁-1/3*ε̄₂₂)
        s̄₂₂ = 2G*(-1/3*ε̄₁₁+2/3*ε̄₂₂)
        s̄₁₂ = 2G*ε̄₁₂
        s̄₃₃ = 2G*(-1/3*ε̄₁₁-1/3*ε̄₂₂)
        u₁ = 0.
        u₂ = 0.
        ε₁₁ = 0.
        ε₂₂ = 0.
        ε₁₂ = 0.
        for (i,xᵢ) in enumerate(apu.𝓒)
            
            u₁ += Nᵘ[i]*xᵢ.d₁
            u₂ += Nᵘ[i]*xᵢ.d₂
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += B₂[i]*xᵢ.d₁ + B₁[i]*xᵢ.d₂
        end
        p = K*(ε₁₁+ε₂₂)
        s₁₁ = 2G*(2/3*ε₁₁-1/3*ε₂₂)
        s₂₂ = 2G*(-1/3*ε₁₁+2/3*ε₂₂)
        s₁₂ = 2G*ε₁₂
        s₃₃ = 2G*(-1/3*ε₁₁-1/3*ε₂₂)
        ΔW² += (3*(p-p̄)^2/2/K+(s₁₁-s̄₁₁)^2/4/G+(s₂₂-s̄₂₂)^2/4/G+(s₃₃-s̄₃₃)^2/4/G+(s₁₂-s̄₁₂)^2/4/G)*𝑤
        W̄² += (3*p̄^2/2/K + s̄₁₁^2/4/G + s̄₂₂^2/4/G + s̄₃₃^2/4/G + s̄₁₂^2/4/G)*𝑤
        # ΔW² += (3*(p-p̄)*(ε₁₁-ε̄₁₁+ε₂₂-ε̄₂₂)/2+(s₁₁-s̄₁₁)*(2/3*(ε₁₁-ε̄₁₁)-1/3*(ε₂₂-ε̄₂₂))/2+(s₂₂-s̄₂₂)*(-1/3*(ε₁₁-ε̄₁₁)+2/3*(ε₂₂-ε̄₂₂))/2+(s₃₃-s̄₃₃)*(-1/3*(ε₁₁-ε̄₁₁)-1/3*(ε₂₂-ε̄₂₂))/2+(s₁₂-s̄₁₂)*(ε₁₂-ε̄₁₂)/2)*𝑤
        # W̄² += (3*p̄*(ε̄₁₁+ε̄₂₂)/2+s̄₁₁*(2/3*ε̄₁₁-1/3*ε̄₂₂)/2+s̄₂₂*(-1/3*ε̄₁₁+2/3*ε̄₂₂)/2+s̄₃₃*(-1/3*ε̄₁₁-1/3*ε̄₂₂)/2+s̄₁₂*ε̄₁₂/2)*𝑤
        Δu² += ((u₁ - ū₁)^2 + (u₂ - ū₂)^2)*𝑤
        ū² += (ū₁^2 + ū₂^2)*𝑤
        Δp² += ((p - p̄)^2 )*𝑤
        p̄² += (p̄^2)*𝑤
        # ΔW²_dil += (3*(p-p̄)*(ε₁₁-ε̄₁₁+ε₂₂-ε̄₂₂)/2)*𝑤
        ΔW²_dil +=(3*(p-p̄)^2/2/K)*𝑤
        # ΔW²_dev += ((s₁₁-s̄₁₁)*(2/3*(ε₁₁-ε̄₁₁)-1/3*(ε₂₂-ε̄₂₂))/2+(s₂₂-s̄₂₂)*(-1/3*(ε₁₁-ε̄₁₁)+2/3*(ε₂₂-ε̄₂₂))/2+(s₃₃-s̄₃₃)*(-1/3*(ε₁₁-ε̄₁₁)-1/3*(ε₂₂-ε̄₂₂))/2+(s₁₂-s̄₁₂)*(ε₁₂-ε̄₁₂)/2)*𝑤
        ΔW²_dev +=((s₁₁-s̄₁₁)^2/4/G+(s₂₂-s̄₂₂)^2/4/G+(s₃₃-s̄₃₃)^2/4/G+(s₁₂-s̄₁₂)^2/4/G)*𝑤
        W̄²_dil +=(3*p̄^2/2/K)*𝑤
        W̄²_dev +=(s̄₁₁^2/4/G + s̄₂₂^2/4/G + s̄₃₃^2/4/G + s̄₁₂^2/4/G)*𝑤
        # p²  +=p^2*𝑤
        # println( p̄²)
        # println( p²)
        

    end
    return ΔW², W̄², Δu², ū²,ΔW²_dil, ΔW²_dev, W̄²_dil,W̄²_dev,Δp²,p̄²
end
function (op::Operator{:Hₑ_up_mix})(apu::T,app::S) where {T<:AbstractElement,S<:AbstractElement}
    ΔW²= 0
    W̄² = 0
    ΔW²_dil= 0.0
    ΔW²_dev= 0.0
    W̄²_dil = 0.0
    W̄²_dev = 0.0
    Δu²= 0
    Δp²= 0
    ū² = 0
    p̄² =0
    p² =0
    E = op.E
    ν = op.ν
    K=E/3/(1-2ν)
    G=E/2/(1+ν)
    𝓖u = apu.𝓖
    𝓖p = app.𝓖
    for (ξu,ξp) in zip(𝓖u,𝓖p)
    𝑤 = ξu.𝑤
        Nᵖ= ξp[:𝝭]
        Nᵘ= ξu[:𝝭]
        B₁ = ξu[:∂𝝭∂x]
        B₂ = ξu[:∂𝝭∂y] 
        ū₁ = ξu.u
        ū₂ = ξu.v
        ∂ū₁∂x = ξu.∂u∂x
        ∂ū₁∂y = ξu.∂u∂y
        ∂ū₂∂x = ξu.∂v∂x
        ∂ū₂∂y = ξu.∂v∂y
        
        ε̄₁₁ = ∂ū₁∂x
        ε̄₂₂ = ∂ū₂∂y
        ε̄₁₂ = ∂ū₁∂y + ∂ū₂∂x
        σ̄₁₁ = K*(ε̄₁₁+ε̄₂₂)+2G*(2/3*ε̄₁₁-1/3*ε̄₂₂)
        σ̄₂₂ = K*(ε̄₁₁+ε̄₂₂)+2G*(-1/3*ε̄₁₁+2/3*ε̄₂₂)
        σ̄₁₂ = 2G*ε̄₁₂
        σ̄₃₃ = K*(ε̄₁₁+ε̄₂₂)+2G*(-1/3*ε̄₁₁-1/3*ε̄₂₂)

        p̄   = K*(ε̄₁₁+ε̄₂₂)
        s̄₁₁ = 2G*(2/3*ε̄₁₁-1/3*ε̄₂₂)
        s̄₂₂ = 2G*(-1/3*ε̄₁₁+2/3*ε̄₂₂)
        s̄₁₂ = 2G*ε̄₁₂
        s̄₃₃ = 2G*(-1/3*ε̄₁₁-1/3*ε̄₂₂)
        p  = 0.
        u₁ = 0.
        u₂ = 0.
        ε₁₁ = 0.
        ε₂₂ = 0.
        ε₁₂ = 0.
        for (i,xᵢ) in enumerate(apu.𝓒)
            
            u₁ += Nᵘ[i]*xᵢ.d₁
            u₂ += Nᵘ[i]*xᵢ.d₂
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += B₂[i]*xᵢ.d₁ + B₁[i]*xᵢ.d₂
        end
        for (j,xⱼ) in enumerate(app.𝓒)
            p  += Nᵖ[j]*xⱼ.q
        end
        s₁₁ = 2G*(2/3*ε₁₁-1/3*ε₂₂)
        s₂₂ = 2G*(-1/3*ε₁₁+2/3*ε₂₂)
        s₁₂ = 2G*ε₁₂
        s₃₃ = 2G*(-1/3*ε₁₁-1/3*ε₂₂)
        ΔW² += (3*(p-p̄)^2/2/K+(s₁₁-s̄₁₁)^2/4/G+(s₂₂-s̄₂₂)^2/4/G+(s₃₃-s̄₃₃)^2/4/G+(s₁₂-s̄₁₂)^2/4/G)*𝑤
        W̄² += (3*p̄^2/2/K + s̄₁₁^2/4/G + s̄₂₂^2/4/G + s̄₃₃^2/4/G + s̄₁₂^2/4/G)*𝑤
        # ΔW² += (3*(p-p̄)*(ε₁₁-ε̄₁₁+ε₂₂-ε̄₂₂)/2+(s₁₁-s̄₁₁)*(2/3*(ε₁₁-ε̄₁₁)-1/3*(ε₂₂-ε̄₂₂))/2+(s₂₂-s̄₂₂)*(-1/3*(ε₁₁-ε̄₁₁)+2/3*(ε₂₂-ε̄₂₂))/2+(s₃₃-s̄₃₃)*(-1/3*(ε₁₁-ε̄₁₁)-1/3*(ε₂₂-ε̄₂₂))/2+(s₁₂-s̄₁₂)*(ε₁₂-ε̄₁₂)/2)*𝑤
        # W̄² += (3*p̄*(ε̄₁₁+ε̄₂₂)/2+s̄₁₁*(2/3*ε̄₁₁-1/3*ε̄₂₂)/2+s̄₂₂*(-1/3*ε̄₁₁+2/3*ε̄₂₂)/2+s̄₃₃*(-1/3*ε̄₁₁-1/3*ε̄₂₂)/2+s̄₁₂*ε̄₁₂/2)*𝑤
        Δu² += ((u₁ - ū₁)^2 + (u₂ - ū₂)^2)*𝑤
        ū² += (ū₁^2 + ū₂^2)*𝑤
        Δp² += ((p - p̄)^2 )*𝑤
        p̄² += (p̄^2)*𝑤
        # ΔW²_dil += (3*(p-p̄)*(ε₁₁-ε̄₁₁+ε₂₂-ε̄₂₂)/2)*𝑤
        ΔW²_dil +=(3*(p-p̄)^2/2/K)*𝑤
        # ΔW²_dev += ((s₁₁-s̄₁₁)*(2/3*(ε₁₁-ε̄₁₁)-1/3*(ε₂₂-ε̄₂₂))/2+(s₂₂-s̄₂₂)*(-1/3*(ε₁₁-ε̄₁₁)+2/3*(ε₂₂-ε̄₂₂))/2+(s₃₃-s̄₃₃)*(-1/3*(ε₁₁-ε̄₁₁)-1/3*(ε₂₂-ε̄₂₂))/2+(s₁₂-s̄₁₂)*(ε₁₂-ε̄₁₂)/2)*𝑤
        ΔW²_dev +=((s₁₁-s̄₁₁)^2/4/G+(s₂₂-s̄₂₂)^2/4/G+(s₃₃-s̄₃₃)^2/4/G+(s₁₂-s̄₁₂)^2/4/G)*𝑤
        W̄²_dil +=(3*p̄^2/2/K)*𝑤
        W̄²_dev +=(s̄₁₁^2/4/G + s̄₂₂^2/4/G + s̄₃₃^2/4/G + s̄₁₂^2/4/G)*𝑤
        # p²  +=p^2*𝑤
        # println( p̄²)
        # println( p²)
        

    end
    return ΔW², W̄², Δu², ū²,ΔW²_dil, ΔW²_dev, W̄²_dil,W̄²_dev,Δp²,p̄²
end

function (op::Operator{:Hₑ_up_fem})(apsu::Vector{T}) where {T<:AbstractElement,S<:AbstractElement}
    HₑNorm_ΔW²= 0.0
    HₑNorm_W̄² = 0.0
    L₂Norm_Δp²= 0.0
    L₂Norm_p̄² = 0.0
    L₂Norm_Δu²= 0.0
    L₂Norm_ū² = 0.0
    HₑNorm_ΔW²_dil= 0.0
    HₑNorm_ΔW²_dev= 0.0
    HₑNorm_W̄²_dil = 0.0
    HₑNorm_W̄²_dev = 0.0
    for apu in apsu
    # for apu in apsu
        # for app in apsp
            ΔW², W̄², Δu², ū²,ΔW²_dil, ΔW²_dev, W̄²_dil,W̄²_dev, Δp²,p̄²  = op(apu)
            HₑNorm_ΔW² += ΔW²
            HₑNorm_W̄²  += W̄²
            L₂Norm_Δu² += Δu²
            L₂Norm_ū²  += ū²
            L₂Norm_Δp² += Δp²
            L₂Norm_p̄²  += p̄²
            
            HₑNorm_ΔW²_dil += ΔW²_dil
            HₑNorm_ΔW²_dev += ΔW²_dev
            HₑNorm_W̄²_dil  += W̄²_dil
            HₑNorm_W̄²_dev  += W̄²_dev
        # end
    end
    return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5,(HₑNorm_ΔW²_dil/HₑNorm_W̄²_dil)^0.5,(HₑNorm_ΔW²_dev/HₑNorm_W̄²_dev)^0.5,(L₂Norm_Δp²/L₂Norm_p̄²)^0.5
    # return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end    
function (op::Operator{:Hₑ_up_mix})(apsu::Vector{T},apsp::Vector{S}) where {T<:AbstractElement,S<:AbstractElement}
    HₑNorm_ΔW²= 0.0
    HₑNorm_W̄² = 0.0
    L₂Norm_Δp²= 0.0
    L₂Norm_p̄² = 0.0
    L₂Norm_Δu²= 0.0
    L₂Norm_ū² = 0.0
    HₑNorm_ΔW²_dil= 0.0
    HₑNorm_ΔW²_dev= 0.0
    HₑNorm_W̄²_dil = 0.0
    HₑNorm_W̄²_dev = 0.0
    for (apu,app) in zip(apsu,apsp)
    # for apu in apsu
        # for app in apsp
            ΔW², W̄², Δu², ū²,ΔW²_dil, ΔW²_dev, W̄²_dil,W̄²_dev, Δp²,p̄²  = op(apu,app)
            HₑNorm_ΔW² += ΔW²
            HₑNorm_W̄²  += W̄²
            L₂Norm_Δu² += Δu²
            L₂Norm_ū²  += ū²
            L₂Norm_Δp² += Δp²
            L₂Norm_p̄²  += p̄²
            
            HₑNorm_ΔW²_dil += ΔW²_dil
            HₑNorm_ΔW²_dev += ΔW²_dev
            HₑNorm_W̄²_dil  += W̄²_dil
            HₑNorm_W̄²_dev  += W̄²_dev
        # end
    end
    return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5,(HₑNorm_ΔW²_dil/HₑNorm_W̄²_dil)^0.5,(HₑNorm_ΔW²_dev/HₑNorm_W̄²_dev)^0.5,(L₂Norm_Δp²/L₂Norm_p̄²)^0.5
    # return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end    

function (op::Operator{:Hₑ_up_mix_Q4P1})(apu::T,app::S) where {T<:AbstractElement,S<:AbstractElement}
    ΔW²= 0
    W̄² = 0
    ΔW²_dil= 0.0
    ΔW²_dev= 0.0
    W̄²_dil = 0.0
    W̄²_dev = 0.0
    Δu²= 0
    Δp²= 0
    ū² = 0
    p̄² =0
    p² =0
    E = op.E
    ν = op.ν
    p = op.p
    K=E/3/(1-2ν)
    G=E/2/(1+ν)
    𝓖u = apu.𝓖
    𝓖p = app.𝓖
    
    for (ξu,ξp) in zip(𝓖u,𝓖p)
    𝑤 = ξu.𝑤
        Nᵖ= ξp[:𝝭]
        Nᵘ= ξu[:𝝭]
        B₁ = ξu[:∂𝝭∂x]
        B₂ = ξu[:∂𝝭∂y] 
        ū₁ = ξu.u
        ū₂ = ξu.v
        ∂ū₁∂x = ξu.∂u∂x
        ∂ū₁∂y = ξu.∂u∂y
        ∂ū₂∂x = ξu.∂v∂x
        ∂ū₂∂y = ξu.∂v∂y
        
        ε̄₁₁ = ∂ū₁∂x
        ε̄₂₂ = ∂ū₂∂y
        ε̄₁₂ = ∂ū₁∂y + ∂ū₂∂x
        σ̄₁₁ = K*(ε̄₁₁+ε̄₂₂)+2G*(2/3*ε̄₁₁-1/3*ε̄₂₂)
        σ̄₂₂ = K*(ε̄₁₁+ε̄₂₂)+2G*(-1/3*ε̄₁₁+2/3*ε̄₂₂)
        σ̄₁₂ = 2G*ε̄₁₂
        σ̄₃₃ = K*(ε̄₁₁+ε̄₂₂)+2G*(-1/3*ε̄₁₁-1/3*ε̄₂₂)

        p̄   = K*(ε̄₁₁+ε̄₂₂)
        s̄₁₁ = 2G*(2/3*ε̄₁₁-1/3*ε̄₂₂)
        s̄₂₂ = 2G*(-1/3*ε̄₁₁+2/3*ε̄₂₂)
        s̄₁₂ = 2G*ε̄₁₂
        s̄₃₃ = 2G*(-1/3*ε̄₁₁-1/3*ε̄₂₂)
        p  = 0.
        u₁ = 0.
        u₂ = 0.
        ε₁₁ = 0.
        ε₂₂ = 0.
        ε₁₂ = 0.
        for (i,xᵢ) in enumerate(apu.𝓒)
            
            u₁ += Nᵘ[i]*xᵢ.d₁
            u₂ += Nᵘ[i]*xᵢ.d₂
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += B₂[i]*xᵢ.d₁ + B₁[i]*xᵢ.d₂
        end
        
        p  = p
        
        s₁₁ = 2G*(2/3*ε₁₁-1/3*ε₂₂)
        s₂₂ = 2G*(-1/3*ε₁₁+2/3*ε₂₂)
        s₁₂ = 2G*ε₁₂
        s₃₃ = 2G*(-1/3*ε₁₁-1/3*ε₂₂)
        ΔW² += (3*(p-p̄)^2/2/K+(s₁₁-s̄₁₁)^2/4/G+(s₂₂-s̄₂₂)^2/4/G+(s₃₃-s̄₃₃)^2/4/G+(s₁₂-s̄₁₂)^2/4/G)*𝑤
        W̄² += (3*p̄^2/2/K + s̄₁₁^2/4/G + s̄₂₂^2/4/G + s̄₃₃^2/4/G + s̄₁₂^2/4/G)*𝑤
        # ΔW² += (3*(p-p̄)*(ε₁₁-ε̄₁₁+ε₂₂-ε̄₂₂)/2+(s₁₁-s̄₁₁)*(2/3*(ε₁₁-ε̄₁₁)-1/3*(ε₂₂-ε̄₂₂))/2+(s₂₂-s̄₂₂)*(-1/3*(ε₁₁-ε̄₁₁)+2/3*(ε₂₂-ε̄₂₂))/2+(s₃₃-s̄₃₃)*(-1/3*(ε₁₁-ε̄₁₁)-1/3*(ε₂₂-ε̄₂₂))/2+(s₁₂-s̄₁₂)*(ε₁₂-ε̄₁₂)/2)*𝑤
        # W̄² += (3*p̄*(ε̄₁₁+ε̄₂₂)/2+s̄₁₁*(2/3*ε̄₁₁-1/3*ε̄₂₂)/2+s̄₂₂*(-1/3*ε̄₁₁+2/3*ε̄₂₂)/2+s̄₃₃*(-1/3*ε̄₁₁-1/3*ε̄₂₂)/2+s̄₁₂*ε̄₁₂/2)*𝑤
        Δu² += ((u₁ - ū₁)^2 + (u₂ - ū₂)^2)*𝑤
        ū² += (ū₁^2 + ū₂^2)*𝑤
        Δp² += ((p - p̄)^2 )*𝑤
        p̄² += (p̄^2)*𝑤
        # ΔW²_dil += (3*(p-p̄)*(ε₁₁-ε̄₁₁+ε₂₂-ε̄₂₂)/2)*𝑤
        ΔW²_dil +=(3*(p-p̄)^2/2/K)*𝑤
        # ΔW²_dev += ((s₁₁-s̄₁₁)*(2/3*(ε₁₁-ε̄₁₁)-1/3*(ε₂₂-ε̄₂₂))/2+(s₂₂-s̄₂₂)*(-1/3*(ε₁₁-ε̄₁₁)+2/3*(ε₂₂-ε̄₂₂))/2+(s₃₃-s̄₃₃)*(-1/3*(ε₁₁-ε̄₁₁)-1/3*(ε₂₂-ε̄₂₂))/2+(s₁₂-s̄₁₂)*(ε₁₂-ε̄₁₂)/2)*𝑤
        ΔW²_dev +=((s₁₁-s̄₁₁)^2/4/G+(s₂₂-s̄₂₂)^2/4/G+(s₃₃-s̄₃₃)^2/4/G+(s₁₂-s̄₁₂)^2/4/G)*𝑤
        W̄²_dil +=(3*p̄^2/2/K)*𝑤
        W̄²_dev +=(s̄₁₁^2/4/G + s̄₂₂^2/4/G + s̄₃₃^2/4/G + s̄₁₂^2/4/G)*𝑤
        # p²  +=p^2*𝑤
        # println( p̄²)
        # println( p²)
        

    end
    return ΔW², W̄², Δu², ū²,ΔW²_dil, ΔW²_dev, W̄²_dil,W̄²_dev,Δp²,p̄²
end

function (op::Operator{:Hₑ_up_mix_Q4P1})(apsu::Vector{T},apsp::Vector{S}) where {T<:AbstractElement,S<:AbstractElement}
    HₑNorm_ΔW²= 0.0
    HₑNorm_W̄² = 0.0
    L₂Norm_Δp²= 0.0
    L₂Norm_p̄² = 0.0
    L₂Norm_Δu²= 0.0
    L₂Norm_ū² = 0.0
    HₑNorm_ΔW²_dil= 0.0
    HₑNorm_ΔW²_dev= 0.0
    HₑNorm_W̄²_dil = 0.0
    HₑNorm_W̄²_dev = 0.0
    for (apu,app) in zip(apsu,apsp)
    # for apu in apsu
        # for app in apsp
            ΔW², W̄², Δu², ū²,ΔW²_dil, ΔW²_dev, W̄²_dil,W̄²_dev, Δp²,p̄²  = op(apu,app)
            HₑNorm_ΔW² += ΔW²
            HₑNorm_W̄²  += W̄²
            L₂Norm_Δu² += Δu²
            L₂Norm_ū²  += ū²
            L₂Norm_Δp² += Δp²
            L₂Norm_p̄²  += p̄²
            
            HₑNorm_ΔW²_dil += ΔW²_dil
            HₑNorm_ΔW²_dev += ΔW²_dev
            HₑNorm_W̄²_dil  += W̄²_dil
            HₑNorm_W̄²_dev  += W̄²_dev
        # end
    end
    return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5,(HₑNorm_ΔW²_dil/HₑNorm_W̄²_dil)^0.5,(HₑNorm_ΔW²_dev/HₑNorm_W̄²_dev)^0.5,(L₂Norm_Δp²/L₂Norm_p̄²)^0.5
    # return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end    
function (op::Operator{:Hₑ_Incompressible})(ap::T) where T<:AbstractElement
    ΔW²= 0
    W̄² = 0
    Δu²= 0
    ū² = 0
    E = op.E
    ν = op.ν
    Cᵈ = E/(1+ν)
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
        ε̄₁₂ = 0.5*(∂ū₁∂y + ∂ū₂∂x)
        σ̄₁₁ = Cᵈ*(2*ε̄₁₁ -   ε̄₂₂)/3
        σ̄₂₂ = Cᵈ*(- ε̄₁₁ + 2*ε̄₂₂)/3
        σ̄₁₂ = Cᵈ*ε̄₁₂
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
            ε₁₂ += 0.5*(B₂[i]*xᵢ.d₁ + B₁[i]*xᵢ.d₂)
        end
        σ₁₁ = Cᵈ*(2*ε₁₁ -   ε₂₂)/3
        σ₂₂ = Cᵈ*(- ε₁₁ + 2*ε₂₂)/3
        σ₁₂ = Cᵈ*ε₁₂
        ΔW² += 0.5*((σ₁₁-σ̄₁₁)*(ε₁₁-ε̄₁₁) + (σ₂₂-σ̄₂₂)*(ε₂₂-ε̄₂₂) + (σ₁₂-σ̄₁₂)*(ε₁₂-ε̄₁₂))*𝑤
        W̄² += 0.5*(σ̄₁₁*ε̄₁₁ + σ̄₂₂*ε̄₂₂ + σ̄₁₂*ε̄₁₂)*𝑤
        Δu² += ((u₁ - ū₁)^2 + (u₂ - ū₂)^2)*𝑤
        ū² += (ū₁^2 + ū₂^2)*𝑤
    end
    return ΔW², W̄², Δu², ū²
end

function (op::Operator{:Hₑ_Incompressible})(aps::Vector{T}) where T<:AbstractElement
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