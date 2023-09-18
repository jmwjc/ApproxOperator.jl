"""
1D Plasticity
"""
function (op::Operator{:∫vₓσdx})(ap::T,k::AbstractMatrix{Float64},fint::AbstractVector) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = op.E
    K = op.K
    σy = op.σy
    for ξ in 𝓖
        B = ξ[:∂𝝭∂x]
        σₙ = ξ.σₙ
        αₙ = ξ.αₙ
        εᵖₙ = ξ.εᵖₙ
        𝑤 = ξ.𝑤
        Δεₙ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            Δεₙ += B[i]*xᵢ.Δd
        end
        # predict phase
        σᵗʳ = σₙ+E*Δεₙ
        fᵗʳ = abs(σᵗʳ) - (σy+K*αₙ)
        if fᵗʳ > 1E-13
            Δγ = fᵗʳ/(E+K)
            ξ.σₙ = σᵗʳ - Δγ*E*sign(σᵗʳ)
            ξ.εᵖₙ = εᵖₙ + Δγ*sign(σᵗʳ)
            ξ.αₙ = αₙ + Δγ
            Eₜ = (E*K)/(E+K)
        else
            ξ.σₙ = σᵗʳ
            Eₜ = E
        end
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += B[i]*Eₜ*B[j]*𝑤
            end
            fint[I] += B[i]*ξ.σₙ*𝑤
        end
    end
end

"""
Tresca
"""
function (op::Operator{:∫vᵢσdΩ_tresca})(ap::T;k::AbstractMatrix{Float64},fint::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    λ = op.λ
    μ = op.μ
    c = op.c
    Cᵢᵢᵢᵢ = λ + 2.0*μ
    Cᵢᵢⱼⱼ = λ
    Cᵢⱼᵢⱼ = μ
    
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        σ₁₁ = ξ.σ₁₁
        σ₂₂ = ξ.σ₂₂
        σ₁₂ = ξ.σ₁₂
        εᵖ₁₁ = ξ.εᵖ₁₁
        εᵖ₂₂ = ξ.εᵖ₂₂
        εᵖ₁₂ = ξ.εᵖ₁₂
        𝑤 = ξ.𝑤
        Δε₁₁ = 0.0
        Δε₂₂ = 0.0
        Δε₁₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            Δε₁₁ += B₁[i]*xᵢ.Δd₁
            Δε₂₂ += B₂[i]*xᵢ.Δd₂
            Δε₁₂ += B₁[i]*xᵢ.Δd₂ + B₂[i]*xᵢ.Δd₁
        end 

        # predict phase
        σ₁₁ᵗʳ = σ₁₁ + Cᵢᵢᵢᵢ*Δε₁₁+Cᵢᵢⱼⱼ*Δε₂₂
        σ₂₂ᵗʳ = σ₂₂ + Cᵢᵢᵢᵢ*Δε₂₂+Cᵢᵢⱼⱼ*Δε₁₁
        σ₁₂ᵗʳ = σ₁₂ + Cᵢⱼᵢⱼ*Δε₁₂
        σ₁,σ₂,n₁,n₂ = getσₙ(σ₁₁ᵗʳ,σ₂₂ᵗʳ,σ₁₂ᵗʳ)
        fᵗʳ = σ₁ - σ₂ - 2.0*c

        if fᵗʳ > 1E-13
            Δγ = fᵗʳ/4.0/μ
           
            ξ.σ₁₁ = σ₁₁ᵗʳ - Δγ*2.0*μ*(n₁[1]*n₁[1] - n₂[1]*n₂[1])
            ξ.σ₂₂ = σ₂₂ᵗʳ - Δγ*2.0*μ*(n₁[2]*n₁[2] - n₂[2]*n₂[2])
            ξ.σ₁₂ = σ₁₂ᵗʳ - Δγ*2.0*μ*(n₁[1]*n₁[2] - n₂[1]*n₂[2])
            ξ.εᵖ₁₁ = εᵖ₁₁ + Δγ*(n₁[1]*n₁[1] - n₂[1]*n₂[1])
            ξ.εᵖ₂₂ = εᵖ₂₂ + Δγ*(n₁[2]*n₁[2] - n₂[2]*n₂[2])
            ξ.εᵖ₁₂ = εᵖ₁₂ + Δγ*(n₁[1]*n₁[2] - n₂[1]*n₂[2])

            Cᵗ₁₁₁₁ = Cᵢᵢᵢᵢ - μ*(n₁[1]*n₁[1] - n₂[1]*n₂[1])^2
            Cᵗ₂₂₂₂ = Cᵢᵢᵢᵢ - μ*(n₁[2]*n₁[2] - n₂[2]*n₂[2])^2
            Cᵗ₁₁₂₂ = Cᵢᵢⱼⱼ - μ*(n₁[1]*n₁[1] - n₂[1]*n₂[1])*(n₁[2]*n₁[2] - n₂[2]*n₂[2])
            Cᵗ₁₂₁₂ = Cᵢⱼᵢⱼ - 2.0*μ*(n₁[1]*n₁[2] - n₂[1]*n₂[2])^2
            Cᵗ₁₁₁₂ = - μ*(n₁[1]*n₁[1] - n₂[1]*n₂[1])*(n₁[1]*n₁[2] - n₂[1]*n₂[2])
            Cᵗ₂₂₁₂ = - μ*(n₁[2]*n₁[2] - n₂[2]*n₂[2])*(n₁[1]*n₁[2] - n₂[1]*n₂[2])
             #println(Cᵗ₁₁₁₁)
             #println(Cᵗ₂₂₂₂)
             #println(Cᵗ₁₁₂₂)
             #println(Cᵗ₁₂₁₂)
             #println(Cᵗ₁₁₁₂)
             #println(Cᵗ₂₂₁₂)
        else
            ξ.σ₁₁ = σ₁₁ᵗʳ
            ξ.σ₂₂ = σ₂₂ᵗʳ
            ξ.σ₁₂ = σ₁₂ᵗʳ
            Cᵗ₁₁₁₁ = Cᵢᵢᵢᵢ
            Cᵗ₂₂₂₂ = Cᵢᵢᵢᵢ
            Cᵗ₁₁₂₂ = Cᵢᵢⱼⱼ
            Cᵗ₁₂₁₂ = Cᵢⱼᵢⱼ
            Cᵗ₁₁₁₂ = 0.0
            Cᵗ₂₂₁₂ = 0.0
        end
       
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += (Cᵗ₁₁₁₁*B₁[i]*B₁[j] + Cᵗ₁₂₁₂*B₂[i]*B₂[j] + Cᵗ₁₁₁₂*(B₁[i]*B₂[j]+B₂[i]*B₁[j]))*𝑤
                k[2*I-1,2*J]   += (Cᵗ₁₁₂₂*B₁[i]*B₂[j] + Cᵗ₁₂₁₂*B₂[i]*B₁[j] + Cᵗ₁₁₁₂*B₁[i]*B₁[j] + Cᵗ₂₂₁₂*B₂[i]*B₂[j])*𝑤
                k[2*I,2*J-1]   += (Cᵗ₁₁₂₂*B₂[i]*B₁[j] + Cᵗ₁₂₁₂*B₁[i]*B₂[j] + Cᵗ₁₁₁₂*B₁[i]*B₁[j] + Cᵗ₂₂₁₂*B₂[i]*B₂[j])*𝑤
                k[2*I,2*J]     += (Cᵗ₂₂₂₂*B₂[i]*B₂[j] + Cᵗ₁₂₁₂*B₁[i]*B₁[j] + Cᵗ₂₂₁₂*(B₁[i]*B₂[j]+B₂[i]*B₁[j]))*𝑤
     
            end
            fint[2*I-1] += B₁[i]*ξ.σ₁₁*𝑤 + B₂[i]*ξ.σ₁₂*𝑤
            fint[2*I]   += B₁[i]*ξ.σ₁₂*𝑤 + B₂[i]*ξ.σ₂₂*𝑤
        end
    end
end    

"""
morh-coulbom
"""
function (op::Operator{:∫vᵢσdΩ_mohr_coulomb})(ap::T;k::AbstractMatrix{Float64},fint::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    λ = op.λ
    μ = op.μ
    c = op.c
    𝜙 = op.𝜙
    tol = op.tol 
    Cᵢᵢᵢᵢ = λ + 2.0*μ
    Cᵢᵢⱼⱼ = λ
    Cᵢⱼᵢⱼ = μ
    
    for ξ in enumerate(𝓖)
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        σ₁₁ = ξ.σ₁₁
        σ₂₂ = ξ.σ₂₂
        σ₃₃ = ξ.σ₃₃
        σ₁₂ = ξ.σ₁₂
        εᵖ₁₁ = ξ.εᵖ₁₁
        εᵖ₂₂ = ξ.εᵖ₂₂
        εᵖ₁₂ = ξ.εᵖ₁₂
        𝑤 = ξ.𝑤
        Δε₁₁ = 0.0
        Δε₂₂ = 0.0
        Δε₁₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            Δε₁₁ += B₁[i]*xᵢ.Δd₁
            Δε₂₂ += B₂[i]*xᵢ.Δd₂
            Δε₁₂ += B₁[i]*xᵢ.Δd₂ + B₂[i]*xᵢ.Δd₁
        end 
        # predict phase
        σ₁₁ᵗʳ = σ₁₁ + Cᵢᵢᵢᵢ*Δε₁₁+Cᵢᵢⱼⱼ*Δε₂₂
        σ₂₂ᵗʳ = σ₂₂ + Cᵢᵢᵢᵢ*Δε₂₂+Cᵢᵢⱼⱼ*Δε₁₁
        σ₁₂ᵗʳ = σ₁₂ + Cᵢⱼᵢⱼ*Δε₁₂
        σ₁,σ₂,n₁,n₂ = getσₙ(σ₁₁ᵗʳ,σ₂₂ᵗʳ,σ₁₂ᵗʳ)
        fᵗʳ = σ₁-σ₂ + (σ₁+σ₂) * sin(𝜙) - 2.0*c*cos(𝜙)
        if fᵗʳ > tol
            sin𝜙 = sin(𝜙)
            sin²𝜙 = sin(𝜙)*sin(𝜙)
            sin²𝜙psin𝜙 = sin²𝜙 +sin(𝜙) 
            sin²𝜙dsin𝜙 = sin²𝜙 -sin(𝜙) 
            ∂fC∂f=(4.0*sin²𝜙*λ+4.0*(1.0+sin²𝜙)*μ)
            Δγ = fᵗʳ/∂fC∂f
           
            ξ.σ₁₁ = σ₁₁ᵗʳ - Δγ*(2*sin𝜙*λ+2*((sin𝜙+1)*n₁[1]*n₁[1]+(sin𝜙-1)*n₂[1]*n₂[1])*μ)#上标1对应之前下标1
            ξ.σ₂₂ = σ₂₂ᵗʳ - Δγ*(2*sin𝜙*λ+2*((sin𝜙+1)*n₁[2]*n₁[2]+(sin𝜙-1)*n₂[2]*n₂[2])*μ)
            ξ.σ₁₂ = σ₁₂ᵗʳ - Δγ*2*((sin𝜙+1)*n₁[1]*n₁[2]+(sin𝜙-1)*n₂[1]*n₂[2])*μ
            ξ.εᵖ₁₁ = εᵖ₁₁ + Δγ*((sin𝜙+1)*n₁[1]*n₁[1]+(sin𝜙-1)n₂[1]*n₂[1])
            ξ.εᵖ₂₂ = εᵖ₂₂ + Δγ*((sin𝜙+1)*n₁[2]*n₁[2]+(sin𝜙-1)n₂[2]*n₂[2])
            ξ.εᵖ₁₂ = εᵖ₁₂ + Δγ*(((sin𝜙+1)*n₁[1]*n₁[2]+(sin𝜙-1)n₂[1]*n₂[2]))

            Cᵗ₁₁₁₁=Cᵢᵢᵢᵢ-(4.0*sin²𝜙*λ^2 + 8.0*sin²𝜙psin𝜙*λ*μ*(n₁[1])^2 + 8.0*(sin²𝜙dsin𝜙)*λ*μ*(n₂[1])^2-8.0*cos(𝜙)*cos(𝜙)*μ^2*n₁[1]^2*n₂[1]^2+4.0*(sin𝜙+1)^2*μ^2*(n₁[1])^4+4.0*(sin𝜙-1)^2*μ^2*(n₂[1])^4)/∂fC∂f
            Cᵗ₂₂₂₂=Cᵢᵢᵢᵢ-(4.0*sin²𝜙*λ^2 + 8.0*sin²𝜙psin𝜙*λ*μ*(n₁[2])^2 + 8.0*(sin²𝜙dsin𝜙)*λ*μ*(n₂[2])^2-8.0*cos(𝜙)*cos(𝜙)*μ^2*n₁[2]^2*n₂[2]^2+4.0*(sin𝜙+1)^2*μ^2*(n₁[2])^4+4.0*(sin𝜙-1)^2*μ^2*(n₂[2])^4)/∂fC∂f
            Cᵗ₁₁₂₂=Cᵢᵢⱼⱼ-(4.0*sin²𝜙*λ^2 + 4.0*sin²𝜙psin𝜙*λ*μ*((n₁[1])^2+(n₁[2])^2)+4.0*(sin²𝜙dsin𝜙)*λ*μ*((n₂[1])^2+(n₂[2])^2)+4.0*(sin𝜙+1)^2*μ^2*(n₁[1])^2*(n₁[2])^2-4.0*cos(𝜙)*cos(𝜙)*μ^2*((n₂[1])^2*(n₁[2])^2+(n₁[1])^2*(n₂[2])^2)+4.0*(sin𝜙-1)^2*μ^2*(n₂[1])^2*(n₂[2])^2)/∂fC∂f
            Cᵗ₁₂₁₂=Cᵢⱼᵢⱼ-2.0(-8.0*cos(𝜙)*cos(𝜙)*μ^2*n₁[1]*n₁[2]*n₂[1]*n₂[2]+4.0*(sin𝜙+1)^2*μ^2*(n₁[1]*n₁[2])^2+4.0*(sin𝜙-1)^2*μ^2*(n₂[1]*n₂[2])^2)/∂fC∂f
            Cᵗ₁₁₁₂ = - (4.0*sin²𝜙psin𝜙*λ*μ*n₁[1]*n₁[2] + 4.0*(sin²𝜙dsin𝜙)*λ*μ*n₂[1]*n₂[2]-4.0*cos(𝜙)*cos(𝜙)*μ^2*n₂[1]^2*n₁[1]*n₁[2]-4.0*cos(𝜙)*cos(𝜙)*μ^2*n₁[1]^2*n₂[1]*n₂[2]+4.0*(sin𝜙+1)^2*μ^2*(n₁[1])^3*n₁[2]+4.0*(sin𝜙-1)^2*μ^2*(n₂[1])^3*n₂[2])/∂fC∂f
            Cᵗ₂₂₁₂ = - (4.0*sin²𝜙psin𝜙*λ*μ*n₁[2]*n₁[1] + 4.0*(sin²𝜙dsin𝜙)*λ*μ*n₂[2]*n₂[1]-4.0*cos(𝜙)*cos(𝜙)*μ^2*n₂[2]^2*n₁[2]*n₁[1]-4.0*cos(𝜙)*cos(𝜙)*μ^2*n₁[2]^2*n₂[2]*n₂[1]+4.0*(sin𝜙+1)^2*μ^2*(n₁[2])^3*n₁[1]+4.0*(sin𝜙-1)^2*μ^2*(n₂[2])^3*n₂[1])/∂fC∂f
            #println(Cᵗ₁₁₁₁)
            #println(Cᵗ₂₂₂₂)
            #println(Cᵗ₁₁₂₂)
            #println(Cᵗ₁₂₁₂)
            #println(Cᵗ₁₁₁₂)
            #println(Cᵗ₂₂₁₂)
        else
            ξ.σ₁₁ = σ₁₁ᵗʳ
            ξ.σ₂₂ = σ₂₂ᵗʳ
            ξ.σ₁₂ = σ₁₂ᵗʳ
            Cᵗ₁₁₁₁ = Cᵢᵢᵢᵢ
            Cᵗ₂₂₂₂ = Cᵢᵢᵢᵢ
            Cᵗ₁₁₂₂ = Cᵢᵢⱼⱼ
            Cᵗ₁₂₁₂ = Cᵢⱼᵢⱼ
            Cᵗ₁₁₁₂ = 0.0
            Cᵗ₂₂₁₂ = 0.0
        end
       # if isnan(σ₁₁)
       #    
       #       σ₁₁ = NaN
       #       error("程序终止：σ₁₁的值为NaN")
       # end
        # println(Cᵗ₁₁₁₁)
        # println(Cᵗ₂₂₂₂)
        # println(Cᵗ₁₁₂₂)
        # println(Cᵗ₁₂₁₂)
       
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += (Cᵗ₁₁₁₁*B₁[i]*B₁[j] + Cᵗ₁₂₁₂*B₂[i]*B₂[j] + Cᵗ₁₁₁₂*(B₁[i]*B₂[j]+B₂[i]*B₁[j]))*𝑤
                k[2*I-1,2*J]   += (Cᵗ₁₁₂₂*B₁[i]*B₂[j] + Cᵗ₁₂₁₂*B₂[i]*B₁[j] + Cᵗ₁₁₁₂*B₁[i]*B₁[j] + Cᵗ₂₂₁₂*B₂[i]*B₂[j])*𝑤
                k[2*I,2*J-1]   += (Cᵗ₁₁₂₂*B₂[i]*B₁[j] + Cᵗ₁₂₁₂*B₁[i]*B₂[j] + Cᵗ₁₁₁₂*B₁[i]*B₁[j] + Cᵗ₂₂₁₂*B₂[i]*B₂[j])*𝑤
                k[2*I,2*J]     += (Cᵗ₂₂₂₂*B₂[i]*B₂[j] + Cᵗ₁₂₁₂*B₁[i]*B₁[j] + Cᵗ₂₂₁₂*(B₁[i]*B₂[j]+B₂[i]*B₁[j]))*𝑤
     
            end
            fint[2*I-1] += B₁[i]*ξ.σ₁₁*𝑤 + B₂[i]*ξ.σ₁₂*𝑤
            fint[2*I]   += B₁[i]*ξ.σ₁₂*𝑤 + B₂[i]*ξ.σ₂₂*𝑤
    
        end
    end
end    


"""
frictional-contact
"""
function (op::Operator{:∫vᵢσdΩ_frictional_contact})(ap::T;k::AbstractMatrix{Float64},fint::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    μ̄ = op.μ̄ 
    E = op.E
    ν = op.ν
    λ = E*ν/(1.0+ν)/(1.0-2.0*ν)
    μ = 0.5*E/(1.0+ν)
    η = 1e-6
    nodes = op.nodes
    
    for ξ in 𝓖
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        σ₁₁ = ξ.σ₁₁
        σ₂₂ = ξ.σ₂₂
        σ₁₂ = ξ.σ₁₂
        𝑤   = ξ.𝑤
        ε₁₁ = 0.0
        ε₂₂ = 0.0
        ε₁₂ = 0.0
        Δε₁₁ = 0.0
        Δε₂₂ = 0.0
        Δε₁₂ = 0.0
        v = 0.0
        ∂v∂x = 0.0
        ∂v∂y = 0.0
      
        for (i,xᵢ) in enumerate(𝓒)
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += B₁[i]*xᵢ.d₂ + B₂[i]*xᵢ.d₁
            Δε₁₁ += B₁[i]*xᵢ.Δd₁
            Δε₂₂ += B₂[i]*xᵢ.Δd₂
            Δε₁₂ += B₁[i]*xᵢ.Δd₂ + B₂[i]*xᵢ.Δd₁
            v += N[i]*xᵢ.v
            ∂v∂x += B₁[i]*xᵢ.v
            ∂v∂y += B₂[i]*xᵢ.v
        end 

        #norm∇v = (∂v∂x^2+∂v∂y^2)^0.5
        #n₁ = ∂v∂x/norm∇v
        #n₂ = ∂v∂y/norm∇v
        #s₁ = n₂
        #s₂ = -n₁
        Γtmp = Set{Int}()
        Γdam = Set{Int}()
        Γfinal = Set{Int}()
        #step1 找到的所有v>0.98的点
        for (i,xᵢ) in enumerate(𝓒)
            if xᵢ.v > 0.98   
                #step2 并将其存入Γtmp中
                push!(Γtmp, i)
            end
            if xᵢ.v > 0.05   
                push!(Γdam, i)
            end

        end
    
       
        while !isempty(Γtmp)
         #step3 找到Γtmp中最大v最大的点

            max_value = -Inf
            max_node = nothing
            
            for i in Γtmp
                if v > max_value
                    max_value = v
                    max_node = nodes[i]
                end
            end
          
            for node in Γtmp
               
                distance = norm(nodes[max_node] - nodes[node])
                #step4 
                if distance <= l
                    deleteat!(Γtmp, findall(x -> x == node, Γtmp))
                    push!(Γfinal, nodes[max_node])
                end
            end
          
            Γtmp = setdiff(Γtmp, Γfinal)
        end
        #stpe5 将Γfinal中的点进行排序
        sorted_nodes = sort(collect(Γfinal), by = node -> nodes[node][1])  
        crack_path = Vector{Tuple{Float64, Float64}}()
    
        normals = []
        slips = []
        #计算每一个线段的法向量和切向量
        for i in 1:length(sorted_nodes)-1
            node1 = nodes[sorted_nodes[i]]
            node2 = nodes[sorted_nodes[i + 1]]
            segment = (node1, node2)
    
            s = normalize!(node2 - node1)
    

            n = (-s[2], s[1])
            push!(normals, n)

            push!(crack_path, segment)

        end
        for node in Γdam
            min_distance = Inf
            for i in 1:length(crack_path) - 1
                segment = crack_path[i]
                distance = pointdistfromseg(segment, node)
                if distance < min_distance
                    min_distance = distance
                    nearest_segment = segment
                end
            end
            node.n = nearest_segment.n
            node.s = nearest_segment.s
            n₁ = n[1]
            n₂ = n[2]
            s₁ = n₂
            s₂ = -n₁

        end  
        if v< 0.05
            n₁ = 0.0
            n₂ = 0.0
            s₁ = 0.0
            s₂ = 0.0
            
        end

        # predict phase
        εₙ = ε₁₁*n₁*n₁ + ε₂₂*n₂*n₂ + ε₁₂*n₁*n₂

        σ₁₁ᵗʳ = (λ+2μ)*ε₁₁ + λ*ε₂₂
        σ₂₂ᵗʳ = λ*ε₁₁ + (λ+2μ)*ε₂₂
        σ₁₂ᵗʳ = μ*ε₁₂

        τ = σ₁₁ᵗʳ*n₁*s₁ + σ₂₂ᵗʳ*n₂*s₂ + σ₁₂ᵗʳ*(n₁*s₂+n₂*s₁)
        σ = σ₁₁ᵗʳ*n₁*n₁ + σ₂₂ᵗʳ*n₂*n₂ + 2.0*σ₁₂ᵗʳ*n₁*n₂
        𝑝 = -σ

        𝑓  = abs(τ) - μ̄ *𝑝
        if abs(1.0-v) < 2e-2
             ξ.σ₁₁ = σ₁₁ᵗʳ
             ξ.σ₂₂ = σ₂₂ᵗʳ
             ξ.σ₁₂ = σ₁₂ᵗʳ
             Cᵗ₁₁₁₁ = λ+2μ
             Cᵗ₂₂₂₂ = λ+2μ
             Cᵗ₁₁₂₂ = λ 
             Cᵗ₂₂₁₁ = λ
             Cᵗ₁₂₁₂ = μ 
             Cᵗ₁₁₁₂ = 0.0
             Cᵗ₂₂₁₂ = 0.0
             Cᵗ₁₂₁₁ = 0.0
             Cᵗ₁₂₂₂ = 0.0
        else
            if εₙ > 0.0
                ξ.σ₁₁ = (v^2+η)*σ₁₁ᵗʳ
                ξ.σ₂₂ = (v^2+η)*σ₂₂ᵗʳ
                ξ.σ₁₂ = (v^2+η)*σ₁₂ᵗʳ

                Cᵗ₁₁₁₁ = (v^2+η)*(λ+2μ)
                Cᵗ₂₂₂₂ = (v^2+η)*(λ+2μ)
                Cᵗ₁₁₂₂ = (v^2+η)*λ 
                Cᵗ₂₂₁₁ = (v^2+η)*λ
                Cᵗ₁₂₁₂ = (v^2+η)*μ 
                Cᵗ₁₁₁₂ = 0.0
                Cᵗ₂₂₁₂ = 0.0
                Cᵗ₁₂₁₁ = 0.0
                Cᵗ₁₂₂₂ = 0.0
            else
                if 𝑓 < 0.0
                    ξ.σ₁₁ = σ₁₁ᵗʳ
                    ξ.σ₂₂ = σ₂₂ᵗʳ
                    ξ.σ₁₂ = σ₁₂ᵗʳ
    
                    Cᵗ₁₁₁₁ = λ+2μ
                    Cᵗ₂₂₂₂ = λ+2μ
                    Cᵗ₁₁₂₂ = λ 
                    Cᵗ₂₂₁₁ = λ
                    Cᵗ₁₂₁₂ = μ 
                    Cᵗ₁₁₁₂ = 0.0
                    Cᵗ₂₂₁₂ = 0.0
                    Cᵗ₁₂₁₁ = 0.0
                    Cᵗ₁₂₂₂ = 0.0

                else
                    signτ = sign(τ)
                    Cᶠ₁₁₁₁ = -μ̄*signτ*(λ+2.0*μ*n₁*n₁)*2.0*n₁*s₁
                    Cᶠ₂₂₂₂ = -μ̄*signτ*(λ+2.0*μ*n₂*n₂)*2.0*n₂*s₂
                    Cᶠ₁₁₂₂ = -μ̄*signτ*(λ+2.0*μ*n₂*n₂)*2.0*n₁*s₁
                    Cᶠ₂₂₁₁ = -μ̄*signτ*(λ+2.0*μ*n₁*n₁)*2.0*n₂*s₂
                    Cᶠ₁₁₁₂ = -μ̄*signτ*2.0*μ*n₁*n₂*2.0*n₁*s₁
                    Cᶠ₁₂₁₁ = -μ̄*signτ*(λ+2.0*μ*n₁*n₁)*(n₁*s₂+n₂*s₁)
                    Cᶠ₂₂₁₂ = -μ̄*signτ*2.0*μ*n₁*n₂*2.0*n₂*s₂ 
                    Cᶠ₁₂₂₂ = -μ̄*signτ*(λ+2.0*μ*n₂*n₂)*(n₁*s₂+n₂*s₁)
                    Cᶠ₁₂₁₂ = -μ̄*signτ*2.0*μ*n₁*n₂*(n₁*s₂+n₂*s₁)
                
                    Cʳ₁₁₁₁ =  μ*(n₁*s₁+n₁*s₁)*(n₁*s₁+n₁*s₁)
                    Cʳ₂₂₂₂ =  μ*(n₂*s₂+n₂*s₂)*(n₂*s₂+n₂*s₂)
                    Cʳ₁₁₂₂ =  μ*(n₂*s₂+n₂*s₂)*(n₁*s₁+n₁*s₁)
                    Cʳ₂₂₁₁ =  μ*(n₁*s₁+n₁*s₁)*(n₂*s₂+n₂*s₂)
                    Cʳ₁₁₁₂ =  μ*(n₁*s₂+n₂*s₁)*(n₁*s₁+n₁*s₁)
                    Cʳ₁₂₁₁ =  μ*(n₁*s₁+n₁*s₁)*(n₁*s₂+n₂*s₁)
                    Cʳ₂₂₁₂ =  μ*(n₁*s₂+n₂*s₁)*(n₂*s₂+n₂*s₂)
                    Cʳ₁₂₂₂ =  μ*(n₂*s₂+n₂*s₂)*(n₁*s₂+n₂*s₁)
                    Cʳ₁₂₁₂ =  μ*(n₁*s₂+n₂*s₁)*(n₁*s₂+n₂*s₁)

                    ξ.σ₁₁ = σ₁₁ᵗʳ + (1.0-v^2-η)*(μ̄*𝑝*signτ - τ)*2*n₁*s₁
                    ξ.σ₂₂ = σ₂₂ᵗʳ + (1.0-v^2-η)*(μ̄*𝑝*signτ - τ)*2*n₂*s₂
                    ξ.σ₁₂ = σ₁₂ᵗʳ + (1.0-v^2-η)*(μ̄*𝑝*signτ - τ)*(n₁*s₂+n₂*s₁)

                    Cᵗ₁₁₁₁ = λ+2μ + (1.0-v^2-η)*(Cᶠ₁₁₁₁-Cʳ₁₁₁₁)
                    Cᵗ₂₂₂₂ = λ+2μ + (1.0-v^2-η)*(Cᶠ₂₂₂₂-Cʳ₂₂₂₂)
                    Cᵗ₁₁₂₂ = λ + (1.0-v^2-η)*(Cᶠ₁₁₂₂-Cʳ₁₁₂₂)
                    Cᵗ₂₂₁₁ = λ + (1.0-v^2-η)*(Cᶠ₂₂₁₁-Cʳ₂₂₁₁)
                    Cᵗ₁₂₁₂ = μ + (1.0-v^2-η)*(Cᶠ₁₂₁₂-Cʳ₁₂₁₂)
                    Cᵗ₁₁₁₂ = (1.0-v^2-η)*(Cᶠ₁₁₁₂-Cʳ₁₁₁₂)
                    Cᵗ₁₂₁₁ = (1.0-v^2-η)*(Cᶠ₁₂₁₁-Cʳ₁₂₁₁)
                    Cᵗ₂₂₁₂ = (1.0-v^2-η)*(Cᶠ₂₂₁₂-Cʳ₂₂₁₂)
                    Cᵗ₁₂₂₂ = (1.0-v^2-η)*(Cᶠ₁₂₂₂-Cʳ₁₂₂₂)
                end
            end
        end

        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
             J = xⱼ.𝐼
             k[2*I-1,2*J-1] += (Cᵗ₁₁₁₁*B₁[i]*B₁[j] + Cᵗ₁₁₁₂*B₁[i]*B₂[j] + Cᵗ₁₂₁₁*B₂[i]*B₁[j] + Cᵗ₁₂₁₂*B₂[i]*B₂[j])*𝑤
             k[2*I-1,2*J]   += (Cᵗ₁₁₁₂*B₁[i]*B₁[j] + Cᵗ₁₁₂₂*B₁[i]*B₂[j] + Cᵗ₁₂₁₂*B₂[i]*B₁[j] + Cᵗ₁₂₂₂*B₂[i]*B₂[j])*𝑤
             k[2*I,2*J-1]   += (Cᵗ₁₂₁₁*B₁[i]*B₁[j] + Cᵗ₁₂₁₂*B₁[i]*B₂[j] + Cᵗ₂₂₁₁*B₂[i]*B₁[j] + Cᵗ₂₂₁₂*B₂[i]*B₂[j])*𝑤
             k[2*I,2*J]     += (Cᵗ₁₂₁₂*B₁[i]*B₁[j] + Cᵗ₁₂₂₂*B₁[i]*B₂[j] + Cᵗ₂₂₁₂*B₂[i]*B₁[j] + Cᵗ₂₂₂₂*B₂[i]*B₂[j])*𝑤
           end
           fint[2*I-1] += (B₁[i]*ξ.σ₁₁ + B₂[i]*ξ.σ₁₂)*𝑤
           fint[2*I]   += (B₁[i]*ξ.σ₁₂ + B₂[i]*ξ.σ₂₂)*𝑤
       end   
    end
end 
function calculate_n_and_s_vectors(nodes, v, l)
        Γtmp = Set{Int}()
        Γdam = Set{Int}()
        Γfinal = Set{Int}()
        #step1 找到的所有v>0.98的点
        for (i,xᵢ) in enumerate(𝓒)
            if xᵢ.v > 0.98   
                #step2 并将其存入Γtmp中
                push!(Γtmp, i)
            end
            if xᵢ.v > 0.05   
                push!(Γdam, i)
            end
        end
    
        
        while !isempty(Γtmp)
         #step3 找到rtmp中最大v最大的点
            max_v, max_node = findmax(v[collect(Γtmp)])
         
            for node in Γtmp
               
                distance = norm(nodes[max_node] - nodes[node])
                #step4 
                if distance <= l
                    deleteat!(Γtmp, findall(x -> x == node, Γtmp))
                    push!(Γfinal, nodes[max_node])
                end
            end
          
            Γtmp = setdiff(Γtmp, Γfinal)
        end
        #stpe5 将Γfinal中的点进行排序
        sorted_nodes = sort(collect(Γfinal), by = node -> nodes[node][1])  
        crack_path = Vector{Tuple{Float64, Float64}}()
    
        normals = []
        slips = []
        #计算每一个线段的法向量和切向量
        for i in 1:length(sorted_nodes)-1
            node1 = nodes[sorted_nodes[i]]
            node2 = nodes[sorted_nodes[i + 1]]
            segment = (node1, node2)
    
            s = normalize!(node2 - node1)
    

            n = (-s[2], s[1])
            push!(normals, n)

            push!(crack_path, segment)

        end
        for node in Γdam
            min_distance = Inf
            for i in 1:length(crack_path) - 1
                distance = pointdistfromseg(segment, node)
                if distance < min_distance
                    min_distance = distance
                    nearest_segment = segment
                end
            end
            node.n = nearest_segment.n
            node.s = nearest_segment.s

        end  
        return normals, slips
    end
"""
frictional-contact2
"""
function (op::Operator{:∫vᵢσdΩ_frictional_contact_2})(ap::T;k::AbstractMatrix{Float64},fint::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    μ̄ = op.μ̄ 
    η = op.η
    E = op.E
    ν = op.ν
    λ = E*ν/(1.0+ν)/(1.0-2.0*ν)
    μ = 0.5*E/(1.0+ν)
    tol = op.tol 
    Cᵢᵢᵢᵢ = E/(1-ν^2)
    Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
    Cᵢⱼᵢⱼ = E/2/(1+ν)
    
    for ξ in 𝓖
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        σ₁₁ = ξ.σ₁₁
        σ₂₂ = ξ.σ₂₂
        σ₁₂ = ξ.σ₁₂
        𝑤   = ξ.𝑤
        ε₁₁ = 0.0
        ε₂₂ = 0.0
        ε₁₂ = 0.0
        Δε₁₁ = 0.0
        Δε₂₂ = 0.0
        Δε₁₂ = 0.0
        v = 0.0
        ∂v∂x = 0.0
        ∂v∂y = 0.0
      
        for (i,xᵢ) in enumerate(𝓒)
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += B₁[i]*xᵢ.d₂ + B₂[i]*xᵢ.d₁
            Δε₁₁ += B₁[i]*xᵢ.Δd₁
            Δε₂₂ += B₂[i]*xᵢ.Δd₂
            Δε₁₂ += B₁[i]*xᵢ.Δd₂ + B₂[i]*xᵢ.Δd₁
            v += N[i]*xᵢ.v
            ∂v∂x += B₁[i]*xᵢ.v
            ∂v∂y += B₂[i]*xᵢ.v
        end 
        #println(ε₁₁)
        norm∇v = (∂v∂x^2+∂v∂y^2)^0.5#n = ∇v/norm(∇v)
        n₁ = (∂v∂x)/(norm∇v)
        n₂ = (∂v∂y)/(norm∇v)
        s₁ = n₂
        s₂ = -n₁


       # println(666,∂v∂x,666,∂v∂y,666,norm∇v)
        
        # predict phase
        εₙ = ε₁₁*n₁*n₁ + ε₂₂*n₂*n₂ + ε₁₂*n₁*n₂#2.所以这里就算一个就行了

        σ₁₁ᵗʳ = σ₁₁ + Cᵢᵢᵢᵢ*Δε₁₁ + Cᵢᵢⱼⱼ*Δε₂₂
        σ₂₂ᵗʳ = σ₂₂ + Cᵢᵢᵢᵢ*Δε₂₂ + Cᵢᵢⱼⱼ*Δε₁₁
        σ₁₂ᵗʳ = σ₁₂ + Cᵢⱼᵢⱼ*Δε₁₂
        
        τ = σ₁₁ᵗʳ*n₁*s₁ + σ₂₂ᵗʳ*n₂*s₂ + σ₁₂ᵗʳ*n₁*s₂ + σ₁₂ᵗʳ*n₂*s₁
        σ = σ₁₁ᵗʳ*n₁*n₁ + σ₂₂ᵗʳ*n₂*n₂ + σ₁₂ᵗʳ*n₁*n₂

        

        𝑓  = abs(τ) + μ̄  * σ
        #println(sign(τ))
        #println(𝑓 )
       if (1-v) < tol  #v=1  
          ξ.σ₁₁ = σ₁₁ᵗʳ
          ξ.σ₂₂ = σ₂₂ᵗʳ
          ξ.σ₁₂ = σ₁₂ᵗʳ
          Cᵗ₁₁₁₁ = Cᵢᵢᵢᵢ
          Cᵗ₂₂₂₂ = Cᵢᵢᵢᵢ
          Cᵗ₁₁₂₂ = Cᵢᵢⱼⱼ 
          Cᵗ₂₂₁₁ = Cᵢᵢⱼⱼ
          Cᵗ₁₂₁₂ = Cᵢⱼᵢⱼ 
          Cᵗ₁₁₁₂ = 0.0
          Cᵗ₂₂₁₂ = 0.0
          Cᵗ₁₂₁₁ = 0.0
          Cᵗ₁₂₂₂ = 0.0
        else
            if εₙ > 0.0  #open
               σ₁₁ᵗʳ = σ₁₁ + (v^2+η)*(Cᵢᵢᵢᵢ*Δε₁₁ + Cᵢᵢⱼⱼ*Δε₂₂)
               σ₂₂ᵗʳ = σ₂₂ + (v^2+η)*(Cᵢᵢᵢᵢ*Δε₂₂ + Cᵢᵢⱼⱼ*Δε₁₁)
               σ₁₂ᵗʳ = σ₁₂ + (v^2+η)*Cᵢⱼᵢⱼ*Δε₁₂
        
               ξ.σ₁₁ = σ₁₁ᵗʳ
               ξ.σ₂₂ = σ₂₂ᵗʳ
               ξ.σ₁₂ = σ₁₂ᵗʳ

               Cᵗ₁₁₁₁ = (v^2+η)*Cᵢᵢᵢᵢ
               Cᵗ₂₂₂₂ = (v^2+η)*Cᵢᵢᵢᵢ
               Cᵗ₁₁₂₂ = (v^2+η)*Cᵢᵢⱼⱼ 
               Cᵗ₂₂₁₁ = (v^2+η)*Cᵢᵢⱼⱼ
               Cᵗ₁₂₁₂ = (v^2+η)*Cᵢⱼᵢⱼ 
               Cᵗ₁₁₁₂ = 0.0
               Cᵗ₂₂₁₂ = 0.0
               Cᵗ₁₂₁₁ = 0.0
               Cᵗ₁₂₂₂ = 0.0
            else
                if 𝑓 < tol

                   ξ.σ₁₁ = σ₁₁ᵗʳ
                   ξ.σ₂₂ = σ₂₂ᵗʳ
                   ξ.σ₁₂ = σ₁₂ᵗʳ
                   Cᵗ₁₁₁₁ = Cᵢᵢᵢᵢ
                   Cᵗ₂₂₂₂ = Cᵢᵢᵢᵢ
                   Cᵗ₁₁₂₂ = Cᵢᵢⱼⱼ 
                   Cᵗ₂₂₁₁ = Cᵢᵢⱼⱼ
                   Cᵗ₁₂₁₂ = Cᵢⱼᵢⱼ 
                   Cᵗ₁₁₁₂ = 0.0             
                   Cᵗ₂₂₁₂ = 0.0
                   Cᵗ₁₂₁₁ = 0.0
                   Cᵗ₁₂₂₂ = 0.0
                else
                  ∂fC∂f = λ*μ̄ ^2.0 + 2.0*(0.5+μ̄ ^2.0)*μ
                  Δγ = fᵗʳ/∂fC∂f
                  #σᵢⱼ⁽ⁿ⁺¹⁾ = σᵢⱼᵗʳ - Δγ*Cᵢⱼₖₗ*∂fᵗʳ/∂σₖₗ  
                  #σᵢⱼ      = σᵢⱼᵗʳ - Δγ*(sign(τ)*μ*(nᵢ*sⱼ +nⱼ*sᵢ) + μ̄ *(2*μ*nᵢ*nⱼ))
                  ξ.σ₁₁ = σ₁₁ᵗʳ - (1-v^2+η)*Δγ*(sign(τ)*μ*2.0*n₁*s₁ + μ̄ *(λ+2*μ*n₁*n₁)) 
                  ξ.σ₂₂ = σ₂₂ᵗʳ - (1-v^2+η)*Δγ*(sign(τ)*μ*2.0*n₂*s₂ + μ̄ *(λ+2*μ*n₂*n₂)) 
                  ξ.σ₁₂ = σ₁₂ᵗʳ - (1-v^2+η)*Δγ*(sign(τ)*μ*(n₁*s₂+n₂*s₁) + μ̄ *(2*μ*n₁*n₂)) 
     
                  #Cᵢⱼₖₗᵗ = Cᵢⱼₖₗ - (4.0*μ^2.0*μ̄ ^2.0*nᵢ*nⱼ*nₖ*nₗ + λ^2.0*δᵢⱼ*δₖₗ*μ̄ ^2.0 + sign(τ)^2.0*μ^2.0*(nᵢ*sⱼ+nⱼ*sᵢ)*(nₖ*sₗ+nₗ*sₖ) + λ*δₖₗ*μ̄ *μ*sign(τ)*(nᵢ*sⱼ+nⱼ*sᵢ)+2.0*μ*λ*μ̄ ^2.0*δᵢⱼ*nᵢ*nⱼ + μ*λ*sign(τ)*μ̄ *δᵢⱼ*(nₖ*sₗ+nₗ*sₖ) + 2.0*μ^2.0*sign(τ)*μ̄ *(nₖ*sₗ+nₗ*sₖ)*nᵢ*nⱼ + 2.0*μ*μ̄ ^2.0*λ*δᵢⱼ*nₖ*nₗ + 2.0*μ^2.0*sign(τ)*μ̄ *nₖ*nₗ*(nᵢ*sⱼ+nⱼ*sᵢ))   

                  Cᵗ₁₁₁₁ = Cᵢᵢᵢᵢ - (1-v^2+η)*(4.0*μ^2*μ̄ ^2*n₁^4 + λ^2*μ̄ ^2 + 4.0*sign(τ)^2.0*μ^2.0*(n₁*s₁)^2 + 4.0*λ*μ̄ *μ*sign(τ)*n₁*s₁ + 4.0*μ*λ*μ̄ ^2*n₁^2 + 8.0*μ^2*sign(τ)*μ̄ *n₁*s₁*n₁^2)/∂fC∂f
                  Cᵗ₂₂₂₂ = Cᵢᵢᵢᵢ - (1-v^2+η)*(4.0*μ^2*μ̄ ^2*n₂^4 + λ^2*μ̄ ^2 + 4.0*sign(τ)^2.0*μ^2.0*(n₂*s₂)^2 + 4.0*λ*μ̄ *μ*sign(τ)*n₂*s₂ + 4.0*μ*λ*μ̄ ^2*n₂^2 + 8.0*μ^2*sign(τ)*μ̄ *n₂*s₂*n₂^2)/∂fC∂f
                  Cᵗ₁₁₂₂ = Cᵢᵢⱼⱼ - (1-v^2+η)*(4.0*μ^2.0*μ̄ ^2.0*n₁^2*n₂^2 + λ^2*μ̄ ^2 + 4.0*sign(τ)^2.0*μ^2.0*n₁*s₁*n₂*s₂ + 2.0*λ*μ̄ *μ*sign(τ)*n₁*s₁ +  4.0*μ*λ*μ̄ ^2*n₁^2 + 4.0*μ^2.0*sign(τ)*μ̄ *(n₂*s₂)*n₁^2 +2.0*λ*μ̄ *μ*sign(τ)*n₂*s₂ + 4.0*μ*λ*μ̄ ^2*n₂^2 + 4.0*μ^2.0*sign(τ)*μ̄ *(n₁*s₁)*n₂^2)/∂fC∂f
                  Cᵗ₂₂₁₁ = Cᵢᵢⱼⱼ - (1-v^2+η)*(4.0*μ^2.0*μ̄ ^2.0*n₁^2*n₂^2 + λ^2*μ̄ ^2 + 4.0*sign(τ)^2.0*μ^2.0*n₁*s₁*n₂*s₂ + 2.0*λ*μ̄ *μ*sign(τ)*n₁*s₁ +  4.0*μ*λ*μ̄ ^2*n₁^2 + 4.0*μ^2.0*sign(τ)*μ̄ *(n₂*s₂)*n₁^2 +2.0*λ*μ̄ *μ*sign(τ)*n₂*s₂ + 4.0*μ*λ*μ̄ ^2*n₂^2 + 4.0*μ^2.0*sign(τ)*μ̄ *(n₁*s₁)*n₂^2)/∂fC∂f
                  Cᵗ₁₂₁₂ = Cᵢⱼᵢⱼ - (1-v^2+η)*(4.0*μ^2.0*μ̄ ^2.0*n₁^2*n₂^2 + sign(τ)^2.0*μ^2.0*(n₁*s₂+n₂*s₁)^2 + 4.0*μ^2.0*sign(τ)*μ̄ *(n₁*s₂+n₂*s₁)*n₁*n₂)/∂fC∂f
                  Cᵗ₁₁₁₂ = 0.0   - (1-v^2+η)*(4.0*μ^2.0*μ̄ ^2.0*n₂*n₁^3 + 2.0*sign(τ)^2*μ^2*(n₁*s₁)*(n₁*s₂+n₂*s₁) + 2.0*λ*μ̄ *μ*sign(τ)*(n₁*s₂+n₂*s₁) + 4.0*μ*λ*μ̄ ^2*n₁*n₂+ 4.0*μ^2.0*sign(τ)*μ̄ *n₁*s₁*n₁*n₂+ 4.0*μ^2.0*sign(τ)*μ̄ *n₁*s₂*n₁*n₁)/∂fC∂f
                  Cᵗ₁₂₁₁ = 0.0   - (1-v^2+η)*(4.0*μ^2.0*μ̄ ^2.0*n₂*n₁^3 + 2.0*sign(τ)^2*μ^2*(n₁*s₁)*(n₁*s₂+n₂*s₁) + 2.0*λ*μ̄ *μ*sign(τ)*(n₁*s₂+n₂*s₁) + 4.0*μ*λ*μ̄ ^2*n₁*n₂+ 4.0*μ^2.0*sign(τ)*μ̄ *n₁*s₁*n₁*n₂+ 4.0*μ^2.0*sign(τ)*μ̄ *n₁*s₂*n₁*n₁)/∂fC∂f
                  Cᵗ₂₂₁₂ = 0.0   - (1-v^2+η)*(4.0*μ^2.0*μ̄ ^2.0*n₁*n₂^3 + 2.0*sign(τ)^2*μ^2*(n₂*s₂)*(n₂*s₁+n₁*s₂) + 2.0*λ*μ̄ *μ*sign(τ)*(n₂*s₁+n₁*s₂) + 4.0*μ*λ*μ̄ ^2*n₂*n₁+ 4.0*μ^2.0*sign(τ)*μ̄ *n₂*s₂*n₂*n₁+ 4.0*μ^2.0*sign(τ)*μ̄ *n₂*s₁*n₂*n₂)/∂fC∂f
                  Cᵗ₁₂₂₂ = 0.0   - (1-v^2+η)*(4.0*μ^2.0*μ̄ ^2.0*n₁*n₂^3 + 2.0*sign(τ)^2*μ^2*(n₂*s₂)*(n₂*s₁+n₁*s₂) + 2.0*λ*μ̄ *μ*sign(τ)*(n₂*s₁+n₁*s₂) + 4.0*μ*λ*μ̄ ^2*n₂*n₁+ 4.0*μ^2.0*sign(τ)*μ̄ *n₂*s₂*n₂*n₁+ 4.0*μ^2.0*sign(τ)*μ̄ *n₂*s₁*n₂*n₂)/∂fC∂f
               end  
            end    
        end   

       for (i,xᵢ) in enumerate(𝓒)
           I = xᵢ.𝐼
           for (j,xⱼ) in enumerate(𝓒)
             J = xⱼ.𝐼
             k[2*I-1,2*J-1] += (Cᵗ₁₁₁₁*B₁[i]*B₁[j] + Cᵗ₁₂₁₂*B₂[i]*B₂[j] + Cᵗ₁₁₁₂*B₁[i]*B₂[j] + Cᵗ₁₂₁₁*B₂[i]*B₁[j])*𝑤
             k[2*I-1,2*J]   += (Cᵗ₁₁₂₂*B₁[i]*B₂[j] + Cᵗ₁₂₁₂*B₂[i]*B₁[j] + Cᵗ₁₁₁₂*B₁[i]*B₁[j] + Cᵗ₁₂₂₂*B₂[i]*B₂[j])*𝑤
             k[2*I,2*J-1]   += (Cᵗ₂₂₁₁*B₂[i]*B₁[j] + Cᵗ₁₂₁₂*B₁[i]*B₂[j] + Cᵗ₁₂₁₁*B₁[i]*B₁[j] + Cᵗ₂₂₁₂*B₂[i]*B₂[j])*𝑤
             k[2*I,2*J]     += (Cᵗ₂₂₂₂*B₂[i]*B₂[j] + Cᵗ₁₂₁₂*B₁[i]*B₁[j] + Cᵗ₁₂₂₂*B₁[i]*B₂[j] + Cᵗ₂₂₁₂*B₂[i]*B₁[j])*𝑤
           end
           fint[2*I-1] += (B₁[i]*ξ.σ₁₁ + B₂[i]*ξ.σ₁₂)*𝑤
           fint[2*I]   += (B₁[i]*ξ.σ₁₂ + B₂[i]*ξ.σ₂₂)*𝑤
       end   
    end
end 

"""
mc_phasefield
"""
function (op::Operator{:∫vᵢσdΩ_mc_phasefield})(ap::T;k::AbstractMatrix{Float64},fint::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    η = op.η
    tol = op.tol
    E = op.E
    ν = op.ν
    λ = E*ν/(1.0+ν)/(1.0-2.0*ν)
    μ = 0.5*E/(1.0+ν)
    c = 10
    𝜙 = π/3.0
    tol = op.tol 
    η = op.η
    Cᵢᵢᵢᵢ = E/(1-ν^2)
    Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
    Cᵢⱼᵢⱼ = E/2/(1+ν)
    
    for  ξ in 𝓖
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        σ₁₁ = ξ.σ₁₁
        σ₂₂ = ξ.σ₂₂
        σ₁₂ = ξ.σ₁₂
        εᵖ₁₁ = ξ.εᵖ₁₁
        εᵖ₂₂ = ξ.εᵖ₂₂
        εᵖ₁₂ = ξ.εᵖ₁₂
        𝑤 = ξ.𝑤
        Δε₁₁ = 0.0
        Δε₂₂ = 0.0
        Δε₁₂ = 0.0
        ε₁₁ = 0.0
        ε₂₂ = 0.0
        ε₁₂ = 0.0
        ∂v∂x = 0.0
        ∂v∂y = 0.0
        v    = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += B₁[i]*xᵢ.d₂ + B₂[i]*xᵢ.d₁
            Δε₁₁ += B₁[i]*xᵢ.Δd₁
            Δε₂₂ += B₂[i]*xᵢ.Δd₂
            Δε₁₂ += B₁[i]*xᵢ.Δd₂ + B₂[i]*xᵢ.Δd₁
            v += N[i]*xᵢ.v
            ∂v∂x += B₁[i]*xᵢ.v
            ∂v∂y += B₂[i]*xᵢ.v
        end 
        norm∇v = (∂v∂x^2+∂v∂y^2)^0.5#n = ∇v/norm(∇v)
        n₁ = (∂v∂x)/(norm∇v)
        n₂ = (∂v∂y)/(norm∇v)


        # predict phase
        εₙ = ε₁₁*n₁*n₁ + ε₂₂*n₂*n₂ + ε₁₂*n₁*n₂

        σ₁₁ᵗʳ = σ₁₁ + Cᵢᵢᵢᵢ*Δε₁₁ + Cᵢᵢⱼⱼ*Δε₂₂
        σ₂₂ᵗʳ = σ₂₂ + Cᵢᵢᵢᵢ*Δε₂₂ + Cᵢᵢⱼⱼ*Δε₁₁
        σ₁₂ᵗʳ = σ₁₂ + Cᵢⱼᵢⱼ*Δε₁₂

        σ₁,σ₂,n₁,n₂ = getσₙ(σ₁₁ᵗʳ,σ₂₂ᵗʳ,σ₁₂ᵗʳ)
        fᵗʳ = σ₁-σ₂ + (σ₁+σ₂) * sin(𝜙) - 2.0*c*cos(𝜙)


 
        #if fᵗʳ > tol
        #    
        #  error("fᵗʳ > tol!")
        #end
        if (1-v) < tol 
            ξ.σ₁₁ = σ₁₁ᵗʳ
            ξ.σ₂₂ = σ₂₂ᵗʳ
            ξ.σ₁₂ = σ₁₂ᵗʳ
            Cᵗ₁₁₁₁ = Cᵢᵢᵢᵢ 
            Cᵗ₂₂₂₂ = Cᵢᵢᵢᵢ 
            Cᵗ₁₁₂₂ = Cᵢᵢⱼⱼ  
            Cᵗ₁₂₁₂ = Cᵢⱼᵢⱼ 
            Cᵗ₁₁₁₂ = 0.0
            Cᵗ₂₂₁₂ = 0.0
        else   
            if εₙ > 0.0
               σ₁₁ᵗʳ = σ₁₁ + (v^2+η)*(Cᵢᵢᵢᵢ*Δε₁₁ + Cᵢᵢⱼⱼ*Δε₂₂)
               σ₂₂ᵗʳ = σ₂₂ + (v^2+η)*(Cᵢᵢᵢᵢ*Δε₂₂ + Cᵢᵢⱼⱼ*Δε₁₁)
               σ₁₂ᵗʳ = σ₁₂ + (v^2+η)*Cᵢⱼᵢⱼ*Δε₁₂
               ξ.σ₁₁ = σ₁₁ᵗʳ
               ξ.σ₂₂ = σ₂₂ᵗʳ
               ξ.σ₁₂ = σ₁₂ᵗʳ
               Cᵗ₁₁₁₁ = (v^2+η)*Cᵢᵢᵢᵢ
               Cᵗ₂₂₂₂ = (v^2+η)*Cᵢᵢᵢᵢ
               Cᵗ₁₁₂₂ = (v^2+η)*Cᵢᵢⱼⱼ 
               Cᵗ₁₂₁₂ = (v^2+η)*Cᵢⱼᵢⱼ 
               Cᵗ₁₁₁₂ = 0.0
               Cᵗ₂₂₁₂ = 0.0
            else
                if fᵗʳ > tol
                    sin𝜙 = sin(𝜙)
                    sin²𝜙 = sin(𝜙)*sin(𝜙)
                    sin²𝜙psin𝜙 = sin²𝜙 +sin(𝜙) 
                    sin²𝜙dsin𝜙 = sin²𝜙 -sin(𝜙) 
                    ∂fC∂f=(4.0*sin²𝜙*λ+4.0*(1.0+sin²𝜙)*μ)
                    Δγ = fᵗʳ/∂fC∂f
                   
                    ξ.σ₁₁ = σ₁₁ᵗʳ - (1-v^2+η)*Δγ*(2*sin𝜙*λ+2*((sin𝜙+1)*n₁[1]*n₁[1]+(sin𝜙-1)*n₂[1]*n₂[1])*μ)#上标1对应之前下标1
                    ξ.σ₂₂ = σ₂₂ᵗʳ - (1-v^2+η)*Δγ*(2*sin𝜙*λ+2*((sin𝜙+1)*n₁[2]*n₁[2]+(sin𝜙-1)*n₂[2]*n₂[2])*μ)
                    ξ.σ₁₂ = σ₁₂ᵗʳ - (1-v^2+η)*Δγ*2*((sin𝜙+1)*n₁[1]*n₁[2]+(sin𝜙-1)*n₂[1]*n₂[2])*μ
                    ξ.εᵖ₁₁ = εᵖ₁₁ + (1-v^2+η)*Δγ*((sin𝜙+1)*n₁[1]*n₁[1]+(sin𝜙-1)n₂[1]*n₂[1])
                    ξ.εᵖ₂₂ = εᵖ₂₂ + (1-v^2+η)*Δγ*((sin𝜙+1)*n₁[2]*n₁[2]+(sin𝜙-1)n₂[2]*n₂[2])
                    ξ.εᵖ₁₂ = εᵖ₁₂ + (1-v^2+η)*Δγ*(((sin𝜙+1)*n₁[1]*n₁[2]+(sin𝜙-1)n₂[1]*n₂[2]))
    
                    Cᵗ₁₁₁₁=Cᵢᵢᵢᵢ-(1-v^2+η)*(4.0*sin²𝜙*λ^2 + 8.0*sin²𝜙psin𝜙*λ*μ*(n₁[1])^2 + 8.0*(sin²𝜙dsin𝜙)*λ*μ*(n₂[1])^2-8.0*cos(𝜙)*cos(𝜙)*μ^2*n₁[1]^2*n₂[1]^2+4.0*(sin𝜙+1)^2*μ^2*(n₁[1])^4+4.0*(sin𝜙-1)^2*μ^2*(n₂[1])^4)/∂fC∂f
                    Cᵗ₂₂₂₂=Cᵢᵢᵢᵢ-(1-v^2+η)*(4.0*sin²𝜙*λ^2 + 8.0*sin²𝜙psin𝜙*λ*μ*(n₁[2])^2 + 8.0*(sin²𝜙dsin𝜙)*λ*μ*(n₂[2])^2-8.0*cos(𝜙)*cos(𝜙)*μ^2*n₁[2]^2*n₂[2]^2+4.0*(sin𝜙+1)^2*μ^2*(n₁[2])^4+4.0*(sin𝜙-1)^2*μ^2*(n₂[2])^4)/∂fC∂f
                    Cᵗ₁₁₂₂=Cᵢᵢⱼⱼ-(1-v^2+η)*(4.0*sin²𝜙*λ^2 + 4.0*sin²𝜙psin𝜙*λ*μ*((n₁[1])^2+(n₁[2])^2)+4.0*(sin²𝜙dsin𝜙)*λ*μ*((n₂[1])^2+(n₂[2])^2)+4.0*(sin𝜙+1)^2*μ^2*(n₁[1])^2*(n₁[2])^2-4.0*cos(𝜙)*cos(𝜙)*μ^2*((n₂[1])^2*(n₁[2])^2+(n₁[1])^2*(n₂[2])^2)+4.0*(sin𝜙-1)^2*μ^2*(n₂[1])^2*(n₂[2])^2)/∂fC∂f
                    Cᵗ₁₂₁₂=Cᵢⱼᵢⱼ-2.0*(1-v^2+η)*(-8.0*cos(𝜙)*cos(𝜙)*μ^2*n₁[1]*n₁[2]*n₂[1]*n₂[2]+4.0*(sin𝜙+1)^2*μ^2*(n₁[1]*n₁[2])^2+4.0*(sin𝜙-1)^2*μ^2*(n₂[1]*n₂[2])^2)/∂fC∂f
                    Cᵗ₁₁₁₂ = - (1-v^2+η)*(4.0*sin²𝜙psin𝜙*λ*μ*n₁[1]*n₁[2] + 4.0*(sin²𝜙dsin𝜙)*λ*μ*n₂[1]*n₂[2]-4.0*cos(𝜙)*cos(𝜙)*μ^2*n₂[1]^2*n₁[1]*n₁[2]-4.0*cos(𝜙)*cos(𝜙)*μ^2*n₁[1]^2*n₂[1]*n₂[2]+4.0*(sin𝜙+1)^2*μ^2*(n₁[1])^3*n₁[2]+4.0*(sin𝜙-1)^2*μ^2*(n₂[1])^3*n₂[2])/∂fC∂f
                    Cᵗ₂₂₁₂ = - (1-v^2+η)*(4.0*sin²𝜙psin𝜙*λ*μ*n₁[2]*n₁[1] + 4.0*(sin²𝜙dsin𝜙)*λ*μ*n₂[2]*n₂[1]-4.0*cos(𝜙)*cos(𝜙)*μ^2*n₂[2]^2*n₁[2]*n₁[1]-4.0*cos(𝜙)*cos(𝜙)*μ^2*n₁[2]^2*n₂[2]*n₂[1]+4.0*(sin𝜙+1)^2*μ^2*(n₁[2])^3*n₁[1]+4.0*(sin𝜙-1)^2*μ^2*(n₂[2])^3*n₂[1])/∂fC∂f
                else
                    ξ.σ₁₁ = σ₁₁ᵗʳ
                    ξ.σ₂₂ = σ₂₂ᵗʳ
                    ξ.σ₁₂ = σ₁₂ᵗʳ
                    Cᵗ₁₁₁₁ = Cᵢᵢᵢᵢ
                    Cᵗ₂₂₂₂ = Cᵢᵢᵢᵢ
                    Cᵗ₁₁₂₂ = Cᵢᵢⱼⱼ
                    Cᵗ₁₂₁₂ = Cᵢⱼᵢⱼ
                    Cᵗ₁₁₁₂ = 0.0
                    Cᵗ₂₂₁₂ = 0.0
                end             
            end
        end    
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += (Cᵗ₁₁₁₁*B₁[i]*B₁[j] + Cᵗ₁₂₁₂*B₂[i]*B₂[j] + Cᵗ₁₁₁₂*(B₁[i]*B₂[j]+B₂[i]*B₁[j]))*𝑤
                k[2*I-1,2*J]   += (Cᵗ₁₁₂₂*B₁[i]*B₂[j] + Cᵗ₁₂₁₂*B₂[i]*B₁[j] + Cᵗ₁₁₁₂*B₁[i]*B₁[j] + Cᵗ₂₂₁₂*B₂[i]*B₂[j])*𝑤
                k[2*I,2*J-1]   += (Cᵗ₁₁₂₂*B₂[i]*B₁[j] + Cᵗ₁₂₁₂*B₁[i]*B₂[j] + Cᵗ₁₁₁₂*B₁[i]*B₁[j] + Cᵗ₂₂₁₂*B₂[i]*B₂[j])*𝑤
                k[2*I,2*J]     += (Cᵗ₂₂₂₂*B₂[i]*B₂[j] + Cᵗ₁₂₁₂*B₁[i]*B₁[j] + Cᵗ₂₂₁₂*(B₁[i]*B₂[j]+B₂[i]*B₁[j]))*𝑤
     
            end
            fint[2*I-1] += B₁[i]*ξ.σ₁₁*𝑤 + B₂[i]*ξ.σ₁₂*𝑤
            fint[2*I]   += B₁[i]*ξ.σ₁₂*𝑤 + B₂[i]*ξ.σ₂₂*𝑤
    
        end
    end
end    
