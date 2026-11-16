module KirchhoffShell

using ..ApproxOperator: AbstractElement

"""
    ∫εᵢⱼNᵢⱼκᵢⱼMᵢⱼdΩ(ap::T;k::AbstractMatrix{Float64}, E, ν, h)

Kirchhoff shell stiffness, membrane plus bending,
`∫ (εᵢⱼNᵢⱼ + κᵢⱼMᵢⱼ) dΩ`.

The shell metric `a¹¹, a¹², a²²` and the Christoffel symbols come from the
surface parameterisation and are read from the quadrature point, so this works
on curved surfaces without a flat projection.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set∇²𝝭!` — reads `ξ[:𝝭]`, the first and second derivatives, the metric data
`ξ.a¹¹`, `ξ.a¹²`, `ξ.a²²`, `ξ.Γ¹₁₁`, ... and `ξ.𝑤`.
"""
function ∫εᵢⱼNᵢⱼκᵢⱼMᵢⱼdΩ(ap::T, k::AbstractMatrix{Float64}; E::Real, ν::Real, h::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        Dᵐ = E*h/(1-ν^2)
        Dᵇ = E*h^3/12/(1-ν^2)
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₁₁ = ξ[:∂²𝝭∂x²]
        B₁₂ = ξ[:∂²𝝭∂x∂y]
        B₂₂ = ξ[:∂²𝝭∂y²]
        𝑤 = ξ.𝑤
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
        𝐃 = [                  a¹¹*a¹¹ ν*a¹¹*a²² + (1-ν)*a¹²*a¹²                            a¹¹*a¹²;
        ν*a¹¹*a²² + (1-ν)*a¹²*a¹²                   a²²*a²²                            a²²*a¹²;
        a¹¹*a¹²                   a²²*a¹² 0.5*((1-ν)*a¹¹*a²² + (1+ν)*a¹²*a¹²)]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            𝐁ᵐᵢ = [          𝒂₁₍₁₎*B₁[i]           𝒂₁₍₂₎*B₁[i]           𝒂₁₍₃₎*B₁[i];
            𝒂₂₍₁₎*B₂[i]           𝒂₂₍₂₎*B₂[i]           𝒂₂₍₃₎*B₂[i];
            𝒂₁₍₁₎*B₂[i]+𝒂₂₍₁₎*B₁[i] 𝒂₁₍₂₎*B₂[i]+𝒂₂₍₂₎*B₁[i] 𝒂₁₍₃₎*B₂[i]+𝒂₂₍₃₎*B₁[i]]
            Bᵇ₁₁ = Γ¹₁₁*B₁[i]+Γ²₁₁*B₂[i]-B₁₁[i]
            Bᵇ₂₂ = Γ¹₂₂*B₁[i]+Γ²₂₂*B₂[i]-B₂₂[i]
            Bᵇ₁₂ = Γ¹₁₂*B₁[i]+Γ²₁₂*B₂[i]-B₁₂[i]
            𝐁ᵇᵢ = [  Bᵇ₁₁*𝒂₃₍₁₎   Bᵇ₁₁*𝒂₃₍₂₎   Bᵇ₁₁*𝒂₃₍₃₎;
            Bᵇ₂₂*𝒂₃₍₁₎   Bᵇ₂₂*𝒂₃₍₂₎   Bᵇ₂₂*𝒂₃₍₃₎;
            2*Bᵇ₁₂*𝒂₃₍₁₎ 2*Bᵇ₁₂*𝒂₃₍₂₎ 2*Bᵇ₁₂*𝒂₃₍₃₎]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                𝐁ᵐⱼ = [          𝒂₁₍₁₎*B₁[j]           𝒂₁₍₂₎*B₁[j]           𝒂₁₍₃₎*B₁[j];
                𝒂₂₍₁₎*B₂[j]           𝒂₂₍₂₎*B₂[j]           𝒂₂₍₃₎*B₂[j];
                𝒂₁₍₁₎*B₂[j]+𝒂₂₍₁₎*B₁[j] 𝒂₁₍₂₎*B₂[j]+𝒂₂₍₂₎*B₁[j] 𝒂₁₍₃₎*B₂[j]+𝒂₂₍₃₎*B₁[j]]
                Bᵇ₁₁ = Γ¹₁₁*B₁[j]+Γ²₁₁*B₂[j]-B₁₁[j]
                Bᵇ₂₂ = Γ¹₂₂*B₁[j]+Γ²₂₂*B₂[j]-B₂₂[j]
                Bᵇ₁₂ = Γ¹₁₂*B₁[j]+Γ²₁₂*B₂[j]-B₁₂[j]
                𝐁ᵇⱼ = [  Bᵇ₁₁*𝒂₃₍₁₎   Bᵇ₁₁*𝒂₃₍₂₎   Bᵇ₁₁*𝒂₃₍₃₎;
                Bᵇ₂₂*𝒂₃₍₁₎   Bᵇ₂₂*𝒂₃₍₂₎   Bᵇ₂₂*𝒂₃₍₃₎;
                2*Bᵇ₁₂*𝒂₃₍₁₎ 2*Bᵇ₁₂*𝒂₃₍₂₎ 2*Bᵇ₁₂*𝒂₃₍₃₎]

                kᵐ = Dᵐ*𝐁ᵐᵢ'*𝐃*𝐁ᵐⱼ
                kᵇ = Dᵇ*𝐁ᵇᵢ'*𝐃*𝐁ᵇⱼ
                k[3*I-2,3*J-2] += (kᵐ[1,1] + kᵇ[1,1])*𝑤
                k[3*I-2,3*J-1] += (kᵐ[1,2] + kᵇ[1,2])*𝑤
                k[3*I-2,3*J]   += (kᵐ[1,3] + kᵇ[1,3])*𝑤
                k[3*I-1,3*J-2] += (kᵐ[2,1] + kᵇ[2,1])*𝑤
                k[3*I-1,3*J-1] += (kᵐ[2,2] + kᵇ[2,2])*𝑤
                k[3*I-1,3*J]   += (kᵐ[2,3] + kᵇ[2,3])*𝑤
                k[3*I,3*J-2]   += (kᵐ[3,1] + kᵇ[3,1])*𝑤
                k[3*I,3*J-1]   += (kᵐ[3,2] + kᵇ[3,2])*𝑤
                k[3*I,3*J]     += (kᵐ[3,3] + kᵇ[3,3])*𝑤
            end
        end
    end
end

"""
    ∫εᵢⱼNᵢⱼdΩ(ap::T;k::AbstractMatrix{Float64}, E, ν, h)

Membrane (stretching) part of the Kirchhoff shell,
`∫ εᵢⱼ(u) Nᵢⱼ(v) dΩ` with `N = Dᵐ C ε`, `Dᵐ = Eh/(1-ν²)`.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set∇𝝭!` — reads `ξ[:𝝭]`, the first derivatives and the surface metric.
"""
function ∫εᵢⱼNᵢⱼdΩ(ap::T, k::AbstractMatrix{Float64}; E::Real, ν::Real, h::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        Dᵐ = E*h/(1-ν^2)
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
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
        𝐃 = [                  a¹¹*a¹¹ ν*a¹¹*a²² + (1-ν)*a¹²*a¹²                            a¹¹*a¹²;
        ν*a¹¹*a²² + (1-ν)*a¹²*a¹²                   a²²*a²²                            a²²*a¹²;
        a¹¹*a¹²                   a²²*a¹² 0.5*((1-ν)*a¹¹*a²² + (1+ν)*a¹²*a¹²)]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            𝐁ᵐᵢ = [          𝒂₁₍₁₎*B₁[i]           𝒂₁₍₂₎*B₁[i]           𝒂₁₍₃₎*B₁[i];
            𝒂₂₍₁₎*B₂[i]           𝒂₂₍₂₎*B₂[i]           𝒂₂₍₃₎*B₂[i];
            𝒂₁₍₁₎*B₂[i]+𝒂₂₍₁₎*B₁[i] 𝒂₁₍₂₎*B₂[i]+𝒂₂₍₂₎*B₁[i] 𝒂₁₍₃₎*B₂[i]+𝒂₂₍₃₎*B₁[i]]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                𝐁ᵐⱼ = [          𝒂₁₍₁₎*B₁[j]           𝒂₁₍₂₎*B₁[j]           𝒂₁₍₃₎*B₁[j];
                𝒂₂₍₁₎*B₂[j]           𝒂₂₍₂₎*B₂[j]           𝒂₂₍₃₎*B₂[j];
                𝒂₁₍₁₎*B₂[j]+𝒂₂₍₁₎*B₁[j] 𝒂₁₍₂₎*B₂[j]+𝒂₂₍₂₎*B₁[j] 𝒂₁₍₃₎*B₂[j]+𝒂₂₍₃₎*B₁[j]]

                kᵐ = Dᵐ*𝐁ᵐᵢ'*𝐃*𝐁ᵐⱼ
                k[3*I-2,3*J-2] += kᵐ[1,1]*𝑤
                k[3*I-2,3*J-1] += kᵐ[1,2]*𝑤
                k[3*I-2,3*J]   += kᵐ[1,3]*𝑤
                k[3*I-1,3*J-2] += kᵐ[2,1]*𝑤
                k[3*I-1,3*J-1] += kᵐ[2,2]*𝑤
                k[3*I-1,3*J]   += kᵐ[2,3]*𝑤
                k[3*I,3*J-2]   += kᵐ[3,1]*𝑤
                k[3*I,3*J-1]   += kᵐ[3,2]*𝑤
                k[3*I,3*J]     += kᵐ[3,3]*𝑤
            end
        end
    end
end

"""
    ∫κᵢⱼMᵢⱼdΩ(ap::T;k::AbstractMatrix{Float64}, E, ν, h)

Bending part of the Kirchhoff shell,
`∫ κᵢⱼ(θ) Mᵢⱼ(θ̄) dΩ` with `M = Dᵇ C κ`, `Dᵇ = Eh³/(12(1-ν²))`.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set∇²𝝭!` — reads the second derivatives and the surface metric.
"""
function ∫κᵢⱼMᵢⱼdΩ(ap::T, k::AbstractMatrix{Float64}; E::Real, ν::Real, h::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        Dᵇ = E*h^3/12/(1-ν^2)
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₁₁ = ξ[:∂²𝝭∂x²]
        B₁₂ = ξ[:∂²𝝭∂x∂y]
        B₂₂ = ξ[:∂²𝝭∂y²]
        𝑤 = ξ.𝑤
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
        𝐃 = [                  a¹¹*a¹¹ ν*a¹¹*a²² + (1-ν)*a¹²*a¹²                            a¹¹*a¹²;
        ν*a¹¹*a²² + (1-ν)*a¹²*a¹²                   a²²*a²²                            a²²*a¹²;
        a¹¹*a¹²                   a²²*a¹² 0.5*((1-ν)*a¹¹*a²² + (1+ν)*a¹²*a¹²)]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            Bᵇ₁₁ = Γ¹₁₁*B₁[i]+Γ²₁₁*B₂[i]-B₁₁[i]
            Bᵇ₂₂ = Γ¹₂₂*B₁[i]+Γ²₂₂*B₂[i]-B₂₂[i]
            Bᵇ₁₂ = Γ¹₁₂*B₁[i]+Γ²₁₂*B₂[i]-B₁₂[i]
            𝐁ᵇᵢ = [  Bᵇ₁₁*𝒂₃₍₁₎   Bᵇ₁₁*𝒂₃₍₂₎   Bᵇ₁₁*𝒂₃₍₃₎;
            Bᵇ₂₂*𝒂₃₍₁₎   Bᵇ₂₂*𝒂₃₍₂₎   Bᵇ₂₂*𝒂₃₍₃₎;
            2*Bᵇ₁₂*𝒂₃₍₁₎ 2*Bᵇ₁₂*𝒂₃₍₂₎ 2*Bᵇ₁₂*𝒂₃₍₃₎]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                Bᵇ₁₁ = Γ¹₁₁*B₁[j]+Γ²₁₁*B₂[j]-B₁₁[j]
                Bᵇ₂₂ = Γ¹₂₂*B₁[j]+Γ²₂₂*B₂[j]-B₂₂[j]
                Bᵇ₁₂ = Γ¹₁₂*B₁[j]+Γ²₁₂*B₂[j]-B₁₂[j]
                𝐁ᵇⱼ = [  Bᵇ₁₁*𝒂₃₍₁₎   Bᵇ₁₁*𝒂₃₍₂₎   Bᵇ₁₁*𝒂₃₍₃₎;
                Bᵇ₂₂*𝒂₃₍₁₎   Bᵇ₂₂*𝒂₃₍₂₎   Bᵇ₂₂*𝒂₃₍₃₎;
                2*Bᵇ₁₂*𝒂₃₍₁₎ 2*Bᵇ₁₂*𝒂₃₍₂₎ 2*Bᵇ₁₂*𝒂₃₍₃₎]

                kᵇ = Dᵇ*𝐁ᵇᵢ'*𝐃*𝐁ᵇⱼ
                k[3*I-2,3*J-2] += kᵇ[1,1]*𝑤
                k[3*I-2,3*J-1] += kᵇ[1,2]*𝑤
                k[3*I-2,3*J]   += kᵇ[1,3]*𝑤
                k[3*I-1,3*J-2] += kᵇ[2,1]*𝑤
                k[3*I-1,3*J-1] += kᵇ[2,2]*𝑤
                k[3*I-1,3*J]   += kᵇ[2,3]*𝑤
                k[3*I,3*J-2]   += kᵇ[3,1]*𝑤
                k[3*I,3*J-1]   += kᵇ[3,2]*𝑤
                k[3*I,3*J]     += kᵇ[3,3]*𝑤
            end
        end
    end
end

"""
    ∫δθθdΓ(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}, E, ν, h)

Boundary term coupling the rotation trace, `∫ δθ θ dΓ`.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.g`, `ξ.𝑤`.
"""
function ∫δθθdΓ(ap::T, k::AbstractMatrix{Float64}, f::AbstractVector{Float64}; E::Real, ν::Real, h::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        n¹ = ξ.n₁
        n² = ξ.n₂
        𝒂₃₍₁₎ = ξ.𝒂₃₍₁₎
        𝒂₃₍₂₎ = ξ.𝒂₃₍₂₎
        𝒂₃₍₃₎ = ξ.𝒂₃₍₃₎
        θ = ξ.θ
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            𝐁ₙᵢ = [(B₁[i]*n¹ + B₂[i]*n²)*𝒂₃₍₁₎;(B₁[i]*n¹ + B₂[i]*n²)*𝒂₃₍₂₎;(B₁[i]*n¹ + B₂[i]*n²)*𝒂₃₍₃₎]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                𝐁ₙⱼ = [(B₁[j]*n¹ + B₂[j]*n²)*𝒂₃₍₁₎;(B₁[j]*n¹ + B₂[j]*n²)*𝒂₃₍₂₎;(B₁[j]*n¹ + B₂[j]*n²)*𝒂₃₍₃₎]
                k[3*I-2,3*J-2] += α*𝐁ₙᵢ[1]*𝐁ₙⱼ[1]*𝑤
                k[3*I-2,3*J-1] += α*𝐁ₙᵢ[1]*𝐁ₙⱼ[2]*𝑤
                k[3*I-2,3*J]   += α*𝐁ₙᵢ[1]*𝐁ₙⱼ[3]*𝑤
                k[3*I-1,3*J-2] += α*𝐁ₙᵢ[2]*𝐁ₙⱼ[1]*𝑤
                k[3*I-1,3*J-1] += α*𝐁ₙᵢ[2]*𝐁ₙⱼ[2]*𝑤
                k[3*I-1,3*J]   += α*𝐁ₙᵢ[2]*𝐁ₙⱼ[3]*𝑤
                k[3*I,3*J-2]   += α*𝐁ₙᵢ[3]*𝐁ₙⱼ[1]*𝑤
                k[3*I,3*J-1]   += α*𝐁ₙᵢ[3]*𝐁ₙⱼ[2]*𝑤
                k[3*I,3*J]     += α*𝐁ₙᵢ[3]*𝐁ₙⱼ[3]*𝑤
            end
            f[3*I-2] += α*𝐁ₙᵢ[1]*θ*𝑤
            f[3*I-1] += α*𝐁ₙᵢ[2]*θ*𝑤
            f[3*I]   += α*𝐁ₙᵢ[3]*θ*𝑤
        end
    end
end

"""
    ∫δεCεdΩ(ap::T;k::AbstractMatrix{Float64}, E, ν, h)

Membrane energy in strain form, `∫ δε : C ε dΩ`.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set∇𝝭!` — reads `ξ[:𝝭]`, the first derivatives and the surface metric.
"""
function ∫δεCεdΩ(ap::T, k::AbstractMatrix{Float64}; E::Real, ν::Real, h::Real) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        Dᵐ = E*h
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        a¹¹ = ξ.a¹¹
        a¹² = ξ.a¹²
        a²² = ξ.a²²
        𝒂₁₍₁₎ = ξ.𝒂₁₍₁₎
        𝒂₁₍₂₎ = ξ.𝒂₁₍₂₎
        𝒂₁₍₃₎ = ξ.𝒂₁₍₃₎
        𝒂₂₍₁₎ = ξ.𝒂₂₍₁₎
        𝒂₂₍₂₎ = ξ.𝒂₂₍₂₎
        𝒂₂₍₃₎ = ξ.𝒂₂₍₃₎

        𝐃 = [                  a¹¹*a¹¹ ν*a¹¹*a²² + (1-ν)*a¹²*a¹²                            a¹¹*a¹²;
        ν*a¹¹*a²² + (1-ν)*a¹²*a¹²                   a²²*a²²                            a²²*a¹²;
        a¹¹*a¹²                   a²²*a¹² 0.5*((1-ν)*a¹¹*a²² + (1+ν)*a¹²*a¹²)]

        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            𝐁ᵐᵢ = [𝒂₁₍₁₎*N[i] 𝒂₁₍₂₎*N[i] 𝒂₁₍₃₎*N[i]        0.0       0.0         0.0;
            0.0        0.0        0.0 𝒂₂₍₁₎*N[i] 𝒂₂₍₂₎*N[i] 𝒂₂₍₃₎*N[i];
            𝒂₂₍₁₎*N[i] 𝒂₂₍₂₎*N[i] 𝒂₂₍₃₎*N[i] 𝒂₁₍₁₎*N[i] 𝒂₁₍₂₎*N[i] 𝒂₁₍₃₎*N[i]]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                𝐁ᵐⱼ = [𝒂₁₍₁₎*N[j] 𝒂₁₍₂₎*N[j] 𝒂₁₍₃₎*N[j]        0.0       0.0         0.0;
                0.0        0.0        0.0 𝒂₂₍₁₎*N[j] 𝒂₂₍₂₎*N[j] 𝒂₂₍₃₎*N[j];
                𝒂₂₍₁₎*N[j] 𝒂₂₍₂₎*N[j] 𝒂₂₍₃₎*N[j] 𝒂₁₍₁₎*N[j] 𝒂₁₍₂₎*N[j] 𝒂₁₍₃₎*N[j]]
                kᵐ = Dᵐ*𝐁ᵐᵢ'*𝐃*𝐁ᵐⱼ
                k[6*I-5,6*J-5] -= kᵐ[1,1]*𝑤
                k[6*I-4,6*J-5] -= kᵐ[2,1]*𝑤
                k[6*I-3,6*J-5] -= kᵐ[3,1]*𝑤
                k[6*I-2,6*J-5] -= kᵐ[4,1]*𝑤
                k[6*I-1,6*J-5] -= kᵐ[5,1]*𝑤
                k[6*I,6*J-5]   -= kᵐ[6,1]*𝑤
                k[6*I-5,6*J-4] -= kᵐ[1,2]*𝑤
                k[6*I-4,6*J-4] -= kᵐ[2,2]*𝑤
                k[6*I-3,6*J-4] -= kᵐ[3,2]*𝑤
                k[6*I-2,6*J-4] -= kᵐ[4,2]*𝑤
                k[6*I-1,6*J-4] -= kᵐ[5,2]*𝑤
                k[6*I,6*J-4]   -= kᵐ[6,2]*𝑤
                k[6*I-5,6*J-3] -= kᵐ[1,3]*𝑤
                k[6*I-4,6*J-3] -= kᵐ[2,3]*𝑤
                k[6*I-3,6*J-3] -= kᵐ[3,3]*𝑤
                k[6*I-2,6*J-3] -= kᵐ[4,3]*𝑤
                k[6*I-1,6*J-3] -= kᵐ[5,3]*𝑤
                k[6*I,6*J-3]   -= kᵐ[6,3]*𝑤
                k[6*I-5,6*J-2] -= kᵐ[1,4]*𝑤
                k[6*I-4,6*J-2] -= kᵐ[2,4]*𝑤
                k[6*I-3,6*J-2] -= kᵐ[3,4]*𝑤
                k[6*I-2,6*J-2] -= kᵐ[4,4]*𝑤
                k[6*I-1,6*J-2] -= kᵐ[5,4]*𝑤
                k[6*I,6*J-2]   -= kᵐ[6,4]*𝑤
                k[6*I-5,6*J-1] -= kᵐ[1,5]*𝑤
                k[6*I-4,6*J-1] -= kᵐ[2,5]*𝑤
                k[6*I-3,6*J-1] -= kᵐ[3,5]*𝑤
                k[6*I-2,6*J-1] -= kᵐ[4,5]*𝑤
                k[6*I-1,6*J-1] -= kᵐ[5,5]*𝑤
                k[6*I,6*J-1]   -= kᵐ[6,5]*𝑤
                k[6*I-5,6*J]   -= kᵐ[1,6]*𝑤
                k[6*I-4,6*J]   -= kᵐ[2,6]*𝑤
                k[6*I-3,6*J]   -= kᵐ[3,6]*𝑤
                k[6*I-2,6*J]   -= kᵐ[4,6]*𝑤
                k[6*I-1,6*J]   -= kᵐ[5,6]*𝑤
                k[6*I,6*J]     -= kᵐ[6,6]*𝑤
            end
        end
    end
end

"""
    ∫δNεdΩ(ap::T;k::AbstractMatrix{Float64}, E, ν, h)

Membrane coupling `∫ δN : ε dΩ`, the stress–strain dual of
[`∫δεCεdΩ`](@ref).

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set∇𝝭!` — reads `ξ[:𝝭]`, the first derivatives and the surface metric.
"""
function ∫δNεdΩ(ap::T, k::AbstractMatrix{Float64}; E::Real, ν::Real, h::Real) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[6*I-5,6*J-5] += N[i]*N[j]*𝑤
                k[6*I-4,6*J-4] += N[i]*N[j]*𝑤
                k[6*I-3,6*J-3] += N[i]*N[j]*𝑤
                k[6*I-2,6*J-2] += N[i]*N[j]*𝑤
                k[6*I-1,6*J-1] += N[i]*N[j]*𝑤
                k[6*I,6*J]     += N[i]*N[j]*𝑤
            end
        end
    end
end

"""
    ∫δκCκdΩ(ap::T;k::AbstractMatrix{Float64}, E, ν, h)

Bending energy in curvature form, `∫ δκ : C κ dΩ`.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set∇²𝝭!` — reads the second derivatives and the surface metric.
"""
function ∫δκCκdΩ(ap::T, k::AbstractMatrix{Float64}; E::Real, ν::Real, h::Real) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        Dᵇ = E*h^3/12
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
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

        𝐃 = [                  a¹¹*a¹¹ ν*a¹¹*a²² + (1-ν)*a¹²*a¹²                            a¹¹*a¹²;
        ν*a¹¹*a²² + (1-ν)*a¹²*a¹²                   a²²*a²²                            a²²*a¹²;
        a¹¹*a¹²                   a²²*a¹² 0.5*((1-ν)*a¹¹*a²² + (1+ν)*a¹²*a¹²)]

        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            𝐁ᵇᵢ = [𝒂₃₍₁₎*N[i] 𝒂₃₍₂₎*N[i] 𝒂₃₍₃₎*N[i] 0.0 0.0 0.0 0.0 0.0 0.0;
            0.0 0.0 0.0 𝒂₃₍₁₎*N[i] 𝒂₃₍₂₎*N[i] 𝒂₃₍₃₎*N[i] 0.0 0.0 0.0;
            0.0 0.0 0.0 0.0 0.0 0.0 𝒂₃₍₁₎*N[i] 𝒂₃₍₂₎*N[i] 𝒂₃₍₃₎*N[i]]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                𝐁ᵇⱼ = [𝒂₃₍₁₎*N[j] 𝒂₃₍₂₎*N[j] 𝒂₃₍₃₎*N[j] 0.0 0.0 0.0 0.0 0.0 0.0;
                0.0 0.0 0.0 𝒂₃₍₁₎*N[j] 𝒂₃₍₂₎*N[j] 𝒂₃₍₃₎*N[j] 0.0 0.0 0.0;
                0.0 0.0 0.0 0.0 0.0 0.0 𝒂₃₍₁₎*N[j] 𝒂₃₍₂₎*N[j] 𝒂₃₍₃₎*N[j]]

                kᵇ = Dᵇ*𝐁ᵇᵢ'*𝐃*𝐁ᵇⱼ
                k[9*I-8,9*J-8] -= kᵇ[1,1]*𝑤
                k[9*I-7,9*J-8] -= kᵇ[2,1]*𝑤
                k[9*I-6,9*J-8] -= kᵇ[3,1]*𝑤
                k[9*I-5,9*J-8] -= kᵇ[4,1]*𝑤
                k[9*I-4,9*J-8] -= kᵇ[5,1]*𝑤
                k[9*I-3,9*J-8] -= kᵇ[6,1]*𝑤
                k[9*I-2,9*J-8] -= kᵇ[7,1]*𝑤
                k[9*I-1,9*J-8] -= kᵇ[8,1]*𝑤
                k[9*I,9*J-8]   -= kᵇ[9,1]*𝑤
                k[9*I-8,9*J-7] -= kᵇ[1,2]*𝑤
                k[9*I-7,9*J-7] -= kᵇ[2,2]*𝑤
                k[9*I-6,9*J-7] -= kᵇ[3,2]*𝑤
                k[9*I-5,9*J-7] -= kᵇ[4,2]*𝑤
                k[9*I-4,9*J-7] -= kᵇ[5,2]*𝑤
                k[9*I-3,9*J-7] -= kᵇ[6,2]*𝑤
                k[9*I-2,9*J-7] -= kᵇ[7,2]*𝑤
                k[9*I-1,9*J-7] -= kᵇ[8,2]*𝑤
                k[9*I,9*J-7]   -= kᵇ[9,2]*𝑤
                k[9*I-8,9*J-6] -= kᵇ[1,3]*𝑤
                k[9*I-7,9*J-6] -= kᵇ[2,3]*𝑤
                k[9*I-6,9*J-6] -= kᵇ[3,3]*𝑤
                k[9*I-5,9*J-6] -= kᵇ[4,3]*𝑤
                k[9*I-4,9*J-6] -= kᵇ[5,3]*𝑤
                k[9*I-3,9*J-6] -= kᵇ[6,3]*𝑤
                k[9*I-2,9*J-6] -= kᵇ[7,3]*𝑤
                k[9*I-1,9*J-6] -= kᵇ[8,3]*𝑤
                k[9*I,9*J-6]   -= kᵇ[9,3]*𝑤
                k[9*I-8,9*J-5] -= kᵇ[1,4]*𝑤
                k[9*I-7,9*J-5] -= kᵇ[2,4]*𝑤
                k[9*I-6,9*J-5] -= kᵇ[3,4]*𝑤
                k[9*I-5,9*J-5] -= kᵇ[4,4]*𝑤
                k[9*I-4,9*J-5] -= kᵇ[5,4]*𝑤
                k[9*I-3,9*J-5] -= kᵇ[6,4]*𝑤
                k[9*I-2,9*J-5] -= kᵇ[7,4]*𝑤
                k[9*I-1,9*J-5] -= kᵇ[8,4]*𝑤
                k[9*I,9*J-5]   -= kᵇ[9,4]*𝑤
                k[9*I-8,9*J-4] -= kᵇ[1,5]*𝑤
                k[9*I-7,9*J-4] -= kᵇ[2,5]*𝑤
                k[9*I-6,9*J-4] -= kᵇ[3,5]*𝑤
                k[9*I-5,9*J-4] -= kᵇ[4,5]*𝑤
                k[9*I-4,9*J-4] -= kᵇ[5,5]*𝑤
                k[9*I-3,9*J-4] -= kᵇ[6,5]*𝑤
                k[9*I-2,9*J-4] -= kᵇ[7,5]*𝑤
                k[9*I-1,9*J-4] -= kᵇ[8,5]*𝑤
                k[9*I,9*J-4]   -= kᵇ[9,5]*𝑤
                k[9*I-8,9*J-3] -= kᵇ[1,6]*𝑤
                k[9*I-7,9*J-3] -= kᵇ[2,6]*𝑤
                k[9*I-6,9*J-3] -= kᵇ[3,6]*𝑤
                k[9*I-5,9*J-3] -= kᵇ[4,6]*𝑤
                k[9*I-4,9*J-3] -= kᵇ[5,6]*𝑤
                k[9*I-3,9*J-3] -= kᵇ[6,6]*𝑤
                k[9*I-2,9*J-3] -= kᵇ[7,6]*𝑤
                k[9*I-1,9*J-3] -= kᵇ[8,6]*𝑤
                k[9*I,9*J-3]   -= kᵇ[9,6]*𝑤
                k[9*I-8,9*J-2] -= kᵇ[1,7]*𝑤
                k[9*I-7,9*J-2] -= kᵇ[2,7]*𝑤
                k[9*I-6,9*J-2] -= kᵇ[3,7]*𝑤
                k[9*I-5,9*J-2] -= kᵇ[4,7]*𝑤
                k[9*I-4,9*J-2] -= kᵇ[5,7]*𝑤
                k[9*I-3,9*J-2] -= kᵇ[6,7]*𝑤
                k[9*I-2,9*J-2] -= kᵇ[7,7]*𝑤
                k[9*I-1,9*J-2] -= kᵇ[8,7]*𝑤
                k[9*I,9*J-2]   -= kᵇ[9,7]*𝑤
                k[9*I-8,9*J-1] -= kᵇ[1,8]*𝑤
                k[9*I-7,9*J-1] -= kᵇ[2,8]*𝑤
                k[9*I-6,9*J-1] -= kᵇ[3,8]*𝑤
                k[9*I-5,9*J-1] -= kᵇ[4,8]*𝑤
                k[9*I-4,9*J-1] -= kᵇ[5,8]*𝑤
                k[9*I-3,9*J-1] -= kᵇ[6,8]*𝑤
                k[9*I-2,9*J-1] -= kᵇ[7,8]*𝑤
                k[9*I-1,9*J-1] -= kᵇ[8,8]*𝑤
                k[9*I,9*J-1]   -= kᵇ[9,8]*𝑤
                k[9*I-8,9*J]   -= kᵇ[1,9]*𝑤
                k[9*I-7,9*J]   -= kᵇ[2,9]*𝑤
                k[9*I-6,9*J]   -= kᵇ[3,9]*𝑤
                k[9*I-5,9*J]   -= kᵇ[4,9]*𝑤
                k[9*I-4,9*J]   -= kᵇ[5,9]*𝑤
                k[9*I-3,9*J]   -= kᵇ[6,9]*𝑤
                k[9*I-2,9*J]   -= kᵇ[7,9]*𝑤
                k[9*I-1,9*J]   -= kᵇ[8,9]*𝑤
                k[9*I,9*J]     -= kᵇ[9,9]*𝑤
            end
        end
    end
end

"""
    ∫δMκdΩ(ap::T;k::AbstractMatrix{Float64}, E, ν, h)

Bending coupling `∫ δM : κ dΩ`.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set∇²𝝭!` — reads the second derivatives and the surface metric.
"""
function ∫δMκdΩ(ap::T, k::AbstractMatrix{Float64}; E::Real, ν::Real, h::Real) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[9*I-8,9*J-8] += N[i]*N[j]*𝑤
                k[9*I-7,9*J-7] += N[i]*N[j]*𝑤
                k[9*I-6,9*J-6] += N[i]*N[j]*𝑤
                k[9*I-5,9*J-5] += N[i]*N[j]*𝑤
                k[9*I-4,9*J-4] += N[i]*N[j]*𝑤
                k[9*I-3,9*J-3] += N[i]*N[j]*𝑤
                k[9*I-2,9*J-2] += N[i]*N[j]*𝑤
                k[9*I-1,9*J-1] += N[i]*N[j]*𝑤
                k[9*I,9*J]     += N[i]*N[j]*𝑤
            end
        end
    end
end

"""
    ∫𝒏𝑵𝒈dΓ_Nitsche(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}, E, ν, h)

Nitsche weak boundary condition on the shell, `∫ 𝒏·𝑵 𝑔 dΓ`.

The Nitsche variant adds the consistency and penalty terms that make the
imposition of a Dirichlet condition stable without extra unknowns; the `_HR`
sibling is the Hellinger–Reissner form of the same term.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.n₁`, `ξ.n₂`, `ξ.n₃`, `ξ.g`, `ξ.𝑤`.
"""
function ∫𝒏𝑵𝒈dΓ_Nitsche(ap::T, k::AbstractMatrix{Float64}, f::AbstractVector{Float64}; E::Real, ν::Real, h::Real) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        Dᵐ = E*h/(1-ν^2)
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        a¹¹ = ξ.a¹¹
        a¹² = ξ.a¹²
        a²² = ξ.a²²
        𝒂₁₍₁₎ = ξ.𝒂₁₍₁₎
        𝒂₁₍₂₎ = ξ.𝒂₁₍₂₎
        𝒂₁₍₃₎ = ξ.𝒂₁₍₃₎
        𝒂₂₍₁₎ = ξ.𝒂₂₍₁₎
        𝒂₂₍₂₎ = ξ.𝒂₂₍₂₎
        𝒂₂₍₃₎ = ξ.𝒂₂₍₃₎

        C¹¹¹¹ = a¹¹*a¹¹*Dᵐ
        C²²²² = a²²*a²²*Dᵐ
        C¹¹²² = (ν*a¹¹*a²² + (1-ν)*a¹²*a¹²)*Dᵐ
        C¹¹¹² = a¹¹*a¹²*Dᵐ
        C²²¹² = a²²*a¹²*Dᵐ
        C¹²¹² = (0.5*((1-ν)*a¹¹*a²² + (1+ν)*a¹²*a¹²))*Dᵐ

        n₁ = ξ.n₁
        n₂ = ξ.n₂

        n₁₁ = ξ.n₁₁
        n₁₂ = ξ.n₁₂
        n₁₃ = ξ.n₁₃
        n₂₂ = ξ.n₂₂
        n₂₃ = ξ.n₂₃
        n₃₃ = ξ.n₃₃
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        g₃ = ξ.g₃

        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            ε₁₁₍₁₎ = 𝒂₁₍₁₎*B₁[i]
            ε₁₁₍₂₎ = 𝒂₁₍₂₎*B₁[i]
            ε₁₁₍₃₎ = 𝒂₁₍₃₎*B₁[i]
            ε₂₂₍₁₎ = 𝒂₂₍₁₎*B₂[i]
            ε₂₂₍₂₎ = 𝒂₂₍₂₎*B₂[i]
            ε₂₂₍₃₎ = 𝒂₂₍₃₎*B₂[i]
            ε₁₂₍₁₎ = 𝒂₁₍₁₎*B₂[i] + 𝒂₂₍₁₎*B₁[i]
            ε₁₂₍₂₎ = 𝒂₁₍₂₎*B₂[i] + 𝒂₂₍₂₎*B₁[i]
            ε₁₂₍₃₎ = 𝒂₁₍₃₎*B₂[i] + 𝒂₂₍₃₎*B₁[i]

            N¹¹₍₁₎ = C¹¹¹¹*ε₁₁₍₁₎ + C¹¹²²*ε₂₂₍₁₎ + C¹¹¹²*ε₁₂₍₁₎
            N¹¹₍₂₎ = C¹¹¹¹*ε₁₁₍₂₎ + C¹¹²²*ε₂₂₍₂₎ + C¹¹¹²*ε₁₂₍₂₎
            N¹¹₍₃₎ = C¹¹¹¹*ε₁₁₍₃₎ + C¹¹²²*ε₂₂₍₃₎ + C¹¹¹²*ε₁₂₍₃₎
            N²²₍₁₎ = C¹¹²²*ε₁₁₍₁₎ + C²²²²*ε₂₂₍₁₎ + C²²¹²*ε₁₂₍₁₎
            N²²₍₂₎ = C¹¹²²*ε₁₁₍₂₎ + C²²²²*ε₂₂₍₂₎ + C²²¹²*ε₁₂₍₂₎
            N²²₍₃₎ = C¹¹²²*ε₁₁₍₃₎ + C²²²²*ε₂₂₍₃₎ + C²²¹²*ε₁₂₍₃₎
            N¹²₍₁₎ = C¹¹¹²*ε₁₁₍₁₎ + C²²¹²*ε₂₂₍₁₎ + C¹²¹²*ε₁₂₍₁₎
            N¹²₍₂₎ = C¹¹¹²*ε₁₁₍₂₎ + C²²¹²*ε₂₂₍₂₎ + C¹²¹²*ε₁₂₍₂₎
            N¹²₍₃₎ = C¹¹¹²*ε₁₁₍₃₎ + C²²¹²*ε₂₂₍₃₎ + C¹²¹²*ε₁₂₍₃₎

            Bᵢ₍₁₁₎ = N¹¹₍₁₎*𝒂₁₍₁₎*n₁ + N¹²₍₁₎*(𝒂₁₍₁₎*n₂+𝒂₂₍₁₎*n₁) + N²²₍₁₎*𝒂₂₍₁₎*n₂
            Bᵢ₍₂₁₎ = N¹¹₍₂₎*𝒂₁₍₁₎*n₁ + N¹²₍₂₎*(𝒂₁₍₁₎*n₂+𝒂₂₍₁₎*n₁) + N²²₍₂₎*𝒂₂₍₁₎*n₂
            Bᵢ₍₃₁₎ = N¹¹₍₃₎*𝒂₁₍₁₎*n₁ + N¹²₍₃₎*(𝒂₁₍₁₎*n₂+𝒂₂₍₁₎*n₁) + N²²₍₃₎*𝒂₂₍₁₎*n₂
            Bᵢ₍₁₂₎ = N¹¹₍₁₎*𝒂₁₍₂₎*n₁ + N¹²₍₁₎*(𝒂₁₍₂₎*n₂+𝒂₂₍₂₎*n₁) + N²²₍₁₎*𝒂₂₍₂₎*n₂
            Bᵢ₍₂₂₎ = N¹¹₍₂₎*𝒂₁₍₂₎*n₁ + N¹²₍₂₎*(𝒂₁₍₂₎*n₂+𝒂₂₍₂₎*n₁) + N²²₍₂₎*𝒂₂₍₂₎*n₂
            Bᵢ₍₃₂₎ = N¹¹₍₃₎*𝒂₁₍₂₎*n₁ + N¹²₍₃₎*(𝒂₁₍₂₎*n₂+𝒂₂₍₂₎*n₁) + N²²₍₃₎*𝒂₂₍₂₎*n₂
            Bᵢ₍₁₃₎ = N¹¹₍₁₎*𝒂₁₍₃₎*n₁ + N¹²₍₁₎*(𝒂₁₍₃₎*n₂+𝒂₂₍₃₎*n₁) + N²²₍₁₎*𝒂₂₍₃₎*n₂
            Bᵢ₍₂₃₎ = N¹¹₍₂₎*𝒂₁₍₃₎*n₁ + N¹²₍₂₎*(𝒂₁₍₃₎*n₂+𝒂₂₍₃₎*n₁) + N²²₍₂₎*𝒂₂₍₃₎*n₂
            Bᵢ₍₃₃₎ = N¹¹₍₃₎*𝒂₁₍₃₎*n₁ + N¹²₍₃₎*(𝒂₁₍₃₎*n₂+𝒂₂₍₃₎*n₁) + N²²₍₃₎*𝒂₂₍₃₎*n₂

            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼

                ε₁₁₍₁₎ = 𝒂₁₍₁₎*B₁[j]
                ε₁₁₍₂₎ = 𝒂₁₍₂₎*B₁[j]
                ε₁₁₍₃₎ = 𝒂₁₍₃₎*B₁[j]
                ε₂₂₍₁₎ = 𝒂₂₍₁₎*B₂[j]
                ε₂₂₍₂₎ = 𝒂₂₍₂₎*B₂[j]
                ε₂₂₍₃₎ = 𝒂₂₍₃₎*B₂[j]
                ε₁₂₍₁₎ = 𝒂₁₍₁₎*B₂[j] + 𝒂₂₍₁₎*B₁[j]
                ε₁₂₍₂₎ = 𝒂₁₍₂₎*B₂[j] + 𝒂₂₍₂₎*B₁[j]
                ε₁₂₍₃₎ = 𝒂₁₍₃₎*B₂[j] + 𝒂₂₍₃₎*B₁[j]

                N¹¹₍₁₎ = C¹¹¹¹*ε₁₁₍₁₎ + C¹¹²²*ε₂₂₍₁₎ + C¹¹¹²*ε₁₂₍₁₎
                N¹¹₍₂₎ = C¹¹¹¹*ε₁₁₍₂₎ + C¹¹²²*ε₂₂₍₂₎ + C¹¹¹²*ε₁₂₍₂₎
                N¹¹₍₃₎ = C¹¹¹¹*ε₁₁₍₃₎ + C¹¹²²*ε₂₂₍₃₎ + C¹¹¹²*ε₁₂₍₃₎
                N²²₍₁₎ = C¹¹²²*ε₁₁₍₁₎ + C²²²²*ε₂₂₍₁₎ + C²²¹²*ε₁₂₍₁₎
                N²²₍₂₎ = C¹¹²²*ε₁₁₍₂₎ + C²²²²*ε₂₂₍₂₎ + C²²¹²*ε₁₂₍₂₎
                N²²₍₃₎ = C¹¹²²*ε₁₁₍₃₎ + C²²²²*ε₂₂₍₃₎ + C²²¹²*ε₁₂₍₃₎
                N¹²₍₁₎ = C¹¹¹²*ε₁₁₍₁₎ + C²²¹²*ε₂₂₍₁₎ + C¹²¹²*ε₁₂₍₁₎
                N¹²₍₂₎ = C¹¹¹²*ε₁₁₍₂₎ + C²²¹²*ε₂₂₍₂₎ + C¹²¹²*ε₁₂₍₂₎
                N¹²₍₃₎ = C¹¹¹²*ε₁₁₍₃₎ + C²²¹²*ε₂₂₍₃₎ + C¹²¹²*ε₁₂₍₃₎

                Bⱼ₍₁₁₎ = N¹¹₍₁₎*𝒂₁₍₁₎*n₁ + N¹²₍₁₎*(𝒂₁₍₁₎*n₂+𝒂₂₍₁₎*n₁) + N²²₍₁₎*𝒂₂₍₁₎*n₂
                Bⱼ₍₂₁₎ = N¹¹₍₁₎*𝒂₁₍₂₎*n₁ + N¹²₍₁₎*(𝒂₁₍₂₎*n₂+𝒂₂₍₂₎*n₁) + N²²₍₁₎*𝒂₂₍₂₎*n₂
                Bⱼ₍₃₁₎ = N¹¹₍₁₎*𝒂₁₍₃₎*n₁ + N¹²₍₁₎*(𝒂₁₍₃₎*n₂+𝒂₂₍₃₎*n₁) + N²²₍₁₎*𝒂₂₍₃₎*n₂
                Bⱼ₍₁₂₎ = N¹¹₍₂₎*𝒂₁₍₁₎*n₁ + N¹²₍₂₎*(𝒂₁₍₁₎*n₂+𝒂₂₍₁₎*n₁) + N²²₍₂₎*𝒂₂₍₁₎*n₂
                Bⱼ₍₂₂₎ = N¹¹₍₂₎*𝒂₁₍₂₎*n₁ + N¹²₍₂₎*(𝒂₁₍₂₎*n₂+𝒂₂₍₂₎*n₁) + N²²₍₂₎*𝒂₂₍₂₎*n₂
                Bⱼ₍₃₂₎ = N¹¹₍₂₎*𝒂₁₍₃₎*n₁ + N¹²₍₂₎*(𝒂₁₍₃₎*n₂+𝒂₂₍₃₎*n₁) + N²²₍₂₎*𝒂₂₍₃₎*n₂
                Bⱼ₍₁₃₎ = N¹¹₍₃₎*𝒂₁₍₁₎*n₁ + N¹²₍₃₎*(𝒂₁₍₁₎*n₂+𝒂₂₍₁₎*n₁) + N²²₍₃₎*𝒂₂₍₁₎*n₂
                Bⱼ₍₂₃₎ = N¹¹₍₃₎*𝒂₁₍₂₎*n₁ + N¹²₍₃₎*(𝒂₁₍₂₎*n₂+𝒂₂₍₂₎*n₁) + N²²₍₃₎*𝒂₂₍₂₎*n₂
                Bⱼ₍₃₃₎ = N¹¹₍₃₎*𝒂₁₍₃₎*n₁ + N¹²₍₃₎*(𝒂₁₍₃₎*n₂+𝒂₂₍₃₎*n₁) + N²²₍₃₎*𝒂₂₍₃₎*n₂

                k[3*I-2,3*J-2] -= ((Bᵢ₍₁₁₎*n₁₁+Bᵢ₍₁₂₎*n₁₂+Bᵢ₍₁₃₎*n₁₃)*N[j] + N[i]*(n₁₁*Bⱼ₍₁₁₎+n₁₂*Bⱼ₍₂₁₎+n₁₃*Bⱼ₍₃₁₎))*𝑤
                k[3*I-1,3*J-2] -= ((Bᵢ₍₂₁₎*n₁₁+Bᵢ₍₂₂₎*n₁₂+Bᵢ₍₂₃₎*n₁₃)*N[j] + N[i]*(n₁₂*Bⱼ₍₁₁₎+n₂₂*Bⱼ₍₂₁₎+n₂₃*Bⱼ₍₃₁₎))*𝑤
                k[3*I,3*J-2]   -= ((Bᵢ₍₃₁₎*n₁₁+Bᵢ₍₃₂₎*n₁₂+Bᵢ₍₃₃₎*n₁₃)*N[j] + N[i]*(n₁₃*Bⱼ₍₁₁₎+n₂₃*Bⱼ₍₂₁₎+n₃₃*Bⱼ₍₃₁₎))*𝑤
                k[3*I-2,3*J-1] -= ((Bᵢ₍₁₁₎*n₁₂+Bᵢ₍₁₂₎*n₂₂+Bᵢ₍₁₃₎*n₂₃)*N[j] + N[i]*(n₁₁*Bⱼ₍₁₂₎+n₁₂*Bⱼ₍₂₂₎+n₁₃*Bⱼ₍₃₂₎))*𝑤
                k[3*I-1,3*J-1] -= ((Bᵢ₍₂₁₎*n₁₂+Bᵢ₍₂₂₎*n₂₂+Bᵢ₍₂₃₎*n₂₃)*N[j] + N[i]*(n₁₂*Bⱼ₍₁₂₎+n₂₂*Bⱼ₍₂₂₎+n₂₃*Bⱼ₍₃₂₎))*𝑤
                k[3*I,3*J-1]   -= ((Bᵢ₍₃₁₎*n₁₂+Bᵢ₍₃₂₎*n₂₂+Bᵢ₍₃₃₎*n₂₃)*N[j] + N[i]*(n₁₃*Bⱼ₍₁₂₎+n₂₃*Bⱼ₍₂₂₎+n₃₃*Bⱼ₍₃₂₎))*𝑤
                k[3*I-2,3*J]   -= ((Bᵢ₍₁₁₎*n₁₃+Bᵢ₍₁₂₎*n₂₃+Bᵢ₍₁₃₎*n₃₃)*N[j] + N[i]*(n₁₁*Bⱼ₍₁₃₎+n₁₂*Bⱼ₍₂₃₎+n₁₃*Bⱼ₍₃₃₎))*𝑤
                k[3*I-1,3*J]   -= ((Bᵢ₍₂₁₎*n₁₃+Bᵢ₍₂₂₎*n₂₃+Bᵢ₍₂₃₎*n₃₃)*N[j] + N[i]*(n₁₂*Bⱼ₍₁₃₎+n₂₂*Bⱼ₍₂₃₎+n₂₃*Bⱼ₍₃₃₎))*𝑤
                k[3*I,3*J]     -= ((Bᵢ₍₃₁₎*n₁₃+Bᵢ₍₃₂₎*n₂₃+Bᵢ₍₃₃₎*n₃₃)*N[j] + N[i]*(n₁₃*Bⱼ₍₁₃₎+n₂₃*Bⱼ₍₂₃₎+n₃₃*Bⱼ₍₃₃₎))*𝑤

            end

            f[3*I-2] -= ((Bᵢ₍₁₁₎*n₁₁+Bᵢ₍₁₂₎*n₁₂+Bᵢ₍₁₃₎*n₁₃)*g₁
            + (Bᵢ₍₁₁₎*n₁₂+Bᵢ₍₁₂₎*n₂₂+Bᵢ₍₁₃₎*n₂₃)*g₂
            + (Bᵢ₍₁₁₎*n₁₃+Bᵢ₍₁₂₎*n₂₃+Bᵢ₍₁₃₎*n₃₃)*g₃)*𝑤
            f[3*I-1] -= ((Bᵢ₍₂₁₎*n₁₁+Bᵢ₍₂₂₎*n₁₂+Bᵢ₍₂₃₎*n₁₃)*g₁
            + (Bᵢ₍₂₁₎*n₁₂+Bᵢ₍₂₂₎*n₂₂+Bᵢ₍₂₃₎*n₂₃)*g₂
            + (Bᵢ₍₂₁₎*n₁₃+Bᵢ₍₂₂₎*n₂₃+Bᵢ₍₂₃₎*n₃₃)*g₃)*𝑤
            f[3*I]   -= ((Bᵢ₍₃₁₎*n₁₁+Bᵢ₍₃₂₎*n₁₂+Bᵢ₍₃₃₎*n₁₃)*g₁
            + (Bᵢ₍₃₁₎*n₁₂+Bᵢ₍₃₂₎*n₂₂+Bᵢ₍₃₃₎*n₂₃)*g₂
            + (Bᵢ₍₃₁₎*n₁₃+Bᵢ₍₃₂₎*n₂₃+Bᵢ₍₃₃₎*n₃₃)*g₃)*𝑤
        end
    end
end

"""
    ∫𝒏𝑵𝒈dΓ_HR(a::T,b::S;k::AbstractMatrix{Float64},f::AbstractVector{Float64})

Overload or Nitsche/Hellinger–Reissner sibling of the same shell term;
see the documented operators with a similar name.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set∇²𝝭!` — reads the shell fields described in the sibling operator.
"""
function ∫𝒏𝑵𝒈dΓ_HR(a::T,b::S, k::AbstractMatrix{Float64}, f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵃ = a.𝓒;𝓖ᵃ = a.𝓖
    𝓒ᵇ = b.𝓒;𝓖ᵇ = b.𝓖
    for (ξᵃ,ξᵇ) in zip(𝓖ᵃ,𝓖ᵇ)
        𝑤 = ξᵃ.𝑤
        N = ξᵃ[:𝝭]
        N̄ = ξᵇ[:𝝭]
        n₁ = ξᵃ.n₁
        n₂ = ξᵃ.n₂

        n₁₁ = ξᵇ.n₁₁
        n₁₂ = ξᵇ.n₁₂
        n₁₃ = ξᵇ.n₁₃
        n₂₂ = ξᵇ.n₂₂
        n₂₃ = ξᵇ.n₂₃
        n₃₃ = ξᵇ.n₃₃
        g₁ = ξᵇ.g₁
        g₂ = ξᵇ.g₂
        g₃ = ξᵇ.g₃
        for (i,xᵢ) in enumerate(𝓒ᵃ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵇ)
                J = xⱼ.𝐼
                k[6*I-5,3*J-2] += n₁*N[i]*n₁₁*N̄[j]*𝑤
                k[6*I-4,3*J-2] += n₁*N[i]*n₁₂*N̄[j]*𝑤
                k[6*I-3,3*J-2] += n₁*N[i]*n₁₃*N̄[j]*𝑤
                k[6*I-2,3*J-2] += n₂*N[i]*n₁₁*N̄[j]*𝑤
                k[6*I-1,3*J-2] += n₂*N[i]*n₁₂*N̄[j]*𝑤
                k[6*I,3*J-2]   += n₂*N[i]*n₁₃*N̄[j]*𝑤
                k[6*I-5,3*J-1] += n₁*N[i]*n₁₂*N̄[j]*𝑤
                k[6*I-4,3*J-1] += n₁*N[i]*n₂₂*N̄[j]*𝑤
                k[6*I-3,3*J-1] += n₁*N[i]*n₂₃*N̄[j]*𝑤
                k[6*I-2,3*J-1] += n₂*N[i]*n₁₂*N̄[j]*𝑤
                k[6*I-1,3*J-1] += n₂*N[i]*n₂₂*N̄[j]*𝑤
                k[6*I,3*J-1]   += n₂*N[i]*n₂₃*N̄[j]*𝑤
                k[6*I-5,3*J]   += n₁*N[i]*n₁₃*N̄[j]*𝑤
                k[6*I-4,3*J]   += n₁*N[i]*n₂₃*N̄[j]*𝑤
                k[6*I-3,3*J]   += n₁*N[i]*n₃₃*N̄[j]*𝑤
                k[6*I-2,3*J]   += n₂*N[i]*n₁₃*N̄[j]*𝑤
                k[6*I-1,3*J]   += n₂*N[i]*n₂₃*N̄[j]*𝑤
                k[6*I,3*J]     += n₂*N[i]*n₃₃*N̄[j]*𝑤
            end
            f[6*I-5] += n₁*N[i]*(n₁₁*g₁+n₁₂*g₂+n₁₃*g₃)*𝑤
            f[6*I-4] += n₁*N[i]*(n₁₂*g₁+n₂₂*g₂+n₂₃*g₃)*𝑤
            f[6*I-3] += n₁*N[i]*(n₁₃*g₁+n₂₃*g₂+n₃₃*g₃)*𝑤
            f[6*I-2] += n₂*N[i]*(n₁₁*g₁+n₁₂*g₂+n₁₃*g₃)*𝑤
            f[6*I-1] += n₂*N[i]*(n₁₂*g₁+n₂₂*g₂+n₂₃*g₃)*𝑤
            f[6*I]   += n₂*N[i]*(n₁₃*g₁+n₂₃*g₂+n₃₃*g₃)*𝑤
        end
    end
end

"""
    ∫𝒏𝑵𝒗dΓ(a::T,b::S;k::AbstractMatrix{Float64})

Overload or Nitsche/Hellinger–Reissner sibling of the same shell term;
see the documented operators with a similar name.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set∇²𝝭!` — reads the shell fields described in the sibling operator.
"""
function ∫𝒏𝑵𝒗dΓ(a::T,b::S, k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵃ = a.𝓒;𝓖ᵃ = a.𝓖
    𝓒ᵇ = b.𝓒;𝓖ᵇ = b.𝓖
    for (ξᵃ,ξᵇ) in zip(𝓖ᵃ,𝓖ᵇ)
        𝑤 = ξᵃ.𝑤
        N = ξᵃ[:𝝭]
        N̄ = ξᵇ[:𝝭]
        n₁ = ξᵃ.n₁
        n₂ = ξᵃ.n₂
        # println(N[1]-1.0)
        # println(N[2]-ξᵃ.x)
        # println(N[3]-ξᵃ.y)
        for (i,xᵢ) in enumerate(𝓒ᵃ)
            I = xᵢ.𝐼
            u1 = 0.0
            u2 = 0.0
            u3 = 0.0
            for (j,xⱼ) in enumerate(𝓒ᵇ)
                J = xⱼ.𝐼
                k[6*I-5,3*J-2] -= N[i]*N̄[j]*n₁*𝑤
                k[6*I-4,3*J-1] -= N[i]*N̄[j]*n₁*𝑤
                k[6*I-3,3*J]   -= N[i]*N̄[j]*n₁*𝑤
                k[6*I-2,3*J-2] -= N[i]*N̄[j]*n₂*𝑤
                k[6*I-1,3*J-1] -= N[i]*N̄[j]*n₂*𝑤
                k[6*I,3*J]     -= N[i]*N̄[j]*n₂*𝑤
                # u1 += N̄[j]*xⱼ.d₁
                # u2 += N̄[j]*xⱼ.d₂
                # u3 += N̄[j]*xⱼ.d₃
            end
            # println(ξᵃ.x)
            # println(u1)
            # println(u1 - ξᵃ.x)
            # println(u1 - ξᵇ.x)
        end
    end
end

"""
    ∫∇𝑵𝒗dΩ(a::T,b::S;k::AbstractMatrix{Float64})

Divergence coupling `∫ ∇·𝑵 v dΩ` for the shell.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set∇𝝭!` — reads `ξ[:𝝭]`, the first derivatives and the surface metric.
"""
function ∫∇𝑵𝒗dΩ(a::T,b::S, k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵃ = a.𝓒;𝓖ᵃ = a.𝓖
    𝓒ᵇ = b.𝓒;𝓖ᵇ = b.𝓖
    for (ξᵃ,ξᵇ) in zip(𝓖ᵃ,𝓖ᵇ)
        𝑤 = ξᵃ.𝑤
        N = ξᵃ[:𝝭]
        B₁ = ξᵃ[:∂𝝭∂x]
        B₂ = ξᵃ[:∂𝝭∂y]
        N̄ = ξᵇ[:𝝭]

        Γ¹₁₁ = ξᵃ.Γ¹₁₁
        Γ¹₁₂ = ξᵃ.Γ¹₁₂
        Γ¹₂₂ = ξᵃ.Γ¹₂₂
        Γ²₁₁ = ξᵃ.Γ²₁₁
        Γ²₁₂ = ξᵃ.Γ²₁₂
        Γ²₂₂ = ξᵃ.Γ²₂₂
        Γᵞᵧ₁ = Γ¹₁₁ + Γ²₁₂
        Γᵞᵧ₂ = Γ¹₁₂ + Γ²₂₂
        # println(Γᵞᵧ₁)
        # println(Γᵞᵧ₂)

        for (i,xᵢ) in enumerate(𝓒ᵃ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵇ)
                J = xⱼ.𝐼
                k[6*I-5,3*J-2] += (B₁[i] + Γᵞᵧ₁*N[i])*N̄[j]*𝑤
                k[6*I-4,3*J-1] += (B₁[i] + Γᵞᵧ₁*N[i])*N̄[j]*𝑤
                k[6*I-3,3*J]   += (B₁[i] + Γᵞᵧ₁*N[i])*N̄[j]*𝑤
                k[6*I-2,3*J-2] += (B₂[i] + Γᵞᵧ₂*N[i])*N̄[j]*𝑤
                k[6*I-1,3*J-1] += (B₂[i] + Γᵞᵧ₂*N[i])*N̄[j]*𝑤
                k[6*I,3*J]     += (B₂[i] + Γᵞᵧ₂*N[i])*N̄[j]*𝑤
                # k[6*I-5,3*J-2] += (B₁[i])*N̄[j]*𝑤
                # k[6*I-4,3*J-1] += (B₁[i])*N̄[j]*𝑤
                # k[6*I-3,3*J]   += (B₁[i])*N̄[j]*𝑤
                # k[6*I-2,3*J-2] += (B₂[i])*N̄[j]*𝑤
                # k[6*I-1,3*J-1] += (B₂[i])*N̄[j]*𝑤
                # k[6*I,3*J]     += (B₂[i])*N̄[j]*𝑤
            end
        end
    end
end

"""
    ∫∇𝑴𝒏𝒂₃𝒈dΓ_Nitsche(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}, E, ν, h)

Overload or Nitsche/Hellinger–Reissner sibling of the same shell term;
see the documented operators with a similar name.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set∇²𝝭!` — reads the shell fields described in the sibling operator.
"""
function ∫∇𝑴𝒏𝒂₃𝒈dΓ_Nitsche(ap::T, k::AbstractMatrix{Float64}, f::AbstractVector{Float64}; E::Real, ν::Real, h::Real) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        Dᵇ = E*h^3/12/(1-ν^2)
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₁₁ = ξ[:∂²𝝭∂x²]
        B₁₂ = ξ[:∂²𝝭∂x∂y]
        B₂₂ = ξ[:∂²𝝭∂y²]
        B₁₁₁ = ξ[:∂³𝝭∂x³]
        B₁₁₂ = ξ[:∂³𝝭∂x²∂y]
        B₁₂₂ = ξ[:∂³𝝭∂x∂y²]
        B₂₂₂ = ξ[:∂³𝝭∂x³]

        a¹¹ = ξ.a¹¹
        a¹² = ξ.a¹²
        a²² = ξ.a²²
        # b = SymmetricTensor{2,2}((ξ.b₁₁,ξ.b₁₂,ξ.b₂₂))
        𝒂₃₍₁₎ = ξ.𝒂₃₍₁₎
        𝒂₃₍₂₎ = ξ.𝒂₃₍₂₎
        𝒂₃₍₃₎ = ξ.𝒂₃₍₃₎
        # 𝒂ᵃ₍₁₎ = Vec{2}((ξ.𝒂¹₍₁₎,ξ.𝒂²₍₁₎))
        # 𝒂ᵃ₍₂₎ = Vec{2}((ξ.𝒂¹₍₂₎,ξ.𝒂²₍₂₎))
        # 𝒂ᵃ₍₃₎ = Vec{2}((ξ.𝒂¹₍₃₎,ξ.𝒂²₍₃₎))

        # Γᵞᵧₐ = Vec{2}((ξ.Γ¹₁₁+ξ.Γ²₁₂,ξ.Γ¹₁₂+ξ.Γ²₂₂))
        # Γ³₁₂ = Tensor{3,2}((ξ.Γ¹₁₁,ξ.Γ¹₁₂,ξ.Γ¹₁₂,ξ.Γ¹₂₂,ξ.Γ²₁₁,ξ.Γ²₁₂,ξ.Γ²₁₂,ξ.Γ²₂₂))
        # Γ⁴₁₂₃ = Tensor{4,2}((ξ.Γ¹₁₁₁,ξ.Γ¹₁₂₁,ξ.Γ¹₁₂₁,ξ.Γ¹₂₂₁,ξ.Γ¹₁₁₂,ξ.Γ¹₁₂₂,ξ.Γ¹₁₂₂,ξ.Γ¹₂₂₂,ξ.Γ²₁₁₁,ξ.Γ²₁₂₁,ξ.Γ²₁₂₁,ξ.Γ²₂₂₁,ξ.Γ²₁₁₂,ξ.Γ²₁₂₂,ξ.Γ²₁₂₂,ξ.Γ²₂₂₂))

        C¹¹¹¹ = a¹¹*a¹¹*Dᵇ
        C²²²² = a²²*a²²*Dᵇ
        C¹¹²² = (ν*a¹¹*a²² + (1-ν)*a¹²*a¹²)*Dᵇ
        C¹¹¹² = a¹¹*a¹²*Dᵇ
        C²²¹² = a²²*a¹²*Dᵇ
        C¹²¹² = (0.5*((1-ν)*a¹¹*a²² + (1+ν)*a¹²*a¹²))*Dᵇ
        # 𝑪 = SymmetricTensor{4,2}((C¹¹¹¹,C¹¹¹²,C¹¹²²,C¹¹¹²,C¹²¹²,C²²¹²,C¹¹²²,C²²¹²,C²²²²))

        # 𝒏ₐ = Vec{2}((ξ.n₁,ξ.n₂))
        # 𝒔ₐ = Vec{2}((ξ.s₁,ξ.s₂))
        # 𝒔ᵃ = Vec{2}((ξ.s¹,ξ.s²))
        # ∂₂𝒔₁ = Tensor{2,2}((ξ.∂₁s₁,ξ.∂₁s₂,ξ.∂₂s₁,ξ.∂₂s₂))
        # ∂₂𝒏₁ = Tensor{2,2}((ξ.∂₁n₁,ξ.∂₁n₂,ξ.∂₂n₁,ξ.∂₂n₂))

        n₁₁ = ξ.n₁₁
        n₁₂ = ξ.n₁₂
        n₁₃ = ξ.n₁₃
        n₂₂ = ξ.n₂₂
        n₂₃ = ξ.n₂₃
        n₃₃ = ξ.n₃₃
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        g₃ = ξ.g₃
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            # ∇N = Tensor{1,2}((B₁[i],B₂[i]))
            # ∇²N = Tensor{2,2}((B₁₁[i],B₁₂[i],B₁₂[i],B₂₂[i]))
            # ∇³N = Tensor{3,2}((B₁₁₁[i],B₁₁₂[i],B₁₁₂[i],B₁₂₂[i],B₁₁₂[i],B₁₂₂[i],B₁₂₂[i],B₂₂₂[i]))
            # 𝜿₍₁₎ = Tensor{2,2}(((Γ³₁₂⋅∇N - ∇²N).*𝒂₃₍₁₎))
            # 𝜿₍₂₎ = Tensor{2,2}(((Γ³₁₂⋅∇N - ∇²N).*𝒂₃₍₂₎))
            # 𝜿₍₃₎ = Tensor{2,2}(((Γ³₁₂⋅∇N - ∇²N).*𝒂₃₍₃₎))

            # ∇𝜿₍₁₎ = Tensor{3,2}(((Γ⁴₁₂₃[:,:,:,1].*∇N[1]+Γ⁴₁₂₃[:,:,:,2].*∇N[2]).*𝒂₃₍₁₎ - Γ³₁₂⋅∇N⊗(b⋅𝒂ᵃ₍₁₎) + Γ³₁₂⋅∇²N*𝒂₃₍₁₎ + ∇²N⊗(b⋅𝒂ᵃ₍₁₎) - ∇³N.*𝒂₃₍₁₎))
            # ∇𝜿₍₂₎ = Tensor{3,2}(((Γ⁴₁₂₃[:,:,:,1].*∇N[1]+Γ⁴₁₂₃[:,:,:,2].*∇N[2]).*𝒂₃₍₂₎ - Γ³₁₂⋅∇N⊗(b⋅𝒂ᵃ₍₂₎) + Γ³₁₂⋅∇²N*𝒂₃₍₂₎ + ∇²N⊗(b⋅𝒂ᵃ₍₂₎) - ∇³N.*𝒂₃₍₂₎))
            # ∇𝜿₍₃₎ = Tensor{3,2}(((Γ⁴₁₂₃[:,:,:,1].*∇N[1]+Γ⁴₁₂₃[:,:,:,2].*∇N[2]).*𝒂₃₍₃₎ - Γ³₁₂⋅∇N⊗(b⋅𝒂ᵃ₍₃₎) + Γ³₁₂⋅∇²N*𝒂₃₍₃₎ + ∇²N⊗(b⋅𝒂ᵃ₍₃₎) - ∇³N.*𝒂₃₍₃₎))

            # 𝑴₍₁₎ = 𝑪 ⊡ 𝜿₍₁₎
            # 𝑴₍₂₎ = 𝑪 ⊡ 𝜿₍₂₎
            # 𝑴₍₃₎ = 𝑪 ⊡ 𝜿₍₃₎

            # ∇𝑴₍₁₎ = 𝑪 ⊡ ∇𝜿₍₁₎
            # ∇𝑴₍₂₎ = 𝑪 ⊡ ∇𝜿₍₂₎
            # ∇𝑴₍₃₎ = 𝑪 ⊡ ∇𝜿₍₃₎

            # div𝑴₍₁₎ = Vec{2}((∇𝑴₍₁₎[:,1,1]+∇𝑴₍₁₎[:,2,2]))
            # div𝑴₍₂₎ = Vec{2}((∇𝑴₍₂₎[:,1,1]+∇𝑴₍₂₎[:,2,2]))
            # div𝑴₍₃₎ = Vec{2}((∇𝑴₍₃₎[:,1,1]+∇𝑴₍₃₎[:,2,2]))

            # Qᵢ₍₁₁₎ = div𝑴₍₁₎⋅𝒏ₐ*𝒂₃₍₁₎ + 𝒏ₐ⋅𝑴₍₁₎⋅Γᵞᵧₐ*𝒂₃₍₁₎ - 𝒏ₐ⋅𝑴₍₁₎⋅b⋅𝒂ᵃ₍₁₎
            # Qᵢ₍₁₂₎ = div𝑴₍₁₎⋅𝒏ₐ*𝒂₃₍₂₎ + 𝒏ₐ⋅𝑴₍₁₎⋅Γᵞᵧₐ*𝒂₃₍₂₎ - 𝒏ₐ⋅𝑴₍₁₎⋅b⋅𝒂ᵃ₍₂₎
            # Qᵢ₍₁₃₎ = div𝑴₍₁₎⋅𝒏ₐ*𝒂₃₍₃₎ + 𝒏ₐ⋅𝑴₍₁₎⋅Γᵞᵧₐ*𝒂₃₍₂₎ - 𝒏ₐ⋅𝑴₍₁₎⋅b⋅𝒂ᵃ₍₂₎
            # Qᵢ₍₂₁₎ = div𝑴₍₂₎⋅𝒏ₐ*𝒂₃₍₁₎ + 𝒏ₐ⋅𝑴₍₂₎⋅Γᵞᵧₐ*𝒂₃₍₁₎ - 𝒏ₐ⋅𝑴₍₂₎⋅b⋅𝒂ᵃ₍₁₎
            # Qᵢ₍₂₂₎ = div𝑴₍₂₎⋅𝒏ₐ*𝒂₃₍₂₎ + 𝒏ₐ⋅𝑴₍₂₎⋅Γᵞᵧₐ*𝒂₃₍₂₎ - 𝒏ₐ⋅𝑴₍₂₎⋅b⋅𝒂ᵃ₍₂₎
            # Qᵢ₍₂₃₎ = div𝑴₍₂₎⋅𝒏ₐ*𝒂₃₍₃₎ + 𝒏ₐ⋅𝑴₍₂₎⋅Γᵞᵧₐ*𝒂₃₍₂₎ - 𝒏ₐ⋅𝑴₍₂₎⋅b⋅𝒂ᵃ₍₂₎
            # Qᵢ₍₃₁₎ = div𝑴₍₃₎⋅𝒏ₐ*𝒂₃₍₁₎ + 𝒏ₐ⋅𝑴₍₃₎⋅Γᵞᵧₐ*𝒂₃₍₁₎ - 𝒏ₐ⋅𝑴₍₃₎⋅b⋅𝒂ᵃ₍₁₎
            # Qᵢ₍₃₂₎ = div𝑴₍₃₎⋅𝒏ₐ*𝒂₃₍₂₎ + 𝒏ₐ⋅𝑴₍₃₎⋅Γᵞᵧₐ*𝒂₃₍₂₎ - 𝒏ₐ⋅𝑴₍₃₎⋅b⋅𝒂ᵃ₍₂₎
            # Qᵢ₍₃₃₎ = div𝑴₍₃₎⋅𝒏ₐ*𝒂₃₍₃₎ + 𝒏ₐ⋅𝑴₍₃₎⋅Γᵞᵧₐ*𝒂₃₍₂₎ - 𝒏ₐ⋅𝑴₍₃₎⋅b⋅𝒂ᵃ₍₂₎

            # Mᵢ₍₁₁₎ = 𝒏ₐ⋅(𝒔ₐ⋅∇𝑴₍₁₎)⋅𝒔ᵃ*𝒂₃₍₁₎ + 𝒏ₐ⋅𝑴₍₁₎⋅∂₂𝒔₁⋅𝒔ᵃ*𝒂₃₍₁₎ + 𝒔ₐ⋅𝑴₍₁₎⋅∂₂𝒏₁⋅𝒔ᵃ*𝒂₃₍₁₎ - 𝒔ₐ⋅𝑴₍₁₎⋅𝒏ₐ*(𝒔ᵃ⋅b⋅𝒂ᵃ₍₁₎)
            # Mᵢ₍₁₂₎ = 𝒏ₐ⋅(𝒔ₐ⋅∇𝑴₍₁₎)⋅𝒔ᵃ*𝒂₃₍₂₎ + 𝒏ₐ⋅𝑴₍₁₎⋅∂₂𝒔₁⋅𝒔ᵃ*𝒂₃₍₂₎ + 𝒔ₐ⋅𝑴₍₁₎⋅∂₂𝒏₁⋅𝒔ᵃ*𝒂₃₍₂₎ - 𝒔ₐ⋅𝑴₍₁₎⋅𝒏ₐ*(𝒔ᵃ⋅b⋅𝒂ᵃ₍₂₎)
            # Mᵢ₍₁₃₎ = 𝒏ₐ⋅(𝒔ₐ⋅∇𝑴₍₁₎)⋅𝒔ᵃ*𝒂₃₍₃₎ + 𝒏ₐ⋅𝑴₍₁₎⋅∂₂𝒔₁⋅𝒔ᵃ*𝒂₃₍₃₎ + 𝒔ₐ⋅𝑴₍₁₎⋅∂₂𝒏₁⋅𝒔ᵃ*𝒂₃₍₃₎ - 𝒔ₐ⋅𝑴₍₁₎⋅𝒏ₐ*(𝒔ᵃ⋅b⋅𝒂ᵃ₍₃₎)
            # Mᵢ₍₂₁₎ = 𝒏ₐ⋅(𝒔ₐ⋅∇𝑴₍₂₎)⋅𝒔ᵃ*𝒂₃₍₁₎ + 𝒏ₐ⋅𝑴₍₂₎⋅∂₂𝒔₁⋅𝒔ᵃ*𝒂₃₍₁₎ + 𝒔ₐ⋅𝑴₍₂₎⋅∂₂𝒏₁⋅𝒔ᵃ*𝒂₃₍₁₎ - 𝒔ₐ⋅𝑴₍₂₎⋅𝒏ₐ*(𝒔ᵃ⋅b⋅𝒂ᵃ₍₁₎)
            # Mᵢ₍₂₂₎ = 𝒏ₐ⋅(𝒔ₐ⋅∇𝑴₍₂₎)⋅𝒔ᵃ*𝒂₃₍₂₎ + 𝒏ₐ⋅𝑴₍₂₎⋅∂₂𝒔₁⋅𝒔ᵃ*𝒂₃₍₂₎ + 𝒔ₐ⋅𝑴₍₂₎⋅∂₂𝒏₁⋅𝒔ᵃ*𝒂₃₍₂₎ - 𝒔ₐ⋅𝑴₍₂₎⋅𝒏ₐ*(𝒔ᵃ⋅b⋅𝒂ᵃ₍₂₎)
            # Mᵢ₍₂₃₎ = 𝒏ₐ⋅(𝒔ₐ⋅∇𝑴₍₂₎)⋅𝒔ᵃ*𝒂₃₍₃₎ + 𝒏ₐ⋅𝑴₍₂₎⋅∂₂𝒔₁⋅𝒔ᵃ*𝒂₃₍₃₎ + 𝒔ₐ⋅𝑴₍₂₎⋅∂₂𝒏₁⋅𝒔ᵃ*𝒂₃₍₃₎ - 𝒔ₐ⋅𝑴₍₂₎⋅𝒏ₐ*(𝒔ᵃ⋅b⋅𝒂ᵃ₍₃₎)
            # Mᵢ₍₃₁₎ = 𝒏ₐ⋅(𝒔ₐ⋅∇𝑴₍₃₎)⋅𝒔ᵃ*𝒂₃₍₁₎ + 𝒏ₐ⋅𝑴₍₃₎⋅∂₂𝒔₁⋅𝒔ᵃ*𝒂₃₍₁₎ + 𝒔ₐ⋅𝑴₍₃₎⋅∂₂𝒏₁⋅𝒔ᵃ*𝒂₃₍₁₎ - 𝒔ₐ⋅𝑴₍₃₎⋅𝒏ₐ*(𝒔ᵃ⋅b⋅𝒂ᵃ₍₁₎)
            # Mᵢ₍₃₂₎ = 𝒏ₐ⋅(𝒔ₐ⋅∇𝑴₍₃₎)⋅𝒔ᵃ*𝒂₃₍₂₎ + 𝒏ₐ⋅𝑴₍₃₎⋅∂₂𝒔₁⋅𝒔ᵃ*𝒂₃₍₂₎ + 𝒔ₐ⋅𝑴₍₃₎⋅∂₂𝒏₁⋅𝒔ᵃ*𝒂₃₍₂₎ - 𝒔ₐ⋅𝑴₍₃₎⋅𝒏ₐ*(𝒔ᵃ⋅b⋅𝒂ᵃ₍₂₎)
            # Mᵢ₍₃₃₎ = 𝒏ₐ⋅(𝒔ₐ⋅∇𝑴₍₃₎)⋅𝒔ᵃ*𝒂₃₍₃₎ + 𝒏ₐ⋅𝑴₍₃₎⋅∂₂𝒔₁⋅𝒔ᵃ*𝒂₃₍₃₎ + 𝒔ₐ⋅𝑴₍₃₎⋅∂₂𝒏₁⋅𝒔ᵃ*𝒂₃₍₃₎ - 𝒔ₐ⋅𝑴₍₃₎⋅𝒏ₐ*(𝒔ᵃ⋅b⋅𝒂ᵃ₍₃₎)

            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                # ∇N = Tensor{1,2}((B₁[j],B₂[j]))
                # ∇²N = Tensor{2,2}((B₁₁[j],B₁₂[j],B₁₂[j],B₂₂[j]))
                # ∇³N = Tensor{3,2}((B₁₁₁[j],B₁₁₂[j],B₁₁₂[j],B₁₂₂[j],B₁₁₂[j],B₁₂₂[j],B₁₂₂[j],B₂₂₂[j]))
                # 𝜿₍₁₎ = Tensor{2,2}(((Γ³₁₂⋅∇N - ∇²N).*𝒂₃₍₁₎))
                # 𝜿₍₂₎ = Tensor{2,2}(((Γ³₁₂⋅∇N - ∇²N).*𝒂₃₍₂₎))
                # 𝜿₍₃₎ = Tensor{2,2}(((Γ³₁₂⋅∇N - ∇²N).*𝒂₃₍₃₎))

                # ∇𝜿₍₁₎ = Tensor{3,2}(((Γ⁴₁₂₃[:,:,:,1].*∇N[1]+Γ⁴₁₂₃[:,:,:,2].*∇N[2]).*𝒂₃₍₁₎ - Γ³₁₂⋅∇N⊗(b⋅𝒂ᵃ₍₁₎) + Γ³₁₂⋅∇²N*𝒂₃₍₁₎ + ∇²N⊗(b⋅𝒂ᵃ₍₁₎) - ∇³N.*𝒂₃₍₁₎))
                # ∇𝜿₍₂₎ = Tensor{3,2}(((Γ⁴₁₂₃[:,:,:,1].*∇N[1]+Γ⁴₁₂₃[:,:,:,2].*∇N[2]).*𝒂₃₍₂₎ - Γ³₁₂⋅∇N⊗(b⋅𝒂ᵃ₍₂₎) + Γ³₁₂⋅∇²N*𝒂₃₍₂₎ + ∇²N⊗(b⋅𝒂ᵃ₍₂₎) - ∇³N.*𝒂₃₍₂₎))
                # ∇𝜿₍₃₎ = Tensor{3,2}(((Γ⁴₁₂₃[:,:,:,1].*∇N[1]+Γ⁴₁₂₃[:,:,:,2].*∇N[2]).*𝒂₃₍₃₎ - Γ³₁₂⋅∇N⊗(b⋅𝒂ᵃ₍₃₎) + Γ³₁₂⋅∇²N*𝒂₃₍₃₎ + ∇²N⊗(b⋅𝒂ᵃ₍₃₎) - ∇³N.*𝒂₃₍₃₎))

                # 𝑴₍₁₎ = 𝑪 ⊡ 𝜿₍₁₎
                # 𝑴₍₂₎ = 𝑪 ⊡ 𝜿₍₂₎
                # 𝑴₍₃₎ = 𝑪 ⊡ 𝜿₍₃₎

                # ∇𝑴₍₁₎ = 𝑪 ⊡ ∇𝜿₍₁₎
                # ∇𝑴₍₂₎ = 𝑪 ⊡ ∇𝜿₍₂₎
                # ∇𝑴₍₃₎ = 𝑪 ⊡ ∇𝜿₍₃₎

                # div𝑴₍₁₎ = Vec{2}((∇𝑴₍₁₎[:,1,1]+∇𝑴₍₁₎[:,2,2]))
                # div𝑴₍₂₎ = Vec{2}((∇𝑴₍₂₎[:,1,1]+∇𝑴₍₂₎[:,2,2]))
                # div𝑴₍₃₎ = Vec{2}((∇𝑴₍₃₎[:,1,1]+∇𝑴₍₃₎[:,2,2]))

                # Qⱼ₍₁₁₎ = div𝑴₍₁₎⋅𝒏ₐ*𝒂₃₍₁₎ + 𝒏ₐ⋅𝑴₍₁₎⋅Γᵞᵧₐ*𝒂₃₍₁₎ - 𝒏ₐ⋅𝑴₍₁₎⋅b⋅𝒂ᵃ₍₁₎
                # Qⱼ₍₁₂₎ = div𝑴₍₂₎⋅𝒏ₐ*𝒂₃₍₁₎ + 𝒏ₐ⋅𝑴₍₂₎⋅Γᵞᵧₐ*𝒂₃₍₁₎ - 𝒏ₐ⋅𝑴₍₂₎⋅b⋅𝒂ᵃ₍₁₎
                # Qⱼ₍₁₃₎ = div𝑴₍₃₎⋅𝒏ₐ*𝒂₃₍₁₎ + 𝒏ₐ⋅𝑴₍₃₎⋅Γᵞᵧₐ*𝒂₃₍₁₎ - 𝒏ₐ⋅𝑴₍₃₎⋅b⋅𝒂ᵃ₍₁₎
                # Qⱼ₍₂₁₎ = div𝑴₍₁₎⋅𝒏ₐ*𝒂₃₍₂₎ + 𝒏ₐ⋅𝑴₍₁₎⋅Γᵞᵧₐ*𝒂₃₍₂₎ - 𝒏ₐ⋅𝑴₍₁₎⋅b⋅𝒂ᵃ₍₂₎
                # Qⱼ₍₂₂₎ = div𝑴₍₂₎⋅𝒏ₐ*𝒂₃₍₂₎ + 𝒏ₐ⋅𝑴₍₂₎⋅Γᵞᵧₐ*𝒂₃₍₂₎ - 𝒏ₐ⋅𝑴₍₂₎⋅b⋅𝒂ᵃ₍₂₎
                # Qⱼ₍₂₃₎ = div𝑴₍₃₎⋅𝒏ₐ*𝒂₃₍₂₎ + 𝒏ₐ⋅𝑴₍₃₎⋅Γᵞᵧₐ*𝒂₃₍₂₎ - 𝒏ₐ⋅𝑴₍₃₎⋅b⋅𝒂ᵃ₍₂₎
                # Qⱼ₍₃₁₎ = div𝑴₍₁₎⋅𝒏ₐ*𝒂₃₍₃₎ + 𝒏ₐ⋅𝑴₍₁₎⋅Γᵞᵧₐ*𝒂₃₍₃₎ - 𝒏ₐ⋅𝑴₍₁₎⋅b⋅𝒂ᵃ₍₃₎
                # Qⱼ₍₃₂₎ = div𝑴₍₂₎⋅𝒏ₐ*𝒂₃₍₃₎ + 𝒏ₐ⋅𝑴₍₂₎⋅Γᵞᵧₐ*𝒂₃₍₃₎ - 𝒏ₐ⋅𝑴₍₂₎⋅b⋅𝒂ᵃ₍₃₎
                # Qⱼ₍₃₃₎ = div𝑴₍₃₎⋅𝒏ₐ*𝒂₃₍₃₎ + 𝒏ₐ⋅𝑴₍₃₎⋅Γᵞᵧₐ*𝒂₃₍₃₎ - 𝒏ₐ⋅𝑴₍₃₎⋅b⋅𝒂ᵃ₍₃₎

                # Mⱼ₍₁₁₎ = 𝒏ₐ⋅(𝒔ₐ⋅∇𝑴₍₁₎)⋅𝒔ᵃ*𝒂₃₍₁₎ + 𝒏ₐ⋅𝑴₍₁₎⋅∂₂𝒔₁⋅𝒔ᵃ*𝒂₃₍₁₎ + 𝒔ₐ⋅𝑴₍₁₎⋅∂₂𝒏₁⋅𝒔ᵃ*𝒂₃₍₁₎ - 𝒔ₐ⋅𝑴₍₁₎⋅𝒏ₐ*(𝒔ᵃ⋅b⋅𝒂ᵃ₍₁₎)
                # Mⱼ₍₁₂₎ = 𝒏ₐ⋅(𝒔ₐ⋅∇𝑴₍₂₎)⋅𝒔ᵃ*𝒂₃₍₁₎ + 𝒏ₐ⋅𝑴₍₂₎⋅∂₂𝒔₁⋅𝒔ᵃ*𝒂₃₍₁₎ + 𝒔ₐ⋅𝑴₍₂₎⋅∂₂𝒏₁⋅𝒔ᵃ*𝒂₃₍₁₎ - 𝒔ₐ⋅𝑴₍₂₎⋅𝒏ₐ*(𝒔ᵃ⋅b⋅𝒂ᵃ₍₁₎)
                # Mⱼ₍₁₃₎ = 𝒏ₐ⋅(𝒔ₐ⋅∇𝑴₍₃₎)⋅𝒔ᵃ*𝒂₃₍₁₎ + 𝒏ₐ⋅𝑴₍₃₎⋅∂₂𝒔₁⋅𝒔ᵃ*𝒂₃₍₁₎ + 𝒔ₐ⋅𝑴₍₃₎⋅∂₂𝒏₁⋅𝒔ᵃ*𝒂₃₍₁₎ - 𝒔ₐ⋅𝑴₍₃₎⋅𝒏ₐ*(𝒔ᵃ⋅b⋅𝒂ᵃ₍₁₎)
                # Mⱼ₍₂₁₎ = 𝒏ₐ⋅(𝒔ₐ⋅∇𝑴₍₁₎)⋅𝒔ᵃ*𝒂₃₍₂₎ + 𝒏ₐ⋅𝑴₍₁₎⋅∂₂𝒔₁⋅𝒔ᵃ*𝒂₃₍₂₎ + 𝒔ₐ⋅𝑴₍₁₎⋅∂₂𝒏₁⋅𝒔ᵃ*𝒂₃₍₂₎ - 𝒔ₐ⋅𝑴₍₁₎⋅𝒏ₐ*(𝒔ᵃ⋅b⋅𝒂ᵃ₍₂₎)
                # Mⱼ₍₂₂₎ = 𝒏ₐ⋅(𝒔ₐ⋅∇𝑴₍₂₎)⋅𝒔ᵃ*𝒂₃₍₂₎ + 𝒏ₐ⋅𝑴₍₂₎⋅∂₂𝒔₁⋅𝒔ᵃ*𝒂₃₍₂₎ + 𝒔ₐ⋅𝑴₍₂₎⋅∂₂𝒏₁⋅𝒔ᵃ*𝒂₃₍₂₎ - 𝒔ₐ⋅𝑴₍₂₎⋅𝒏ₐ*(𝒔ᵃ⋅b⋅𝒂ᵃ₍₂₎)
                # Mⱼ₍₂₃₎ = 𝒏ₐ⋅(𝒔ₐ⋅∇𝑴₍₃₎)⋅𝒔ᵃ*𝒂₃₍₂₎ + 𝒏ₐ⋅𝑴₍₃₎⋅∂₂𝒔₁⋅𝒔ᵃ*𝒂₃₍₂₎ + 𝒔ₐ⋅𝑴₍₃₎⋅∂₂𝒏₁⋅𝒔ᵃ*𝒂₃₍₂₎ - 𝒔ₐ⋅𝑴₍₃₎⋅𝒏ₐ*(𝒔ᵃ⋅b⋅𝒂ᵃ₍₂₎)
                # Mⱼ₍₃₁₎ = 𝒏ₐ⋅(𝒔ₐ⋅∇𝑴₍₁₎)⋅𝒔ᵃ*𝒂₃₍₃₎ + 𝒏ₐ⋅𝑴₍₁₎⋅∂₂𝒔₁⋅𝒔ᵃ*𝒂₃₍₃₎ + 𝒔ₐ⋅𝑴₍₁₎⋅∂₂𝒏₁⋅𝒔ᵃ*𝒂₃₍₃₎ - 𝒔ₐ⋅𝑴₍₁₎⋅𝒏ₐ*(𝒔ᵃ⋅b⋅𝒂ᵃ₍₃₎)
                # Mⱼ₍₃₂₎ = 𝒏ₐ⋅(𝒔ₐ⋅∇𝑴₍₂₎)⋅𝒔ᵃ*𝒂₃₍₃₎ + 𝒏ₐ⋅𝑴₍₂₎⋅∂₂𝒔₁⋅𝒔ᵃ*𝒂₃₍₃₎ + 𝒔ₐ⋅𝑴₍₂₎⋅∂₂𝒏₁⋅𝒔ᵃ*𝒂₃₍₃₎ - 𝒔ₐ⋅𝑴₍₂₎⋅𝒏ₐ*(𝒔ᵃ⋅b⋅𝒂ᵃ₍₃₎)
                # Mⱼ₍₃₃₎ = 𝒏ₐ⋅(𝒔ₐ⋅∇𝑴₍₃₎)⋅𝒔ᵃ*𝒂₃₍₃₎ + 𝒏ₐ⋅𝑴₍₃₎⋅∂₂𝒔₁⋅𝒔ᵃ*𝒂₃₍₃₎ + 𝒔ₐ⋅𝑴₍₃₎⋅∂₂𝒏₁⋅𝒔ᵃ*𝒂₃₍₃₎ - 𝒔ₐ⋅𝑴₍₃₎⋅𝒏ₐ*(𝒔ᵃ⋅b⋅𝒂ᵃ₍₃₎)

                k[3*I-2,3*J-2] -= (((Qᵢ₍₁₁₎ + Mᵢ₍₁₁₎)*n₁₁+(Qᵢ₍₁₂₎ + Mᵢ₍₁₂₎)*n₁₂+(Qᵢ₍₁₃₎ + Mᵢ₍₁₃₎)*n₁₃)*N[j]
                + N[i]*((Qⱼ₍₁₁₎ + Mⱼ₍₁₁₎)*n₁₁+(Qⱼ₍₂₁₎ + Mⱼ₍₂₁₎)*n₁₂+(Qⱼ₍₃₁₎ + Mⱼ₍₃₁₎)*n₁₃))*𝑤
                k[3*I-1,3*J-2] -= (((Qᵢ₍₂₁₎ + Mᵢ₍₂₁₎)*n₁₁+(Qᵢ₍₂₂₎ + Mᵢ₍₂₂₎)*n₁₂+(Qᵢ₍₂₃₎ + Mᵢ₍₂₃₎)*n₁₃)*N[j]
                + N[i]*((Qⱼ₍₁₁₎ + Mⱼ₍₁₁₎)*n₁₂+(Qⱼ₍₂₁₎ + Mⱼ₍₂₁₎)*n₂₂+(Qⱼ₍₃₁₎ + Mⱼ₍₃₁₎)*n₂₃))*𝑤
                k[3*I,3*J-2]   -= (((Qᵢ₍₃₁₎ + Mᵢ₍₃₁₎)*n₁₁+(Qᵢ₍₃₂₎ + Mᵢ₍₃₂₎)*n₁₂+(Qᵢ₍₃₃₎ + Mᵢ₍₃₃₎)*n₁₃)*N[j]
                + N[i]*((Qⱼ₍₁₁₎ + Mⱼ₍₁₁₎)*n₁₃+(Qⱼ₍₂₁₎ + Mⱼ₍₂₁₎)*n₂₃+(Qⱼ₍₃₁₎ + Mⱼ₍₃₁₎)*n₃₃))*𝑤
                k[3*I-2,3*J-1] -= (((Qᵢ₍₁₁₎ + Mᵢ₍₁₁₎)*n₁₂+(Qᵢ₍₁₂₎ + Mᵢ₍₁₂₎)*n₂₂+(Qᵢ₍₁₃₎ + Mᵢ₍₁₃₎)*n₂₃)*N[j]
                + N[i]*((Qⱼ₍₁₂₎ + Mⱼ₍₁₂₎)*n₁₁+(Qⱼ₍₂₂₎ + Mⱼ₍₂₂₎)*n₁₂+(Qⱼ₍₃₂₎ + Mⱼ₍₃₂₎)*n₁₃))*𝑤
                k[3*I-1,3*J-1] -= (((Qᵢ₍₂₁₎ + Mᵢ₍₂₁₎)*n₁₂+(Qᵢ₍₂₂₎ + Mᵢ₍₂₂₎)*n₂₂+(Qᵢ₍₂₃₎ + Mᵢ₍₂₃₎)*n₂₃)*N[j]
                + N[i]*((Qⱼ₍₁₂₎ + Mⱼ₍₁₂₎)*n₁₂+(Qⱼ₍₂₂₎ + Mⱼ₍₂₂₎)*n₂₂+(Qⱼ₍₃₂₎ + Mⱼ₍₃₂₎)*n₂₃))*𝑤
                k[3*I,3*J-1]   -= (((Qᵢ₍₃₁₎ + Mᵢ₍₃₁₎)*n₁₂+(Qᵢ₍₃₂₎ + Mᵢ₍₃₂₎)*n₂₂+(Qᵢ₍₃₃₎ + Mᵢ₍₃₃₎)*n₂₃)*N[j]
                + N[i]*((Qⱼ₍₁₂₎ + Mⱼ₍₁₂₎)*n₁₃+(Qⱼ₍₂₂₎ + Mⱼ₍₂₂₎)*n₂₃+(Qⱼ₍₃₂₎ + Mⱼ₍₃₂₎)*n₃₃))*𝑤
                k[3*I-2,3*J]   -= (((Qᵢ₍₁₁₎ + Mᵢ₍₁₁₎)*n₁₃+(Qᵢ₍₁₂₎ + Mᵢ₍₁₂₎)*n₂₃+(Qᵢ₍₁₃₎ + Mᵢ₍₁₃₎)*n₃₃)*N[j]
                + N[i]*((Qⱼ₍₁₃₎ + Mⱼ₍₁₃₎)*n₁₁+(Qⱼ₍₂₃₎ + Mⱼ₍₂₃₎)*n₁₂+(Qⱼ₍₃₃₎ + Mⱼ₍₃₃₎)*n₁₃))*𝑤
                k[3*I-1,3*J]   -= (((Qᵢ₍₂₁₎ + Mᵢ₍₂₁₎)*n₁₃+(Qᵢ₍₂₂₎ + Mᵢ₍₂₂₎)*n₂₃+(Qᵢ₍₂₃₎ + Mᵢ₍₂₃₎)*n₃₃)*N[j]
                + N[i]*((Qⱼ₍₁₃₎ + Mⱼ₍₁₃₎)*n₁₂+(Qⱼ₍₂₃₎ + Mⱼ₍₂₃₎)*n₂₂+(Qⱼ₍₃₃₎ + Mⱼ₍₃₃₎)*n₂₃))*𝑤
                k[3*I,3*J]     -= (((Qᵢ₍₃₁₎ + Mᵢ₍₃₁₎)*n₁₃+(Qᵢ₍₃₂₎ + Mᵢ₍₃₂₎)*n₂₃+(Qᵢ₍₃₃₎ + Mᵢ₍₃₃₎)*n₃₃)*N[j]
                + N[i]*((Qⱼ₍₁₃₎ + Mⱼ₍₁₃₎)*n₁₃+(Qⱼ₍₂₃₎ + Mⱼ₍₂₃₎)*n₂₃+(Qⱼ₍₃₃₎ + Mⱼ₍₃₃₎)*n₃₃))*𝑤

            end

            f[3*I-2] -= (((Qᵢ₍₁₁₎ + Mᵢ₍₁₁₎)*n₁₁ + (Qᵢ₍₁₂₎ + Mᵢ₍₁₂₎)*n₁₂ + (Qᵢ₍₁₃₎ + Mᵢ₍₁₃₎)*n₁₃)*g₁
            + ((Qᵢ₍₁₁₎ + Mᵢ₍₁₁₎)*n₁₂ + (Qᵢ₍₁₂₎ + Mᵢ₍₁₂₎)*n₂₂ + (Qᵢ₍₁₃₎ + Mᵢ₍₁₃₎)*n₂₃)*g₂
            + ((Qᵢ₍₁₁₎ + Mᵢ₍₁₁₎)*n₁₃ + (Qᵢ₍₁₂₎ + Mᵢ₍₁₂₎)*n₂₃ + (Qᵢ₍₁₃₎ + Mᵢ₍₁₃₎)*n₃₃)*g₃)*𝑤
            f[3*I-1] -= (((Qᵢ₍₂₁₎ + Mᵢ₍₂₁₎)*n₁₁ + (Qᵢ₍₂₂₎ + Mᵢ₍₂₂₎)*n₁₂ + (Qᵢ₍₂₃₎ + Mᵢ₍₂₃₎)*n₁₃)*g₁
            + ((Qᵢ₍₂₁₎ + Mᵢ₍₂₁₎)*n₁₂ + (Qᵢ₍₂₂₎ + Mᵢ₍₂₂₎)*n₂₂ + (Qᵢ₍₂₃₎ + Mᵢ₍₂₃₎)*n₂₃)*g₂
            + ((Qᵢ₍₂₁₎ + Mᵢ₍₂₁₎)*n₁₃ + (Qᵢ₍₂₂₎ + Mᵢ₍₂₂₎)*n₂₃ + (Qᵢ₍₂₃₎ + Mᵢ₍₂₃₎)*n₃₃)*g₃)*𝑤
            f[3*I]   -= (((Qᵢ₍₃₁₎ + Mᵢ₍₃₁₎)*n₁₁ + (Qᵢ₍₃₂₎ + Mᵢ₍₃₂₎)*n₁₂ + (Qᵢ₍₃₃₎ + Mᵢ₍₃₃₎)*n₁₃)*g₁
            + ((Qᵢ₍₃₁₎ + Mᵢ₍₃₁₎)*n₁₂ + (Qᵢ₍₃₂₎ + Mᵢ₍₃₂₎)*n₂₂ + (Qᵢ₍₃₃₎ + Mᵢ₍₃₃₎)*n₂₃)*g₂
            + ((Qᵢ₍₃₁₎ + Mᵢ₍₃₁₎)*n₁₃ + (Qᵢ₍₃₂₎ + Mᵢ₍₃₂₎)*n₂₃ + (Qᵢ₍₃₃₎ + Mᵢ₍₃₃₎)*n₃₃)*g₃)*𝑤
        end
    end
end

"""
    ∫∇𝑴𝒏𝒂₃𝒈dΓ_HR(a::T,b::S;k::AbstractMatrix{Float64},f::AbstractVector{Float64})

Overload or Nitsche/Hellinger–Reissner sibling of the same shell term;
see the documented operators with a similar name.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set∇²𝝭!` — reads the shell fields described in the sibling operator.
"""
function ∫∇𝑴𝒏𝒂₃𝒈dΓ_HR(a::T,b::S, k::AbstractMatrix{Float64}, f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵃ = a.𝓒;𝓖ᵃ = a.𝓖
    𝓒ᵇ = b.𝓒;𝓖ᵇ = b.𝓖
    for (ξᵃ,ξᵇ) in zip(𝓖ᵃ,𝓖ᵇ)
        𝑤 = ξᵃ.𝑤
        N = ξᵃ[:𝝭]
        B₁ = ξᵃ[:∂𝝭∂x]
        B₂ = ξᵃ[:∂𝝭∂y]
        N̄ = ξᵇ[:𝝭]
        n₁ = ξᵃ.n₁
        n₂ = ξᵃ.n₂
        s₁ = ξᵃ.s₁
        s₂ = ξᵃ.s₂
        s¹ = ξᵃ.s¹
        s² = ξᵃ.s²
        ∂₁n₁ = ξᵃ.∂₁n₁
        ∂₂n₁ = ξᵃ.∂₂n₁
        ∂₁n₂ = ξᵃ.∂₁n₂
        ∂₂n₂ = ξᵃ.∂₂n₂
        ∂₁s₁ = ξᵃ.∂₁s₁
        ∂₂s₁ = ξᵃ.∂₂s₁
        ∂₁s₂ = ξᵃ.∂₁s₂
        ∂₂s₂ = ξᵃ.∂₂s₂
        Γ¹₁₁ = ξᵃ.Γ¹₁₁
        Γ¹₁₂ = ξᵃ.Γ¹₁₂
        Γ¹₂₂ = ξᵃ.Γ¹₂₂
        Γ²₁₁ = ξᵃ.Γ²₁₁
        Γ²₁₂ = ξᵃ.Γ²₁₂
        Γ²₂₂ = ξᵃ.Γ²₂₂
        Γᵞᵧ₁n₁ = Γ¹₁₁*n₁ + Γ²₁₂*n₁
        Γᵞᵧ₂n₂ = Γ¹₁₂*n₂ + Γ²₂₂*n₂
        Γᵞᵧ₂n₁ = Γ¹₁₂*n₁ + Γ²₂₂*n₁
        Γᵞᵧ₁n₂ = Γ¹₁₁*n₂ + Γ²₁₂*n₂
        Γᵞ₁₁nᵧ = Γ¹₁₁*n₁ + Γ²₁₁*n₂
        Γᵞ₂₂nᵧ = Γ¹₂₂*n₁ + Γ²₂₂*n₂
        Γᵞ₁₂nᵧ = Γ¹₁₂*n₁ + Γ²₁₂*n₂
        s₁n₁s¹ = s₁*n₁*s¹
        s₁n₁s² = s₁*n₁*s²
        s₂n₂s¹ = s₂*n₂*s¹
        s₂n₂s² = s₂*n₂*s²
        s₁n₂s¹ = s₁*n₂*s¹
        s₁n₂s² = s₁*n₂*s²
        s₂n₁s¹ = s₂*n₁*s¹
        s₂n₁s² = s₂*n₁*s²
        ∂ᵧs₁n₁sᵞ = ∂₁s₁*n₁*s¹ + ∂₂s₁*n₁*s²
        ∂ᵧs₂n₂sᵞ = ∂₁s₂*n₂*s¹ + ∂₂s₂*n₂*s²
        ∂ᵧs₁n₂sᵞ = ∂₁s₁*n₂*s¹ + ∂₂s₁*n₂*s²
        ∂ᵧs₂n₁sᵞ = ∂₁s₂*n₁*s¹ + ∂₂s₂*n₁*s²
        s₁∂ᵧn₁sᵞ = s₁*∂₁n₁*s¹ + s₁*∂₂n₁*s²
        s₂∂ᵧn₂sᵞ = s₂*∂₁n₂*s¹ + s₂*∂₂n₂*s²
        s₁∂ᵧn₂sᵞ = s₁*∂₁n₂*s¹ + s₁*∂₂n₂*s²
        s₂∂ᵧn₁sᵞ = s₂*∂₁n₁*s¹ + s₂*∂₂n₁*s²

        n₁₁ = ξᵇ.n₁₁
        n₁₂ = ξᵇ.n₁₂
        n₁₃ = ξᵇ.n₁₃
        n₂₂ = ξᵇ.n₂₂
        n₂₃ = ξᵇ.n₂₃
        n₃₃ = ξᵇ.n₃₃

        g₁ = ξᵇ.g₁
        g₂ = ξᵇ.g₂
        g₃ = ξᵇ.g₃
        for (i,xᵢ) in enumerate(𝓒ᵃ)
            I = xᵢ.𝐼
            B̃₁₁ = B₁[i]*s₁n₁s¹ + B₂[i]*s₁n₁s² + N[i]*(∂ᵧs₁n₁sᵞ + s₁∂ᵧn₁sᵞ) + B₁[i]*n₁ + N[i]*Γᵞ₁₁nᵧ + N[i]*Γᵞᵧ₁n₁
            B̃₂₂ = B₁[i]*s₂n₂s¹ + B₂[i]*s₂n₂s² + N[i]*(∂ᵧs₂n₂sᵞ + s₂∂ᵧn₂sᵞ) + B₂[i]*n₂ + N[i]*Γᵞ₂₂nᵧ + N[i]*Γᵞᵧ₂n₂
            B̃₁₂ = B₁[i]*(s₁n₂s¹ +s₂n₁s¹) + B₂[i]*(s₁n₂s² + s₂n₁s²) + N[i]*(∂ᵧs₁n₂sᵞ + ∂ᵧs₂n₁sᵞ + s₁∂ᵧn₂sᵞ + s₂∂ᵧn₁sᵞ) + B₁[i]*n₂ + B₂[i]*n₁ + 2*N[i]*Γᵞ₁₂nᵧ + N[i]*(Γᵞᵧ₁n₂ + Γᵞᵧ₂n₁)
            for (j,xⱼ) in enumerate(𝓒ᵇ)
                J = xⱼ.𝐼
                k[9*I-8,3*J-2] += B̃₁₁*n₁₁*N̄[j]*𝑤
                k[9*I-7,3*J-2] += B̃₁₁*n₁₂*N̄[j]*𝑤
                k[9*I-6,3*J-2] += B̃₁₁*n₁₃*N̄[j]*𝑤
                k[9*I-5,3*J-2] += B̃₂₂*n₁₁*N̄[j]*𝑤
                k[9*I-4,3*J-2] += B̃₂₂*n₁₂*N̄[j]*𝑤
                k[9*I-3,3*J-2] += B̃₂₂*n₁₃*N̄[j]*𝑤
                k[9*I-2,3*J-2] += B̃₁₂*n₁₁*N̄[j]*𝑤
                k[9*I-1,3*J-2] += B̃₁₂*n₁₂*N̄[j]*𝑤
                k[9*I,3*J-2]   += B̃₁₂*n₁₃*N̄[j]*𝑤
                k[9*I-8,3*J-1] += B̃₁₁*n₁₂*N̄[j]*𝑤
                k[9*I-7,3*J-1] += B̃₁₁*n₂₂*N̄[j]*𝑤
                k[9*I-6,3*J-1] += B̃₁₁*n₂₃*N̄[j]*𝑤
                k[9*I-5,3*J-1] += B̃₂₂*n₁₂*N̄[j]*𝑤
                k[9*I-4,3*J-1] += B̃₂₂*n₂₂*N̄[j]*𝑤
                k[9*I-3,3*J-1] += B̃₂₂*n₂₃*N̄[j]*𝑤
                k[9*I-2,3*J-1] += B̃₁₂*n₁₂*N̄[j]*𝑤
                k[9*I-1,3*J-1] += B̃₁₂*n₂₂*N̄[j]*𝑤
                k[9*I,3*J-1]   += B̃₁₂*n₂₃*N̄[j]*𝑤
                k[9*I-8,3*J]   += B̃₁₁*n₁₃*N̄[j]*𝑤
                k[9*I-7,3*J]   += B̃₁₁*n₂₃*N̄[j]*𝑤
                k[9*I-6,3*J]   += B̃₁₁*n₃₃*N̄[j]*𝑤
                k[9*I-5,3*J]   += B̃₂₂*n₁₃*N̄[j]*𝑤
                k[9*I-4,3*J]   += B̃₂₂*n₂₃*N̄[j]*𝑤
                k[9*I-3,3*J]   += B̃₂₂*n₃₃*N̄[j]*𝑤
                k[9*I-2,3*J]   += B̃₁₂*n₁₃*N̄[j]*𝑤
                k[9*I-1,3*J]   += B̃₁₂*n₂₃*N̄[j]*𝑤
                k[9*I,3*J]     += B̃₁₂*n₃₃*N̄[j]*𝑤
            end
            f[9*I-8] += B̃₁₁*(n₁₁*g₁+n₁₂*g₂+n₁₃*g₃)*𝑤
            f[9*I-7] += B̃₁₁*(n₁₂*g₁+n₂₂*g₂+n₂₃*g₃)*𝑤
            f[9*I-6] += B̃₁₁*(n₁₃*g₁+n₂₃*g₂+n₃₃*g₃)*𝑤
            f[9*I-5] += B̃₂₂*(n₁₁*g₁+n₁₂*g₂+n₁₃*g₃)*𝑤
            f[9*I-4] += B̃₂₂*(n₁₂*g₁+n₂₂*g₂+n₂₃*g₃)*𝑤
            f[9*I-3] += B̃₂₂*(n₁₃*g₁+n₂₃*g₂+n₃₃*g₃)*𝑤
            f[9*I-2] += B̃₁₂*(n₁₁*g₁+n₁₂*g₂+n₁₃*g₃)*𝑤
            f[9*I-1] += B̃₁₂*(n₁₂*g₁+n₂₂*g₂+n₂₃*g₃)*𝑤
            f[9*I]   += B̃₁₂*(n₁₃*g₁+n₂₃*g₂+n₃₃*g₃)*𝑤
        end
    end
end

"""
    ∫∇𝑴𝒏𝒂₃𝒗dΓ(a::T,b::S;k::AbstractMatrix{Float64})

Overload or Nitsche/Hellinger–Reissner sibling of the same shell term;
see the documented operators with a similar name.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set∇²𝝭!` — reads the shell fields described in the sibling operator.
"""
function ∫∇𝑴𝒏𝒂₃𝒗dΓ(a::T,b::S, k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵃ = a.𝓒;𝓖ᵃ = a.𝓖
    𝓒ᵇ = b.𝓒;𝓖ᵇ = b.𝓖
    for (ξᵃ,ξᵇ) in zip(𝓖ᵃ,𝓖ᵇ)
        𝑤 = ξᵃ.𝑤
        N = ξᵃ[:𝝭]
        B₁ = ξᵃ[:∂𝝭∂x]
        B₂ = ξᵃ[:∂𝝭∂y]
        N̄ = ξᵇ[:𝝭]
        n₁ = ξᵃ.n₁
        n₂ = ξᵃ.n₂
        s₁ = ξᵃ.s₁
        s₂ = ξᵃ.s₂
        s¹ = ξᵃ.s¹
        s² = ξᵃ.s²
        ∂₁n₁ = ξᵃ.∂₁n₁
        ∂₂n₁ = ξᵃ.∂₂n₁
        ∂₁n₂ = ξᵃ.∂₁n₂
        ∂₂n₂ = ξᵃ.∂₂n₂
        ∂₁s₁ = ξᵃ.∂₁s₁
        ∂₂s₁ = ξᵃ.∂₂s₁
        ∂₁s₂ = ξᵃ.∂₁s₂
        ∂₂s₂ = ξᵃ.∂₂s₂
        Γ¹₁₁ = ξᵃ.Γ¹₁₁
        Γ¹₁₂ = ξᵃ.Γ¹₁₂
        Γ¹₂₂ = ξᵃ.Γ¹₂₂
        Γ²₁₁ = ξᵃ.Γ²₁₁
        Γ²₁₂ = ξᵃ.Γ²₁₂
        Γ²₂₂ = ξᵃ.Γ²₂₂
        Γᵞᵧ₁n₁ = Γ¹₁₁*n₁ + Γ²₁₂*n₁
        Γᵞᵧ₂n₂ = Γ¹₁₂*n₂ + Γ²₂₂*n₂
        Γᵞᵧ₂n₁ = Γ¹₁₂*n₁ + Γ²₂₂*n₁
        Γᵞᵧ₁n₂ = Γ¹₁₁*n₂ + Γ²₁₂*n₂
        Γᵞ₁₁nᵧ = Γ¹₁₁*n₁ + Γ²₁₁*n₂
        Γᵞ₂₂nᵧ = Γ¹₂₂*n₁ + Γ²₂₂*n₂
        Γᵞ₁₂nᵧ = Γ¹₁₂*n₁ + Γ²₁₂*n₂
        s₁n₁s¹ = s₁*n₁*s¹
        s₁n₁s² = s₁*n₁*s²
        s₂n₂s¹ = s₂*n₂*s¹
        s₂n₂s² = s₂*n₂*s²
        s₁n₂s¹ = s₁*n₂*s¹
        s₁n₂s² = s₁*n₂*s²
        s₂n₁s¹ = s₂*n₁*s¹
        s₂n₁s² = s₂*n₁*s²
        ∂ᵧs₁n₁sᵞ = ∂₁s₁*n₁*s¹ + ∂₂s₁*n₁*s²
        ∂ᵧs₂n₂sᵞ = ∂₁s₂*n₂*s¹ + ∂₂s₂*n₂*s²
        ∂ᵧs₁n₂sᵞ = ∂₁s₁*n₂*s¹ + ∂₂s₁*n₂*s²
        ∂ᵧs₂n₁sᵞ = ∂₁s₂*n₁*s¹ + ∂₂s₂*n₁*s²
        s₁∂ᵧn₁sᵞ = s₁*∂₁n₁*s¹ + s₁*∂₂n₁*s²
        s₂∂ᵧn₂sᵞ = s₂*∂₁n₂*s¹ + s₂*∂₂n₂*s²
        s₁∂ᵧn₂sᵞ = s₁*∂₁n₂*s¹ + s₁*∂₂n₂*s²
        s₂∂ᵧn₁sᵞ = s₂*∂₁n₁*s¹ + s₂*∂₂n₁*s²


        for (i,xᵢ) in enumerate(𝓒ᵃ)
            I = xᵢ.𝐼
            # B̃₁₁ = 0.
            # B̃₂₂ = 0.
            B̃₁₁ = B₁[i]*s₁n₁s¹ + B₂[i]*s₁n₁s² + N[i]*(∂ᵧs₁n₁sᵞ + s₁∂ᵧn₁sᵞ) + B₁[i]*n₁ + N[i]*Γᵞ₁₁nᵧ + N[i]*Γᵞᵧ₁n₁
            B̃₂₂ = B₁[i]*s₂n₂s¹ + B₂[i]*s₂n₂s² + N[i]*(∂ᵧs₂n₂sᵞ + s₂∂ᵧn₂sᵞ) + B₂[i]*n₂ + N[i]*Γᵞ₂₂nᵧ + N[i]*Γᵞᵧ₂n₂
            B̃₁₂ = B₁[i]*(s₁n₂s¹ +s₂n₁s¹) + B₂[i]*(s₁n₂s² + s₂n₁s²) + N[i]*(∂ᵧs₁n₂sᵞ + ∂ᵧs₂n₁sᵞ + s₁∂ᵧn₂sᵞ + s₂∂ᵧn₁sᵞ) + B₁[i]*n₂ + B₂[i]*n₁ + 2*N[i]*Γᵞ₁₂nᵧ + N[i]*(Γᵞᵧ₁n₂ + Γᵞᵧ₂n₁)

            for (j,xⱼ) in enumerate(𝓒ᵇ)
                J = xⱼ.𝐼
                k[9*I-8,3*J-2] -= B̃₁₁*N̄[j]*𝑤
                k[9*I-7,3*J-1] -= B̃₁₁*N̄[j]*𝑤
                k[9*I-6,3*J]   -= B̃₁₁*N̄[j]*𝑤
                k[9*I-5,3*J-2] -= B̃₂₂*N̄[j]*𝑤
                k[9*I-4,3*J-1] -= B̃₂₂*N̄[j]*𝑤
                k[9*I-3,3*J]   -= B̃₂₂*N̄[j]*𝑤
                k[9*I-2,3*J-2] -= B̃₁₂*N̄[j]*𝑤
                k[9*I-1,3*J-1] -= B̃₁₂*N̄[j]*𝑤
                k[9*I,3*J]     -= B̃₁₂*N̄[j]*𝑤
            end
        end
    end
end

"""
    ∫MₙₙθₙdΓ_Nitsche(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}, E, ν, h)

Overload or Nitsche/Hellinger–Reissner sibling of the same shell term;
see the documented operators with a similar name.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set∇²𝝭!` — reads the shell fields described in the sibling operator.
"""
function ∫MₙₙθₙdΓ_Nitsche(ap::T, k::AbstractMatrix{Float64}, f::AbstractVector{Float64}; E::Real, ν::Real, h::Real) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        Dᵇ = E*h^3/12/(1-ν^2)
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₁₁ = ξ[:∂²𝝭∂x²]
        B₁₂ = ξ[:∂²𝝭∂x∂y]
        B₂₂ = ξ[:∂²𝝭∂y²]

        n₁ = ξ.n₁
        n₂ = ξ.n₂
        n¹ = ξ.n¹
        n² = ξ.n²
        a¹¹ = ξ.a¹¹
        a¹² = ξ.a¹²
        a²² = ξ.a²²

        Γ¹₁₁ = ξ.Γ¹₁₁
        Γ¹₁₂ = ξ.Γ¹₁₂
        Γ¹₂₂ = ξ.Γ¹₂₂
        Γ²₁₁ = ξ.Γ²₁₁
        Γ²₁₂ = ξ.Γ²₁₂
        Γ²₂₂ = ξ.Γ²₂₂

        𝒂₃₍₁₎ = ξ.𝒂₃₍₁₎
        𝒂₃₍₂₎ = ξ.𝒂₃₍₂₎
        𝒂₃₍₃₎ = ξ.𝒂₃₍₃₎

        C¹¹¹¹ = a¹¹*a¹¹*Dᵇ
        C²²²² = a²²*a²²*Dᵇ
        C¹¹²² = (ν*a¹¹*a²² + (1-ν)*a¹²*a¹²)*Dᵇ
        C¹¹¹² = a¹¹*a¹²*Dᵇ
        C²²¹² = a²²*a¹²*Dᵇ
        C¹²¹² = (0.5*((1-ν)*a¹¹*a²² + (1+ν)*a¹²*a¹²))*Dᵇ

        θ = ξ.θ
        M11 = 0.
        M22 = 0.
        M12 = 0.
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼

            κ₁₁₍₁₎ =   (Γ¹₁₁*B₁[i] + Γ²₁₁*B₂[i] - B₁₁[i])*𝒂₃₍₁₎
            κ₁₁₍₂₎ =   (Γ¹₁₁*B₁[i] + Γ²₁₁*B₂[i] - B₁₁[i])*𝒂₃₍₂₎
            κ₁₁₍₃₎ =   (Γ¹₁₁*B₁[i] + Γ²₁₁*B₂[i] - B₁₁[i])*𝒂₃₍₃₎
            κ₂₂₍₁₎ =   (Γ¹₂₂*B₁[i] + Γ²₂₂*B₂[i] - B₂₂[i])*𝒂₃₍₁₎
            κ₂₂₍₂₎ =   (Γ¹₂₂*B₁[i] + Γ²₂₂*B₂[i] - B₂₂[i])*𝒂₃₍₂₎
            κ₂₂₍₃₎ =   (Γ¹₂₂*B₁[i] + Γ²₂₂*B₂[i] - B₂₂[i])*𝒂₃₍₃₎
            κ₁₂₍₁₎ = 2*(Γ¹₁₂*B₁[i] + Γ²₁₂*B₂[i] - B₁₂[i])*𝒂₃₍₁₎
            κ₁₂₍₂₎ = 2*(Γ¹₁₂*B₁[i] + Γ²₁₂*B₂[i] - B₁₂[i])*𝒂₃₍₂₎
            κ₁₂₍₃₎ = 2*(Γ¹₁₂*B₁[i] + Γ²₁₂*B₂[i] - B₁₂[i])*𝒂₃₍₃₎

            M¹¹₍₁₎ = C¹¹¹¹*κ₁₁₍₁₎ + C¹¹²²*κ₂₂₍₁₎ + C¹¹¹²*κ₁₂₍₁₎
            M¹¹₍₂₎ = C¹¹¹¹*κ₁₁₍₂₎ + C¹¹²²*κ₂₂₍₂₎ + C¹¹¹²*κ₁₂₍₂₎
            M¹¹₍₃₎ = C¹¹¹¹*κ₁₁₍₃₎ + C¹¹²²*κ₂₂₍₃₎ + C¹¹¹²*κ₁₂₍₃₎
            M²²₍₁₎ = C¹¹²²*κ₁₁₍₁₎ + C²²²²*κ₂₂₍₁₎ + C²²¹²*κ₁₂₍₁₎
            M²²₍₂₎ = C¹¹²²*κ₁₁₍₂₎ + C²²²²*κ₂₂₍₂₎ + C²²¹²*κ₁₂₍₂₎
            M²²₍₃₎ = C¹¹²²*κ₁₁₍₃₎ + C²²²²*κ₂₂₍₃₎ + C²²¹²*κ₁₂₍₃₎
            M¹²₍₁₎ = C¹¹¹²*κ₁₁₍₁₎ + C²²¹²*κ₂₂₍₁₎ + C¹²¹²*κ₁₂₍₁₎
            M¹²₍₂₎ = C¹¹¹²*κ₁₁₍₂₎ + C²²¹²*κ₂₂₍₂₎ + C¹²¹²*κ₁₂₍₂₎
            M¹²₍₃₎ = C¹¹¹²*κ₁₁₍₃₎ + C²²¹²*κ₂₂₍₃₎ + C¹²¹²*κ₁₂₍₃₎

            M₍₁₎ᵢ = M¹¹₍₁₎*n₁*n₁ + M²²₍₁₎*n₂*n₂ + 2*M¹²₍₁₎*n₁*n₂
            M₍₂₎ᵢ = M¹¹₍₂₎*n₁*n₁ + M²²₍₂₎*n₂*n₂ + 2*M¹²₍₂₎*n₁*n₂
            M₍₃₎ᵢ = M¹¹₍₃₎*n₁*n₁ + M²²₍₃₎*n₂*n₂ + 2*M¹²₍₃₎*n₁*n₂

            θ₍₁₎ᵢ = (B₁[i]*n¹+B₂[i]*n²)*𝒂₃₍₁₎
            θ₍₂₎ᵢ = (B₁[i]*n¹+B₂[i]*n²)*𝒂₃₍₂₎
            θ₍₃₎ᵢ = (B₁[i]*n¹+B₂[i]*n²)*𝒂₃₍₃₎

            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼

                κ₁₁₍₁₎ =   (Γ¹₁₁*B₁[j] + Γ²₁₁*B₂[j] - B₁₁[j])*𝒂₃₍₁₎
                κ₁₁₍₂₎ =   (Γ¹₁₁*B₁[j] + Γ²₁₁*B₂[j] - B₁₁[j])*𝒂₃₍₂₎
                κ₁₁₍₃₎ =   (Γ¹₁₁*B₁[j] + Γ²₁₁*B₂[j] - B₁₁[j])*𝒂₃₍₃₎
                κ₂₂₍₁₎ =   (Γ¹₂₂*B₁[j] + Γ²₂₂*B₂[j] - B₂₂[j])*𝒂₃₍₁₎
                κ₂₂₍₂₎ =   (Γ¹₂₂*B₁[j] + Γ²₂₂*B₂[j] - B₂₂[j])*𝒂₃₍₂₎
                κ₂₂₍₃₎ =   (Γ¹₂₂*B₁[j] + Γ²₂₂*B₂[j] - B₂₂[j])*𝒂₃₍₃₎
                κ₁₂₍₁₎ = 2*(Γ¹₁₂*B₁[j] + Γ²₁₂*B₂[j] - B₁₂[j])*𝒂₃₍₁₎
                κ₁₂₍₂₎ = 2*(Γ¹₁₂*B₁[j] + Γ²₁₂*B₂[j] - B₁₂[j])*𝒂₃₍₂₎
                κ₁₂₍₃₎ = 2*(Γ¹₁₂*B₁[j] + Γ²₁₂*B₂[j] - B₁₂[j])*𝒂₃₍₃₎

                M¹¹₍₁₎ = C¹¹¹¹*κ₁₁₍₁₎ + C¹¹²²*κ₂₂₍₁₎ + C¹¹¹²*κ₁₂₍₁₎
                M¹¹₍₂₎ = C¹¹¹¹*κ₁₁₍₂₎ + C¹¹²²*κ₂₂₍₂₎ + C¹¹¹²*κ₁₂₍₂₎
                M¹¹₍₃₎ = C¹¹¹¹*κ₁₁₍₃₎ + C¹¹²²*κ₂₂₍₃₎ + C¹¹¹²*κ₁₂₍₃₎
                M²²₍₁₎ = C¹¹²²*κ₁₁₍₁₎ + C²²²²*κ₂₂₍₁₎ + C²²¹²*κ₁₂₍₁₎
                M²²₍₂₎ = C¹¹²²*κ₁₁₍₂₎ + C²²²²*κ₂₂₍₂₎ + C²²¹²*κ₁₂₍₂₎
                M²²₍₃₎ = C¹¹²²*κ₁₁₍₃₎ + C²²²²*κ₂₂₍₃₎ + C²²¹²*κ₁₂₍₃₎
                M¹²₍₁₎ = C¹¹¹²*κ₁₁₍₁₎ + C²²¹²*κ₂₂₍₁₎ + C¹²¹²*κ₁₂₍₁₎
                M¹²₍₂₎ = C¹¹¹²*κ₁₁₍₂₎ + C²²¹²*κ₂₂₍₂₎ + C¹²¹²*κ₁₂₍₂₎
                M¹²₍₃₎ = C¹¹¹²*κ₁₁₍₃₎ + C²²¹²*κ₂₂₍₃₎ + C¹²¹²*κ₁₂₍₃₎

                M₍₁₎ⱼ = M¹¹₍₁₎*n₁*n₁ + M²²₍₁₎*n₂*n₂ + 2*M¹²₍₁₎*n₁*n₂
                M₍₂₎ⱼ = M¹¹₍₂₎*n₁*n₁ + M²²₍₂₎*n₂*n₂ + 2*M¹²₍₂₎*n₁*n₂
                M₍₃₎ⱼ = M¹¹₍₃₎*n₁*n₁ + M²²₍₃₎*n₂*n₂ + 2*M¹²₍₃₎*n₁*n₂

                θ₍₁₎ⱼ = (B₁[j]*n¹+B₂[j]*n²)*𝒂₃₍₁₎
                θ₍₂₎ⱼ = (B₁[j]*n¹+B₂[j]*n²)*𝒂₃₍₂₎
                θ₍₃₎ⱼ = (B₁[j]*n¹+B₂[j]*n²)*𝒂₃₍₃₎

                k[3*I-2,3*J-2] += (M₍₁₎ᵢ*θ₍₁₎ⱼ + θ₍₁₎ᵢ*M₍₁₎ⱼ)*𝑤
                k[3*I-1,3*J-2] += (M₍₂₎ᵢ*θ₍₁₎ⱼ + θ₍₂₎ᵢ*M₍₁₎ⱼ)*𝑤
                k[3*I,3*J-2]   += (M₍₃₎ᵢ*θ₍₁₎ⱼ + θ₍₃₎ᵢ*M₍₁₎ⱼ)*𝑤
                k[3*I-2,3*J-1] += (M₍₁₎ᵢ*θ₍₂₎ⱼ + θ₍₁₎ᵢ*M₍₂₎ⱼ)*𝑤
                k[3*I-1,3*J-1] += (M₍₂₎ᵢ*θ₍₂₎ⱼ + θ₍₂₎ᵢ*M₍₂₎ⱼ)*𝑤
                k[3*I,3*J-1]   += (M₍₃₎ᵢ*θ₍₂₎ⱼ + θ₍₃₎ᵢ*M₍₂₎ⱼ)*𝑤
                k[3*I-2,3*J]   += (M₍₁₎ᵢ*θ₍₃₎ⱼ + θ₍₁₎ᵢ*M₍₃₎ⱼ)*𝑤
                k[3*I-1,3*J]   += (M₍₂₎ᵢ*θ₍₃₎ⱼ + θ₍₂₎ᵢ*M₍₃₎ⱼ)*𝑤
                k[3*I,3*J]     += (M₍₃₎ᵢ*θ₍₃₎ⱼ + θ₍₃₎ᵢ*M₍₃₎ⱼ)*𝑤

            end
            f[3*I-2] += M₍₁₎ᵢ*θ*𝑤
            f[3*I-1] += M₍₂₎ᵢ*θ*𝑤
            f[3*I]   += M₍₃₎ᵢ*θ*𝑤
        end
    end
end

"""
    ∫MₙₙθₙdΓ_HR(a::T,b::S;k::AbstractMatrix{Float64},f::AbstractVector{Float64})

Overload or Nitsche/Hellinger–Reissner sibling of the same shell term;
see the documented operators with a similar name.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set∇²𝝭!` — reads the shell fields described in the sibling operator.
"""
function ∫MₙₙθₙdΓ_HR(a::T,b::S, k::AbstractMatrix{Float64}, f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵃ = a.𝓒;𝓖ᵃ = a.𝓖
    𝓒ᵇ = b.𝓒;𝓖ᵇ = b.𝓖
    for (ξᵃ,ξᵇ) in zip(𝓖ᵃ,𝓖ᵇ)
        𝑤 = ξᵃ.𝑤
        N = ξᵃ[:𝝭]
        B̄₁ = ξᵇ[:∂𝝭∂x]
        B̄₂ = ξᵇ[:∂𝝭∂y]
        n₁ = ξᵃ.n₁
        n₂ = ξᵃ.n₂
        n¹ = ξᵃ.n¹
        n² = ξᵃ.n²
        𝒂₃₍₁₎ = ξᵃ.𝒂₃₍₁₎
        𝒂₃₍₂₎ = ξᵃ.𝒂₃₍₂₎
        𝒂₃₍₃₎ = ξᵃ.𝒂₃₍₃₎
        ∂ₙu₁ = ξᵇ.∂ₙu₁
        ∂ₙu₂ = ξᵇ.∂ₙu₂
        ∂ₙu₃ = ξᵇ.∂ₙu₃

        θ = ξᵇ.θ

        for (i,xᵢ) in enumerate(𝓒ᵃ)
            I = xᵢ.𝐼
            N₁₁ = N[i]*n₁*n₁
            N₂₂ = N[i]*n₂*n₂
            N₁₂ = N[i]*2*n₁*n₂
            for (j,xⱼ) in enumerate(𝓒ᵇ)
                J = xⱼ.𝐼
                B̄ₙ = n¹*B̄₁[j] + n²*B̄₂[j]
                # k[9*I-8,3*J-2] -= N₁₁*𝒂₃₍₁₎*𝒂₃₍₁₎*B̄ₙ*𝑤
                # k[9*I-7,3*J-2] -= N₁₁*𝒂₃₍₂₎*𝒂₃₍₁₎*B̄ₙ*𝑤
                # k[9*I-6,3*J-2] -= N₁₁*𝒂₃₍₃₎*𝒂₃₍₁₎*B̄ₙ*𝑤
                # k[9*I-5,3*J-2] -= N₂₂*𝒂₃₍₁₎*𝒂₃₍₁₎*B̄ₙ*𝑤
                # k[9*I-4,3*J-2] -= N₂₂*𝒂₃₍₂₎*𝒂₃₍₁₎*B̄ₙ*𝑤
                # k[9*I-3,3*J-2] -= N₂₂*𝒂₃₍₃₎*𝒂₃₍₁₎*B̄ₙ*𝑤
                # k[9*I-2,3*J-2] -= N₁₂*𝒂₃₍₁₎*𝒂₃₍₁₎*B̄ₙ*𝑤
                # k[9*I-1,3*J-2] -= N₁₂*𝒂₃₍₂₎*𝒂₃₍₁₎*B̄ₙ*𝑤
                # k[9*I,3*J-2]   -= N₁₂*𝒂₃₍₃₎*𝒂₃₍₁₎*B̄ₙ*𝑤
                # k[9*I-8,3*J-1] -= N₁₁*𝒂₃₍₁₎*𝒂₃₍₂₎*B̄ₙ*𝑤
                # k[9*I-7,3*J-1] -= N₁₁*𝒂₃₍₂₎*𝒂₃₍₂₎*B̄ₙ*𝑤
                # k[9*I-6,3*J-1] -= N₁₁*𝒂₃₍₃₎*𝒂₃₍₂₎*B̄ₙ*𝑤
                # k[9*I-5,3*J-1] -= N₂₂*𝒂₃₍₁₎*𝒂₃₍₂₎*B̄ₙ*𝑤
                # k[9*I-4,3*J-1] -= N₂₂*𝒂₃₍₂₎*𝒂₃₍₂₎*B̄ₙ*𝑤
                # k[9*I-3,3*J-1] -= N₂₂*𝒂₃₍₃₎*𝒂₃₍₂₎*B̄ₙ*𝑤
                # k[9*I-2,3*J-1] -= N₁₂*𝒂₃₍₁₎*𝒂₃₍₂₎*B̄ₙ*𝑤
                # k[9*I-1,3*J-1] -= N₁₂*𝒂₃₍₂₎*𝒂₃₍₂₎*B̄ₙ*𝑤
                # k[9*I,3*J-1]   -= N₁₂*𝒂₃₍₃₎*𝒂₃₍₂₎*B̄ₙ*𝑤
                # k[9*I-8,3*J]   -= N₁₁*𝒂₃₍₁₎*𝒂₃₍₃₎*B̄ₙ*𝑤
                # k[9*I-7,3*J]   -= N₁₁*𝒂₃₍₂₎*𝒂₃₍₃₎*B̄ₙ*𝑤
                # k[9*I-6,3*J]   -= N₁₁*𝒂₃₍₃₎*𝒂₃₍₃₎*B̄ₙ*𝑤
                # k[9*I-5,3*J]   -= N₂₂*𝒂₃₍₁₎*𝒂₃₍₃₎*B̄ₙ*𝑤
                # k[9*I-4,3*J]   -= N₂₂*𝒂₃₍₂₎*𝒂₃₍₃₎*B̄ₙ*𝑤
                # k[9*I-3,3*J]   -= N₂₂*𝒂₃₍₃₎*𝒂₃₍₃₎*B̄ₙ*𝑤
                # k[9*I-2,3*J]   -= N₁₂*𝒂₃₍₁₎*𝒂₃₍₃₎*B̄ₙ*𝑤
                # k[9*I-1,3*J]   -= N₁₂*𝒂₃₍₂₎*𝒂₃₍₃₎*B̄ₙ*𝑤
                # k[9*I,3*J]     -= N₁₂*𝒂₃₍₃₎*𝒂₃₍₃₎*B̄ₙ*𝑤
                k[9*I-8,3*J-2] -= N₁₁*B̄ₙ*𝑤
                k[9*I-7,3*J-1] -= N₁₁*B̄ₙ*𝑤
                k[9*I-6,3*J]   -= N₁₁*B̄ₙ*𝑤
                k[9*I-5,3*J-2] -= N₂₂*B̄ₙ*𝑤
                k[9*I-4,3*J-1] -= N₂₂*B̄ₙ*𝑤
                k[9*I-3,3*J]   -= N₂₂*B̄ₙ*𝑤
                k[9*I-2,3*J-2] -= N₁₂*B̄ₙ*𝑤
                k[9*I-1,3*J-1] -= N₁₂*B̄ₙ*𝑤
                k[9*I,3*J]     -= N₁₂*B̄ₙ*𝑤
            end
            f[9*I-8] -= N₁₁*𝒂₃₍₁₎*θ*𝑤
            f[9*I-7] -= N₁₁*𝒂₃₍₂₎*θ*𝑤
            f[9*I-6] -= N₁₁*𝒂₃₍₃₎*θ*𝑤
            f[9*I-5] -= N₂₂*𝒂₃₍₁₎*θ*𝑤
            f[9*I-4] -= N₂₂*𝒂₃₍₂₎*θ*𝑤
            f[9*I-3] -= N₂₂*𝒂₃₍₃₎*θ*𝑤
            f[9*I-2] -= N₁₂*𝒂₃₍₁₎*θ*𝑤
            f[9*I-1] -= N₁₂*𝒂₃₍₂₎*θ*𝑤
            f[9*I]   -= N₁₂*𝒂₃₍₃₎*θ*𝑤
            # f[9*I-8] -= N₁₁*∂ₙu₁*𝑤
            # f[9*I-7] -= N₁₁*∂ₙu₂*𝑤
            # f[9*I-6] -= N₁₁*∂ₙu₃*𝑤
            # f[9*I-5] -= N₂₂*∂ₙu₁*𝑤
            # f[9*I-4] -= N₂₂*∂ₙu₂*𝑤
            # f[9*I-3] -= N₂₂*∂ₙu₃*𝑤
            # f[9*I-2] -= N₁₂*∂ₙu₁*𝑤
            # f[9*I-1] -= N₁₂*∂ₙu₂*𝑤
            # f[9*I]   -= N₁₂*∂ₙu₃*𝑤
        end
    end
end

"""
    ∫𝑴ₙₙ𝜽ₙdΓ(a::T,b::S;k::AbstractMatrix{Float64})

Overload or Nitsche/Hellinger–Reissner sibling of the same shell term;
see the documented operators with a similar name.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set∇²𝝭!` — reads the shell fields described in the sibling operator.
"""
function ∫𝑴ₙₙ𝜽ₙdΓ(a::T,b::S, k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵃ = a.𝓒;𝓖ᵃ = a.𝓖
    𝓒ᵇ = b.𝓒;𝓖ᵇ = b.𝓖
    for (ξᵃ,ξᵇ) in zip(𝓖ᵃ,𝓖ᵇ)
        𝑤 = ξᵃ.𝑤
        N = ξᵃ[:𝝭]
        N̄ = ξᵇ[:𝝭]
        B̄₁ = ξᵇ[:∂𝝭∂x]
        B̄₂ = ξᵇ[:∂𝝭∂y]
        n₁ = ξᵃ.n₁
        n₂ = ξᵃ.n₂
        n¹ = ξᵃ.n¹
        n² = ξᵃ.n²
        𝒂₃₍₁₎ = ξᵃ.𝒂₃₍₁₎
        𝒂₃₍₂₎ = ξᵃ.𝒂₃₍₂₎
        𝒂₃₍₃₎ = ξᵃ.𝒂₃₍₃₎
        for (i,xᵢ) in enumerate(𝓒ᵃ)
            I = xᵢ.𝐼
            N₁₁ = N[i]*n₁*n₁
            N₂₂ = N[i]*n₂*n₂
            N₁₂ = N[i]*2*n₁*n₂
            for (j,xⱼ) in enumerate(𝓒ᵇ)
                J = xⱼ.𝐼
                B̄ₙ = n¹*B̄₁[j] + n²*B̄₂[j]
                # k[9*I-8,3*J-2] += N₁₁*𝒂₃₍₁₎*𝒂₃₍₁₎*B̄ₙ*𝑤
                # k[9*I-7,3*J-2] += N₁₁*𝒂₃₍₂₎*𝒂₃₍₁₎*B̄ₙ*𝑤
                # k[9*I-6,3*J-2] += N₁₁*𝒂₃₍₃₎*𝒂₃₍₁₎*B̄ₙ*𝑤
                # k[9*I-5,3*J-2] += N₂₂*𝒂₃₍₁₎*𝒂₃₍₁₎*B̄ₙ*𝑤
                # k[9*I-4,3*J-2] += N₂₂*𝒂₃₍₂₎*𝒂₃₍₁₎*B̄ₙ*𝑤
                # k[9*I-3,3*J-2] += N₂₂*𝒂₃₍₃₎*𝒂₃₍₁₎*B̄ₙ*𝑤
                # k[9*I-2,3*J-2] += N₁₂*𝒂₃₍₁₎*𝒂₃₍₁₎*B̄ₙ*𝑤
                # k[9*I-1,3*J-2] += N₁₂*𝒂₃₍₂₎*𝒂₃₍₁₎*B̄ₙ*𝑤
                # k[9*I,3*J-2]   += N₁₂*𝒂₃₍₃₎*𝒂₃₍₁₎*B̄ₙ*𝑤
                # k[9*I-8,3*J-1] += N₁₁*𝒂₃₍₁₎*𝒂₃₍₂₎*B̄ₙ*𝑤
                # k[9*I-7,3*J-1] += N₁₁*𝒂₃₍₂₎*𝒂₃₍₂₎*B̄ₙ*𝑤
                # k[9*I-6,3*J-1] += N₁₁*𝒂₃₍₃₎*𝒂₃₍₂₎*B̄ₙ*𝑤
                # k[9*I-5,3*J-1] += N₂₂*𝒂₃₍₁₎*𝒂₃₍₂₎*B̄ₙ*𝑤
                # k[9*I-4,3*J-1] += N₂₂*𝒂₃₍₂₎*𝒂₃₍₂₎*B̄ₙ*𝑤
                # k[9*I-3,3*J-1] += N₂₂*𝒂₃₍₃₎*𝒂₃₍₂₎*B̄ₙ*𝑤
                # k[9*I-2,3*J-1] += N₁₂*𝒂₃₍₁₎*𝒂₃₍₂₎*B̄ₙ*𝑤
                # k[9*I-1,3*J-1] += N₁₂*𝒂₃₍₂₎*𝒂₃₍₂₎*B̄ₙ*𝑤
                # k[9*I,3*J-1]   += N₁₂*𝒂₃₍₃₎*𝒂₃₍₂₎*B̄ₙ*𝑤
                # k[9*I-8,3*J]   += N₁₁*𝒂₃₍₁₎*𝒂₃₍₃₎*B̄ₙ*𝑤
                # k[9*I-7,3*J]   += N₁₁*𝒂₃₍₂₎*𝒂₃₍₃₎*B̄ₙ*𝑤
                # k[9*I-6,3*J]   += N₁₁*𝒂₃₍₃₎*𝒂₃₍₃₎*B̄ₙ*𝑤
                # k[9*I-5,3*J]   += N₂₂*𝒂₃₍₁₎*𝒂₃₍₃₎*B̄ₙ*𝑤
                # k[9*I-4,3*J]   += N₂₂*𝒂₃₍₂₎*𝒂₃₍₃₎*B̄ₙ*𝑤
                # k[9*I-3,3*J]   += N₂₂*𝒂₃₍₃₎*𝒂₃₍₃₎*B̄ₙ*𝑤
                # k[9*I-2,3*J]   += N₁₂*𝒂₃₍₁₎*𝒂₃₍₃₎*B̄ₙ*𝑤
                # k[9*I-1,3*J]   += N₁₂*𝒂₃₍₂₎*𝒂₃₍₃₎*B̄ₙ*𝑤
                # k[9*I,3*J]     += N₁₂*𝒂₃₍₃₎*𝒂₃₍₃₎*B̄ₙ*𝑤
                k[9*I-8,3*J-2] += N₁₁*B̄ₙ*𝑤
                k[9*I-7,3*J-1] += N₁₁*B̄ₙ*𝑤
                k[9*I-6,3*J]   += N₁₁*B̄ₙ*𝑤
                k[9*I-5,3*J-2] += N₂₂*B̄ₙ*𝑤
                k[9*I-4,3*J-1] += N₂₂*B̄ₙ*𝑤
                k[9*I-3,3*J]   += N₂₂*B̄ₙ*𝑤
                k[9*I-2,3*J-2] += N₁₂*B̄ₙ*𝑤
                k[9*I-1,3*J-1] += N₁₂*B̄ₙ*𝑤
                k[9*I,3*J]     += N₁₂*B̄ₙ*𝑤
            end
        end
    end
end
"""
    ΔMₙₛ𝒂₃𝒈_Nitsche(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}, E, ν, h)

Overload or Nitsche/Hellinger–Reissner sibling of the same shell term;
see the documented operators with a similar name.


    ∫ (ΔMₙₛ 𝒂₃ · 𝒈) dΓ  — corner force in the Nitsche form.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set∇²𝝭!` — reads the shell fields described in the sibling operator.
"""
function ΔMₙₛ𝒂₃𝒈_Nitsche(ap::T, k::AbstractMatrix{Float64}, f::AbstractVector{Float64}; E::Real, ν::Real, h::Real) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        Dᵇ = E*h^3/12/(1-ν^2)
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₁₁ = ξ[:∂²𝝭∂x²]
        B₁₂ = ξ[:∂²𝝭∂x∂y]
        B₂₂ = ξ[:∂²𝝭∂y²]

        n₁ = ξ.n₁
        n₂ = ξ.n₂
        s₁ = ξ.s₁
        s₂ = ξ.s₂
        a¹¹ = ξ.a¹¹
        a¹² = ξ.a¹²
        a²² = ξ.a²²
        Δ = ξ.Δ

        Γ¹₁₁ = ξ.Γ¹₁₁
        Γ¹₁₂ = ξ.Γ¹₁₂
        Γ¹₂₂ = ξ.Γ¹₂₂
        Γ²₁₁ = ξ.Γ²₁₁
        Γ²₁₂ = ξ.Γ²₁₂
        Γ²₂₂ = ξ.Γ²₂₂

        𝒂₃₍₁₎ = ξ.𝒂₃₍₁₎
        𝒂₃₍₂₎ = ξ.𝒂₃₍₂₎
        𝒂₃₍₃₎ = ξ.𝒂₃₍₃₎

        C¹¹¹¹ = a¹¹*a¹¹*Dᵇ
        C²²²² = a²²*a²²*Dᵇ
        C¹¹²² = (ν*a¹¹*a²² + (1-ν)*a¹²*a¹²)*Dᵇ
        C¹¹¹² = a¹¹*a¹²*Dᵇ
        C²²¹² = a²²*a¹²*Dᵇ
        C¹²¹² = (0.5*((1-ν)*a¹¹*a²² + (1+ν)*a¹²*a¹²))*Dᵇ

        n₁₁ = ξ.n₁₁
        n₁₂ = ξ.n₁₂
        n₁₃ = ξ.n₁₃
        n₂₂ = ξ.n₂₂
        n₂₃ = ξ.n₂₃
        n₃₃ = ξ.n₃₃
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        g₃ = ξ.g₃

        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼

            κ₁₁₍₁₎ =   (Γ¹₁₁*B₁[i] + Γ²₁₁*B₂[i] - B₁₁[i])*𝒂₃₍₁₎
            κ₁₁₍₂₎ =   (Γ¹₁₁*B₁[i] + Γ²₁₁*B₂[i] - B₁₁[i])*𝒂₃₍₂₎
            κ₁₁₍₃₎ =   (Γ¹₁₁*B₁[i] + Γ²₁₁*B₂[i] - B₁₁[i])*𝒂₃₍₃₎
            κ₂₂₍₁₎ =   (Γ¹₂₂*B₁[i] + Γ²₂₂*B₂[i] - B₂₂[i])*𝒂₃₍₁₎
            κ₂₂₍₂₎ =   (Γ¹₂₂*B₁[i] + Γ²₂₂*B₂[i] - B₂₂[i])*𝒂₃₍₂₎
            κ₂₂₍₃₎ =   (Γ¹₂₂*B₁[i] + Γ²₂₂*B₂[i] - B₂₂[i])*𝒂₃₍₃₎
            κ₁₂₍₁₎ = 2*(Γ¹₁₂*B₁[i] + Γ²₁₂*B₂[i] - B₁₂[i])*𝒂₃₍₁₎
            κ₁₂₍₂₎ = 2*(Γ¹₁₂*B₁[i] + Γ²₁₂*B₂[i] - B₁₂[i])*𝒂₃₍₂₎
            κ₁₂₍₃₎ = 2*(Γ¹₁₂*B₁[i] + Γ²₁₂*B₂[i] - B₁₂[i])*𝒂₃₍₃₎

            M¹¹₍₁₎ = C¹¹¹¹*κ₁₁₍₁₎ + C¹¹²²*κ₂₂₍₁₎ + C¹¹¹²*κ₁₂₍₁₎
            M¹¹₍₂₎ = C¹¹¹¹*κ₁₁₍₂₎ + C¹¹²²*κ₂₂₍₂₎ + C¹¹¹²*κ₁₂₍₂₎
            M¹¹₍₃₎ = C¹¹¹¹*κ₁₁₍₃₎ + C¹¹²²*κ₂₂₍₃₎ + C¹¹¹²*κ₁₂₍₃₎
            M²²₍₁₎ = C¹¹²²*κ₁₁₍₁₎ + C²²²²*κ₂₂₍₁₎ + C²²¹²*κ₁₂₍₁₎
            M²²₍₂₎ = C¹¹²²*κ₁₁₍₂₎ + C²²²²*κ₂₂₍₂₎ + C²²¹²*κ₁₂₍₂₎
            M²²₍₃₎ = C¹¹²²*κ₁₁₍₃₎ + C²²²²*κ₂₂₍₃₎ + C²²¹²*κ₁₂₍₃₎
            M¹²₍₁₎ = C¹¹¹²*κ₁₁₍₁₎ + C²²¹²*κ₂₂₍₁₎ + C¹²¹²*κ₁₂₍₁₎
            M¹²₍₂₎ = C¹¹¹²*κ₁₁₍₂₎ + C²²¹²*κ₂₂₍₂₎ + C¹²¹²*κ₁₂₍₂₎
            M¹²₍₃₎ = C¹¹¹²*κ₁₁₍₃₎ + C²²¹²*κ₂₂₍₃₎ + C¹²¹²*κ₁₂₍₃₎

            M₍₁₎ᵢ = M¹¹₍₁₎*s₁*n₁ + M²²₍₁₎*s₂*n₂ + M¹²₍₁₎*(s₁*n₂+s₂*n₁)
            M₍₂₎ᵢ = M¹¹₍₂₎*s₁*n₁ + M²²₍₂₎*s₂*n₂ + M¹²₍₂₎*(s₁*n₂+s₂*n₁)
            M₍₃₎ᵢ = M¹¹₍₃₎*s₁*n₁ + M²²₍₃₎*s₂*n₂ + M¹²₍₃₎*(s₁*n₂+s₂*n₁)

            N₍₁₎ᵢ = N[i]*𝒂₃₍₁₎
            N₍₂₎ᵢ = N[i]*𝒂₃₍₂₎
            N₍₃₎ᵢ = N[i]*𝒂₃₍₃₎

            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼

                κ₁₁₍₁₎ =   (Γ¹₁₁*B₁[j] + Γ²₁₁*B₂[j] - B₁₁[j])*𝒂₃₍₁₎
                κ₁₁₍₂₎ =   (Γ¹₁₁*B₁[j] + Γ²₁₁*B₂[j] - B₁₁[j])*𝒂₃₍₂₎
                κ₁₁₍₃₎ =   (Γ¹₁₁*B₁[j] + Γ²₁₁*B₂[j] - B₁₁[j])*𝒂₃₍₃₎
                κ₂₂₍₁₎ =   (Γ¹₂₂*B₁[j] + Γ²₂₂*B₂[j] - B₂₂[j])*𝒂₃₍₁₎
                κ₂₂₍₂₎ =   (Γ¹₂₂*B₁[j] + Γ²₂₂*B₂[j] - B₂₂[j])*𝒂₃₍₂₎
                κ₂₂₍₃₎ =   (Γ¹₂₂*B₁[j] + Γ²₂₂*B₂[j] - B₂₂[j])*𝒂₃₍₃₎
                κ₁₂₍₁₎ = 2*(Γ¹₁₂*B₁[j] + Γ²₁₂*B₂[j] - B₁₂[j])*𝒂₃₍₁₎
                κ₁₂₍₂₎ = 2*(Γ¹₁₂*B₁[j] + Γ²₁₂*B₂[j] - B₁₂[j])*𝒂₃₍₂₎
                κ₁₂₍₃₎ = 2*(Γ¹₁₂*B₁[j] + Γ²₁₂*B₂[j] - B₁₂[j])*𝒂₃₍₃₎

                M¹¹₍₁₎ = C¹¹¹¹*κ₁₁₍₁₎ + C¹¹²²*κ₂₂₍₁₎ + C¹¹¹²*κ₁₂₍₁₎
                M¹¹₍₂₎ = C¹¹¹¹*κ₁₁₍₂₎ + C¹¹²²*κ₂₂₍₂₎ + C¹¹¹²*κ₁₂₍₂₎
                M¹¹₍₃₎ = C¹¹¹¹*κ₁₁₍₃₎ + C¹¹²²*κ₂₂₍₃₎ + C¹¹¹²*κ₁₂₍₃₎
                M²²₍₁₎ = C¹¹²²*κ₁₁₍₁₎ + C²²²²*κ₂₂₍₁₎ + C²²¹²*κ₁₂₍₁₎
                M²²₍₂₎ = C¹¹²²*κ₁₁₍₂₎ + C²²²²*κ₂₂₍₂₎ + C²²¹²*κ₁₂₍₂₎
                M²²₍₃₎ = C¹¹²²*κ₁₁₍₃₎ + C²²²²*κ₂₂₍₃₎ + C²²¹²*κ₁₂₍₃₎
                M¹²₍₁₎ = C¹¹¹²*κ₁₁₍₁₎ + C²²¹²*κ₂₂₍₁₎ + C¹²¹²*κ₁₂₍₁₎
                M¹²₍₂₎ = C¹¹¹²*κ₁₁₍₂₎ + C²²¹²*κ₂₂₍₂₎ + C¹²¹²*κ₁₂₍₂₎
                M¹²₍₃₎ = C¹¹¹²*κ₁₁₍₃₎ + C²²¹²*κ₂₂₍₃₎ + C¹²¹²*κ₁₂₍₃₎

                M₍₁₎ⱼ = M¹¹₍₁₎*s₁*n₁ + M²²₍₁₎*s₂*n₂ + M¹²₍₁₎*(s₁*n₂+s₂*n₁)
                M₍₂₎ⱼ = M¹¹₍₂₎*s₁*n₁ + M²²₍₂₎*s₂*n₂ + M¹²₍₂₎*(s₁*n₂+s₂*n₁)
                M₍₃₎ⱼ = M¹¹₍₃₎*s₁*n₁ + M²²₍₃₎*s₂*n₂ + M¹²₍₃₎*(s₁*n₂+s₂*n₁)

                N₍₁₎ⱼ = N[j]*𝒂₃₍₁₎
                N₍₂₎ⱼ = N[j]*𝒂₃₍₂₎
                N₍₃₎ⱼ = N[j]*𝒂₃₍₃₎

                k[3*I-2,3*J-2] += (M₍₁₎ᵢ*N₍₁₎ⱼ + N₍₁₎ᵢ*M₍₁₎ⱼ)*Δ
                k[3*I-1,3*J-2] += (M₍₂₎ᵢ*N₍₁₎ⱼ + N₍₂₎ᵢ*M₍₁₎ⱼ)*Δ
                k[3*I,3*J-2]   += (M₍₃₎ᵢ*N₍₁₎ⱼ + N₍₃₎ᵢ*M₍₁₎ⱼ)*Δ
                k[3*I-2,3*J-1] += (M₍₁₎ᵢ*N₍₂₎ⱼ + N₍₁₎ᵢ*M₍₂₎ⱼ)*Δ
                k[3*I-1,3*J-1] += (M₍₂₎ᵢ*N₍₂₎ⱼ + N₍₂₎ᵢ*M₍₂₎ⱼ)*Δ
                k[3*I,3*J-1]   += (M₍₃₎ᵢ*N₍₂₎ⱼ + N₍₃₎ᵢ*M₍₂₎ⱼ)*Δ
                k[3*I-2,3*J]   += (M₍₁₎ᵢ*N₍₃₎ⱼ + N₍₁₎ᵢ*M₍₃₎ⱼ)*Δ
                k[3*I-1,3*J]   += (M₍₂₎ᵢ*N₍₃₎ⱼ + N₍₂₎ᵢ*M₍₃₎ⱼ)*Δ
                k[3*I,3*J]     += (M₍₃₎ᵢ*N₍₃₎ⱼ + N₍₃₎ᵢ*M₍₃₎ⱼ)*Δ
            end
            f[3*I-2] += M₍₁₎ᵢ*(𝒂₃₍₁₎*g₁+𝒂₃₍₂₎*g₂+𝒂₃₍₃₎*g₃)*Δ
            f[3*I-1] += M₍₂₎ᵢ*(𝒂₃₍₁₎*g₁+𝒂₃₍₂₎*g₂+𝒂₃₍₃₎*g₃)*Δ
            f[3*I]   += M₍₃₎ᵢ*(𝒂₃₍₁₎*g₁+𝒂₃₍₂₎*g₂+𝒂₃₍₃₎*g₃)*Δ
        end
    end
end

"""
    ΔMₙₛ𝒂₃𝒈_HR(a::T,b::S;k::AbstractMatrix{Float64},f::AbstractVector{Float64})

Overload or Nitsche/Hellinger–Reissner sibling of the same shell term;
see the documented operators with a similar name.


    ∫ (ΔMₙₛ 𝒂₃ · 𝒈) dΓ  — corner force in the Hellinger–Reissner form.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set∇²𝝭!` — reads the shell fields described in the sibling operator.
"""
function ΔMₙₛ𝒂₃𝒈_HR(a::T,b::S, k::AbstractMatrix{Float64}, f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵃ = a.𝓒;𝓖ᵃ = a.𝓖
    𝓒ᵇ = b.𝓒;𝓖ᵇ = b.𝓖
    for (ξᵃ,ξᵇ) in zip(𝓖ᵃ,𝓖ᵇ)
        𝑤 = ξᵃ.𝑤
        N = ξᵃ[:𝝭]
        N̄ = ξᵇ[:𝝭]
        n₁ = ξᵃ.n₁
        n₂ = ξᵃ.n₂
        s₁ = ξᵃ.s₁
        s₂ = ξᵃ.s₂
        Δ = ξᵃ.Δ

        n₁₁ = ξᵇ.n₁₁
        n₁₂ = ξᵇ.n₁₂
        n₁₃ = ξᵇ.n₁₃
        n₂₂ = ξᵇ.n₂₂
        n₂₃ = ξᵇ.n₂₃
        n₃₃ = ξᵇ.n₃₃

        g₁ = ξᵇ.g₁
        g₂ = ξᵇ.g₂
        g₃ = ξᵇ.g₃
        for (i,xᵢ) in enumerate(𝓒ᵃ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵇ)
                J = xⱼ.𝐼
                k[9*I-8,3*J-2] -= Δ*N[i]*n₁₁*N̄[j]*s₁*n₁
                k[9*I-7,3*J-2] -= Δ*N[i]*n₁₂*N̄[j]*s₁*n₁
                k[9*I-6,3*J-2] -= Δ*N[i]*n₁₃*N̄[j]*s₁*n₁
                k[9*I-5,3*J-2] -= Δ*N[i]*n₁₁*N̄[j]*s₂*n₂
                k[9*I-4,3*J-2] -= Δ*N[i]*n₁₂*N̄[j]*s₂*n₂
                k[9*I-3,3*J-2] -= Δ*N[i]*n₁₃*N̄[j]*s₂*n₂
                k[9*I-2,3*J-2] -= Δ*N[i]*n₁₁*N̄[j]*(s₁*n₂+s₂*n₁)
                k[9*I-1,3*J-2] -= Δ*N[i]*n₁₂*N̄[j]*(s₁*n₂+s₂*n₁)
                k[9*I,3*J-2]   -= Δ*N[i]*n₁₃*N̄[j]*(s₁*n₂+s₂*n₁)
                k[9*I-8,3*J-1] -= Δ*N[i]*n₁₂*N̄[j]*s₁*n₁
                k[9*I-7,3*J-1] -= Δ*N[i]*n₂₂*N̄[j]*s₁*n₁
                k[9*I-6,3*J-1] -= Δ*N[i]*n₂₃*N̄[j]*s₁*n₁
                k[9*I-5,3*J-1] -= Δ*N[i]*n₁₂*N̄[j]*s₂*n₂
                k[9*I-4,3*J-1] -= Δ*N[i]*n₂₂*N̄[j]*s₂*n₂
                k[9*I-3,3*J-1] -= Δ*N[i]*n₂₃*N̄[j]*s₂*n₂
                k[9*I-2,3*J-1] -= Δ*N[i]*n₁₂*N̄[j]*(s₁*n₂+s₂*n₁)
                k[9*I-1,3*J-1] -= Δ*N[i]*n₂₂*N̄[j]*(s₁*n₂+s₂*n₁)
                k[9*I,3*J-1]   -= Δ*N[i]*n₂₃*N̄[j]*(s₁*n₂+s₂*n₁)
                k[9*I-8,3*J]   -= Δ*N[i]*n₁₃*N̄[j]*s₁*n₁
                k[9*I-7,3*J]   -= Δ*N[i]*n₂₃*N̄[j]*s₁*n₁
                k[9*I-6,3*J]   -= Δ*N[i]*n₃₃*N̄[j]*s₁*n₁
                k[9*I-5,3*J]   -= Δ*N[i]*n₁₃*N̄[j]*s₂*n₂
                k[9*I-4,3*J]   -= Δ*N[i]*n₂₃*N̄[j]*s₂*n₂
                k[9*I-3,3*J]   -= Δ*N[i]*n₃₃*N̄[j]*s₂*n₂
                k[9*I-2,3*J]   -= Δ*N[i]*n₁₃*N̄[j]*(s₁*n₂+s₂*n₁)
                k[9*I-1,3*J]   -= Δ*N[i]*n₂₃*N̄[j]*(s₁*n₂+s₂*n₁)
                k[9*I,3*J]     -= Δ*N[i]*n₃₃*N̄[j]*(s₁*n₂+s₂*n₁)
            end
            f[9*I-8] -= Δ*N[i]*(n₁₁*g₁+n₁₂*g₂+n₁₃*g₃)*s₁*n₁
            f[9*I-7] -= Δ*N[i]*(n₁₂*g₁+n₂₂*g₂+n₂₃*g₃)*s₁*n₁
            f[9*I-6] -= Δ*N[i]*(n₁₃*g₁+n₂₃*g₂+n₃₃*g₃)*s₁*n₁
            f[9*I-5] -= Δ*N[i]*(n₁₁*g₁+n₁₂*g₂+n₁₃*g₃)*s₂*n₂
            f[9*I-4] -= Δ*N[i]*(n₁₂*g₁+n₂₂*g₂+n₂₃*g₃)*s₂*n₂
            f[9*I-3] -= Δ*N[i]*(n₁₃*g₁+n₂₃*g₂+n₃₃*g₃)*s₂*n₂
            f[9*I-2] -= Δ*N[i]*(n₁₁*g₁+n₁₂*g₂+n₁₃*g₃)*(s₁*n₂+s₂*n₁)
            f[9*I-1] -= Δ*N[i]*(n₁₂*g₁+n₂₂*g₂+n₂₃*g₃)*(s₁*n₂+s₂*n₁)
            f[9*I]   -= Δ*N[i]*(n₁₃*g₁+n₂₃*g₂+n₃₃*g₃)*(s₁*n₂+s₂*n₁)
        end
    end
end

"""
    ΔMₙₛ𝒂₃𝒗(a::T,b::S;k::AbstractMatrix{Float64})

Overload or Nitsche/Hellinger–Reissner sibling of the same shell term;
see the documented operators with a similar name.


    ∫ (ΔMₙₛ 𝒂₃ · 𝒗) dΓ  — the trial-side counterpart of `ΔMₙₛ𝒂₃𝒈`, used to
    assemble the left-hand side of the corner condition.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set∇²𝝭!` — reads the shell fields described in the sibling operator.
"""
function ΔMₙₛ𝒂₃𝒗(a::T,b::S, k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵃ = a.𝓒;𝓖ᵃ = a.𝓖
    𝓒ᵇ = b.𝓒;𝓖ᵇ = b.𝓖
    for (ξᵃ,ξᵇ) in zip(𝓖ᵃ,𝓖ᵇ)
        𝑤 = ξᵃ.𝑤
        N = ξᵃ[:𝝭]
        N̄ = ξᵇ[:𝝭]
        n₁ = ξᵃ.n₁
        n₂ = ξᵃ.n₂
        s₁ = ξᵃ.s₁
        s₂ = ξᵃ.s₂
        Δ = ξᵃ.Δ
        for (i,xᵢ) in enumerate(𝓒ᵃ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵇ)
                J = xⱼ.𝐼
                k[9*I-8,3*J-2] += Δ*N[i]*N̄[j]*s₁*n₁
                k[9*I-7,3*J-1] += Δ*N[i]*N̄[j]*s₁*n₁
                k[9*I-6,3*J]   += Δ*N[i]*N̄[j]*s₁*n₁
                k[9*I-5,3*J-2] += Δ*N[i]*N̄[j]*s₂*n₂
                k[9*I-4,3*J-1] += Δ*N[i]*N̄[j]*s₂*n₂
                k[9*I-3,3*J]   += Δ*N[i]*N̄[j]*s₂*n₂
                k[9*I-2,3*J-2] += Δ*N[i]*N̄[j]*(s₁*n₂+s₂*n₁)
                k[9*I-1,3*J-1] += Δ*N[i]*N̄[j]*(s₁*n₂+s₂*n₁)
                k[9*I,3*J]     += Δ*N[i]*N̄[j]*(s₁*n₂+s₂*n₁)
            end
        end
    end
end

"""
    ∫∇𝑴∇𝒂₃𝒗dΩ(a::T,b::S;k::AbstractMatrix{Float64})

Bending-moment/deflection coupling in gradient form.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set∇²𝝭!` — reads the second derivatives and the surface metric.
"""
function ∫∇𝑴∇𝒂₃𝒗dΩ(a::T,b::S, k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵃ = a.𝓒;𝓖ᵃ = a.𝓖
    𝓒ᵇ = b.𝓒;𝓖ᵇ = b.𝓖
    for (ξᵃ,ξᵇ) in zip(𝓖ᵃ,𝓖ᵇ)
        𝑤 = ξᵃ.𝑤
        N = ξᵃ[:𝝭]
        B₁ = ξᵃ[:∂𝝭∂x]
        B₂ = ξᵃ[:∂𝝭∂y]
        B₁₁ = ξᵃ[:∂²𝝭∂x²]
        B₁₂ = ξᵃ[:∂²𝝭∂x∂y]
        B₂₂ = ξᵃ[:∂²𝝭∂y²]
        N̄ = ξᵇ[:𝝭]
        Γ¹₁₁ = ξᵃ.Γ¹₁₁
        Γ¹₁₂ = ξᵃ.Γ¹₁₂
        Γ¹₂₂ = ξᵃ.Γ¹₂₂
        Γ²₁₁ = ξᵃ.Γ²₁₁
        Γ²₁₂ = ξᵃ.Γ²₁₂
        Γ²₂₂ = ξᵃ.Γ²₂₂
        Γ¹₁₁₁ = ξᵃ.Γ¹₁₁₁
        Γ¹₁₁₂ = ξᵃ.Γ¹₁₁₂
        Γ¹₁₂₁ = ξᵃ.Γ¹₁₂₁
        Γ¹₁₂₂ = ξᵃ.Γ¹₁₂₂
        Γ¹₂₂₁ = ξᵃ.Γ¹₂₂₁
        Γ¹₂₂₂ = ξᵃ.Γ¹₂₂₂
        Γ²₁₁₁ = ξᵃ.Γ²₁₁₁
        Γ²₁₁₂ = ξᵃ.Γ²₁₁₂
        Γ²₁₂₁ = ξᵃ.Γ²₁₂₁
        Γ²₁₂₂ = ξᵃ.Γ²₁₂₂
        Γ²₂₂₁ = ξᵃ.Γ²₂₂₁
        Γ²₂₂₂ = ξᵃ.Γ²₂₂₂

        Γᵞᵧ₁₁ = Γ¹₁₁₁+Γ²₁₂₁
        Γᵞᵧ₂₂ = Γ¹₁₂₂+Γ²₂₂₂
        Γᵞᵧ₁₂ = Γ¹₁₁₂+Γ²₁₂₂
        Γᵞᵧ₂₁ = Γ¹₁₂₁+Γ²₂₂₁
        Γᵞ₁₁ᵧ = Γ¹₁₁₁+Γ²₁₁₂
        Γᵞ₂₂ᵧ = Γ¹₂₂₁+Γ²₂₂₂
        Γᵞ₁₂ᵧ = Γ¹₁₂₁+Γ²₁₂₂

        Γᵞᵧ₁ = Γ¹₁₁+Γ²₁₂
        Γᵞᵧ₂ = Γ¹₁₂+Γ²₂₂
        Γᵞᵧ₁Γᵝᵦ₁ = Γᵞᵧ₁^2
        Γᵞᵧ₂Γᵝᵦ₂ = Γᵞᵧ₂^2
        Γᵞᵧ₁Γᵝᵦ₂ = Γᵞᵧ₁*Γᵞᵧ₂
        Γᵞᵧ₂Γᵝᵦ₁ = Γᵞᵧ₁*Γᵞᵧ₂
        ΓᵞᵧᵦΓᵝ₁₁ = Γᵞᵧ₁*Γ¹₁₁ + Γᵞᵧ₂*Γ²₁₁
        ΓᵞᵧᵦΓᵝ₂₂ = Γᵞᵧ₁*Γ¹₂₂ + Γᵞᵧ₂*Γ²₂₂
        ΓᵞᵧᵦΓᵝ₁₂ = Γᵞᵧ₁*Γ¹₁₂ + Γᵞᵧ₂*Γ²₁₂

        for (i,xᵢ) in enumerate(𝓒ᵃ)
            I = xᵢ.𝐼
            B̃₁₁ = B₁₁[i] + 2*Γᵞᵧ₁*B₁[i] + Γ¹₁₁*B₁[i] + Γ²₁₁*B₂[i] + (Γᵞᵧ₁₁+Γᵞ₁₁ᵧ+Γᵞᵧ₁Γᵝᵦ₁+ΓᵞᵧᵦΓᵝ₁₁)*N[i]
            B̃₂₂ = B₂₂[i] + 2*Γᵞᵧ₂*B₂[i] + Γ¹₂₂*B₁[i] + Γ²₂₂*B₂[i] + (Γᵞᵧ₂₂+Γᵞ₂₂ᵧ+Γᵞᵧ₂Γᵝᵦ₂+ΓᵞᵧᵦΓᵝ₂₂)*N[i]
            B̃₁₂ = 2*B₁₂[i] + 2*(Γᵞᵧ₁*B₂[i]+Γᵞᵧ₂*B₁[i]) + 2*(Γ¹₁₂*B₁[i] + Γ²₁₂*B₂[i]) + (Γᵞᵧ₁₂+Γᵞᵧ₂₁+2*Γᵞ₁₂ᵧ+Γᵞᵧ₁Γᵝᵦ₂+Γᵞᵧ₂Γᵝᵦ₁+2*ΓᵞᵧᵦΓᵝ₁₂)*N[i]
            # B̃₁₁ = Γᵞ₁₁ᵧ*N[i] + Γ¹₁₁*B₁[i] + Γ²₁₁*B₂[i] + ΓᵞᵧᵦΓᵝ₁₁*N[i] + Γᵞᵧ₁₁*N[i] + B₁₁[i] + 2*Γᵞᵧ₁*B₁[i] + Γᵞᵧ₁Γᵝᵦ₁*N[i]
            # B̃₁₁ = Γᵞ₁₁ᵧ*N[i] + Γ¹₁₁*B₁[i] + Γ²₁₁*B₂[i] + ΓᵞᵧᵦΓᵝ₁₁*N[i] + Γᵞᵧ₁₁*N[i] + B₁₁[i] + 2*Γᵞᵧ₁*B₁[i] + Γᵞᵧ₁Γᵝᵦ₁*N[i]
            # B̃₁₁ = Γᵞ₁₁ᵧ*N[i] + Γ¹₁₁*B₁[i] + Γ²₁₁*B₂[i] + ΓᵞᵧᵦΓᵝ₁₁*N[i] + Γᵞᵧ₁₁*N[i] + B₁₁[i] + 2*Γᵞᵧ₁*B₁[i] + Γᵞᵧ₁Γᵝᵦ₁*N[i]
            for (j,xⱼ) in enumerate(𝓒ᵇ)
                J = xⱼ.𝐼
                k[9*I-8,3*J-2] += B̃₁₁*N̄[j]*𝑤
                k[9*I-7,3*J-1] += B̃₁₁*N̄[j]*𝑤
                k[9*I-6,3*J]   += B̃₁₁*N̄[j]*𝑤
                k[9*I-5,3*J-2] += B̃₂₂*N̄[j]*𝑤
                k[9*I-4,3*J-1] += B̃₂₂*N̄[j]*𝑤
                k[9*I-3,3*J]   += B̃₂₂*N̄[j]*𝑤
                k[9*I-2,3*J-2] += B̃₁₂*N̄[j]*𝑤
                k[9*I-1,3*J-1] += B̃₁₂*N̄[j]*𝑤
                k[9*I,3*J]     += B̃₁₂*N̄[j]*𝑤
            end
        end
    end
end

"""
    ScordelisLoRoof_𝐴(ap::T)

Geometry helper for the Scordelis–Lo roof benchmark: the
physical coordinates of a quadrature point on the shell surface.


    𝒙(ξ,η) = ( L ξ,  R sin(φ),  R cos(φ) ),   φ = (1 - η) π/6

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
reads `ξ.ξ`, `ξ.η` only.
"""
function ScordelisLoRoof_𝐴(ap::T) where T<:AbstractElement
    𝓒 = ap.𝓒
    ξ = ap.𝓖[1]
    w = 0.0
    N = ξ[:𝝭]
    for (i,xᵢ) in enumerate(𝓒)
        w += N[i]*xᵢ.d₃
    end
    return w
end

"""
    SphericalShell_𝐴(ap::T)

Geometry helper for the pinched spherical shell benchmark:
physical coordinates on the sphere.


    𝒙(ξ,η) = R ( sinθ cosφ,  sinθ sinφ,  cosθ ),   θ = η π/2,  φ = 2ξ π

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
reads `ξ.ξ`, `ξ.η` only.
"""
function SphericalShell_𝐴(ap::T) where T<:AbstractElement
    𝓒 = ap.𝓒
    ξ = ap.𝓖[1]
    w = 0.0
    N = ξ[:𝝭]
    for (i,xᵢ) in enumerate(𝓒)
        w += N[i]*xᵢ.d₁
    end
    return w
end
"""
    SphericalShell_𝐵(ap::T)

Geometry helper for the pinched spherical shell benchmark:
the second coordinate of the paramaterisation.


    θ(η) = η π/2,   the polar angle of the spherical parameterisation.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
reads `ξ.ξ`, `ξ.η` only.
"""
function SphericalShell_𝐵(ap::T) where T<:AbstractElement
    𝓒 = ap.𝓒
    ξ = ap.𝓖[1]
    w = 0.0
    N = ξ[:𝝭]
    for (i,xᵢ) in enumerate(𝓒)
        w += N[i]*xᵢ.d₂
    end
    return w
end

"""
    ∫𝒂ᵢvᵢgᵢdΓ(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}, E, ν, h)

Overload or Nitsche/Hellinger–Reissner sibling of the same shell term;
see the documented operators with a similar name.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set∇²𝝭!` — reads the shell fields described in the sibling operator.
"""
function ∫𝒂ᵢvᵢgᵢdΓ(ap::T, k::AbstractMatrix{Float64}, f::AbstractVector{Float64}; E::Real, ν::Real, h::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        𝒂₁₍₁₎ = ξ.𝒂₁₍₁₎
        𝒂₁₍₃₎ = ξ.𝒂₁₍₃₎
        𝒂₂₍₁₎ = ξ.𝒂₂₍₁₎
        𝒂₂₍₃₎ = ξ.𝒂₂₍₃₎
        𝒂₃₍₁₎ = ξ.𝒂₃₍₁₎
        𝒂₃₍₃₎ = ξ.𝒂₃₍₃₎
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        g₃ = ξ.g₃
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += N[i]*(α₁*𝒂₁₍₁₎^2+α₁*𝒂₂₍₁₎^2+α₂*𝒂₃₍₁₎^2)*N[j]*𝑤
                k[3*I-2,3*J-1] += 0
                k[3*I-2,3*J]   += N[i]*(α₁*𝒂₁₍₁₎*𝒂₁₍₃₎+α₁*𝒂₂₍₁₎*𝒂₂₍₃₎+α₂*𝒂₃₍₁₎*𝒂₃₍₃₎)*N[j]*𝑤
                k[3*I-1,3*J-2] += 0
                k[3*I-1,3*J-1] += 0
                k[3*I-1,3*J]   += 0
                k[3*I,3*J-2]   += N[i]*(α₁*𝒂₁₍₁₎*𝒂₁₍₃₎+α₁*𝒂₂₍₁₎*𝒂₂₍₃₎+α₂*𝒂₃₍₁₎*𝒂₃₍₃₎)*N[j]*𝑤
                k[3*I,3*J-1]   += 0
                k[3*I,3*J]     += N[i]*(α₁*𝒂₁₍₃₎^2+α₂*𝒂₃₍₃₎^2)*N[j]*𝑤
            end
            f[3*I-2] += N[i]*(α₁*𝒂₁₍₁₎*(𝒂₁₍₁₎*g₁+𝒂₁₍₃₎*g₃)+α₁*𝒂₂₍₁₎*(𝒂₂₍₁₎*g₁+𝒂₂₍₃₎*g₃)+α₂*𝒂₃₍₁₎*(𝒂₃₍₁₎*g₁+𝒂₃₍₃₎*g₃))*𝑤
            f[3*I-1] += 0
            f[3*I]   += N[i]*(α₁*𝒂₁₍₃₎*(𝒂₁₍₁₎*g₁+𝒂₁₍₃₎*g₃)+α₁*𝒂₂₍₃₎*(𝒂₂₍₁₎*g₁+𝒂₂₍₃₎*g₃)+α₂*𝒂₃₍₃₎*(𝒂₃₍₁₎*g₁+𝒂₃₍₃₎*g₃))*𝑤
        end
    end
end

"""
    Hₑ_ThinShell(ap::T)

Relative energy-norm error of a Kirchhoff shell solution,
`‖e‖_E / ‖u‖_E` over membrane plus bending energy.

The exact field is supplied with `prescribe!`.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set∇²𝝭!` — reads `ξ[:𝝭]`, the first and second derivatives, the surface
metric and the exact field.
"""
function Hₑ_ThinShell(ap::T) where T<:AbstractElement
    ΔW²= 0
    W̄² = 0
    Δu²= 0
    ū² = 0
    E = ap.𝓖[1].E
    ν = ap.𝓖[1].ν
    h = ap.𝓖[1].h
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
        # 𝑪 = @SArray [                  a¹¹*a¹¹ ν*a¹¹*a²² + (1-ν)*a¹²*a¹²                            a¹¹*a¹²;
        #              ν*a¹¹*a²² + (1-ν)*a¹²*a¹²                   a²²*a²²                            a²²*a¹²;
        #                                a¹¹*a¹²                   a²²*a¹² 0.5*((1-ν)*a¹¹*a²² + (1+ν)*a¹²*a¹²)]
        # 𝜺̄ = @SArray [ε̄₁₁,ε̄₂₂,2*ε̄₁₂]
        # 𝜿̄ = @SArray [κ̄₁₁,κ̄₂₂,2*κ̄₁₂]
        # u₁ = 0.
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
        # 𝜺 = @SArray [ε₁₁,ε₂₂,ε₁₂]
        # 𝜿 = @SArray [κ₁₁,κ₂₂,κ₁₂]
        ΔW² += 0.5*(Dᵐ*((𝜺-𝜺̄)'*𝑪*(𝜺-𝜺̄)) + Dᵇ*((𝜿-𝜿̄)'*𝑪*(𝜿-𝜿̄)))*𝑤
        W̄² += 0.5*(Dᵐ*(𝜺̄'*𝑪*𝜺̄) + Dᵇ*(𝜿̄'*𝑪*𝜿̄))*𝑤
        Δu² += ((u₁ - ū₁)^2 + (u₂ - ū₂)^2 + (u₃ - ū₃)^2)*𝑤
        ū² += (ū₁^2 + ū₂^2 + ū₃^2)*𝑤
    end
    return ΔW², W̄², Δu², ū²
end

"""
    Hₑ_ThinShell(aps::Vector{T})

Relative energy-norm error of a Kirchhoff shell solution,
`‖e‖_E / ‖u‖_E` over membrane plus bending energy.

The exact field is supplied with `prescribe!`.

# Keyword arguments
Material constants are keywords where the operator needs them; see the
signature above.

# Element requirements
`set∇²𝝭!` — reads `ξ[:𝝭]`, the first and second derivatives, the surface
metric and the exact field.
"""
function Hₑ_ThinShell(aps::Vector{T}) where T<:AbstractElement
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

end