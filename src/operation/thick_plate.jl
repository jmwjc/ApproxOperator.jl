module MindlinPlate 
    
using ..ApproxOperator: AbstractElement, Element

"""
    ∫κMγQdΩ(ap, k; E, ν, h)

Reissner–Mindlin plate bilinear form,
`∫ (κ(w,φ)·M(κ) + γ(w,φ)·Q(γ)) dΩ`.

Uses the shear correction factor 5/6. The bending and shear stiffnesses are
`D = Eh³/(12(1-ν²))` and `Dˢ = 5hE/(12(1+ν))`.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio
- `h` : plate thickness

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫κMγQdΩ(ap::T,k::AbstractMatrix; E::Real, ν::Real, h::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        Dᵇᵢᵢᵢᵢ = E*h^3/12/(1-ν^2)
        Dᵇᵢᵢⱼⱼ = E*ν*h^3/12/(1-ν^2)
        Dᵇᵢⱼᵢⱼ = E*h^3/24/(1+ν)
        Dˢ =  5/6*h*E/2/(1+ν)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += ( Dˢ*B₁[i]*B₁[j] + Dˢ*B₂[i]*B₂[j])*𝑤
                k[3*I-2,3*J-1] += (-Dˢ*B₁[i]*N[j])*𝑤
                k[3*I-2,3*J]   += (-Dˢ*B₂[i]*N[j])*𝑤
                k[3*I-1,3*J-2] += (-Dˢ*N[i]*B₁[j])*𝑤
                k[3*I-1,3*J-1] += (-Dᵇᵢᵢᵢᵢ*B₁[i]*B₁[j] - Dᵇᵢⱼᵢⱼ*B₂[i]*B₂[j] + Dˢ*N[i]*N[j])*𝑤
                k[3*I-1,3*J]   += (-Dᵇᵢᵢⱼⱼ*B₁[i]*B₂[j] - Dᵇᵢⱼᵢⱼ*B₂[i]*B₁[j])*𝑤
                k[3*I,3*J-2]   += (-Dˢ*N[i]*B₂[j])*𝑤
                k[3*I,3*J-1]   += (-Dᵇᵢᵢⱼⱼ*B₂[i]*B₁[j] - Dᵇᵢⱼᵢⱼ*B₁[i]*B₂[j])*𝑤
                k[3*I,3*J]     += (-Dᵇᵢᵢᵢᵢ*B₂[i]*B₂[j] - Dᵇᵢⱼᵢⱼ*B₁[i]*B₁[j] + Dˢ*N[i]*N[j])*𝑤
            end
        end
    end
end

"""
    ∫κκdΩ(ap, k; E, ν, h)

Bending-only block `∫ κ(δφ)·κ(φ) dΩ` of the Mindlin plate, carrying the
bending stiffness.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio
- `h` : plate thickness

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫κκdΩ(ap::T,k::AbstractMatrix; E::Real, ν::Real, h::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    Dᵢᵢᵢᵢ = E/12/(1-ν^2)
    Dᵢᵢⱼⱼ = E*ν/12/(1-ν^2)
    Dᵢⱼᵢⱼ = E/24/(1+ν)
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += (Dᵢᵢᵢᵢ*B₁[i]*B₁[j] + Dᵢⱼᵢⱼ*B₂[i]*B₂[j])*𝑤
                k[2*I-1,2*J]   += (Dᵢᵢⱼⱼ*B₁[i]*B₂[j] + Dᵢⱼᵢⱼ*B₂[i]*B₁[j])*𝑤
                k[2*I,2*J-1]   += (Dᵢᵢⱼⱼ*B₂[i]*B₁[j] + Dᵢⱼᵢⱼ*B₁[i]*B₂[j])*𝑤
                k[2*I,2*J]     += (Dᵢᵢᵢᵢ*B₂[i]*B₂[j] + Dᵢⱼᵢⱼ*B₁[i]*B₁[j])*𝑤
            end
        end
    end
end

"""
    ∫wwdΩ(ap::T,k::AbstractMatrix)

Bending stiffness of the Mindlin plate in deflection form,
`∫ κ(δw)·κ(w) dΩ` with `κ = ∇²w`.

# Keyword arguments
See the signature above.

# Element requirements
`set∇²𝝭!` — reads `ξ[:∂²𝝭∂x²]`, `ξ[:∂²𝝭∂x∂y]`, `ξ[:∂²𝝭∂y²]`, `ξ.𝑤`.
"""
function ∫wwdΩ(ap::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += N[i]*N[j]*𝑤
            end
        end
    end
end

"""
    ∫∇w∇wdΩ(ap, k; E, ν, h)

Shear block `∫ ∇w·∇w dΩ` of the Mindlin plate. Together with
[`∫φφdΩ`](@ref) and [`∫φwdΩ`](@ref) it forms the shear contribution.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio
- `h` : plate thickness

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫∇w∇wdΩ(ap::T,k::AbstractMatrix; E::Real, ν::Real, h::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    D = 5/6/h^2*E/2/(1+ν)
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += D*(B₁[i]*B₁[j] + B₂[i]*B₂[j])*𝑤
            end
        end
    end
end

"""
    ∫φφdΩ(ap, k; E, ν, h)

Rotation–rotation shear block `∫ φ·φ dΩ` of the Mindlin plate.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio
- `h` : plate thickness

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫φφdΩ(ap::T,k::AbstractMatrix; E::Real, ν::Real, h::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    D = 5/6/h^2*E/2/(1+ν)
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += D*N[i]*N[j]*𝑤
                k[2*I,2*J]     += D*N[i]*N[j]*𝑤
            end
        end
    end
end

"""
    ∫φwdΩ(ap, k; E, ν, h)
    ∫φwdΩ(a₁, a₂, k; E, ν, h)

Rotation–deflection coupling `∫ φ·∇w dΩ` of the Mindlin plate.

The two-element-set overload takes the rotation and deflection approximations
separately.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio
- `h` : plate thickness

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫φwdΩ(ap::T,k::AbstractMatrix; E::Real, ν::Real, h::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    D = 5/6/h^2*E/2/(1+ν)
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        N = ξ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,J] -= D*N[i]*B₁[j]*𝑤
                k[2*I,J]   -= D*N[i]*B₂[j]*𝑤
            end
        end
    end
end

"""
    ∫φwdΩ(a₁::T,a₂::S,k::AbstractMatrix; E, ν, h)

Overload of `∫φwdΩ` — see the documented sibling with the same name.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.𝑤`.
"""
function ∫φwdΩ(a₁::T,a₂::S,k::AbstractMatrix; E::Real, ν::Real, h::Real) where {T,S<:AbstractElement}
    𝓒₁ = a₁.𝓒; 𝓖₁ = a₁.𝓖
    𝓒₂ = a₂.𝓒; 𝓖₂ = a₂.𝓖
    D = 5/6/h^2*E/2/(1+ν)
    for (ξ₁,ξ₂) in zip(𝓖₁,𝓖₂)
        𝑤 = ξ₁.𝑤
        B₁ = ξ₁[:∂𝝭∂x]
        B₂ = ξ₁[:∂𝝭∂y]
        N = ξ₂[:𝝭]
        for (i,xᵢ) in enumerate(𝓒₁)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒₂)
                J = xⱼ.𝐼
                k[2*I-1,J] -= D*N[i]*B₁[j]*𝑤
                k[2*I,J]   -= D*N[i]*B₂[j]*𝑤
            end
        end
    end
end

"""
    ∫wqdΩ(ap::T,f::AbstractVector)

Transverse load vector, `f[I] += ∫ q Nᵢ dΩ`. The load `q` is prescribed
with `prescribe!(el, :q => f)`.

# Keyword arguments
See the signature above.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.q`, `ξ.𝑤`.
"""
function ∫wqdΩ(ap::T,f::AbstractVector) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        q = ξ.q
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] += N[i]*q*𝑤
        end
    end
end

"""
    ∫φmdΩ(ap::T,f::AbstractVector)

Distributed couple load, `f[I] += ∫ m Nᵢ dΩ`.

# Keyword arguments
See the signature above.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.m`, `ξ.𝑤`.
"""
function ∫φmdΩ(ap::T,f::AbstractVector) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        m₁ = ξ.m₁
        m₂ = ξ.m₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += N[i]*m₁*𝑤
            f[2*I]   += N[i]*m₂*𝑤
        end
    end
end

"""
    ∫wVdΓ(ap::T,f::AbstractVector)

Kirchhoff shear boundary load, `f[I] += ∫ V Nᵢ dΓ`.

# Keyword arguments
See the signature above.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.V`, `ξ.𝑤`.
"""
function ∫wVdΓ(ap::T,f::AbstractVector) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        V = ξ.V
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] += N[i]*V*𝑤
        end
    end
end

"""
    ∫φMdΓ(ap::T,f::AbstractVector)

Boundary moment load, `f[I] += ∫ M Nᵢ dΓ`.

# Keyword arguments
See the signature above.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.M`, `ξ.𝑤`.
"""
function ∫φMdΓ(ap::T,f::AbstractVector) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        M₁ = ξ.M₁
        M₂ = ξ.M₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] -= N[i]*M₁*𝑤
            f[2*I]   -= N[i]*M₂*𝑤
        end
    end
end

"""
    ∫Q∇wdΩ(ap::T,k::AbstractMatrix)

Mixed shear force coupling, `∫ Q·∇w dΩ`. The two-element-set overload
takes the shear and deflection approximations separately.

# Keyword arguments
See the signature above.

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ[:𝝭]`, `ξ.𝑤`.
"""
function ∫Q∇wdΩ(ap::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,J] += N[i]*B₁[j]*𝑤
                k[2*I,J]   += N[i]*B₂[j]*𝑤
            end
        end
    end
end

"""
    ∫Q∇wdΩ(a::T,b::S,k::AbstractMatrix)

Mixed shear force coupling, `∫ Q·∇w dΩ`. The two-element-set overload
takes the shear and deflection approximations separately.

# Keyword arguments
See the signature above.

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ[:𝝭]`, `ξ.𝑤`.
"""
function ∫Q∇wdΩ(a::T,b::S,k::AbstractMatrix) where {T<:AbstractElement,S<:AbstractElement}
    𝓒₁ = a.𝓒; 𝓖₁ = a.𝓖
    𝓒₂ = b.𝓒; 𝓖₂ = b.𝓖
    for (ξ₁,ξ₂) in zip(𝓖₁,𝓖₂)
        𝑤 = ξ₁.𝑤
        N = ξ₁[:𝝭]
        B̄₁ = ξ₂[:∂𝝭∂x]
        B̄₂ = ξ₂[:∂𝝭∂y]
        for (i,xᵢ) in enumerate(𝓒₁)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒₂)
                J = xⱼ.𝐼
                k[2*I-1,J] += N[i]*B̄₁[j]*𝑤
                k[2*I,J]   += N[i]*B̄₂[j]*𝑤
            end
        end
    end
end

"""
    ∫∇QwdΩ(ap::T,k::AbstractMatrix)

Mixed shear force coupling in gradient form, `∫ ∇Q · w dΩ`.
The transpose of [`∫Q∇wdΩ`](@ref) for the rotated variable ordering.

# Keyword arguments
See the signature above.

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ[:𝝭]`, `ξ.𝑤`.
"""
function ∫∇QwdΩ(ap::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,J] -= B₁[i]*N[j]*𝑤
                k[2*I,J]   -= B₂[i]*N[j]*𝑤
            end
        end
    end
end

"""
    ∫∇QwdΩ(a::T,b::S,k::AbstractMatrix)

Mixed shear force coupling in gradient form, `∫ ∇Q · w dΩ`.
The transpose of [`∫Q∇wdΩ`](@ref) for the rotated variable ordering.

# Keyword arguments
See the signature above.

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ[:𝝭]`, `ξ.𝑤`.
"""
function ∫∇QwdΩ(a::T,b::S,k::AbstractMatrix) where {T<:AbstractElement,S<:AbstractElement}
    𝓒₁ = a.𝓒; 𝓖₁ = a.𝓖
    𝓒₂ = b.𝓒; 𝓖₂ = b.𝓖
    for (ξ₁,ξ₂) in zip(𝓖₁,𝓖₂)
        𝑤 = ξ₁.𝑤
        B₁ = ξ₁[:∂𝝭∂x]
        B₂ = ξ₁[:∂𝝭∂y]
        N̄ = ξ₂[:𝝭]
        for (i,xᵢ) in enumerate(𝓒₁)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒₂)
                J = xⱼ.𝐼
                k[2*I-1,J] -= B₁[i]*N̄[j]*𝑤
                k[2*I,J]   -= B₂[i]*N̄[j]*𝑤
            end
        end
    end
end

"""
    ∫QφdΩ(ap::T,k::AbstractMatrix)

Mixed shear–rotation coupling, `∫ Q·φ dΩ`.

# Keyword arguments
See the signature above.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.𝑤`.
"""
function ∫QφdΩ(ap::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] -= N[i]*N[j]*𝑤
                k[2*I,2*J]     -= N[i]*N[j]*𝑤
            end
        end
    end
end

"""
    ∫QφdΩ(a::T,b::S,k::AbstractMatrix)

Mixed shear–rotation coupling, `∫ Q·φ dΩ`.

# Keyword arguments
See the signature above.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.𝑤`.
"""
function ∫QφdΩ(a::T,b::S,k::AbstractMatrix) where {T<:AbstractElement,S<:AbstractElement}
    𝓒₁ = a.𝓒; 𝓖₁ = a.𝓖
    𝓒₂ = b.𝓒; 𝓖₂ = b.𝓖
    for (ξ₁,ξ₂) in zip(𝓖₁,𝓖₂)
        𝑤 = ξ₁.𝑤
        N = ξ₁[:𝝭]
        N̄ = ξ₂[:𝝭]
        for (i,xᵢ) in enumerate(𝓒₁)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒₂)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] -= N[i]*N̄[j]*𝑤
                k[2*I,2*J]     -= N[i]*N̄[j]*𝑤
            end
        end
    end
end

"""
    ∫QQdΩ(ap, k; E, ν, h)

Shear-force mass term `∫ Q·Q dΩ`, used in mixed and stabilised Mindlin
formulations.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio
- `h` : plate thickness

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫QQdΩ(ap::T,k::AbstractMatrix; E::Real, ν::Real, h::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    D = 5/6/h^2*E/2/(1+ν)
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] -= N[i]*N[j]/D*𝑤
                k[2*I,2*J]     -= N[i]*N[j]/D*𝑤
            end
        end
    end
end

"""
    ∫QwdΓ(ap::T,k::AbstractMatrix)

Mixed shear force on the boundary, `∫ (Q·n) w dΓ`; the overloads differ in
whether the boundary value and the two element sets are supplied.

# Keyword arguments
See the signature above.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.n₁`, `ξ.n₂`, `ξ.𝑤`.
"""
function ∫QwdΓ(ap::T,k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,J] += n₁*N[i]*N[j]*𝑤
                k[2*I,J]   += n₂*N[i]*N[j]*𝑤
            end
        end
    end
end

"""
    ∫QwdΓ(a::T,b::S,k::AbstractMatrix,f::AbstractVector)

Mixed shear force on the boundary, `∫ (Q·n) w dΓ`; the overloads differ in
whether the boundary value and the two element sets are supplied.

# Keyword arguments
See the signature above.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.n₁`, `ξ.n₂`, `ξ.𝑤`.
"""
function ∫QwdΓ(a::T,b::S,k::AbstractMatrix,f::AbstractVector) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵛ = a.𝓒; 𝓖ᵛ= a.𝓖
    𝓒ʷ = b.𝓒; 𝓖ʷ= b.𝓖
    for (ξᵛ,ξʷ) in zip(𝓖ᵛ,𝓖ʷ)
        𝑤 = ξᵛ.𝑤
        Nᵛ = ξᵛ[:𝝭]
        Nʷ = ξʷ[:𝝭]
        n₁ = ξᵛ.n₁
        n₂ = ξᵛ.n₂
        g = ξʷ.g
        for (i,xᵢ) in enumerate(𝓒ᵛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ʷ)
                J = xⱼ.𝐼
                k[2*I-1,J] -= n₁*Nᵛ[i]*Nʷ[j]*𝑤
                k[2*I,J]   -= n₂*Nᵛ[i]*Nʷ[j]*𝑤
            end
            f[2*I-1] -= n₁*Nᵛ[i]*g*𝑤
            f[2*I]   -= n₂*Nᵛ[i]*g*𝑤
        end
    end
end

"""
    ∫QwdΓ(a::T,b::S,k::AbstractMatrix)

Mixed shear force on the boundary, `∫ (Q·n) w dΓ`; the overloads differ in
whether the boundary value and the two element sets are supplied.

# Keyword arguments
See the signature above.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.n₁`, `ξ.n₂`, `ξ.𝑤`.
"""
function ∫QwdΓ(a::T,b::S,k::AbstractMatrix) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵛ = a.𝓒; 𝓖ᵛ= a.𝓖
    𝓒ʷ = b.𝓒; 𝓖ʷ= b.𝓖
    for (ξᵛ,ξʷ) in zip(𝓖ᵛ,𝓖ʷ)
        𝑤 = ξᵛ.𝑤
        Nᵛ = ξᵛ[:𝝭]
        Nʷ = ξʷ[:𝝭]
        n₁ = ξᵛ.n₁
        n₂ = ξᵛ.n₂
        for (i,xᵢ) in enumerate(𝓒ᵛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ʷ)
                J = xⱼ.𝐼
                k[2*I-1,J] += n₁*Nᵛ[i]*Nʷ[j]*𝑤
                k[2*I,J]   += n₂*Nᵛ[i]*Nʷ[j]*𝑤
            end
        end
    end
end

"""
    ∫QwdΓ(ap::T,k::AbstractMatrix,f::AbstractVector)

Mixed shear force on the boundary, `∫ (Q·n) w dΓ`; the overloads differ in
whether the boundary value and the two element sets are supplied.

# Keyword arguments
See the signature above.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.n₁`, `ξ.n₂`, `ξ.𝑤`.
"""
function ∫QwdΓ(ap::T,k::AbstractMatrix,f::AbstractVector) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,J] -= n₁*N[i]*N[j]*𝑤
                k[2*I,J]   -= n₂*N[i]*N[j]*𝑤
            end
            f[2*I-1] -= n₁*N[i]*g*𝑤
            f[2*I]   -= n₂*N[i]*g*𝑤
        end
    end
end

"""
    ∫MMdΩ(ap, k; E, ν, h)

Bending-moment mass term `∫ M·M dΩ`, the dual counterpart of
[`∫QQdΩ`](@ref) in mixed plate formulations.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio
- `h` : plate thickness

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫MMdΩ(ap::T,k::AbstractMatrix{Float64}; E::Real, ν::Real, h::Real) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        C⁻¹ᵢᵢᵢᵢ = 1/E*12
        C⁻¹ᵢᵢⱼⱼ = -ν/E*12
        C⁻¹ᵢⱼᵢⱼ = 2*(1+ν)/E*12
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] -= N[i]*C⁻¹ᵢᵢᵢᵢ*N[j]*𝑤
                k[3*I-2,3*J-1] -= N[i]*C⁻¹ᵢᵢⱼⱼ*N[j]*𝑤
                k[3*I-1,3*J-2] -= N[i]*C⁻¹ᵢᵢⱼⱼ*N[j]*𝑤
                k[3*I-1,3*J-1] -= N[i]*C⁻¹ᵢᵢᵢᵢ*N[j]*𝑤
                k[3*I,3*J]     -= N[i]*C⁻¹ᵢⱼᵢⱼ*N[j]*𝑤
            end
        end
    end
end

"""
    ∫MφdΓ(aₘ::T,aᵩ::S,k::AbstractMatrix{Float64})

Bending moment on the boundary, `∫ (M·n) φ dΓ`.

# Keyword arguments
See the signature above.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.n₁`, `ξ.n₂`, `ξ.𝑤`.
"""
function ∫MφdΓ(aₘ::T,aᵩ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₘ = aₘ.𝓒;𝓖ₘ = aₘ.𝓖
    𝓒ᵩ = aᵩ.𝓒;𝓖ᵩ = aᵩ.𝓖
    for (ξₘ,ξᵩ) in zip(𝓖ₘ,𝓖ᵩ)
        𝑤 = ξₘ.𝑤
        N = ξₘ[:𝝭]
        N̄ = ξᵩ[:𝝭]
        n₁ = ξᵩ.n₁
        n₂ = ξᵩ.n₂
        for (i,xᵢ) in enumerate(𝓒ₘ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵩ)
                J = xⱼ.𝐼
                k[3*I-2,2*J-1] -= N[i]*n₁*N̄[j]*𝑤
                k[3*I-1,2*J]   -= N[i]*n₂*N̄[j]*𝑤
                k[3*I,2*J-1]   -= N[i]*n₂*N̄[j]*𝑤
                k[3*I,2*J]     -= N[i]*n₁*N̄[j]*𝑤
            end
        end
    end
end

"""
    ∫MφdΓ(aₘ::T,aᵩ::S,k::AbstractMatrix{Float64},f::AbstractVector{Float64})

Bending moment on the boundary, `∫ (M·n) φ dΓ`.

# Keyword arguments
See the signature above.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.n₁`, `ξ.n₂`, `ξ.𝑤`.
"""
function ∫MφdΓ(aₘ::T,aᵩ::S,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₘ = aₘ.𝓒;𝓖ₘ = aₘ.𝓖
    𝓒ᵩ = aᵩ.𝓒;𝓖ᵩ = aᵩ.𝓖
    for (ξₘ,ξᵩ) in zip(𝓖ₘ,𝓖ᵩ)
        𝑤 = ξₘ.𝑤
        N = ξₘ[:𝝭]
        N̄ = ξᵩ[:𝝭]
        n₁ = ξᵩ.n₁
        n₂ = ξᵩ.n₂
        n₁₁ = ξᵩ.n₁₁
        n₁₂ = ξᵩ.n₁₂
        n₂₂ = ξᵩ.n₂₂
        g₁ = ξᵩ.g₁
        g₂ = ξᵩ.g₂
        for (i,xᵢ) in enumerate(𝓒ₘ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵩ)
                J = xⱼ.𝐼
                k[3*I-2,2*J-1] += N[i]*n₁*n₁₁*N̄[j]*𝑤
                k[3*I-2,2*J]   += N[i]*n₁*n₁₂*N̄[j]*𝑤
                k[3*I-1,2*J-1] += N[i]*n₂*n₁₂*N̄[j]*𝑤
                k[3*I-1,2*J]   += N[i]*n₂*n₂₂*N̄[j]*𝑤
                k[3*I,2*J-1]   += N[i]*(n₁*n₁₂ + n₂*n₁₁)*N̄[j]*𝑤
                k[3*I,2*J]     += N[i]*(n₁*n₂₂ + n₂*n₁₂)*N̄[j]*𝑤
            end
            f[3*I-2] += N[i]*(n₁*n₁₁*g₁ + n₁*n₁₂*g₂)*𝑤
            f[3*I-1] += N[i]*(n₂*n₁₂*g₁ + n₂*n₂₂*g₂)*𝑤
            f[3*I]   += N[i]*((n₁*n₁₂+n₂*n₁₁)*g₁ + (n₁*n₂₂+n₂*n₁₂)*g₂)*𝑤 
        end
    end
end

"""
    ∫M∇φdΩ(aₘ::T,aᵩ::S,k::AbstractMatrix{Float64})

Mixed bending moment coupling, `∫ M:∇φ dΩ`.

# Keyword arguments
See the signature above.

# Element requirements
`set∇𝝭!` — reads `ξ[:𝝭]`, `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫M∇φdΩ(aₘ::T,aᵩ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₘ = aₘ.𝓒;𝓖ₘ = aₘ.𝓖
    𝓒ᵩ = aᵩ.𝓒;𝓖ᵩ = aᵩ.𝓖
    for (ξₘ,ξᵩ) in zip(𝓖ₘ,𝓖ᵩ)
        𝑤 = ξₘ.𝑤
        N = ξₘ[:𝝭]
        B₁ = ξᵩ[:∂𝝭∂x]
        B₂ = ξᵩ[:∂𝝭∂y]
        for (i,xᵢ) in enumerate(𝓒ₘ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵩ)
                J = xⱼ.𝐼
                k[3*I-2,2*J-1] -= N[i]*B₁[j]*𝑤
                k[3*I-1,2*J]   -= N[i]*B₂[j]*𝑤
                k[3*I,2*J-1]   -= N[i]*B₂[j]*𝑤
                k[3*I,2*J]     -= N[i]*B₁[j]*𝑤
            end
        end
    end
end

"""
    ∫∇MφdΩ(aₘ::T,aᵩ::S,k::AbstractMatrix{Float64})

Mixed bending moment coupling in gradient form, `∫ ∇·M · φ dΩ`.

# Keyword arguments
See the signature above.

# Element requirements
`set∇𝝭!` — reads `ξ[:𝝭]`, `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫∇MφdΩ(aₘ::T,aᵩ::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₘ = aₘ.𝓒;𝓖ₘ = aₘ.𝓖
    𝓒ᵩ = aᵩ.𝓒;𝓖ᵩ = aᵩ.𝓖
    for (ξₘ,ξᵩ) in zip(𝓖ₘ,𝓖ᵩ)
        𝑤 = ξₘ.𝑤
        B₁ = ξₘ[:∂𝝭∂x]
        B₂ = ξₘ[:∂𝝭∂y]
        N = ξᵩ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒ₘ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵩ)
                J = xⱼ.𝐼
                k[3*I-2,2*J-1] += B₁[i]*N[j]*𝑤
                k[3*I-1,2*J]   += B₂[i]*N[j]*𝑤
                k[3*I,2*J-1]   += B₂[i]*N[j]*𝑤
                k[3*I,2*J]     += B₁[i]*N[j]*𝑤
            end
        end
    end
end

"""
    ∫αwwdΓ(ap, k[, f]; α)

Penalty term on the deflection, `α ∫ w w̄ dΓ`, for a weak simply-supported or
clamped edge.

# Keyword arguments
- `α` : penalty coefficient

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.g`, `ξ.𝑤`.
"""
function ∫αwwdΓ(ap::T,k::AbstractMatrix,f::AbstractVector; α::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += α*N[i]*N[j]*𝑤
            end
            f[I] += α*N[i]*g*𝑤
        end
    end
end

"""
    ∫αwwdΓ(ap::T,k::AbstractMatrix; α)

Overload of `∫αwwdΓ` — see the documented sibling with the same name.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.𝑤`.
"""
function ∫αwwdΓ(ap::T,k::AbstractMatrix; α::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += α*N[i]*N[j]*𝑤
            end
        end
    end
end

"""
    ∫αφφdΓ(ap, k[, f]; α)

Penalty term on the rotation, `α ∫ φ φ̄ dΓ`, for a weakly clamped edge.

# Keyword arguments
- `α` : penalty coefficient

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.g`, `ξ.𝑤`.
"""
function ∫αφφdΓ(ap::T,k::AbstractMatrix,f::AbstractVector; α::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        n₁₁ = ξ.n₁₁
        n₁₂ = ξ.n₁₂
        n₂₂ = ξ.n₂₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += α*N[i]*n₁₁*N[j]*𝑤
                k[2*I,2*J-1]   += α*N[i]*n₁₂*N[j]*𝑤
                k[2*I-1,2*J]   += α*N[i]*n₁₂*N[j]*𝑤
                k[2*I,2*J]     += α*N[i]*n₂₂*N[j]*𝑤
            end
            f[2*I-1] += α*N[i]*(n₁₁*g₁+n₁₂*g₂)*𝑤
            f[2*I]   += α*N[i]*(n₁₂*g₁+n₂₂*g₂)*𝑤
        end
    end
end

"""
    ∫αφφdΓ(ap::T,k::AbstractMatrix; α)

Overload of `∫αφφdΓ` — see the documented sibling with the same name.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.𝑤`.
"""
function ∫αφφdΓ(ap::T,k::AbstractMatrix; α::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        n₁₁ = ξ.n₁₁
        n₁₂ = ξ.n₁₂
        n₂₂ = ξ.n₂₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += α*N[i]*n₁₁*N[j]*𝑤
                k[2*I,2*J-1]   += α*N[i]*n₁₂*N[j]*𝑤
                k[2*I-1,2*J]   += α*N[i]*n₁₂*N[j]*𝑤
                k[2*I,2*J]     += α*N[i]*n₂₂*N[j]*𝑤
            end
        end
    end
end

"""
    ∫∇w∇wdΩ_MITC(ap::Element{:Quad4}, k; E, ν, h)

Shear block of the **MITC4** element: the transverse shear is sampled at the
mid-side points instead of the Gauss points, which removes shear locking.

Requires a `Quad4` element; the sampling points are taken from the element's
edge DOFs.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio
- `h` : plate thickness

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫∇w∇wdΩ_MITC(ap::Element{:Quad4},k::AbstractMatrix; E::Real, ν::Real, h::Real)
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    D = 5/6/h^2*E/2/(1+ν)

    x₁,x₂,x₃,x₄ = (xᵢ.x for xᵢ in 𝓒)
    y₁,y₂,y₃,y₄ = (xᵢ.y for xᵢ in 𝓒)
    ∂Nᵃ∂ξ = (-0.5, 0.5, 0.0, 0.0)
    ∂Nᵇ∂η = ( 0.0,-0.5, 0.5, 0.0)
    ∂Nᶜ∂ξ = ( 0.0, 0.0, 0.5,-0.5)
    ∂Nᵈ∂η = (-0.5, 0.0, 0.0, 0.5)
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        ξ̄ = ξ.ξ
        η̄ = ξ.η
        ∂N₁∂ξ = - 0.25*(1-η̄)
        ∂N₂∂ξ =   0.25*(1-η̄)
        ∂N₃∂ξ =   0.25*(1+η̄)
        ∂N₄∂ξ = - 0.25*(1+η̄)
        ∂N₁∂η = - 0.25*(1-ξ̄)
        ∂N₂∂η = - 0.25*(1+ξ̄)
        ∂N₃∂η =   0.25*(1+ξ̄)
        ∂N₄∂η =   0.25*(1-ξ̄)
        ∂x∂ξ = ∂N₁∂ξ*x₁ + ∂N₂∂ξ*x₂ + ∂N₃∂ξ*x₃ + ∂N₄∂ξ*x₄
        ∂x∂η = ∂N₁∂η*x₁ + ∂N₂∂η*x₂ + ∂N₃∂η*x₃ + ∂N₄∂η*x₄
        ∂y∂ξ = ∂N₁∂ξ*y₁ + ∂N₂∂ξ*y₂ + ∂N₃∂ξ*y₃ + ∂N₄∂ξ*y₄
        ∂y∂η = ∂N₁∂η*y₁ + ∂N₂∂η*y₂ + ∂N₃∂η*y₃ + ∂N₄∂η*y₄
        detJ = ∂x∂ξ*∂y∂η - ∂x∂η*∂y∂ξ
        ∂ξ∂x =   ∂y∂η/detJ
        ∂η∂x = - ∂y∂ξ/detJ
        ∂ξ∂y = - ∂x∂η/detJ
        ∂η∂y =   ∂x∂ξ/detJ
        ∂N̄∂ξ = 0.5*(1-η̄).*∂Nᵃ∂ξ .+ 0.5*(1+η̄).*∂Nᶜ∂ξ
        ∂N̄∂η = 0.5*(1-ξ̄).*∂Nᵈ∂η .+ 0.5*(1+ξ̄).*∂Nᵇ∂η
        B̄₁ = ∂ξ∂x.*∂N̄∂ξ .+ ∂η∂x.*∂N̄∂η
        B̄₂ = ∂ξ∂y.*∂N̄∂ξ .+ ∂η∂y.*∂N̄∂η
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += D*(B̄₁[i]*B̄₁[j] + B̄₂[i]*B̄₂[j])*𝑤
            end
        end
    end
end


"""
    ∫φφdΩ_MITC(ap::Element{:Quad4}, k; E, ν, h)

Rotation–rotation shear block of the **MITC4** element (mid-side sampling).

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio
- `h` : plate thickness

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫φφdΩ_MITC(ap::Element{:Quad4},k::AbstractMatrix; E::Real, ν::Real, h::Real)
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    D = 5/6/h^2*E/2/(1+ν)
    Nᵃ = (0.5,0.5,0.0,0.0)
    Nᵇ = (0.0,0.5,0.5,0.0)
    Nᶜ = (0.0,0.0,0.5,0.5)
    Nᵈ = (0.5,0.0,0.0,0.5)

    for ξ in 𝓖
        𝑤 = ξ.𝑤
        ξ̄ = ξ.ξ
        η̄ = ξ.η
        N̄₁ = 0.5*(1-η̄).*Nᵃ .+ 0.5*(1+η̄).*Nᶜ
        N̄₂ = 0.5*(1-ξ̄).*Nᵈ .+ 0.5*(1+ξ̄).*Nᵇ
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += D*N̄₁[i]*N̄₁[j]*𝑤
                k[2*I,2*J]     += D*N̄₂[i]*N̄₂[j]*𝑤
            end
        end
    end
end

"""
    ∫φwdΩ_MITC(ap::Element{:Quad4}, k; E, ν, h)

Rotation–deflection shear coupling of the **MITC4** element (mid-side
sampling).

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio
- `h` : plate thickness

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫φwdΩ_MITC(ap::Element{:Quad4},k::AbstractMatrix; E::Real, ν::Real, h::Real)
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    D = 5/6/h^2*E/2/(1+ν)

    x₁,x₂,x₃,x₄ = (xᵢ.x for xᵢ in 𝓒)
    y₁,y₂,y₃,y₄ = (xᵢ.y for xᵢ in 𝓒)
    Nᵃ = (0.5,0.5,0.0,0.0)
    Nᵇ = (0.0,0.5,0.5,0.0)
    Nᶜ = (0.0,0.0,0.5,0.5)
    Nᵈ = (0.5,0.0,0.0,0.5)
    ∂Nᵃ∂ξ = (-0.5, 0.5, 0.0, 0.0)
    ∂Nᵇ∂η = ( 0.0,-0.5, 0.5, 0.0)
    ∂Nᶜ∂ξ = ( 0.0, 0.0, 0.5,-0.5)
    ∂Nᵈ∂η = (-0.5, 0.0, 0.0, 0.5)
 
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        ξ̄ = ξ.ξ
        η̄ = ξ.η
        N̄₁ = 0.5*(1-η̄).*Nᵃ .+ 0.5*(1+η̄).*Nᶜ
        N̄₂ = 0.5*(1-ξ̄).*Nᵈ .+ 0.5*(1+ξ̄).*Nᵇ
        ∂N₁∂ξ = - 0.25*(1-η̄)
        ∂N₂∂ξ =   0.25*(1-η̄)
        ∂N₃∂ξ =   0.25*(1+η̄)
        ∂N₄∂ξ = - 0.25*(1+η̄)
        ∂N₁∂η = - 0.25*(1-ξ̄)
        ∂N₂∂η = - 0.25*(1+ξ̄)
        ∂N₃∂η =   0.25*(1+ξ̄)
        ∂N₄∂η =   0.25*(1-ξ̄)
        ∂x∂ξ = ∂N₁∂ξ*x₁ + ∂N₂∂ξ*x₂ + ∂N₃∂ξ*x₃ + ∂N₄∂ξ*x₄
        ∂x∂η = ∂N₁∂η*x₁ + ∂N₂∂η*x₂ + ∂N₃∂η*x₃ + ∂N₄∂η*x₄
        ∂y∂ξ = ∂N₁∂ξ*y₁ + ∂N₂∂ξ*y₂ + ∂N₃∂ξ*y₃ + ∂N₄∂ξ*y₄
        ∂y∂η = ∂N₁∂η*y₁ + ∂N₂∂η*y₂ + ∂N₃∂η*y₃ + ∂N₄∂η*y₄
        detJ = ∂x∂ξ*∂y∂η - ∂x∂η*∂y∂ξ
        ∂ξ∂x =   ∂y∂η/detJ
        ∂η∂x = - ∂y∂ξ/detJ
        ∂ξ∂y = - ∂x∂η/detJ
        ∂η∂y =   ∂x∂ξ/detJ
        ∂N̄∂ξ = 0.5*(1-η̄).*∂Nᵃ∂ξ .+ 0.5*(1+η̄).*∂Nᶜ∂ξ
        ∂N̄∂η = 0.5*(1-ξ̄).*∂Nᵈ∂η .+ 0.5*(1+ξ̄).*∂Nᵇ∂η
        B̄₁ = ∂ξ∂x.*∂N̄∂ξ .+ ∂η∂x.*∂N̄∂η
        B̄₂ = ∂ξ∂y.*∂N̄∂ξ .+ ∂η∂y.*∂N̄∂η
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,J] -= D*N̄₁[i]*B̄₁[j]*𝑤
                k[2*I,J]   -= D*N̄₂[i]*B̄₂[j]*𝑤
            end
        end
    end
end

"""
    ∫φφdΩ_DSG(ap::Element{:Tri3}, k; E, ν, h)

Rotation–rotation block of the **DSG** (Discrete Shear Gap) triangle.
The shear gap is integrated exactly on the element boundary, which avoids
shear locking on `Tri3`.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio
- `h` : plate thickness

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫φφdΩ_DSG(ap::Element{:Tri3},k::AbstractMatrix; E::Real, ν::Real, h::Real)
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    D = 5/6/h^2*E/2/(1+ν)
    x₁,x₂,x₃ = (xᵢ.x for xᵢ in 𝓒)
    y₁,y₂,y₃ = (xᵢ.y for xᵢ in 𝓒)
    𝐽 = ap.𝐽
    𝑤 = 𝐽/2
    a = x₂-x₁
    b = y₂-y₁
    c = y₃-y₁
    d = x₃-x₁
    B₁₁ = 0.5/𝐽.*(  𝐽, a*c,-b*d)
    B₁₂ = 0.5/𝐽.*(0.0, b*c,-b*c)
    B₂₁ = 0.5/𝐽.*(0.0,-a*d, a*d)
    B₂₂ = 0.5/𝐽.*(  𝐽,-b*d, a*c)

    for (i,xᵢ) in enumerate(𝓒)
        I = xᵢ.𝐼
        for (j,xⱼ) in enumerate(𝓒)
            J = xⱼ.𝐼
            k[2*I-1,2*J-1] += D*(B₁₁[i]*B₁₁[j]+B₂₁[i]*B₂₁[j])*𝑤
            k[2*I-1,2*J]   += D*(B₁₁[i]*B₁₂[j]+B₂₁[i]*B₂₂[j])*𝑤
            k[2*I,2*J-1]   += D*(B₁₂[i]*B₁₁[j]+B₂₂[i]*B₂₁[j])*𝑤
            k[2*I,2*J]     += D*(B₁₂[i]*B₁₂[j]+B₂₂[i]*B₂₂[j])*𝑤
        end
    end
end

"""
    ∫φwdΩ_DSG(ap::Element{:Tri3}, k; E, ν, h)

Rotation–deflection coupling of the **DSG** triangle.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio
- `h` : plate thickness

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫φwdΩ_DSG(ap::Element{:Tri3},k::AbstractMatrix; E::Real, ν::Real, h::Real)
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    D = 5/6/h^2*E/2/(1+ν)
    x₁,x₂,x₃ = (xᵢ.x for xᵢ in 𝓒)
    y₁,y₂,y₃ = (xᵢ.y for xᵢ in 𝓒)
    𝐽 = ap.𝐽
    a = x₂-x₁
    b = y₂-y₁
    c = y₃-y₁
    d = x₃-x₁
    Bᵠ₁₁ = 0.5/𝐽.*(  𝐽, a*c,-b*d)
    Bᵠ₁₂ = 0.5/𝐽.*(0.0, b*c,-b*c)
    Bᵠ₂₁ = 0.5/𝐽.*(0.0,-a*d, a*d)
    Bᵠ₂₂ = 0.5/𝐽.*(  𝐽,-b*d, a*c)
 
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,J] -= D*(Bᵠ₁₁[i]*B₁[j]+Bᵠ₂₁[i]*B₂[j])*𝑤
                k[2*I,J]   -= D*(Bᵠ₁₂[i]*B₁[j]+Bᵠ₂₂[i]*B₂[j])*𝑤
            end
        end
    end
end

"""
    ∫∇wσ∇wdΩ(ap, k; h)

Deflection-gradient block of a `u/p` or `u/σ` mixed formulation,
`∫ ∇w·σ·∇w dΩ`, where `σ` comes from its own element set.

# Keyword arguments
- `h` : plate thickness

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫∇wσ∇wdΩ(ap::T, k::AbstractMatrix; h::Real) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        σ₁₁ = ξ.σ₁₁
        σ₂₂ = ξ.σ₂₂
        σ₁₂ = ξ.σ₁₂
        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j, xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I, J] += (
                    σ₁₁ * B₁[i] * B₁[j] +
                    σ₂₂ * B₂[i] * B₂[j] +
                    σ₁₂ * (B₁[i] * B₂[j] + B₂[i] * B₁[j])
                )/h^2 * 𝑤
            end
        end
    end
end

"""
    ∫∇φσ∇φdΩ(ap, k; h)

Rotation-gradient block of a mixed formulation, `∫ ∇φ·σ·∇φ dΩ`.

# Keyword arguments
- `h` : plate thickness

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫∇φσ∇φdΩ(ap::T, k::AbstractMatrix; h::Real) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        σ₁₁ = ξ.σ₁₁
        σ₂₂ = ξ.σ₂₂
        σ₁₂ = ξ.σ₁₂
        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j, xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1, 2*J-1] += (
                    σ₁₁ * B₁[i] * B₁[j] +
                    σ₂₂ * B₂[i] * B₂[j] +
                    σ₁₂ * (B₁[i] * B₂[j] + B₂[i] * B₁[j])
                )/12 * 𝑤
                k[2*I, 2*J]     += (
                    σ₁₁ * B₁[i] * B₁[j] +
                    σ₂₂ * B₂[i] * B₂[j] +
                    σ₁₂ * (B₁[i] * B₂[j] + B₂[i] * B₁[j])
                )/12 * 𝑤
            end
        end
    end
end

"""
    ∫ρwwdΩ(ap, k; ρ, h)

Deflection mass matrix of the Mindlin plate, `∫ ρ w w / h² dΩ`.

# Keyword arguments
- `ρ` : mass density
- `h` : plate thickness

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.𝑤`.
"""
function ∫ρwwdΩ(ap::T, k::AbstractMatrix; ρ::Real, h::Real) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j, xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I, J] += ρ*N[i]*N[j]/h^2*𝑤
            end
        end
    end
end

"""
    ∫ρφφdΩ(ap, k; ρ)

Rotary-inertia mass matrix, `∫ ρ φ·φ dΩ`.

# Keyword arguments
- `ρ` : mass density

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.𝑤`.
"""
function ∫ρφφdΩ(ap::T, k::AbstractMatrix; ρ::Real) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j, xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2I-1, 2J-1] += ρ*N[i]*N[j]/12*𝑤
                k[2I,   2J]   += ρ*N[i]*N[j]/12*𝑤
            end
        end
    end
end

"""
    L₂Q(ap::T)

Relative L₂ error of the shear force, `‖Q-Q̄‖/‖Q̄‖`.

# Keyword arguments
See the signature above.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.Q₁`, `ξ.Q₂`, `ξ.𝑤`.
"""
function L₂Q(ap::T) where T<:AbstractElement
    ΔQ²= BigFloat(0.0)
    Q̄² = BigFloat(0.0)
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        Q̄₁ = ξ.Q₁
        Q̄₂ = ξ.Q₂
        Q₁ = 0.0
        Q₂ = 0.0
        for (i,xᵢ) in enumerate(ap.𝓒)
            Q₁ += N[i]*xᵢ.q₁
            Q₂ += N[i]*xᵢ.q₂
        end
        ΔQ² +=((Q₁ - Q̄₁)^2 + (Q₂ - Q̄₂)^2)*𝑤
        Q̄²  += (Q̄₁^2 + Q̄₂^2)*𝑤
    end
    return ΔQ², Q̄²
end

"""
    L₂Q(aps::Vector{T})

Relative L₂ error of the shear force, `‖Q-Q̄‖/‖Q̄‖`.

# Keyword arguments
See the signature above.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.Q₁`, `ξ.Q₂`, `ξ.𝑤`.
"""
function L₂Q(aps::Vector{T}) where T<:AbstractElement
    L₂Norm_ΔQ²= BigFloat(0.0)
    L₂Norm_Q̄² = BigFloat(0.0)
    for ap in aps
        ΔQ², Q̄² = L₂Q(ap)
        L₂Norm_ΔQ² += ΔQ²
        L₂Norm_Q̄²  += Q̄²
    end
    return (L₂Norm_ΔQ²/L₂Norm_Q̄²)^0.5
end

"""
    L₂γ(ap::T)

Relative L₂ error of the transverse shear strain, `‖γ-γ̄‖/‖γ̄‖`.

# Keyword arguments
See the signature above.

# Element requirements
`set∇𝝭!` and `set𝝭!` — reads the gradients, `ξ.γ₁`, `ξ.γ₂`, `ξ.𝑤`.
"""
function L₂γ(ap::T) where T<:AbstractElement
    Δγ²= BigFloat(0.0)
    γ̄² = BigFloat(0.0)
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        γ̄₁ = ξ.γ₁
        γ̄₂ = ξ.γ₂
        γ₁ = 0.0
        γ₂ = 0.0
        for (i,xᵢ) in enumerate(ap.𝓒)
            γ₁ += B₁[i]*xᵢ.d - N[i]*xᵢ.d₁
            γ₂ += B₂[i]*xᵢ.d - N[i]*xᵢ.d₂
        end
        Δγ² += ((γ₁ - γ̄₁)^2 + (γ₂ - γ̄₂)^2)*𝑤
        γ̄²  += (γ̄₁^2 + γ̄₂^2)*𝑤
    end
    return Δγ², γ̄²
end

"""
    L₂γ(aps::Vector{T})

Relative L₂ error of the transverse shear strain, `‖γ-γ̄‖/‖γ̄‖`.

# Keyword arguments
See the signature above.

# Element requirements
`set∇𝝭!` and `set𝝭!` — reads the gradients, `ξ.γ₁`, `ξ.γ₂`, `ξ.𝑤`.
"""
function L₂γ(aps::Vector{T}) where T<:AbstractElement
    L₂Norm_Δγ²= BigFloat(0.0)
    L₂Norm_γ̄² = BigFloat(0.0)
    for ap in aps
        Δγ², γ̄² = L₂γ(ap)
        L₂Norm_Δγ² += Δγ²
        L₂Norm_γ̄²  += γ̄²
    end
    return (L₂Norm_Δγ²/L₂Norm_γ̄²)^0.5
end

"""
    L₂γ(a::T,b::S)

Relative L₂ error of the transverse shear strain, `‖γ-γ̄‖/‖γ̄‖`.

# Keyword arguments
See the signature above.

# Element requirements
`set∇𝝭!` and `set𝝭!` — reads the gradients, `ξ.γ₁`, `ξ.γ₂`, `ξ.𝑤`.
"""
function L₂γ(a::T,b::S) where {T,S<:AbstractElement}
    Δγ²= BigFloat(0.0)
    γ̄² = BigFloat(0.0)
    for (ξʷ,ξᵠ) in zip(a.𝓖,b.𝓖)
        𝑤 = ξʷ.𝑤
        N = ξᵠ[:𝝭]
        B₁ = ξʷ[:∂𝝭∂x]
        B₂ = ξʷ[:∂𝝭∂y]
        γ̄₁ = ξʷ.γ₁
        γ̄₂ = ξʷ.γ₂
        γ₁ = 0.0
        γ₂ = 0.0
        for (i,xᵢ) in enumerate(a.𝓒)
            γ₁ += B₁[i]*xᵢ.d
            γ₂ += B₂[i]*xᵢ.d
        end
        for (i,xᵢ) in enumerate(b.𝓒)
            γ₁ -= N[i]*xᵢ.d₁
            γ₂ -= N[i]*xᵢ.d₂
        end
        Δγ² += ((γ₁ - γ̄₁)^2 + (γ₂ - γ̄₂)^2)*𝑤
        γ̄²  += (γ̄₁^2 + γ̄₂^2)*𝑤
    end
    return Δγ², γ̄²
end

"""
    L₂γ(as::Vector{T},bs::Vector{S})

Relative L₂ error of the transverse shear strain, `‖γ-γ̄‖/‖γ̄‖`.

# Keyword arguments
See the signature above.

# Element requirements
`set∇𝝭!` and `set𝝭!` — reads the gradients, `ξ.γ₁`, `ξ.γ₂`, `ξ.𝑤`.
"""
function L₂γ(as::Vector{T},bs::Vector{S}) where {T,S<:AbstractElement}
    L₂Norm_Δγ²= BigFloat(0.0)
    L₂Norm_γ̄² = BigFloat(0.0)
    for (a,b) in zip(as,bs)
        Δγ², γ̄² = L₂γ(a,b)
        L₂Norm_Δγ² += Δγ²
        L₂Norm_γ̄²  += γ̄²
    end
    return (L₂Norm_Δγ²/L₂Norm_γ̄²)^0.5
end


"""
    L₂φ(ap::T)

Relative L₂ error of the rotation field.


    ‖φʰ - φ̄‖_L₂ / ‖φ̄‖_L₂,   φʰ = Σᵢ Nᵢ φᵢ

# Keyword arguments
See the signature above.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.φ₁`, `ξ.φ₂`, `ξ.𝑤`.
"""
function L₂φ(ap::T) where T<:AbstractElement
    Δφ²= BigFloat(0.0)
    φ̄² = BigFloat(0.0)
    Δφ²= BigFloat(0.0)
    φ̄² = BigFloat(0.0)
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        φ̄₁ = ξ.φ₁
        φ̄₂ = ξ.φ₂
        φ₁ = 0.0
        φ₂ = 0.0
        φ₁ = 0.0
        φ₂ = 0.0
        for (i,xᵢ) in enumerate(ap.𝓒)
            φ₁ += N[i]*xᵢ.d₁
            φ₂ += N[i]*xᵢ.d₂
        end
        Δφ² +=((φ₁ - φ̄₁)^2 + (φ₂ - φ̄₂)^2)*𝑤
        φ̄²  += (φ̄₁^2 + φ̄₂^2)*𝑤
    end
    return Δφ², φ̄²
end

"""
    L₂φ(aps::Vector{T})

Relative L₂ error of the rotation field.


    Element-set overload of `L₂φ(ap)`.

# Keyword arguments
See the signature above.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.φ₁`, `ξ.φ₂`, `ξ.𝑤`.
"""
function L₂φ(aps::Vector{T}) where T<:AbstractElement
    L₂Norm_Δφ²= BigFloat(0.0)
    L₂Norm_φ̄² = BigFloat(0.0)
    for ap in aps
        Δφ², φ̄² = L₂φ(ap)
        L₂Norm_Δφ² += Δφ²
        L₂Norm_φ̄²  += φ̄²
    end
    return (L₂Norm_Δφ²/L₂Norm_φ̄²)^0.5
end

"""
    H₁φ(ap::T)

Relative H¹ error of the rotation field, combining the L₂ and gradient
parts.


    ‖(φʰ - φ̄, ∇φʰ - ∇φ̄)‖ / ‖(φ̄, ∇φ̄)‖

# Keyword arguments
See the signature above.

# Element requirements
`set∇𝝭!` — reads `ξ[:𝝭]`, gradients, `ξ.φ₁`, `ξ.φ₂`, `ξ.𝑤`.
"""
function H₁φ(ap::T) where T<:AbstractElement
    Δ∇φ²= 0.0
    ∇φ̄² = 0.0
    Δφ²= 0.0
    φ̄² = 0.0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        φ̄₁ = ξ.φ₁
        ∂φ̄₁∂x = ξ.∂φ₁∂x
        ∂φ̄₁∂y = ξ.∂φ₁∂y
        φ̄₂ = ξ.φ₂
        ∂φ̄₂∂x = ξ.∂φ₂∂x
        ∂φ̄₂∂y = ξ.∂φ₂∂y
        φ₁ = 0.0
        ∂φ₁∂x = 0.0
        ∂φ₁∂y = 0.0
        φ₂ = 0.0
        ∂φ₂∂x = 0.0
        ∂φ₂∂y = 0.0
        for (i,xᵢ) in enumerate(ap.𝓒)
            φ₁ += N[i]*xᵢ.d₁
            ∂φ₁∂x += B₁[i]*xᵢ.d₁
            ∂φ₁∂y += B₂[i]*xᵢ.d₁
            φ₂ += N[i]*xᵢ.d₂
            ∂φ₂∂x += B₁[i]*xᵢ.d₂
            ∂φ₂∂y += B₂[i]*xᵢ.d₂
        end
        Δ∇φ² += ((∂φ₁∂x - ∂φ̄₁∂x)^2 + (∂φ₁∂y - ∂φ̄₁∂y)^2 + (∂φ₂∂x - ∂φ̄₂∂x)^2 + (∂φ₂∂y - ∂φ̄₂∂y)^2)*𝑤
        ∇φ̄²  += (∂φ̄₁∂x^2 + ∂φ̄₁∂y^2 + ∂φ̄₂∂x^2 + ∂φ̄₂∂y^2)*𝑤
        Δφ² +=((φ₁ - φ̄₁)^2 + (φ₂ - φ̄₂)^2)*𝑤
        φ̄²  += (φ̄₁^2 + φ̄₂^2)*𝑤
    end
    return Δ∇φ², ∇φ̄², Δφ², φ̄²
end

"""
    L₂w(ap::T)

Relative L₂ error of the deflection.


    ‖wʰ - w̄‖_L₂ / ‖w̄‖_L₂,   wʰ = Σᵢ Nᵢ wᵢ

# Keyword arguments
See the signature above.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.w`, `ξ.𝑤`.
"""
function L₂w(ap::T) where T<:AbstractElement
    Δw²= BigFloat(0.0)
    w̄² = BigFloat(0.0)
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        w̄ = ξ.w
        w = 0.0
        for (i,xᵢ) in enumerate(ap.𝓒)
            w += N[i]*xᵢ.d
        end
        Δw² +=(w - w̄)^2*𝑤
        w̄²  += w̄^2*𝑤
    end
    return Δw², w̄²
end

"""
    H₁w(ap::T)

Relative H¹ error of the deflection (L₂ plus the gradient part).


    ‖(wʰ - w̄, ∇wʰ - ∇w̄)‖ / ‖(w̄, ∇w̄)‖

# Keyword arguments
See the signature above.

# Element requirements
`set∇𝝭!` — reads `ξ[:𝝭]`, `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.w`, `ξ.∂w∂x`, ...
"""
function H₁w(ap::T) where T<:AbstractElement
    Δ∇w²= BigFloat(0.0)
    ∇w̄² = BigFloat(0.0)
    Δw²= BigFloat(0.0)
    w̄² = BigFloat(0.0)
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        w̄ = ξ.w
        ∂w̄∂x = ξ.∂w∂x
        ∂w̄∂y = ξ.∂w∂y
        w = 0.0
        ∂w∂x = 0.0
        ∂w∂y = 0.0
        for (i,xᵢ) in enumerate(ap.𝓒)
            w += N[i]*xᵢ.d
            ∂w∂x += B₁[i]*xᵢ.d
            ∂w∂y += B₂[i]*xᵢ.d
        end
        Δw² += (w - w̄)^2*𝑤
        w̄²  += w̄^2*𝑤
        Δ∇w² += ((∂w∂x - ∂w̄∂x)^2 + (∂w∂y - ∂w̄∂y)^2)*𝑤
        ∇w̄² += (∂w̄∂x^2 + ∂w̄∂y^2)*𝑤
    end
    return Δ∇w², ∇w̄², Δw², w̄²
end

"""
    L₂w(aps::Vector{T})

Relative L₂ error of the deflection.


    Element-set overload of `L₂w(ap)`.

# Keyword arguments
See the signature above.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.w`, `ξ.𝑤`.
"""
function L₂w(aps::Vector{T}) where T<:AbstractElement
    L₂Norm_Δw²= BigFloat(0.0)
    L₂Norm_w̄² = BigFloat(0.0)
    for ap in aps
        Δw², w̄² = L₂w(ap)
        L₂Norm_Δw² += Δw²
        L₂Norm_w̄²  += w̄²
    end
    return (L₂Norm_Δw²/L₂Norm_w̄²)^0.5
end

"""
    H₁(aps_w::Vector{Tʷ},aps_φ::Vector{Tᵠ})

Relative H¹ error of the mixed `(w, φ)` solution.


    Relative H¹ error of the pair `(w, φ)`, combining the deflection and
    rotation errors over the two element sets.

# Keyword arguments
See the signature above.

# Element requirements
`set∇𝝭!` — reads the gradients and the exact field.
"""
function H₁(aps_w::Vector{Tʷ},aps_φ::Vector{Tᵠ}) where {Tʷ,Tᵠ<:AbstractElement}
    H₁Norm_Δu²= 0.0
    H₁Norm_ū² = 0.0
    L₂Norm_Δu²= 0.0
    L₂Norm_ū² = 0.0
    for ap in aps_w
        Δ∇w², ∇w̄², Δw², w̄² = H₁w(ap)
        H₁Norm_Δu² += Δ∇w²
        H₁Norm_ū²  += ∇w̄²
        L₂Norm_Δu² += Δw²
        L₂Norm_ū²  += w̄²
    end
    for ap in aps_φ
        Δ∇φ², ∇φ̄², Δφ², φ̄² = H₁φ(ap)
        H₁Norm_Δu² += Δ∇φ²
        H₁Norm_ū²  += ∇φ̄²
        L₂Norm_Δu² += Δφ²
        L₂Norm_ū²  += φ̄²
    end
    return (H₁Norm_Δu²/H₁Norm_ū²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end

"""
    Hₑ(ap; E, ν, h)

Relative energy-norm error of a Mindlin plate solution against an exact field:

    ‖e‖_E / ‖u‖_E,   ‖v‖_E² = ∫ (κ(v)·M(κ) + γ(v)·Q(γ)) dΩ

The exact field is supplied through `prescribe!`.

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio
- `h` : plate thickness

# Element requirements
`set∇𝝭!` and `set𝝭!` — reads `ξ[:𝝭]`, `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`, and the exact field.
"""
function Hₑ(aps::Vector{T}) where T<:AbstractElement
    HₑNorm_Δu² = 0.0
    HₑNorm_ū² = 0.0
    for ap in aps
        Δu², ū² = Hₑ(ap)
        HₑNorm_Δu² += Δu²
        HₑNorm_ū² += ū²
    end
    return (HₑNorm_Δu²/HₑNorm_ū²)^0.5
end

"""
    Hₑ(ap::T; E, ν, h)

Overload of `Hₑ` — see the documented sibling with the same name.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.𝑤`.
"""
function Hₑ(ap::T; E::Real, ν::Real, h::Real) where T<:AbstractElement
    Δu² = BigFloat(0.0)
    ū² = BigFloat(0.0)
    Dᵇᵢᵢᵢᵢ = E*h^3/12/(1-ν^2)
    Dᵇᵢᵢⱼⱼ = E*ν*h^3/12/(1-ν^2)
    Dᵇᵢⱼᵢⱼ = E*h^3/24/(1+ν)
    Dˢ =  5/6*h*E/2/(1+ν)
    for ξ in ap.𝓖
        κ̄₁₁ = -ξ.∂φ₁∂x
        κ̄₂₂ = -ξ.∂φ₂∂y
        κ̄₁₂ = -ξ.∂φ₁∂y-ξ.∂φ₂∂x
        γ̄₁ = ξ.∂w∂x - ξ.φ₁
        γ̄₂ = ξ.∂w∂y - ξ.φ₂
        𝑤 = ξ.𝑤

        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        κ₁₁ = 0.0
        κ₂₂ = 0.0
        κ₁₂ = 0.0
        γ₁ = 0.0
        γ₂ = 0.0
        for (i,xᵢ) in enumerate(ap.𝓒)
            γ₁ += B₁[i]*xᵢ.d - N[i]*xᵢ.d₁
            γ₂ += B₂[i]*xᵢ.d - N[i]*xᵢ.d₂
            κ₁₁ -= B₁[i]*xᵢ.d₁
            κ₂₂ -= B₂[i]*xᵢ.d₂
            κ₁₂ -= B₁[i]*xᵢ.d₂+B₂[i]*xᵢ.d₁
        end
        Δu² += 0.5*(Dᵇᵢᵢᵢᵢ*(κ₁₁-κ̄₁₁)^2+Dᵇᵢᵢᵢᵢ*(κ₂₂-κ̄₂₂)^2+2*Dᵇᵢᵢⱼⱼ*(κ₁₁-κ̄₁₁)*(κ₂₂-κ̄₂₂)+2*Dᵇᵢⱼᵢⱼ*(κ₁₂-κ̄₁₂)^2 + ((γ₁-γ̄₁)^2 + (γ₂-γ̄₂)^2)*Dˢ)*𝑤
        ū² += 0.5*(Dᵇᵢᵢᵢᵢ*κ̄₁₁^2+Dᵇᵢᵢᵢᵢ*κ̄₂₂^2+2*Dᵇᵢᵢⱼⱼ*κ̄₁₁*κ̄₂₂+2*Dᵇᵢⱼᵢⱼ*κ̄₁₂^2 + (γ̄₁^2 + γ̄₂^2)*Dˢ)*𝑤
    end
    return Δu², ū²
end

"""
    Hₑ(aps_w::Vector{Tʷ},aps_φ::Vector{Tᵠ},aps_Q::Vector{Tˢ})

Overload of `Hₑ` — see the documented sibling with the same name.


    Element-set overload of `Hₑ(ap_w, ap_φ, ap_Q)`, summed over the elements.

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.𝑤`.
"""
function Hₑ(aps_w::Vector{Tʷ},aps_φ::Vector{Tᵠ},aps_Q::Vector{Tˢ}) where {Tʷ,Tᵠ,Tˢ<:AbstractElement}
    HₑNorm_Δu² = 0.0
    HₑNorm_ū² = 0.0
    for (ap_w,ap_φ,ap_Q) in zip(aps_w,aps_φ,aps_Q)
        Δu², ū² = Hₑ(ap_w,ap_φ,ap_Q)
        HₑNorm_Δu² += Δu²
        HₑNorm_ū² += ū²
    end
    return (HₑNorm_Δu²/HₑNorm_ū²)^0.5
end

"""
    Hₑ(ap_w::Tʷ,ap_φ::Tᵠ,ap_Q::Tˢ)

Overload of `Hₑ` — see the documented sibling with the same name.


    ‖e‖_E / ‖u‖_E for the mixed `(w, φ, Q)` formulation:

        ‖v‖_E² = ∫ (κ(v)·M(κ) + γ(v)·Q(γ)) dΩ

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.𝑤`.
"""
function Hₑ(ap_w::Tʷ,ap_φ::Tᵠ,ap_Q::Tˢ) where {Tʷ,Tᵠ,Tˢ<:AbstractElement}
    Δu² = BigFloat(0.0)
    ū² = BigFloat(0.0)
    E = ap_φ.E
    ν = ap_φ.ν
    h = ap_φ.h
    Dᵇᵢᵢᵢᵢ = E*h^3/12/(1-ν^2)
    Dᵇᵢᵢⱼⱼ = E*ν*h^3/12/(1-ν^2)
    Dᵇᵢⱼᵢⱼ = E*h^3/24/(1+ν)
    Dˢ =  5/6*h*E/2/(1+ν)
    for (ξʷ,ξᵠ,ξˢ) in zip(ap_w.𝓖,ap_φ.𝓖,ap_Q.𝓖)
        κ̄₁₁ = -ξᵠ.∂φ₁∂x
        κ̄₂₂ = -ξᵠ.∂φ₂∂y
        κ̄₁₂ = -ξᵠ.∂φ₁∂y-ξᵠ.∂φ₂∂x
        γ̄₁ = ξʷ.∂w∂x - ξᵠ.φ₁
        γ̄₂ = ξʷ.∂w∂y - ξᵠ.φ₂
        q̄₁ = Dˢ*γ̄₁
        q̄₂ = Dˢ*γ̄₂
        q̃₁ = ξˢ.Q₁
        q̃₂ = ξˢ.Q₂
        𝑤 = ξᵠ.𝑤
        if abs(q̄₁-q̃₁)>1e3*eps() || abs(q̄₂-q̃₂)>1e3*eps()
            error("Some problem in exact solution! q̄₁: $q̄₁, q̃₁: $q̃₁, q̄₂: $q̄₂, q̃₂: $q̃₂")
        end

        Nʷ = ξʷ[:𝝭]
        Bʷ₁ = ξʷ[:∂𝝭∂x]
        Bʷ₂ = ξʷ[:∂𝝭∂y]
        Nᵠ = ξᵠ[:𝝭]
        Bᵠ₁ = ξᵠ[:∂𝝭∂x]
        Bᵠ₂ = ξᵠ[:∂𝝭∂y]
        Nˢ = ξˢ[:𝝭]
        κ₁₁ = 0.0
        κ₂₂ = 0.0
        κ₁₂ = 0.0
        γ₁ = 0.0
        γ₂ = 0.0
        q₁ = 0.0
        q₂ = 0.0
        t₁₀₀ = 0.0
        t₁₁₀ = 0.0
        t₁₀₁ = 0.0
        for (i,xᵢ) in enumerate(ap_w.𝓒)
            γ₁ += Bʷ₁[i]*xᵢ.d
            γ₂ += Bʷ₂[i]*xᵢ.d
        end
        for (i,xᵢ) in enumerate(ap_φ.𝓒)
            κ₁₁ -= Bᵠ₁[i]*xᵢ.d₁
            κ₂₂ -= Bᵠ₂[i]*xᵢ.d₂
            κ₁₂ -= Bᵠ₁[i]*xᵢ.d₂+Bᵠ₂[i]*xᵢ.d₁
            γ₁ -= Nᵠ[i]*xᵢ.d₁
            γ₂ -= Nᵠ[i]*xᵢ.d₂
            t₁₀₀ += Bᵠ₁[i]
            t₁₁₀ += Bᵠ₁[i]*xᵢ.t
            t₁₀₁ += Bᵠ₂[i]*xᵢ.t
        end
        # if abs(t₁₀₀) > 1e2*eps() || abs(t₁₁₀-1.0) > 1e2*eps() || abs(t₁₀₁) > 1e2*eps()
        #     error("Some problem in consistency condition! t₁₀₀: $t₁₀₀, t₁₁₀: $t₁₁₀, t₁₀₁: $t₁₀₁")
        # end
        for (i,xᵢ) in enumerate(ap_Q.𝓒)
            q₁ += Nˢ[i]*xᵢ.q₁
            q₂ += Nˢ[i]*xᵢ.q₂
        end

        Δu² += 0.5*(Dᵇᵢᵢᵢᵢ*(κ₁₁-κ̄₁₁)^2+Dᵇᵢᵢᵢᵢ*(κ₂₂-κ̄₂₂)^2+2*Dᵇᵢᵢⱼⱼ*(κ₁₁-κ̄₁₁)*(κ₂₂-κ̄₂₂)+2*Dᵇᵢⱼᵢⱼ*(κ₁₂-κ̄₁₂)^2 + ((γ₁-γ̄₁)^2 + (γ₂-γ̄₂)^2)*Dˢ)*𝑤
        ū² += 0.5*(Dᵇᵢᵢᵢᵢ*κ̄₁₁^2+Dᵇᵢᵢᵢᵢ*κ̄₂₂^2+2*Dᵇᵢᵢⱼⱼ*κ̄₁₁*κ̄₂₂+2*Dᵇᵢⱼᵢⱼ*κ̄₁₂^2 + (γ̄₁^2 + γ̄₂^2)*Dˢ)*𝑤
    end
    return Δu², ū²
end

end
