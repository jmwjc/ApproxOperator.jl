module Hamilton
    
using ..ApproxOperator: AbstractElement, Element2D, Element3D

"""
    ∫qmpdΩ(ap, k; ρA)

Mass-like term of a rod, `∫ ρA q p dΩ`.

In a space–time formulation the second coordinate is time; here the
operator acts on one spatial domain slice.

# Keyword arguments
- `ρA` : mass per unit length

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.𝑤`.
"""
function ∫qmpdΩ(ap::T,k::AbstractMatrix{Float64}; ρA::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += N[i]*ρA*N[j]*𝑤
            end    
        end
    end
end

"""
    ∫qkpdΩ(ap, k; EA)

Axial stiffness term of a rod, `∫ EA q p dΩ`.

# Keyword arguments
- `EA` : axial stiffness

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.𝑤`.
"""
function ∫qkpdΩ(ap::T,k::AbstractMatrix{Float64}; EA::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += N[i]*EA*N[j]*𝑤
            end    
        end
    end
end

"""
    ∫∫q̇mṗqkpdxdt(ap, k; m, kᶜ)

Space–time operator `∫∫ (m q̇ ṗ - k q p) dx dt`.

The time direction is the local coordinate `y`, so `q̇ = ξ[:∂𝝭∂y]`.

# Keyword arguments
- `m` : mass per unit length
- `kᶜ` : stiffness coefficient (`kᶜ`, because `k` is the matrix)

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫∫q̇mṗqkpdxdt(ap::T,k::AbstractMatrix{Float64}, kᶜ::Real; m::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B = ξ[:∂𝝭∂x]
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (B[i]*m*B[j] - N[i]*kᶜ*N[j])*𝑤
            end
        end
    end
end

"""
    ∫∫∇q∇pdxdt(ap, k; ρA, EA)
    ∫∫∇q∇pdxdt(a₁, a₂, k; ρA, EA)

Space–time wave operator `∫∫ (ρA q,ₜ p,ₜ - EA q,ₓ p,ₓ) dx dt`.

The two-element-set overload is for mixed formulations, with the
quadrature points of `a₁` and `a₂` zipped.

# Keyword arguments
- `ρA` : mass per unit length
- `EA` : axial stiffness

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫∫∇q∇pdxdt(ap::T,k::AbstractMatrix{Float64}; ρA::Real, EA::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        Bₓ = ξ[:∂𝝭∂x]
        Bₜ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (Bₜ[i]*ρA*Bₜ[j] - Bₓ[i]*EA*Bₓ[j])*𝑤
            end
        end
    end
end

"""
    ∫∫∇q∇pdΩdt(ap, k; c²)

Space–time operator `∫∫ (q,ₜ p,ₜ - c² ∇q·∇p) dΩ dt`.

Taking the speed already squared avoids recomputing `c^2` in the inner loop.

# Keyword arguments
- `c²` : wave speed squared

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫∫∇q∇pdΩdt(ap::T,k::AbstractMatrix{Float64}; c²::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        Bₜ = ξ[:∂𝝭∂z]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (Bₜ[i]*Bₜ[j] - c²*(B₁[i]*B₁[j] + B₂[i]*B₂[j]))*𝑤
            end
        end
    end
end

"""
    ∫∫∇q∇pdxdt(ap, k; ρA, EA)
    ∫∫∇q∇pdxdt(a₁, a₂, k; ρA, EA)

Space–time wave operator `∫∫ (ρA q,ₜ p,ₜ - EA q,ₓ p,ₓ) dx dt`.

The two-element-set overload is for mixed formulations, with the
quadrature points of `a₁` and `a₂` zipped.

# Keyword arguments
- `ρA` : mass per unit length
- `EA` : axial stiffness

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫∫∇q∇pdxdt(a₁::T,a₂::S,k::AbstractMatrix{Float64}; ρA::Real, EA::Real) where {T<:AbstractElement,S<:AbstractElement}
    𝓒₁ = a₁.𝓒; 𝓖₁ = a₁.𝓖
    𝓒₂ = a₂.𝓒; 𝓖₂ = a₂.𝓖
    for (ξ₁,ξ₂) in zip(𝓖₁,𝓖₂)
        B̄ₓ = ξ₁[:∂𝝭∂x]
        B̄ₜ = ξ₁[:∂𝝭∂y]
        Bₓ = ξ₂[:∂𝝭∂x]
        Bₜ = ξ₂[:∂𝝭∂y]
        𝑤 = ξ₂.𝑤
        for (i,xᵢ) in enumerate(𝓒₁)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒₂)
                J = xⱼ.𝐼
                k[I,J] += (-B̄ₜ[i]*ρA*Bₜ[j] + B̄ₓ[i]*EA*Bₓ[j])*𝑤
            end
        end
    end
end

"""
    ∫q∇𝑛pds(ap, k, f; ρA, EA)

Space–time boundary term `∫ q ∂ₙp ds` with the wave coefficients.

`f` is the output vector.

# Keyword arguments
- `ρA` : mass per unit length
- `EA` : axial stiffness

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫q∇𝑛pds(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}; ρA::Real, EA::Real) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        Bₓ = ξ[:∂𝝭∂x]
        Bₜ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        nₓ = ξ.n₁
        nₜ = ξ.n₂
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] -= (-N[i]*ρA*Bₜ[j]*nₜ + N[i]*EA*Bₓ[j]*nₓ)*𝑤
            end
            # f[I] -= (-ρA*Bₜ[i]*nₜ + EA*Bₓ[i]*nₓ)*g*𝑤
        end
    end
end

"""
    ∫∫αqṗdxdt(ap, k; α)

Space–time stabilisation term `α ∫∫ q ṗ dx dt`.

# Keyword arguments
- `α` : stabilisation coefficient

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫∫αqṗdxdt(ap::T,k::AbstractMatrix{Float64}; α::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B = ξ[:∂𝝭∂y]
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += α*(B[i]*N[j] + N[i]*B[j])*𝑤
            end
        end
    end
end


"""
    stabilization_bar_LSG(ap, k; ρA, α, EA)
    stabilization_bar_LSG_Γ(ap, k; ρA, α, EA)

Least-squares Galerkin (LSG) stabilisation for the space–time rod.


    a_h(u,v) = a(u,v) + α Σ_e ∫_Ω ℒu · ℒv dΩ

with `ℒ` the residual of the strong form, `ℒu = ρA u,ₜₜ - EA u,ₓₓ`.

Adds the residual of the strong form, weighted by `α`, to the weak form.
The `_Γ` variant applies the same construction on the boundary.

# Keyword arguments
- `ρA` : mass per unit length
- `α` : stabilisation coefficient
- `EA` : axial stiffness

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function stabilization_bar_LSG(ap::T,k::AbstractMatrix{Float64}; ρA::Real, α::Real, EA::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        Bₓₓ = ξ[:∂²𝝭∂x²]
        Bₜₜ = ξ[:∂²𝝭∂y²]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += α*(ρA*Bₜₜ[i] - EA*Bₓₓ[i])*(ρA*Bₜₜ[j] - EA*Bₓₓ[j])*𝑤
            end
        end
    end
end

"""
    stabilization_bar_LSG_Γ(ap::T,k::AbstractMatrix{Float64}; ρA, α, EA)

LSG stabilisation applied on the boundary of the space–time rod
(the `_Γ` counterpart of `stabilization_bar_LSG`).


    Boundary counterpart of `stabilization_bar_LSG`:

        a_Γ(u,v) = α Σ_e ∫_Γ (ℒu)(ℒv) dΓ

# Element requirements
`set∇𝝭!` — reads `ξ[:𝝭]`, gradients, `ξ.n₁`, `ξ.n₂`, `ξ.𝑤`.
"""
function stabilization_bar_LSG_Γ(ap::T,k::AbstractMatrix{Float64}; ρA::Real, α::Real, EA::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        Bₓ = ξ[:∂𝝭∂x]
        Bₜ = ξ[:∂𝝭∂y]
        nₓ = ξ.n₁
        nₜ = ξ.n₂
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += α*(ρA*Bₜ[i]*nₜ - EA*Bₓ[i]*nₓ)*(ρA*Bₜ[j]*nₜ - EA*Bₓ[j]*nₓ)*𝑤
            end
        end
    end
end

"""
    truncation_error(ap, fₓ, fₜ, fₓₓ, fₜₜ; c)

Truncation-error measures of the space–time discretisation.


    For the exact characteristic direction `ξ = x - c t`, the four defect measures

        fₓ  += BₓᵢBₓⱼ (Δx - cΔt),      fₜ  += BₜᵢBₜⱼ (Δx - cΔt)
        fₓₓ += BₓᵢBₓⱼ (Δx - cΔt)²,     fₜₜ += BₜᵢBₜⱼ (Δx - cΔt)²

Accumulates four defect measures into the supplied vectors, using the exact
characteristic direction `x - c t`.

# Keyword arguments
- `c` : wave speed

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`, `ξ.x`, `ξ.y`.
"""
function truncation_error(ap::T,fₓ::AbstractVector{Float64},fₜ::AbstractVector{Float64},fₓₓ::AbstractVector{Float64},fₜₜ::AbstractVector{Float64}; c::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        Bₓ = ξ[:∂𝝭∂x]
        Bₜ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        x = ξ.x
        t = ξ.y
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                fₓ[I] += Bₓ[i]*Bₓ[j]*(xⱼ.x-xᵢ.x - c*(xⱼ.y-xᵢ.y))*𝑤
                fₜ[I] += Bₜ[i]*Bₜ[j]*(xⱼ.x-xᵢ.x - c*(xⱼ.y-xᵢ.y))*𝑤
                fₓₓ[I] += Bₓ[i]*Bₓ[j]*(xⱼ.x-xᵢ.x - c*(xⱼ.y-xᵢ.y))^2*𝑤
                fₜₜ[I] += Bₜ[i]*Bₜ[j]*(xⱼ.x-xᵢ.x - c*(xⱼ.y-xᵢ.y))^2*𝑤
            end
        end
    end
end

"""
    truncation_error(ap, fₓ, fₜ, fₓₓ, fₜₜ; c)

Truncation-error measures of the space–time discretisation.


    Element-set overload of `truncation_error`, accumulating the same four
    measures into freshly allocated vectors of length `nₚ`.

Accumulates four defect measures into the supplied vectors, using the exact
characteristic direction `x - c t`.

# Keyword arguments
- `c` : wave speed

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`, `ξ.x`, `ξ.y`.
"""
function truncation_error(aps::Vector{T},nₚ::Int) where T<:AbstractElement
    fₓ = zeros(nₚ)
    fₜ = zeros(nₚ)
    fₓₓ = zeros(nₚ)
    fₜₜ = zeros(nₚ)
    for ap in aps
        truncation_error(ap,fₓ,fₜ,fₓₓ,fₜₜ)
    end
    return fₓ,fₜ,fₓₓ,fₜₜ
end

"""
    ∫pudΩ(a₁::T,a₂::S,k::AbstractMatrix{Float64})

Element-set overload of `∫pudΩ`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function ∫pudΩ(a₁::T,a₂::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒₁ = a₁.𝓒; 𝓖₁ = a₁.𝓖
    𝓒₂ = a₂.𝓒; 𝓖₂ = a₂.𝓖
    for (ξ₁,ξ₂) in zip(𝓖₁,𝓖₂)
        Bₜ = ξ₂[:∂𝝭∂y]
        N = ξ₁[:𝝭]
        𝑤 = ξ₂.𝑤
        for (i,xᵢ) in enumerate(𝓒₁)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒₂)
                J = xⱼ.𝐼
                k[I,J] += N[i]*Bₜ[j]*𝑤
            end
        end
    end
end

"""
    ∫ppdΩ(ap, k; ρA)

Mass matrix of the pressure variable, `∫ ρA p p dΩ`.

# Keyword arguments
- `ρA` : mass per unit length

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.𝑤`.
"""
function ∫ppdΩ(ap::T,k::AbstractMatrix{Float64}; ρA::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += -(1/ρA)N[i]*N[j]*𝑤
            end
        end
    end
end

"""
    ∫uudΩ(ap, k; EA)

Stiffness matrix of the displacement variable, `∫ EA u u dΩ`.

# Keyword arguments
- `EA` : axial stiffness

# Element requirements
`set𝝭!` — reads `ξ[:𝝭]`, `ξ.𝑤`.
"""
function ∫uudΩ(ap::T,k::AbstractMatrix{Float64}; EA::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        Bₓ = ξ[:∂𝝭∂x]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += -Bₓ[i]*EA*Bₓ[j]*𝑤
            end
        end
    end
end

"""
    Hₑ(ap; ρA, EA)        # T<:Element2D

Relative energy-norm error of the 2-D space–time rod.

    ‖e‖_E / ‖u‖_E,   ‖v‖_E² = ∫ (ρA v,ₜ² + EA v,ₓ²) dx dt

Supplied with the exact field through `prescribe!`.

# Keyword arguments
- `ρA` : mass per unit length
- `EA` : axial stiffness

# Element requirements
`set∇𝝭!` and `set𝝭!` — reads `ξ[:𝝭]`, `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`, and the
exact field `ξ.u`, `ξ.∂u∂x`, ...
"""
function Hₑ(ap::T; ρA::Real, EA::Real) where T<:Element2D
    Δu²= 0
    Δδu² = 0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        Bₓ = ξ[:∂𝝭∂x]
        Bₜ = ξ[:∂𝝭∂y]
        ∂ūᵢ∂x = ξ.∂u∂x
        ∂ūᵢ∂t = ξ.∂u∂t
        ∂δūᵢ∂x = ξ.∂δu∂x
        ∂δūᵢ∂t = ξ.∂δu∂t
        ∂uᵢ∂x = 0
        ∂uᵢ∂t = 0
        ∂δuᵢ∂x = 0
        ∂δuᵢ∂t = 0
        for (i,xᵢ) in enumerate(ap.𝓒)
            ∂uᵢ∂x += Bₓ[i]*xᵢ.d
            ∂uᵢ∂t += Bₜ[i]*xᵢ.d
            ∂δuᵢ∂x += Bₓ[i]*xᵢ.δd
            ∂δuᵢ∂t += Bₜ[i]*xᵢ.δd
        end
        Δu² += 0.5*abs(ρA*(∂uᵢ∂t - ∂ūᵢ∂t)^2 - EA*(∂uᵢ∂x - ∂ūᵢ∂x)^2)*𝑤
        Δδu² += 0.5*abs(ρA*(∂δuᵢ∂t - ∂δūᵢ∂t)^2 - EA*(∂δuᵢ∂x - ∂δūᵢ∂x)^2)*𝑤
    end
    return Δu², Δδu²
end

"""
    Hₑ(ap; c²)            # T<:Element3D

Relative energy-norm error of the 3-D wave equation.

    ‖v‖_E² = ∫ (v,ₜ² + c² ∇v·∇v) dΩ dt

# Keyword arguments
- `c²` : wave speed squared

# Element requirements
`set∇𝝭!` and `set𝝭!` — reads `ξ[:𝝭]`, the gradients, `ξ.𝑤`, and the exact field.
"""
function Hₑ(ap::T; c²::Real) where T<:Element3D
    Δu²= 0
    Δδu² = 0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        Bₜ = ξ[:∂𝝭∂z]
        ∂ūᵢ∂x = ξ.∂u∂x
        ∂ūᵢ∂y = ξ.∂u∂y
        ∂ūᵢ∂t = ξ.∂u∂t
        ∂δūᵢ∂x = ξ.∂δu∂x
        ∂δūᵢ∂y = ξ.∂δu∂y
        ∂δūᵢ∂t = ξ.∂δu∂t
        ∂uᵢ∂x = 0
        ∂uᵢ∂y = 0
        ∂uᵢ∂t = 0
        ∂δuᵢ∂x = 0
        ∂δuᵢ∂y = 0
        ∂δuᵢ∂t = 0
        for (i,xᵢ) in enumerate(ap.𝓒)
            ∂uᵢ∂x += B₁[i]*xᵢ.d
            ∂uᵢ∂y += B₂[i]*xᵢ.d
            ∂uᵢ∂t += Bₜ[i]*xᵢ.d
            ∂δuᵢ∂x += B₁[i]*xᵢ.δd
            ∂δuᵢ∂y += B₂[i]*xᵢ.δd
            ∂δuᵢ∂t += Bₜ[i]*xᵢ.δd
        end
        Δu² += 0.5*abs((∂uᵢ∂t - ∂ūᵢ∂t)^2 - c²*((∂uᵢ∂x - ∂ūᵢ∂x)^2 + (∂uᵢ∂y - ∂ūᵢ∂y)^2))*𝑤
        Δδu² += 0.5*abs((∂δuᵢ∂t - ∂δūᵢ∂t)^2 - c²*((∂δuᵢ∂x - ∂δūᵢ∂x)^2 + (∂δuᵢ∂y - ∂δūᵢ∂y)^2))*𝑤
    end
    return Δu², Δδu²
end

"""
    Hₑ(aps::Vector{T})

Element-set overload of `Hₑ`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.


    ‖e‖_E / ‖u‖_E,   ‖v‖_E² = Σ_e ∫ (ρA v,ₜ² + EA v,ₓ²) dx dt

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function Hₑ(aps::Vector{T}) where T<:AbstractElement
    HₑNorm_Δu²= 0.
    HₑNorm_Δδu² = 0.
    for ap in aps
        Δu², Δδu² = Hₑ(ap)
        HₑNorm_Δu² += Δu²
        HₑNorm_Δδu² += Δδu²
    end
    return HₑNorm_Δu²^0.5, HₑNorm_Δδu²^0.5
end

end