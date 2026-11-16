
"""
    ∫v²uₓuₓdx(ap::T;k::AbstractMatrix{Float64})

Element-set overload of `∫v²uₓuₓdx`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function ∫v²uₓuₓdx(ap::T, k::AbstractMatrix{Float64}; η::Real = 1.0e-8) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        B = ξ[:∂𝝭∂x]
        v = sum(N[i]*xᵢ.v for (i,xᵢ) in enumerate(𝓒))
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (v^2+η)*B[i]*B[j]*𝑤
            end
        end
    end
end

"""
    ∫vₓvₓvvdx(ap, k, f; l)
    ∫vₓvₓvvdx_hard_device(ap, k, f; l)

1-D phase-field damage equation, AT2 model.

    k[I,J] = ∫ (k_c(2l Bᵢ Bⱼ + Nᵢ Nⱼ/(2l)) + ℋ Nᵢ Nⱼ) w
    f[I]   = ∫ k_c Nᵢ/(2l) w

`l` is the length-scale parameter controlling the diffuse crack width and `k_c`
the fracture toughness (read as `ξ.k`). `ℋ` is the history field driving the
damage, updated by the caller each step. The `_hard_device` variant omits the
history term and uses the current strain directly.

# Keyword arguments
- `l` : phase-field length scale

# Element requirements
`set∇𝝭!` — reads `ξ[:𝝭]`, `ξ[:∂𝝭∂x]`, `ξ.ℋ`, `ξ.k`, `ξ.𝑤`.
"""
function ∫vₓvₓvvdx_hard_device(ap::T, k::AbstractMatrix{Float64}, f::AbstractVector{Float64}; l::Real, η::Real = 1.0e-8, kc::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        B = ξ[:∂𝝭∂x]
        ℋ = ξ.ℋ
        ε = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            ε += B[i]*xᵢ.u
        end
        ℋₜ = max(ℋ,(ε-1)^2)
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (kc*(2*l*B[i]*B[j] + N[i]*N[j]/2/l) + ℋₜ*N[i]*N[j])*𝑤
            end
            f[I] += N[i]*kc/2/l*𝑤
        end
    end
end

"""
    ∫vₓvₓvvdx(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}, l)

Element-set overload of `∫vₓvₓvvdx`: the same integrand, summed over
a vector of elements. See the element-level method for the formula.

# Element requirements
`set∇𝝭!` — reads the quadrature-point fields of the operator body.
"""
function ∫vₓvₓvvdx(ap::T, k::AbstractMatrix{Float64}, f::AbstractVector{Float64}; l::Real, η::Real = 1.0e-8, kc::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
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
        ℋₜ = max(ℋ,(v^2+η)*ε^2)
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (kc*(2*l*B[i]*B[j] + N[i]*N[j]/2/l) + ℋₜ*N[i]*N[j])*𝑤
            end
            f[I] += N[i]*(kc/2/l - η*ℋₜ)*𝑤
        end
    end
end

"""
    UPDATE_PFM_1D(ap::T)

Advance the 1-D phase-field history variable: `ℋ ← max(ℋ, (ε-1)²)` at every
quadrature point. Stateful — call once per load step before assembling.


    ℋ ← max(ℋ, (ε - 1)²)        at every quadrature point

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ.ℋ`, and the nodal displacement.
"""
function UPDATE_PFM_1D(ap::T; η::Real = 1.0e-8) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        B = ξ[:∂𝝭∂x]
        ℋ = ξ.ℋ
        ε = 0.0
        v = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            ε += B[i]*xᵢ.u
            v += N[i]*xᵢ.v
        end
        ξ.ℋ = max(ℋ,(v^2+η)*ε^2)
    end
end
"""
    ∫∫∇v∇vvvdxdy(ap, k, f; l)

2-D phase-field damage equation, AT2 model.

Two-dimensional counterpart of [`∫vₓvₓvvdx`](@ref):

    k[I,J] += (k_c(2l ∇Nᵢ·∇Nⱼ + Nᵢ Nⱼ/(2l)) + ℋ Nᵢ Nⱼ) w

# Keyword arguments
- `l` : phase-field length scale

# Element requirements
`set∇𝝭!` — reads `ξ[:𝝭]`, `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.ℋ`, `ξ.k`, `ξ.𝑤`.
"""
function ∫∫∇v∇vvvdxdy(ap::T, k::AbstractMatrix{Float64}, f::AbstractVector{Float64}; l::Real, η::Real = 1.0e-8, kc::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        ℋ = ξ.ℋ
        ε₁₁ = 0.0
        ε₂₂ = 0.0
        ε₁₂ = 0.0
        σ₁₁ = ξ.σ₁₁
        σ₂₂ = ξ.σ₂₂
        σ₁₂ = ξ.σ₁₂
        for (i,xᵢ) in enumerate(𝓒)
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += B₁[i]*xᵢ.d₂ + B₂[i]*xᵢ.d₁
        end
        ℋₜ = max(ℋ,(ε₁₁*σ₁₁ + ε₂₂*σ₂₂ + ε₁₂*σ₁₂))
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (kc*(2*l*(B₁[i]*B₁[j] + B₂[i]*B₂[j]) + N[i]*N[j]/2/l) + ℋₜ*N[i]*N[j])*𝑤
            end
            f[I] += N[i]*kc/2/l*𝑤
        end
    end
end

"""
    UPDATE_PFM_2D(ap::T)

Advance the 2-D phase-field history variable from the strain energy,
`ℋ ← max(ℋ, ψ⁺(ε))`. Stateful — call once per load step.


    ℋ ← max(ℋ, ψ⁺(ε))           at every quadrature point

# Element requirements
`set∇𝝭!` — reads the gradients, `ξ.ℋ`, and the nodal displacement.
"""
function UPDATE_PFM_2D(ap::T; η::Real = 1.0e-8) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂x]
        σ₁₁ = ξ.σ₁₁
        σ₂₂ = ξ.σ₂₂
        σ₁₂ = ξ.σ₁₂
        v = 0.0
        ε₁₁ = 0.0
        ε₂₂ = 0.0
        ε₁₂ = 0.0
        ℋ = ξ.ℋ
        for (i,xᵢ) in enumerate(𝓒)
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += B₁[i]*xᵢ.d₂ + B₂[i]*xᵢ.d₁
            v += N[i]*xᵢ.v
        end
        ξ.ℋ = max(ℋ,ε₁₁*σ₁₁ + ε₂₂*σ₂₂ + ε₁₂*σ₁₂)
    end
end

"""
    CRACK_NORMAL(aps::Vector{T},nodes::Vector{N},v::Vector{Float64})

Compute the crack normal at each quadrature point from the damage gradient,
`n = -∇v/|∇v|`, and write it back into the quadrature point.


    𝒏 ← -∇v / |∇v|              at every quadrature point

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, and the nodal damage.
"""
function CRACK_NORMAL(aps::Vector{T},nodes::Vector{N},v::Vector{Float64}) where {T<:AbstractElement,N<:Node}
    Γtmp = findall(x->x<0.02,v)
    Γfinal = Int[]
    # 循环，当Γtmp不为空时
    if ~isempty(Γtmp)
        while ~isempty(Γtmp)
            # 找到identity中v[Γtmp]的最小值，并将其索引N1记录
            _,N1 = findmin(identity,v[Γtmp])
            # 将Γtmp[N1]添加到Γfinal中
            push!(Γfinal,Γtmp[N1])
            # 记录Γtmp[N1]的x坐标和y坐标
            x₁ = nodes[Γtmp[N1]].x
            y₁ = nodes[Γtmp[N1]].y
            # 遍历Γtmp中的每一个元素
            for index in Γtmp
                # 记录当前元素对应的节点
                node = nodes[index]
                # 记录当前节点对应的x坐标和y坐标
                x₂ = node.x
                y₂ = node.y
                # 计算当前节点和Γtmp[N1]之间的距离
                Δ = ((x₁-x₂)^2 + (y₁-y₂)^2)^0.5
                # 如果距离小于l，则将当前节点从Γtmp中移除
                Δ < l ? setdiff!(Γtmp,index) : nothing
            end
        end
        sort!(Γfinal, by = i->nodes[i].x)
        n₁ = zeros(length(Γfinal)-1)
        n₂ = zeros(length(Γfinal)-1)
        for i in 1:length(Γfinal)-1
            i₁ = Γfinal[i]
            i₂ = Γfinal[i+1]
            x₁ = nodes[i₁].x
            y₁ = nodes[i₁].y
            x₂ = nodes[i₂].x
            y₂ = nodes[i₂].y
            𝐿 = ((x₁-x₂)^2 + (y₁-y₂)^2)^0.5
            n₁[i] = (y₂-y₁)/𝐿
            n₂[i] = (x₁-x₂)/𝐿
        end
        for ap in aps
            𝓖 = ap.𝓖
            for ξ in 𝓖
                𝐿 = Inf
                x₀ = ξ.x
                y₀ = ξ.y
                for i in 1:length(Γfinal)-1
                    i₁ = Γfinal[i]
                    i₂ = Γfinal[i+1]
                    x₁ = nodes[i₁].x
                    y₁ = nodes[i₁].y
                    x₂ = nodes[i₂].x
                    y₂ = nodes[i₂].y
                    A = y₂-y₁
                    B = x₂-x₁
                    C = y₁*(x₂-x₁) - x₁*(y₂-y₁)
                    𝐿tmp = abs(A*x₀+B*y₀+C)/(A^2+B^2)^0.5
                    if 𝐿tmp < 𝐿
                        ξ.n₁ = n₁[i]
                        ξ.n₂ = n₂[i]
                        𝐿 = 𝐿tmp
                    end
                end
            end
        end
    end
end

"""
    ∫vₓvₓvvdx(aps, k, f; l, η)

Element-set overload of [`∫vₓvₓvvdx`](@ref): the same integrand, summed over
the elements of `aps`.
"""
function ∫vₓvₓvvdx(aps::Vector{T}, k::AbstractMatrix{Float64},
                    f::AbstractVector{Float64}; l::Real, η::Real = 1.0e-8,
                    kc::Real) where T<:AbstractElement
    for ap in aps
        ∫vₓvₓvvdx(ap, k, f; l = l, η = η, kc = kc)
    end
    return k, f
end
