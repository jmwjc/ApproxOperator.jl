"""
Operator

A weak form is a plain Julia function

    form(ap, k[, f]; <material constants>)

assembled over a set of elements with the `Pair` sugar

    (form => elements)(k[, f]; <material constants>)

The `Vector{Pair}` overloads apply a list of operators in sequence, which is how
a global system is assembled from several weak terms:

    [
        ∫κᵢⱼMᵢⱼdΩ       => elements_bending,
        ∫αwwdΓ          => elements_penalty,
    ](k; D = D, ν = 0.3, α = 1e7)

Material constants are keyword arguments so that they appear in the operator's
own signature; anything else the weak form needs (nodal values, boundary data)
is supplied through `prescribe!` and read from the quadrature point `ξ`.
"""

# A weak form that takes a stiffness matrix and, optionally, a load vector.
# Keywords are forwarded verbatim to `form`, so every material constant stays a
# keyword all the way from the call site to the operator body.
function (op::Pair{F,Vector{T}})(k::AbstractMatrix, f::AbstractVector;
                                 kw...) where {F<:Function,T<:AbstractElement}
    form, elms = op
    for elm in elms
        form(elm, k, f; kw...)
    end
    return k, f
end

function (op::Pair{F,Vector{T}})(k::AbstractMatrix; kw...) where {F<:Function,T<:AbstractElement}
    form, elms = op
    for elm in elms
        form(elm, k; kw...)
    end
    return k
end

function (op::Pair{F,Vector{T}})(f::AbstractVector; kw...) where {F<:Function,T<:AbstractElement}
    form, elms = op
    for elm in elms
        form(elm, f; kw...)
    end
    return f
end

# Two-element forms: the two element sets are zipped, as in mixed methods where
# each velocity element is paired with one pressure element.
function (op::Pair{F,Tuple{Vector{T},Vector{S}}})(k::AbstractMatrix, f::AbstractVector;
                                                  kw...) where {F<:Function,T<:AbstractElement,S<:AbstractElement}
    form, elms = op
    for (a,b) in zip(elms...)
        form(a, b, k, f; kw...)
    end
    return k, f
end

function (op::Pair{F,Tuple{Vector{T},Vector{S}}})(k::AbstractMatrix;
                                                  kw...) where {F<:Function,T<:AbstractElement,S<:AbstractElement}
    form, elms = op
    for (a,b) in zip(elms...)
        form(a, b, k; kw...)
    end
    return k
end

function (op::Pair{F,Tuple{Vector{T},Vector{S}}})(f::AbstractVector;
                                                  kw...) where {F<:Function,T<:AbstractElement,S<:AbstractElement}
    form, elms = op
    for (a,b) in zip(elms...)
        form(a, b, f; kw...)
    end
    return f
end

# A list of operators is applied in order, so a global system is assembled by
# one call:  (ops::Vector{Pair})(k; kw...)
function (ops::Vector{Pair{F,Vector{T}}})(k::AbstractMatrix, f::AbstractVector;
                                          kw...) where {F<:Function,T<:AbstractElement}
    for op in ops
        op(k, f; kw...)
    end
end

function (ops::Vector{Pair{F,Vector{T}}})(k::AbstractMatrix;
                                          kw...) where {F<:Function,T<:AbstractElement}
    for op in ops
        op(k; kw...)
    end
end

function (ops::Vector{Pair{F,Vector{T}}})(f::AbstractVector;
                                          kw...) where {F<:Function,T<:AbstractElement}
    for op in ops
        op(f; kw...)
    end
end

function (ops::Vector{Pair{F,Tuple{Vector{T},Vector{S}}}})(k::AbstractMatrix, f::AbstractVector;
                                                           kw...) where {F<:Function,T<:AbstractElement,S<:AbstractElement}
    for op in ops
        op(k, f; kw...)
    end
end

function (ops::Vector{Pair{F,Tuple{Vector{T},Vector{S}}}})(k::AbstractMatrix;
                                                           kw...) where {F<:Function,T<:AbstractElement,S<:AbstractElement}
    for op in ops
        op(k; kw...)
    end
end

function (ops::Vector{Pair{F,Tuple{Vector{T},Vector{S}}}})(f::AbstractVector;
                                                           kw...) where {F<:Function,T<:AbstractElement,S<:AbstractElement}
    for op in ops
        op(f; kw...)
    end
end

function getDOFs(aps::Vector{T}) where T<:AbstractElement
    𝐼 = Set{Int}()
    for ap in aps
        for xᵢ in ap.𝓒
            push!(𝐼,xᵢ.𝐼)
        end
    end
    return 𝐼
end

"""
    MATERIAL_PARAMETERS

Symbols that count as material or model constants, i.e. values the user chooses
and that are the same throughout the domain. Operators in `src/operation/` read
these as keyword arguments; everything else an operator needs (boundary data,
prescribed solution fields, mesh geometry) comes from the quadrature point `ξ`
via `prescribe!` or from `importmsh`.

Deliberately **excluded**:

* load / boundary functions — `:g`, `:t`, `:b`, `:q`, `:g₁`, `:g₂`, `:g₃`
* nodal variables — `:u`, `:u₁`, `:u₂`, `:u₃`, `:v`, `:w`, `:θ`, `:φ`, `:φ₁`, `:φ₂`
* mesh geometry — `:n₁`, `:n₂`, `:s₁`, `:s₂`, `:a¹¹`, `:𝐽`, `:x`, `:y`, `:z`
* quadrature-point state updated in place — `:σ₁₁`, `:σ₂₂`, `:σ₁₂`, `:εᵖ₁₁`, `:αₙ`
* exact-solution data used only by tests — `:p₁`, `:T`, `:∂u∂x`
"""
const MATERIAL_PARAMETERS = Set{Symbol}([
    :E, :ν, :h,            # linear elasticity / plates / shells
    :D,                    # Kirchhoff plate bending stiffness
    :ρ, :ρA, :m,           # inertia
    :μ, :λ, :σy, :K, :tol, # plasticity
    :α,                    # penalty coefficient here; stiffness scaling in beams
    :k, :c, :c²,           # conductivity here; wave number / speed in beams
    :EA, :EI, :kGA, :R,    # beams
    :l,                    # phase-field length scale
])

"""
    ismaterial(name::Symbol)

True when `name` is a material/model constant, i.e. belongs in an operator's
keyword list rather than being read from the quadrature point.
"""
ismaterial(name::Symbol) = name in MATERIAL_PARAMETERS
