"""
Operator
"""

function (ops::Vector{Pair{F,Vector{T}}})(k::AbstractMatrix,f::AbstractVector) where {F<:Function,T<:AbstractElement}
    for op in ops
        op(k,f)
    end
end

function (op::Pair{F,Vector{T}})(k::AbstractMatrix,f::AbstractVector) where {F<:Function,T<:AbstractElement}
    form, elms = op
    for elm in elms
        form(elm,k,f)
    end
    return k,f
end

function (ops::Vector{Pair{F,Vector{T}}})(k::AbstractMatrix) where {F<:Function,T<:AbstractElement}
    for op in ops
        op(k)
    end
end

function (op::Pair{F,Vector{T}})(k::AbstractMatrix) where {F<:Function,T<:AbstractElement}
    form, elms = op
    for elm in elms
        form(elm,k)
    end
    return k
end

function (ops::Vector{Pair{F,Vector{T}}})(f::AbstractVector) where {F<:Function,T<:AbstractElement}
    for op in ops
        op(f)
    end
end

function (op::Pair{F,Vector{T}})(f::AbstractVector) where {F<:Function,T<:AbstractElement}
    form, elms = op
    for elm in elms
        form(elm,f)
    end
    return f
end



function (op::Pair{F,Vector{T}})(f, E_cons_hist, D_hist) where {F<:Function,T<:AbstractElement}
    form, elms = op
        for (e, elm) in enumerate(elms)
            form(e, elm, f, E_cons_hist, D_hist)
        end
    return f, D_hist
end




function (op::Pair{F,Vector{T}})(k,D_hist) where {F<:Function,T<:AbstractElement}
    form, elms = op

        for (e, elm) in enumerate(elms)
            form(e, elm, k, D_hist)
        end

    return k, D_hist
end

function (ops::Vector{Pair{F,Tuple{Vector{T},Vector{S}}}})(k::AbstractMatrix,f::AbstractVector) where {F<:Function,T<:AbstractElement,S<:AbstractElement}
    for op in ops
        op(k,f)
    end
end

function (op::Pair{F,Tuple{Vector{T},Vector{S}}})(k::AbstractMatrix,f::AbstractVector) where {F<:Function,T<:AbstractElement,S<:AbstractElement}
    form, elms = op
    for (a,b) in zip(elms...)
        form(a,b,k,f)
    end
    return k,f
end
function (op::Pair{F,Tuple{Vector{T},Vector{S}}})(k::AbstractMatrix,fₛ::AbstractVector,fᵤ::AbstractVector) where {F<:Function,T<:AbstractElement,S<:AbstractElement}
    form, elms = op
    for (a,b) in zip(elms...)
        form(a,b,k,fₛ,fᵤ)
    end
    return k,fₛ,fᵤ
end
function (ops::Vector{Pair{F,Tuple{Vector{T},Vector{S}}}})(k::AbstractMatrix) where {F<:Function,T<:AbstractElement,S<:AbstractElement}
    for op in ops
        op(k)
    end
end

function (op::Pair{F,Tuple{Vector{T},Vector{S}}})(k::AbstractMatrix) where {F<:Function,T<:AbstractElement,S<:AbstractElement}
    form, elms = op
    for (a,b) in zip(elms...)
        form(a,b,k)
    end
    return k
end

function (ops::Vector{Pair{F,Tuple{Vector{T},Vector{S}}}})(f::AbstractVector) where {F<:Function,T<:AbstractElement,S<:AbstractElement}
    for op in ops
        op(f)
    end
end

function (op::Pair{F,Tuple{Vector{T},Vector{S}}})(f::AbstractVector) where {F<:Function,T<:AbstractElement,S<:AbstractElement}
    form, elms = op
    for (a,b) in zip(elms...)
        form(a,b,f)
    end
    return f
end

function (op::Pair{F,Tuple{Vector{T},Vector{S}}})(f1::AbstractVector,f2::AbstractVector,k1::AbstractMatrix,k2::AbstractMatrix,k3::AbstractMatrix) where {F<:Function,T<:AbstractElement,S<:AbstractElement}
    form, elms = op
    for (a,b) in zip(elms...)
        form(a,b,f1,f2,k1,k2,k3)
    end
    return f1,f2,k1,k2,k3
end



