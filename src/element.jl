"""
Element{T}<:AbstractElement{T}
"""
struct Element{T}<:AbstractElement{T}
    𝓒::Tuple{Int,Int,Vector{Node{(:𝐼,),1}}}
    𝓖::Tuple{Int,Int,Vector{Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}}}
end

get𝓒(ap::T) where T<:AbstractElement = (𝓒[3][𝓒[1]+i] for i in 1:𝓒[2])
get𝓖(ap::T) where T<:AbstractElement = (𝓖[3][𝓖[1]+i] for i in 1:𝓖[2])

function Base.getproperty(ap::T,s::Symbol) where T<:AbstractElement
    𝓖 = getfield(ap,:𝓖)
    ξ = 𝓖[3][𝓖[1]+1]
    return getproperty(ξ,s)
end

function Base.setproperty!(ap::T,s::Symbol) where T<:AbstractElement
    𝓖 = getfield(ap,:𝓖)
    ξ = 𝓖[3][𝓖[1]+1]
    setproperty!(ξ,s)
end