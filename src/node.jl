"""
Returen value (RV)

This is a simple struct
"""
struct RV
    i::Int
    v::Vector{Float64}
end

getindex(r::RV,i::Int) = r.v[r.i+i]
function setindex!(r::RV,x::Float64,i::Int)
    r.v[r.i+i] = x
end

"""
Node
"""
struct Node{T,N}
    index::NTuple{N,Int}
    data::Dict{Symbol,Tuple{Int,Vector{Float64}}}
end

function Base.getproperty(p::Node{T,N},s::Symbol) where {T,N}
    index = getfield(p,:index)
    if s ∈ T
        return index[findfirst(x->x==s,T)]
    else
        i,v = getfield(p,:data)[s]
        return v[index[i]]
    end
end
function Base.setproperty!(p::Node,s::Symbol,x::Float64)
    i,v = getfield(p,:data)[s]
    j = getfield(p,:index)[i]
    v[j] = x
end
function Base.getindex(p::Node,f::Symbol)
    i,v = getfield(p,:data)[f]
    j = getfield(p,:index)[i]
    return RV(j,v)
end

+(a::T,b::S) where {T<:Node,S<:Node} = (a.x+b.x,a.y+b.y,a.z+b.z)
-(a::T,b::S) where {T<:Node,S<:Node} = (a.x-b.x,a.y-b.y,a.z-b.z)

function printinfo(p::Node{T,N},s::Symbol) where {T,N}
    if s∈T
        index = getfield(p,:index)
        i = findfirst(x->x==s,T)
        @printf "Node "*string(s)*" = %i \n" index[i]
        for (key,(j,value)) in getfield(p,:data)
            if i == j
                @printf string(key)*" = %e  " value[index[i]]
            end
        end
    else
        error("No index of $s.")
    end
end