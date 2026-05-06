struct RV
    i::Int
    v::Vector{Float64}
end

getindex(r::RV,i::Int) = r.v[r.i+i]
getindex(r::RV,is::UnitRange{Int}) = [r.v[r.i+i] for i in is]
function setindex!(r::RV,x::Float64,i::Int)
    r.v[r.i+i] = x
end

abstract type AbstractFieldValue{T} end

struct UniformValue{T} <: AbstractFieldValue{T}
    value::T
end

struct PerNodeValue{T} <: AbstractFieldValue{T}
    values::Vector{T}
end

struct LazyValue <: AbstractFieldValue{Float64}
    func::Any
end

const FieldValue = Union{UniformValue{Float64},PerNodeValue{Float64},LazyValue}

struct NodeData
    value::Dict{Symbol,FieldValue}
    index::Dict{Symbol,Int}
end

NodeData(; kwargs...) = begin
    dv = Dict{Symbol,FieldValue}()
    di = Dict{Symbol,Int}()
    for (s, (i, v)) in kwargs
        dv[s] = v
        di[s] = i
    end
    NodeData(dv, di)
end

function Base.getindex(d::NodeData, s::Symbol)
    (d.index[s], d.value[s])
end

function Base.haskey(d::NodeData, s::Symbol)
    haskey(d.value, s)
end

struct Node{T,N}
    index::NamedTuple{T,NTuple{N,Int}}
    data::NodeData
end

const 𝑿ᵢ = Node{(:𝐼,),1}
const 𝑿ₛ = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}

function Base.getproperty(p::Node{T,N},s::Symbol) where {T,N}
    index = getfield(p,:index)
    if s ∈ T
        return index[s]
    end
    dat = getfield(p,:data)
    v = dat.value[s]
    if v isa UniformValue
        return v.value
    elseif v isa PerNodeValue
        i = dat.index[s]
        return v.values[index[i]]
    else
        return getproperty_lazy(p, v)
    end
end

@noinline function getproperty_lazy(p::Node, v::LazyValue)::Float64
    𝒙 = (p.x,p.y,p.z)
    f = v.func
    if applicable(f,𝒙...)
        return Float64(f(𝒙...))
    elseif applicable(f,𝒙...,p.n₁)
        return Float64(f(𝒙...,p.n₁))
    elseif applicable(f,𝒙...,p.n₁,p.n₂)
        return Float64(f(𝒙...,p.n₁,p.n₂))
    else
        return Float64(f(𝒙...,p.n₁,p.n₂,p.n₃))
    end
end

function Base.setproperty!(p::Node,s::Symbol,x::Float64)
    dat = getfield(p,:data)
    v = dat.value[s]
    if v isa PerNodeValue
        i = dat.index[s]
        j = getfield(p,:index)[i]
        v.values[j] = x
    else
        i = dat.index[s]
        nₜ = getfield(p,:index)[end]
        arr = zeros(nₜ)
        if v isa UniformValue
            fill!(arr, v.value)
        end
        arr[getfield(p,:index)[i]] = x
        dat.value[s] = PerNodeValue(arr)
    end
end

function Base.getindex(p::Node,s::Symbol)
    dat = getfield(p,:data)
    v = dat.value[s]
    i = dat.index[s]
    j = getfield(p,:index)[i]
    return RV(j, v.values)
end

+(a::T,b::S) where {T<:Node,S<:Node} = (a.x+b.x,a.y+b.y,a.z+b.z)
-(a::T,b::S) where {T<:Node,S<:Node} = (a.x-b.x,a.y-b.y,a.z-b.z)

push!(ps::Vector{T},svs::Pair{Symbol,S}...) where {T<:Node,S} = push!(ps[1],svs...)

function push!(p::Node{T},svs::Pair{Symbol,Vector{Float64}}...;index::Symbol=:𝐼) where T
    dat = getfield(p,:data)
    for (s,v) in svs
        i = findfirst((x)->x==index,T)
        push!(dat.value,s=>PerNodeValue(v))
        push!(dat.index,s=>i)
    end
end

function Base.getproperty(ps::Vector{N},s::Symbol) where N<:Node
    if s == :ref || isempty(ps)
        return getfield(ps,s)
    end
    dat = getfield(ps[1],:data)
    if !haskey(dat, s)
        return getfield(ps,s)
    end
    v = dat.value[s]
    if v isa PerNodeValue
        return [getproperty(p, s) for p in ps]
    else
        return v
    end
end

function printinfo(p::Node{S}) where S
    index = getfield(p,:index)
    dat = getfield(p,:data)
    print("Node")
    print(index)
    println(":")
    shapes = Symbol[]
    for (name, vs) in dat.value
        n = dat.index[name]
        if n ≠ 0
            s = S[n]
            if s ≠ :𝑠
                i = index[n]
                if vs isa PerNodeValue
                    v = vs.values[i]
                    @printf "  %s(%s = %i): %e\n" string(name) string(s) i v
                elseif vs isa UniformValue
                    @printf "  %s(%s = %i): %e (uniform)\n" string(name) string(s) i vs.value
                else
                    @printf "  %s(%s = %i): (lazy)\n" string(name) string(s) i
                end
            else
                push!(shapes,name)
            end
        end
    end
    return shapes
end

Base.show(io::IO,::MIME"text/plain",p::Node) = printinfo(p)
Base.show(io::IO,p::Node) = printinfo(p)

Base.push!(v::PerNodeValue, x) = push!(v.values, x)
Base.push!(v::PerNodeValue, xs...) = push!(v.values, xs...)
Base.append!(v::PerNodeValue, xs) = append!(v.values, xs)
Base.fill!(v::PerNodeValue, x) = fill!(v.values, x)
