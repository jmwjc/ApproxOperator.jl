"""
Element{T}
"""
struct Element{T} <: AbstractElement
    𝓒::Vector{𝑿ᵢ}
    𝓖::Vector{𝑿ₛ}
end

struct ReproducingKernel{𝑝,𝑠,𝜙}<:AbstractElement
    𝓒::Vector{𝑿ᵢ}
    𝓖::Vector{𝑿ₛ}
end

function Base.getproperty(a::T,s::Symbol) where T<:AbstractElement
    if s∈(:𝓒,:𝓖)
        return getfield(a,s)
    else
        ξ = getfield(a,:𝓖)[1]
        return getproperty(ξ,s)
    end
end

function Base.setproperty!(ap::T,s::Symbol,v::Float64) where T<:AbstractElement
    ξ = getfield(ap,:𝓖)[1]
    setproperty!(ξ,s,v)
end

function count(aps::Vector{T},i::Symbol) where T<:AbstractElement
    index = getfield(aps[end].𝓖[end],:index)
    return i ≠ :𝑠 ? index[i] : index[i]+length(aps[end].𝓒)
end

function Base.push!(aps::Vector{T},ss::Symbol...;index::Symbol=:𝑠) where T<:AbstractElement
    dat = getfield(aps[1].𝓖[1],:data)
    indices = getfield(aps[end].𝓖[end],:index)
    i = findfirst((x)->x==index,keys(indices))
    n = count(aps,index)
    for s in ss
        dat.value[s] = PerNodeValue(zeros(n))
        dat.index[s] = i
    end
end

function Base.push!(aps::Vector{T},svs::Pair{Symbol,Vector{Float64}}...;index::Symbol=:𝑠) where T<:AbstractElement
    dat = getfield(aps[1].𝓖[1],:data)
    indices = getfield(aps[end].𝓖[end],:index)
    i = findfirst((x)->x==index,keys(indices))
    for sv in svs
        s,v = sv
        dat.value[s] = PerNodeValue(v)
        dat.index[s] = i
    end
end

function prescribe!(aps::Vector{T},sf::Pair{Symbol,F};index::Symbol=:𝐺) where {T<:AbstractElement,F<:Function}
    s, f = sf
    dat = getfield((aps[1].𝓖)[1],:data)
    indices = getfield(aps[end].𝓖[end],:index)
    i = findfirst((x)->x==index,keys(indices))
    dat.value[s] = LazyValue(f)
    dat.index[s] = i
end

function prescribe!(aps::Vector{T},sv::Pair{Symbol,Float64};index::Symbol=:𝐺) where T<:AbstractElement
    s,v = sv
    indices = getfield(aps[end].𝓖[end],:index)
    if index ∈ keys(indices)
        dat = getfield((aps[1].𝓖)[1],:data)
        i = findfirst((x)->x==index,keys(indices))
        dat.value[s] = UniformValue(v)
        dat.index[s] = i
    else
        error("prescribe error! Index is not supported.")
    end
end

function prescribe!(aps::Vector{T},sfs::Pair...;index::Symbol=:𝐺) where T<:AbstractElement
    for sf in sfs
        prescribe!(aps,sf;index=index)
    end
end

function Base.show(io::IO,::MIME"text/plain",a::T) where T<:AbstractElement
    𝓒 = a.𝓒
    𝓖 = a.𝓖
    println(T)
    println("𝓒:")
    for (i,p) in enumerate(𝓒)
        print("$i. ")
        printinfo(p)
    end
    println("𝓖:")
    for (j,p) in enumerate(𝓖)
        print("$j. ")
        shapes = printinfo(p)
        @printf "         "
        for shape in shapes
            @printf "%13s" string(shape)
        end
        @printf "\n"
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            @printf "  (𝐼 = %i)" I
            for shape in shapes
                @printf " %e" p[shape][i]
            end
            @printf "\n"
        end
    end
end