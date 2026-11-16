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

struct TRElement{T} <: AbstractElement
    𝓒::Vector{𝑿ⱼ}
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

# Require at least one argument: a variadic-only method makes the empty call
# `push!(aps)` ambiguous with the `Pair` method below, which `Aqua.test_ambiguities`
# flags. Requiring one argument removes the ambiguity without changing any
# existing call.
function Base.push!(aps::Vector{T},s::Symbol,ss::Symbol...;index::Symbol=:𝑠) where T<:AbstractElement
    data = getfield(aps[1].𝓖[1],:data)
    indices = getfield(aps[end].𝓖[end],:index)
    i = findfirst((x)->x==index,keys(indices))
    n = count(aps,index)
    for name in (s, ss...)
        data[name] = (i,zeros(n))
    end
end

function Base.push!(aps::Vector{T},sv::Pair{Symbol,Vector{Float64}},
                   svs::Pair{Symbol,Vector{Float64}}...;index::Symbol=:𝑠) where T<:AbstractElement
    data = getfield(aps[1].𝓖[1],:data)
    indices = getfield(aps[end].𝓖[end],:index)
    i = findfirst((x)->x==index,keys(indices))
    for (name,v) in (sv, svs...)
        data[name] = (i,v)
    end
end

function prescribe!(ξ::Node,sf::Pair{Symbol,F}) where F<:Function
    s,f = sf
    𝒙 = (ξ.x,ξ.y,ξ.z)
    # The normal components are only present when the element was built with
    # `normal = true`, and `n₃` only for 2-D boundaries. Reading a missing
    # field raises `KeyError`, so each candidate arity is guarded by `haskey`
    # rather than relying on `applicable` (whose arguments are evaluated first).
    # A prescribed function may take the coordinates only, or the coordinates
    # plus 1, 2 or 3 normal components. The 3-argument form is tried first, so a
    # plain `(x, y) -> …` must be matched explicitly rather than by relying on
    # `applicable(f, x, y, z)`.
    has(st) = haskey(getfield(ξ,:data), st) || st ∈ keys(getfield(ξ,:index))
    if applicable(f, ξ.x, ξ.y)
        v = f(ξ.x, ξ.y)
    elseif applicable(f, ξ.x, ξ.y, ξ.z)
        v = f(ξ.x, ξ.y, ξ.z)
    elseif has(:n₁) && applicable(f, ξ.x, ξ.y, ξ.z, ξ.n₁)
        v = f(ξ.x, ξ.y, ξ.z, ξ.n₁)
    elseif has(:n₁) && has(:n₂) && applicable(f, ξ.x, ξ.y, ξ.z, ξ.n₁, ξ.n₂)
        v = f(ξ.x, ξ.y, ξ.z, ξ.n₁, ξ.n₂)
    elseif has(:n₁) && has(:n₂) && has(:n₃) &&
           applicable(f, ξ.x, ξ.y, ξ.z, ξ.n₁, ξ.n₂, ξ.n₃)
        v = f(ξ.x, ξ.y, ξ.z, ξ.n₁, ξ.n₂, ξ.n₃)
    else
        error("prescribe!: the function given for :$s accepts none of " *
              "f(x,y), f(x,y,z), f(x,y,z,n₁), f(x,y,z,n₁,n₂), f(x,y,z,n₁,n₂,n₃)")
    end
    setproperty!(ξ,s,v)
end

function prescribe!(aps::Vector{T},sf::Pair{Symbol,F};index::Symbol=:𝐺) where {T<:AbstractElement,F<:Function}
    s,f = sf
    data = getfield((aps[1].𝓖)[1],:data)
    indices = getfield(aps[end].𝓖[end],:index)
    n = count(aps,index)
    i = findfirst((x)->x==index,keys(indices))
    data[s] = (i,zeros(n))
    if index == :𝑔
        𝓖 = aps[1].𝓖
        for ξ in 𝓖
            prescribe!(ξ,sf)
        end
    elseif index == :𝐺
        for ap in aps
            𝓖 = ap.𝓖
            for ξ in 𝓖
                prescribe!(ξ,sf)
            end
        end
    elseif index == :𝐶
        for ap in aps
            ξ = ap.𝓖[1]
            prescribe!(ξ,sf)
        end
    else
        error("prescribe error! Index is not supported.")
    end
end

function prescribe!(aps::Vector{T},sv::Pair{Symbol,Float64};index::Symbol=:𝐺) where T<:AbstractElement
    s,v = sv
    indices = getfield(aps[end].𝓖[end],:index)
    if index ∈ keys(indices)
        data = getfield((aps[1].𝓖)[1],:data)
        n = count(aps,index)
        i = findfirst((x)->x==index,keys(indices))
        data[s] = (i,ones(n).*v)
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