"""
"""
struct Element{T} <: AbstractElement
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

for set𝝭 in (:set𝝭!,:set∇𝝭!,:set∇²𝝭!,:set∇̂³𝝭!)
    @eval begin
        function $set𝝭(a::T) where T<:AbstractElement
            𝓖 = a.𝓖
            for x in 𝓖
                $set𝝭(a,x)
            end
        end
        $set𝝭(as::Vector{T}) where T<:AbstractElement = $set𝝭.(as)
    end
end

function Base.push!(aps::Vector{T},sis::Pair{Symbol, Symbol}...) where T<:AbstractElement
    for si in sis
        s,i = si
        data = getfield(aps[1].𝓖[1],:data)
        index = getfield(aps[end].𝓖[end],:index)
        if i ≠ :𝑠
            j = findfirst((x)->x==i,keys(index))
            data[s] = (j,zeros(index[i]))
        else
            data[s] = (4,zeros(index[:𝑠]+length(aps[end].𝓒)))
        end
    end
end

function Base.push!(aps::Vector{T},svs::Pair{Symbol, Tuple{Int,Vector{Float64}}}...) where T<:AbstractElement
    for sv in svs
        s,v = sv
        data = getfield(aps[1].𝓖[1],:data)
        data[s] = v
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