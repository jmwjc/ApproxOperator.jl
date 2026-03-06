"""
Element{T}
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


function count(aps::Vector{T},i::Symbol) where T<:AbstractElement
    index = getfield(aps[end].𝓖[end],:index)
    return i ≠ :𝑠 ? index[i] : index[i]+length(aps[end].𝓒)
end

function Base.push!(aps::Vector{T},ss::Symbol...;index::Symbol=:𝑠) where T<:AbstractElement
    data = getfield(aps[1].𝓖[1],:data)
    indices = getfield(aps[end].𝓖[end],:index)
    i = findfirst((x)->x==index,keys(indices))
    n = count(aps,index)
    for s in ss
        data[s] = (i,zeros(n))
    end
end

# function Base.push!(aps::Vector{T}, ss::Symbol...; index::Symbol = :𝑠) where {T<:AbstractElement}
#     @assert !isempty(aps) "push!: aps 为空"
#     data = getfield(aps[1].𝓖[1], :data)

#     # 取最后一个积分点的 index（NamedTuple），例如 (𝑔=?, 𝐺=?, 𝐶=?, 𝑠=?)
#     last_idx = getfield(aps[end].𝓖[end], :index)

#     # index 在 NamedTuple 里的位置（你原来的写法保留）
#     i = findfirst(x -> x == index, keys(last_idx))
#     @assert i !== nothing "push!: index=$(index) 不存在于积分点 index 中: $(keys(last_idx))"

#     # 关键：不要再调用 count(aps,index)，直接用最后一个积分点的偏移 + 最后一个 patch 的 nc 计算总长度
#     # 你们的 𝑠 是每个积分点在“形函数数组”里的起始偏移，构造时通常从 0 开始累加 nc
#     s_last = getproperty(last_idx, index)
#     nc_last = getfield(aps[end], :𝓒)[2]  # 𝓒 = (c, nc, patch_nodes)

#     n = s_last + nc_last
#     @assert n ≥ 0 "push!: 计算得到的 n 非法: n=$n (s_last=$s_last, nc_last=$nc_last)"

#     for s in ss
#         data[s] = (i, zeros(n))
#     end
# end


function Base.push!(aps::Vector{T},svs::Pair{Symbol, Vector{Float64}}...;index::Symbol=:𝑠) where T<:AbstractElement
    data = getfield(aps[1].𝓖[1],:data)
    indices = getfield(aps[end].𝓖[end],:index)
    i = findfirst((x)->x==index,keys(indices))
    for sv in svs
        s,v = sv
        data[s] = (i,v)
    end
end

function prescribe!(ξ::Node,sf::Pair{Symbol,F}) where F<:Function
    s,f = sf
    𝒙 = (ξ.x,ξ.y,ξ.z)
    if applicable(f,𝒙...)
        v = f(𝒙...)
    elseif applicable(f,𝒙...,ξ.n₁)
        v = f(𝒙...,ξ.n₁)
    elseif applicable(f,𝒙...,ξ.n₁,ξ.n₂)
        v = f(𝒙...,ξ.n₁,ξ.n₂)
    elseif applicable(f,𝒙...,ξ.n₁,ξ.n₂,ξ.n₃)
        v = f(𝒙...,ξ.n₁,ξ.n₂,ξ.n₃)
    end
    setproperty!(ξ,s,v)
end

function prescribe!(aps::Vector{T},sf::Pair{Symbol,F};index::Symbol=:𝐺) where {T<:AbstractElement,F<:Function}
    s,f = sf
    data = getfield((aps[1].𝓖)[1],:data)
    indices = getfield(aps[end].𝓖[end],:index)
    n = count(aps,index)
    i = findfirst((x)->x==index,keys(indices))
    haskey(data,s) ? nothing : push!(data,s=>(i,zeros(n)))
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