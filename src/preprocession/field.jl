
mutable struct Field{S₁,N₁,S₂,N₂}
    type::DataType
    index::Dict{Symbol,Int}
    p₁::Int
    l₁::Int
    𝓒::Vector{Node{S₁,N₁}}
    data𝓒::Dict{Symbol,Tuple{Int,Vector{Float64}}}
    p₂::Int
    l₂::Int
    𝓖::Vector{Node{S₂,N₂}}
    data𝓖::Dict{Symbol,Tuple{Int,Vector{Float64}}}
end

function Field{S₁,N₁,S₂,N₂}(type::DataType,qtype::Symbol,data𝓒::Dict{Symbol,Tuple{Int,Vector{Float64}}} = Dict{Symbol,Tuple{Int,Vector{Float64}}}();qindex::Int=1) where {S₁,N₁,S₂,N₂}
    index = Dict([s=>0 for s in (S₁...,S₂...)])
    scheme = quadraturerule(qtype)
    p₁ = 0
    l₁ = 0
    𝓒 = Node{S₁,N₁}[]
    p₂ = 0
    l₂ = 0
    𝓖 = Node{S₂,N₂}[]
    data𝓖 = Dict{Symbol,Tuple{Int,Vector{Float64}}}()
    for (s, v) in scheme
        data𝓖[s] = (qindex,v)
    end
    return Field{S₁,N₁,S₂,N₂}(type,index,p₁,l₁,𝓒,data𝓒,p₂,l₂,𝓖,data𝓖)
end

Base.getproperty(f::Field,s::Symbol) = getfield(f,:index)[s]
function Base.setproperty!(f::Field,s::Symbol,I::Int)
    getfield(f,:index)[s] = I
end

function add𝓒!(f::Field{S₁,N₁,S₂,N₂}) where {S₁,N₁,S₂,N₂}
    index = getfield(f,:index)
    l = getfield(f,:l₁)
    𝓢 = getfield(f,:𝓒)
    data = getfield(f,:data𝓒)
    l += 1
    setfield!(f,:l₁,l)
    i = tuple([index[s] for s in S₁]...)
    p = Node{S₁,N₁}(i,data)
    push!(𝓢,p)
end

function add𝓖!(f::Field{S₁,N₁,S₂,N₂}) where {S₁,N₁,S₂,N₂}
    index = getfield(f,:index)
    l = getfield(f,:l₂)
    𝓢 = getfield(f,:𝓖)
    data = getfield(f,:data𝓖)
    l += 1
    setfield!(f,:l₂,l)
    i = tuple([index[s] for s in S₂]...)
    p = Node{S₂,N₂}(i,data)
    push!(𝓢,p)
end

function get𝓒(f::Field)
    p = getfield(f,:p₁)
    l = getfield(f,:l₁)
    𝓢 = getfield(f,:𝓒)
    t = (p,l,𝓢)
    p += l
    l = 0
    setfield!(f,:p₁,p)
    setfield!(f,:l₁,l)
    return t
end

function get𝓖(f::Field)
    p = getfield(f,:p₂)
    l = getfield(f,:l₂)
    𝓢 = getfield(f,:𝓖)
    t = (p,l,𝓢)
    p += l
    l = 0
    setfield!(f,:p₂,p)
    setfield!(f,:l₂,l)
    return t
end

function Base.push!(f::Field{S₁,N₁,S₂,N₂},vs::Pair{Symbol,Symbol}) where {S₁,N₁,S₂,N₂}
    v, s = vs
    l = getfield(f,:index)[s]
    w = zeros(l)
    if s ∈ S₁
        data = getfield(f,:data𝓒)
        i = findfirst(x->x==s,S₁)
        for (s,(j,x)) in data
            if i == j
                resize!(w,length(x))
                break
            end
        end
    else
        data = getfield(f,:data𝓖)
        i = findfirst(x->x==s,S₂)
    end
    data[v] = (i,w)
end

function Base.push!(f::Field{S₁,N₁,S₂,N₂},vs::Pair{Symbol,Tuple{Symbol,Vector{Float64}}}) where {S₁,N₁,S₂,N₂}
    v, (s, t) = vs
    if s ∈ S₁
        data = getfield(f,:data𝓒)
        i = findfirst(x->x==s,S₁)
    elseif s ∈ S₂
        data = getfield(f,:data𝓖)
        i = findfirst(x->x==s,S₂)
    else
        data = getfield(f,:data𝓖)
        i = 0
    end
    data[v] = (i,t)
end

function Base.push!(f::Field,ps::Any...)
    for p in ps
        push!(f,p)
    end
end

function (f::Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠)})(as::Vector{T}) where T<:AbstractGeometry
    if f.type<:AbstractPiecewise
        return setPiecewise(as)
    else
        return setElement(as)
    end
end

function setElement(as::Vector{T}) where T<:AbstractGeometry
    data = getfield(f,:data𝓖)
    weights = data[:w][2]
    if haskey(data,:γ)
        ξ = data[:ξ][2]
        η = data[:η][2]
        γ = data[:γ][2]
        points = zip(ξ,η,γ)
        push!(f,
            :ξ=>(:𝑔,ξ),
            :η=>(:𝑔,η),
            :γ=>(:𝑔,γ),
        )
    elseif haskey(data,:η)
        ξ = data[:ξ][2]
        η = data[:η][2]
        points = zip(ξ,η)
        push!(f,
            :ξ=>(:𝑔,ξ),
            :η=>(:𝑔,η),
        )
    else
        ξ = data[:ξ][2]
        points = ξ
        push!(f,
            :ξ=>(:𝑔,ξ),
        )
    end
    scheme = zip(weights,points)
    ne = length(as)
    ni = length(as[1].i)
    ng = length(weights)
    type = getfield(f,:type)
    elements = type[]
    𝑤 = zeros(ng*ne)
    x = zeros(ng*ne)
    y = zeros(ng*ne)
    z = zeros(ng*ne)
    push!(f,
        :𝑤=>(:𝐺,𝑤),
        :x=>(:𝐺,x),
        :y=>(:𝐺,y),
        :z=>(:𝐺,z),
    )
    for (C,a) in enumerate(as)
        for i in a.i
            f.𝐼 = i
            ApproxOperator.add𝓒!(f)
        end
        for (g,(w,ps)) in enumerate(scheme)
            f.𝑔 = g
            f.𝐺 += 1
            f.𝐶 = C
            ApproxOperator.add𝓖!(f)
            f.𝑠 += ni
            𝑤[f.𝐺] = get𝐽(a,ps...)*w
            x[f.𝐺], y[f.𝐺], z[f.𝐺] = a(ps...)
        end
        𝓒 = ApproxOperator.get𝓒(f)
        𝓖 = ApproxOperator.get𝓖(f)
        push!(elements,type(𝓒,𝓖))
    end
    return elements
end

function setPiecewise(as::Vector{T}) where T<:AbstractGeometry
    data = getfield(f,:data𝓖)
    weights = data[:w][2]
    if haskey(data,:γ)
        ξ = data[:ξ][2]
        η = data[:η][2]
        γ = data[:γ][2]
        points = zip(ξ,η,γ)
        push!(f,
            :ξ=>(:𝑔,ξ),
            :η=>(:𝑔,η),
            :γ=>(:𝑔,γ),
        )
    elseif haskey(data,:η)
        ξ = data[:ξ][2]
        η = data[:η][2]
        points = zip(ξ,η)
        push!(f,
            :ξ=>(:𝑔,ξ),
            :η=>(:𝑔,η),
        )
    else
        ξ = data[:ξ][2]
        points = ξ
        push!(f,
            :ξ=>(:𝑔,ξ),
        )
    end
    scheme = zip(weights,points)
    ne = length(as)
    ni = length(as[1].i)
    ng = length(weights)
    type = getfield(f,:type)
    elements = type[]
    𝑤 = zeros(ng*ne)
    x = zeros(ng*ne)
    y = zeros(ng*ne)
    z = zeros(ng*ne)
    𝐽 = zeros(ng*ne)
    push!(f,
        :𝑤=>(:𝐺,𝑤),
        :x=>(:𝐺,x),
        :y=>(:𝐺,y),
        :z=>(:𝐺,z),
        :𝐽=>(:𝐺,𝐽),
    )
    𝑛𝑝 = get𝑛𝑝(as[1])
    for (C,a) in enumerate(as)
        for i in 1:𝑛𝑝
            f.𝐼 += 1
            ApproxOperator.add𝓒!(f)
        end
        for (g,(w,ps)) in enumerate(scheme)
            f.𝑔 = g
            f.𝐺 += 1
            f.𝐶 = C
            ApproxOperator.add𝓖!(f)
            f.𝑠 += 𝑛𝑝
            𝐽[f.𝐺] = get𝐽(a,ps...)
            𝑤[f.𝐺] = 𝐽[f.𝐺]*w
            x[f.𝐺], y[f.𝐺], z[f.𝐺] = a(ps...)
        end
        𝓒 = ApproxOperator.get𝓒(f)
        𝓖 = ApproxOperator.get𝓖(f)
        push!(elements,type(𝓒,𝓖))
    end
    return elements
end

function (f::Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠)})(as::Vector{T},sp::SpatialPartition) where T<:AbstractGeometry
    data = getfield(f,:data𝓖)
    weights = data[:w][2]
    if haskey(data,:γ)
        ξ = data[:ξ][2]
        η = data[:η][2]
        γ = data[:γ][2]
        points = zip(ξ,η,γ)
        push!(f,
            :ξ=>(:𝑔,ξ),
            :η=>(:𝑔,η),
            :γ=>(:𝑔,γ),
        )
    elseif haskey(data,:η)
        ξ = data[:ξ][2]
        η = data[:η][2]
        points = zip(ξ,η)
        push!(f,
            :ξ=>(:𝑔,ξ),
            :η=>(:𝑔,η),
        )
    else
        ξ = data[:ξ][2]
        points = ξ
        push!(f,
            :ξ=>(:𝑔,ξ),
        )
    end
    scheme = zip(weights,points)
    ne = length(as)
    ng = length(weights)
    push!(f,
        :x=>(:𝐼,sp.x),
        :y=>(:𝐼,sp.y),
        :z=>(:𝐼,sp.z),
    )
    type = getfield(f,:type)
    elements = type[]
    𝑤 = zeros(ng*ne)
    x = zeros(ng*ne)
    y = zeros(ng*ne)
    z = zeros(ng*ne)
    push!(f,
        :𝑤=>(:𝐺,𝑤),
        :x=>(:𝐺,x),
        :y=>(:𝐺,y),
        :z=>(:𝐺,z),
    )
    for (C,a) in enumerate(as)
        indices = Set{Int}()
        for (w,ps) in scheme
            xᵢ,yᵢ,zᵢ = a(ps...)
            union!(indices,sp(xᵢ,yᵢ,zᵢ))
        end
        ni = length(indices)
        for i in indices
            f.𝐼 = i
            ApproxOperator.add𝓒!(f)
        end
        for (g,(w,ps)) in enumerate(scheme)
            f.𝑔 = g
            f.𝐺 += 1
            f.𝐶 = C
            ApproxOperator.add𝓖!(f)
            f.𝑠 += ni
            𝑤[f.𝐺] = get𝐽(a,ps...)*w
            x[f.𝐺], y[f.𝐺], z[f.𝐺] = a(ps...)
        end
        𝓒 = ApproxOperator.get𝓒(f)
        𝓖 = ApproxOperator.get𝓖(f)
        push!(elements,type(𝓒,𝓖))
    end
    return elements
end

function (f::Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠)})(as::Vector{Tuple{T₁,T₂,T₃}}) where {T₁<:AbstractGeometry,T₂<:AbstractGeometry,T₃<:AbstractElement}
    type = getfield(f,:type)
    if type<:AbstractReproducingKernel
        data = getfield(f,:data𝓖)
        weights = data[:w][2]
        ne = length(as)
        ng = length(weights)
        if T₁ ≠ T₂
            if haskey(data,:γ)
                ξ = data[:ξ][2]
                η = data[:η][2]
                γ = data[:γ][2]
                points = zip(ξ,η,γ)
                ξ = zeros(ng*ne)
                η = zeros(ng*ne)
                γ = zeros(ng*ne)
                push!(f,
                    :ξ=>(:𝐺,ξ)
                    :η=>(:𝐺,η)
                    :γ=>(:𝐺,γ)
                )
            elseif haskey(data,:η)
                ξ = data[:ξ][2]
                η = data[:η][2]
                points = zip(ξ,η)
                ξ = zeros(ng*ne)
                η = zeros(ng*ne)
                push!(f,
                    :ξ=>(:𝐺,ξ)
                    :η=>(:𝐺,η)
                )
            else
                points = data[:ξ][2]
                ξ = zeros(ng*ne)
                push!(f,
                    :ξ=>(:𝐺,ξ)
                )
            end
        else
            if haskey(data,:γ)
                ξ = data[:ξ][2]
                η = data[:η][2]
                γ = data[:γ][2]
                points = zip(ξ,η,γ)
                push!(f,
                    :ξ=>(:𝑔,ξ)
                    :η=>(:𝑔,η)
                    :γ=>(:𝑔,γ)
                )
            elseif haskey(data,:η)
                ξ = data[:ξ][2]
                η = data[:η][2]
                points = zip(ξ,η)
                ξ = zeros(ng*ne)
                η = zeros(ng*ne)
                push!(f,
                    :ξ=>(:𝑔,ξ)
                    :η=>(:𝑔,η)
                )
            else
                points = data[:ξ][2]
                ξ = zeros(ng*ne)
                push!(f,
                    :ξ=>(:𝑔,ξ)
                )
            end
        end
        scheme = zip(weights,points)
        elements = type[]
        𝑤 = zeros(ng*ne)
        x = zeros(ng*ne)
        y = zeros(ng*ne)
        z = zeros(ng*ne)
        push!(f,
            :𝑤=>(:𝐺,𝑤),
            :x=>(:𝐺,x),
            :y=>(:𝐺,y),
            :z=>(:𝐺,z),
        )
        for (C,(a,b,c)) in enumerate(as)
            𝓒 = getfield(b,:𝓒)
            𝓖ˢ = getfield(b,:𝓖)
            ni = 𝓒[2]
            for (g,(w,ps)) in enumerate(scheme)
                f.𝑔 = g
                f.𝐺 += 1
                f.𝐶 = C
                ApproxOperator.add𝓖!(f)
                f.𝑠 += ni
                𝑤[f.𝐺] = get𝐽(a,ps...)*w
                if T₁ ≠ T₂
                    (x[f.𝐺],y[f.𝐺],z[f.𝐺]),𝝃 = b(a,ps...)
                    if haskey(data,:γ)
                        ξ[f.𝐺] = 𝝃[1]
                        η[f.𝐺] = 𝝃[2]
                        γ[f.𝐺] = 𝝃[3]
                    elseif haskey(data,:η)
                        ξ[f.𝐺] = 𝝃[1]
                        η[f.𝐺] = 𝝃[2]
                    else
                        ξ[f.𝐺] = 𝝃[1]
                    end
                else
                    x[f.𝐺], y[f.𝐺], z[f.𝐺] = a(ps...)
                end
            end
            𝓖 = ApproxOperator.get𝓖(f)
            push!(elements,type(𝓒,𝓖,𝓖ˢ))
        end
        return elements
    else
        error("Element type is wrong.")
    end
end

function (f::Field{(:𝐼,:𝐽),2,(:𝑔,:𝐺,:𝐶,:𝑠)})(as::Vector{Tri3},𝓑::Vector{Set{Int}})
    data = getfield(f,:data𝓖)
    weights = data[:w][2]
    ξ = data[:ξ][2]
    η = data[:η][2]
    points = zip(ξ,η)
    push!(f,
        :ξ=>(:𝑔,ξ),
        :η=>(:𝑔,η),
    )
    scheme = zip(weights,points)
    ne = length(as)
    ni = length(as[1].i)
    ng = length(weights)
    push!(f,
        :x=>(:𝐽,as[1].x),
        :y=>(:𝐽,as[1].y),
        :z=>(:𝐽,as[1].z),
    )
    type = getfield(f,:type)
    elements = type[]
    𝑤 = zeros(ng*ne)
    x = zeros(ng*ne)
    y = zeros(ng*ne)
    z = zeros(ng*ne)
    push!(f,
        :𝑤=>(:𝐺,𝑤),
        :x=>(:𝐺,x),
        :y=>(:𝐺,y),
        :z=>(:𝐺,z),
    )
    for (C,a) in enumerate(as)
        for i in a.i
            f.𝐼 = findfirst(x->x==Set(setdiff(a.i,i)),𝓑)
            f.𝐽 = i
            ApproxOperator.add𝓒!(f)
        end
        for (g,(w,ps)) in enumerate(scheme)
            f.𝑔 = g
            f.𝐺 += 1
            f.𝐶 = C
            ApproxOperator.add𝓖!(f)
            f.𝑠 += ni
            𝑤[f.𝐺] = get𝐽(a,ps...)*w
            x[f.𝐺], y[f.𝐺], z[f.𝐺] = a(ps...)
        end
        𝓒 = ApproxOperator.get𝓒(f)
        𝓖 = ApproxOperator.get𝓖(f)
        push!(elements,type(𝓒,𝓖))
    end
    return elements
end

function (f::Field{(:𝐼,:𝐽),2,(:𝑔,:𝐺,:𝐶,:𝑠)})(as::Vector{Tuple{Seg2,Tri3}},𝓑::Vector{Set{Int}})
    data = getfield(f,:data𝓖)
    weights = data[:w][2]
    ne = length(as)
    ng = length(weights)
    ξ = data[:ξ][2]
    points = ξ
    ξ = zeros(ne*ng)
    η = zeros(ne*ng)
    n₁ = zeros(ne*ng)
    n₂ = zeros(ne*ng)
    push!(f,
        :ξ=>(:𝐺,ξ),
        :η=>(:𝐺,η),
        :n₁=>(:𝐺,n₁),
        :n₂=>(:𝐺,n₂),
    )
    scheme = zip(weights,points)
    ni = length(as[1][2].i)
    push!(f,
        :x=>(:𝐽,as[1][2].x),
        :y=>(:𝐽,as[1][2].y),
        :z=>(:𝐽,as[1][2].z),
    )
    type = getfield(f,:type)
    elements = type[]
    𝑤 = zeros(ng*ne)
    x = zeros(ng*ne)
    y = zeros(ng*ne)
    z = zeros(ng*ne)
    push!(f,
        :𝑤=>(:𝐺,𝑤),
        :x=>(:𝐺,x),
        :y=>(:𝐺,y),
        :z=>(:𝐺,z),
    )
    for (C,(a,b)) in enumerate(as)
        for i in b.i
            f.𝐼 = findfirst(x->x==Set(setdiff(b.i,i)),𝓑)
            f.𝐽 = i
            ApproxOperator.add𝓒!(f)
        end
        for (g,(w,ps)) in enumerate(scheme)
            f.𝑔 = g
            f.𝐺 += 1
            f.𝐶 = C
            ApproxOperator.add𝓖!(f)
            f.𝑠 += ni
            𝑤[f.𝐺] = get𝐽(a,ps...)*w
            (x[f.𝐺], y[f.𝐺], z[f.𝐺]), (ξ[f.𝐺], η[f.𝐺]) = b(a,ps...)
            n₁[f.𝐺], n₂[f.𝐺] = get𝒏(a)
        end
        𝓒 = ApproxOperator.get𝓒(f)
        𝓖 = ApproxOperator.get𝓖(f)
        push!(elements,type(𝓒,𝓖))
    end
    return elements
end