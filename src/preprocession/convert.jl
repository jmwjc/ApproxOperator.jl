

function Seg2toTri3(seg2::Vector{T},tri3::Vector{S}) where {T,S<:AbstractElement}
    elms = Element{:Tri3}[]
    data = Dict{Symbol,Tuple{Int,Vector{Float64}}}()
    data_seg2 = getfield(seg2[1].𝓖[1],:data)
    nᵢ = length(data_seg2[:w][2])
    data[:w] = data_seg2[:w]
    data[:x] = data_seg2[:x]
    data[:y] = data_seg2[:y]
    data[:z] = data_seg2[:z]
    data[:𝑤] = data_seg2[:𝑤]
    data[:𝐽] = (3,Float64[])
    data[:ξ] = (2,Float64[])
    data[:η] = (2,Float64[])

    data[:n₁] = (3,Float64[])
    data[:n₂] = (3,Float64[])
    G = 0;C = 1;s = 0;
    for elm_seg2 in seg2
        𝓒_seg2 = elm_seg2.𝓒
        𝓖_seg2 = elm_seg2.𝓖
        for elm_tri3 in tri3
            𝓒_tri3 = elm_tri3.𝓒
            indices = indexin(𝓒_seg2,𝓒_tri3)
            if nothing ∉ indices
                x₁ = 𝓒_seg2[1].x
                y₁ = 𝓒_seg2[1].y
                x₂ = 𝓒_seg2[2].x
                y₂ = 𝓒_seg2[2].y
                𝐿 = 2*elm_seg2.𝐽
                push!(data[:𝐽][2],elm_tri3.𝐽)
                push!(data[:n₁][2],(y₂-y₁)/𝐿)
                push!(data[:n₂][2],(x₁-x₂)/𝐿)
                if indices == [2,3]
                    for ξ in 𝓖_seg2
                        push!(data[:ξ][2],0.5*(1-ξ.ξ))
                        push!(data[:η][2],0.5*(1+ξ.ξ))
                    end
                elseif indices == [3,1]
                    for ξ in 𝓖_seg2
                        push!(data[:ξ][2],0.0)
                        push!(data[:η][2],0.5*(1-ξ.ξ))
                    end
                elseif indices == [1,2]
                    for ξ in 𝓖_seg2
                        push!(data[:ξ][2],0.5*(1+ξ.ξ))
                        push!(data[:η][2],0.0)
                    end
                else
                    continue
                end
                𝓖 = [𝑿ₛ((𝑔=g,𝐺=G+g,𝐶=C,𝑠=s+3*(g-1)),data) for g in 1:nᵢ]
                push!(elms, Element{:Tri3}(𝓒_tri3,𝓖))
                G += nᵢ
                C += 1
                s += 3*nᵢ
                break
            end
        end
    end
    return elms
end

function Seg3toTri6(seg3::Vector{T},tri6::Vector{S}) where {T,S<:AbstractElement}
    elms = Element{:Tri6}[]
    data = Dict{Symbol,Tuple{Int,Vector{Float64}}}()
    data_seg3 = getfield(seg3[1].𝓖[1],:data)
    nᵢ = length(data_seg3[:w][2])
    data[:w] = data_seg3[:w]
    data[:x] = data_seg3[:x]
    data[:y] = data_seg3[:y]
    data[:z] = data_seg3[:z]
    data[:𝑤] = data_seg3[:𝑤]
    data[:𝐽] = (3,Float64[])
    data[:ξ] = (2,Float64[])
    data[:η] = (2,Float64[])

    data[:n₁] = (3,Float64[])
    data[:n₂] = (3,Float64[])
    G = 0;C = 1;s = 0;
    for elm_seg3 in seg3
        𝓒_seg3 = elm_seg3.𝓒
        𝓖_seg3 = elm_seg3.𝓖
        for elm_tri6 in tri6
            𝓒_tri6 = elm_tri6.𝓒
            indices = indexin(𝓒_seg3,𝓒_tri6)
            if nothing ∉ indices
                x₁ = 𝓒_seg3[1].x
                y₁ = 𝓒_seg3[1].y
                x₂ = 𝓒_seg3[2].x
                y₂ = 𝓒_seg3[2].y
                𝐿 = 2*elm_seg3.𝐽
                push!(data[:𝐽][2],elm_tri6.𝐽)
                push!(data[:n₁][2],(y₂-y₁)/𝐿)
                push!(data[:n₂][2],(x₁-x₂)/𝐿)
                if indices == [2,3,5]
                    for ξ in 𝓖_seg3
                        push!(data[:ξ][2],0.5*(1-ξ.ξ))
                        push!(data[:η][2],0.5*(1+ξ.ξ))
                    end
                end
                if indices == [3,1,6]
                    for ξ in 𝓖_seg3
                        push!(data[:ξ][2],0.0)
                        push!(data[:η][2],0.5*(1-ξ.ξ))
                    end
                end
                if indices == [1,2,4]
                    for ξ in 𝓖_seg3
                        push!(data[:ξ][2],0.5*(1+ξ.ξ))
                        push!(data[:η][2],0.0)
                    end
                end
                𝓖 = [𝑿ₛ((𝑔=g,𝐺=G+g,𝐶=C,𝑠=s+6*(g-1)),data) for g in 1:nᵢ]
                push!(elms, Element{:Tri6}(𝓒_tri6,𝓖))
                G += nᵢ
                C += 1
                s += 6*nᵢ
            end
        end
    end
    return elms
end

function Seg2toSegHermite(as::Vector{T},nodes::Vector{𝑿ᵢ},edges::Vector{Tuple{Int,Int}}) where T<:AbstractElement
    nₚ = length(getfield(as[1].𝓒[1],:data)[:x][2])
    elms = Element{:SegHermite}[]
    data𝓖 = Dict{Symbol,Tuple{Int,Vector{Float64}}}([
        :𝑤 => getfield(as[1].𝓖[1],:data)[:𝑤],
        :ξ => getfield(as[1].𝓖[1],:data)[:ξ],
        :x => getfield(as[1].𝓖[1],:data)[:x],
        :y => getfield(as[1].𝓖[1],:data)[:y],
        :z => getfield(as[1].𝓖[1],:data)[:z],
        :𝐽 => getfield(as[1].𝓖[1],:data)[:𝐽],
    ])
    s = 0
    for (C,a) in enumerate(as)
        𝓒 = a.𝓒
        𝓖 = a.𝓖
        𝓒_ = [nodes[xᵢ.𝐼] for xᵢ in 𝓒]
        𝓖_ = 𝑿ₛ[]
        ind_edge = indexin([(𝓒[1].𝐼,𝓒[2].𝐼)],edges)[1]
        push!(𝓒_,nodes[nₚ+ind_edge])
        ind_edge = indexin([(𝓒[2].𝐼,𝓒[1].𝐼)],edges)[1]
        push!(𝓒_,nodes[nₚ+ind_edge])
        for ξ in 𝓖
            push!(𝓖_,Node((𝑔=ξ.𝑔,𝐺=ξ.𝐺,𝐶=ξ.𝐶,𝑠=s),data𝓖))
            s += length(𝓒_)
        end
        push!(elms,Element{:SegHermite}(𝓒_,𝓖_))
    end
    return elms
end
function Tri3toTriHermite(as::Vector{T},nodes::Vector{𝑿ᵢ}) where T<:AbstractElement
    elms = Element{:TriHermite}[]
    edges = getTriEdgeIndices(as)
    nds = 𝑿ᵢ[]
    nₚ = length(nodes)
    nₑ = length(as)
    nₗ = length(edges)
    data𝓒 = Dict{Symbol,Tuple{Int,Vector{Float64}}}([
        :x => (1,zeros(nₚ+nₗ+nₑ)),
        :y => (1,zeros(nₚ+nₗ+nₑ)),
        :z => (1,zeros(nₚ+nₗ+nₑ)),
        :s₁ => (1,zeros(nₚ+nₗ+nₑ)),
        :s₂ => (1,zeros(nₚ+nₗ+nₑ)),
    ])
    data𝓖 = Dict{Symbol,Tuple{Int,Vector{Float64}}}([
        :𝑤 => getfield(as[1].𝓖[1],:data)[:𝑤],
        :ξ => getfield(as[1].𝓖[1],:data)[:ξ],
        :η => getfield(as[1].𝓖[1],:data)[:η],
        :x => getfield(as[1].𝓖[1],:data)[:x],
        :y => getfield(as[1].𝓖[1],:data)[:y],
        :z => getfield(as[1].𝓖[1],:data)[:z],
        :𝐽 => getfield(as[1].𝓖[1],:data)[:𝐽],
    ])
    for node in nodes
        𝐼 = node.𝐼
        xᵢ = 𝑿ᵢ((𝐼=𝐼,),data𝓒)
        xᵢ.x = node.x
        xᵢ.y = node.y
        xᵢ.z = node.z
        push!(nds,xᵢ)
    end
    for (i,(𝐼₁,𝐼₂)) in enumerate(edges)
        xᵢ = 𝑿ᵢ((𝐼=nₚ+i,),data𝓒)
        x₁ = nodes[𝐼₁].x
        y₁ = nodes[𝐼₁].y
        x₂ = nodes[𝐼₂].x
        y₂ = nodes[𝐼₂].y
        xᵢ.x = nodes[𝐼₁].x
        xᵢ.y = nodes[𝐼₁].y
        xᵢ.z = nodes[𝐼₁].z
        xᵢ.s₁ = x₂-x₁
        xᵢ.s₂ = y₂-y₁
        push!(nds,xᵢ)
    end
    s = 0
    for (C,a) in enumerate(as)
        𝓒 = a.𝓒
        𝓖 = a.𝓖
        𝓒_ = [nds[xᵢ.𝐼] for xᵢ in 𝓒]
        𝓖_ = 𝑿ₛ[]
        ind_edge = indexin([(𝓒[1].𝐼,𝓒[2].𝐼)],edges)[1]
        push!(𝓒_,nds[nₚ+ind_edge])
        ind_edge = indexin([(𝓒[1].𝐼,𝓒[3].𝐼)],edges)[1]
        push!(𝓒_,nds[nₚ+ind_edge])
        ind_edge = indexin([(𝓒[2].𝐼,𝓒[3].𝐼)],edges)[1]
        push!(𝓒_,nds[nₚ+ind_edge])
        ind_edge = indexin([(𝓒[2].𝐼,𝓒[1].𝐼)],edges)[1]
        push!(𝓒_,nds[nₚ+ind_edge])
        ind_edge = indexin([(𝓒[3].𝐼,𝓒[1].𝐼)],edges)[1]
        push!(𝓒_,nds[nₚ+ind_edge])
        ind_edge = indexin([(𝓒[3].𝐼,𝓒[2].𝐼)],edges)[1]
        push!(𝓒_,nds[nₚ+ind_edge])
        x₁ = 𝓒[1].x
        y₁ = 𝓒[1].y
        x₂ = 𝓒[2].x
        y₂ = 𝓒[2].y
        x₃ = 𝓒[3].x
        y₃ = 𝓒[3].y
        xᵢ = 𝑿ᵢ((𝐼=nₚ+nₗ+C,),data𝓒)
        xᵢ.x = (x₁+x₂+x₃)/3
        xᵢ.y = (y₁+y₂+y₃)/3
        push!(𝓒_,xᵢ)
        push!(nds,xᵢ)
        for ξ in 𝓖
            push!(𝓖_,Node((𝑔=ξ.𝑔,𝐺=ξ.𝐺,𝐶=ξ.𝐶,𝑠=s),data𝓖))
            s += length(𝓒_)
        end
        push!(elms,Element{:TriHermite}(𝓒_,𝓖_))
    end
    return elms, nds, edges
end

function Tri3toTriBell(as::Vector{T},nodes::Vector{𝑿ᵢ}) where T<:AbstractElement
    elms = Element{:TriBell}[]
    nds = 𝑿ᵢ[]
    nₚ = length(nodes)
    data𝓒 = Dict{Symbol,Tuple{Int,Vector{Float64}}}([
        :x => (1,zeros(6*nₚ)),
        :y => (1,zeros(6*nₚ)),
        :z => (1,zeros(6*nₚ)),
    ])
    data𝓖 = Dict{Symbol,Tuple{Int,Vector{Float64}}}([
        :𝑤 => getfield(as[1].𝓖[1],:data)[:𝑤],
        :ξ => getfield(as[1].𝓖[1],:data)[:ξ],
        :η => getfield(as[1].𝓖[1],:data)[:η],
        :x => getfield(as[1].𝓖[1],:data)[:x],
        :y => getfield(as[1].𝓖[1],:data)[:y],
        :z => getfield(as[1].𝓖[1],:data)[:z],
        :𝐽 => getfield(as[1].𝓖[1],:data)[:𝐽],
    ])
    for node in nodes
        𝐼 = node.𝐼
        for i in 1:6
            xᵢ = 𝑿ᵢ((𝐼=6*𝐼-6+i,),data𝓒)
            xᵢ.x = node.x
            xᵢ.y = node.y
            xᵢ.z = node.z
            push!(nds,xᵢ)
        end
    end
    s = 0
    for a in as
        𝓒 = a.𝓒
        𝓖 = a.𝓖
        𝓒_ = 𝑿ᵢ[]
        𝓖_ = 𝑿ₛ[]
        for xᵢ in 𝓒
            𝐼 = xᵢ.𝐼
            append!(𝓒_,nds[6*𝐼-5:6*𝐼])
        end
        for ξ in 𝓖
            push!(𝓖_,Node((𝑔=ξ.𝑔,𝐺=ξ.𝐺,𝐶=ξ.𝐶,𝑠=s),data𝓖))
        end
        s += length(𝓒_)
        push!(elms,Element{:TriBell}(𝓒_,𝓖_))
    end
    return elms, nds
end

function getTriEdgeIndices(as::Vector{T}) where T<:AbstractElement
    indices = Vector{Tuple{Int,Int}}()
    for a in as
        𝓒 = a.𝓒
        push!(indices,(𝓒[1].𝐼,𝓒[2].𝐼))
        push!(indices,(𝓒[1].𝐼,𝓒[3].𝐼))
        push!(indices,(𝓒[2].𝐼,𝓒[3].𝐼))
        push!(indices,(𝓒[2].𝐼,𝓒[1].𝐼))
        push!(indices,(𝓒[3].𝐼,𝓒[1].𝐼))
        push!(indices,(𝓒[3].𝐼,𝓒[2].𝐼))
    end
    return unique!(indices)
end