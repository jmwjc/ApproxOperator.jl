
const Element1D = Union{
    Element{:Seg2},
    Element{:Seg3},
    Element{:SegHermite},
    ReproducingKernel{:Linear1D},
    ReproducingKernel{:Quadratic1D},
    ReproducingKernel{:Cubic1D}
}
const Element2D = Union{
    Element{:Tri3},
    Element{:Tri6},
    Element{:TriBell},
    Element{:TriHermite},
    Element{:Quad4},
    Element{:Quad8},
    ReproducingKernel{:Linear2D},
    ReproducingKernel{:Quadratic2D},
    ReproducingKernel{:Cubic2D}
}
const Element3D = Union{
    Element{:Tet4},
    Element{:Tet10},
    Element{:Hex8},
    ReproducingKernel{:Linear3D},
    ReproducingKernel{:Quadratic3D},
    ReproducingKernel{:Cubic3D}
}

for set𝝭 in (:set𝝭!,:set∇𝝭!,:set∇²𝝭!,:set∇̂³𝝭!)
    @eval begin
        function $set𝝭(a::T) where T<:AbstractElement
            𝓖 = a.𝓖
            for x in 𝓖
                 $set𝝭(a,x)
            end
        end
    end
end

function get𝗠size(T::DataType)
    if T<:ReproducingKernel{:Linear1D}
        return 3
    elseif T<:ReproducingKernel{:Quadratic1D}
        return 6
    elseif T<:ReproducingKernel{:Cubic1D}
        return 10
    elseif T<:ReproducingKernel{:Linear2D}
        return 6
    elseif T<:ReproducingKernel{:Quadratic2D}
        return 21
    elseif T<:ReproducingKernel{:Cubic2D}
        return 55
    elseif T<:ReproducingKernel{:Linear3D}
        return 10
    elseif T<:ReproducingKernel{:Quadratic3D}
        return 55
    elseif T<:ReproducingKernel{:Cubic3D}
        return 220
    else
        error("Element type $T does not have matrix M")
    end
end

function set𝝭!(as::Vector{T}) where T<:AbstractElement
    push!(as,:𝝭)
    type = typeof(as[1])
    if type<:ReproducingKernel
        push!(as,:𝗠=>zeros(get𝗠size(type)))
    end
    set𝝭!.(as)
end

function set∇𝝭!(as::Vector{T}) where T<:AbstractElement
    push!(as,:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z)
    type = typeof(as[1])
    if type<:ReproducingKernel
        n = get𝗠size(type)
        push!(as,:𝗠=>zeros(n),:∂𝗠∂x=>zeros(n),:∂𝗠∂y=>zeros(n),:∂𝗠∂z=>zeros(n))
    end
    set∇𝝭!.(as)
end

function set∇²𝝭!(as::Vector{T}) where T<:AbstractElement
    push!(as,:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z,:∂²𝝭∂x²,:∂²𝝭∂y²,:∂²𝝭∂z²,:∂²𝝭∂x∂y,:∂²𝝭∂x∂z,:∂²𝝭∂y∂z)
    type = typeof(as[1])
    if type<:ReproducingKernel
        n = get𝗠size(type)
        push!(as,:𝗠=>zeros(n),:∂𝗠∂x=>zeros(n),:∂𝗠∂y=>zeros(n),:∂𝗠∂z=>zeros(n),:∂²𝗠∂x²=>zeros(n),:∂²𝗠∂y²=>zeros(n),:∂²𝗠∂z²=>zeros(n),:∂²𝗠∂x∂y=>zeros(n),:∂²𝗠∂x∂z=>zeros(n),:∂²𝗠∂y∂z=>zeros(n))
    end
    set∇²𝝭!.(as)
end

function set∇̂³𝝭!(as::Vector{T}) where T<:AbstractElement
    push!(as,:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z,
             :∂²𝝭∂x²,:∂²𝝭∂y²,:∂²𝝭∂z²,
             :∂²𝝭∂x∂y,:∂²𝝭∂x∂z,:∂²𝝭∂y∂z,
             :∂³𝝭∂x³,:∂³𝝭∂y³,:∂³𝝭∂z³,
             :∂³𝝭∂x²∂y,:∂³𝝭∂x∂y²,
             :∂³𝝭∂x²∂z,:∂³𝝭∂x∂z²,
             :∂³𝝭∂y²∂z,:∂³𝝭∂y∂z²,:∂³𝝭∂x∂y∂z)
    type = typeof(as[1])
    if type<:ReproducingKernel
        n = get𝗠size(type)
        push!(as,:𝗠=>zeros(n),:∂𝗠∂x=>zeros(n),:∂𝗠∂y=>zeros(n),:∂𝗠∂z=>zeros(n),:∂²𝗠∂x²=>zeros(n),:∂²𝗠∂y²=>zeros(n),:∂²𝗠∂z²=>zeros(n),:∂²𝗠∂x∂y=>zeros(n),:∂²𝗠∂x∂z=>zeros(n),:∂²𝗠∂y∂z=>zeros(n),:∂³𝗠∂x³=>zeros(n),:∂³𝗠∂y³=>zeros(n),:∂³𝗠∂z³=>zeros(n),:∂³𝗠∂x²∂y=>zeros(n),:∂³𝗠∂x∂y²=>zeros(n),:∂³𝗠∂x²∂z=>zeros(n),:∂³𝗠∂x∂z²=>zeros(n),:∂³𝗠∂y²∂z=>zeros(n),:∂³𝗠∂y∂z²=>zeros(n),:∂³𝗠∂x∂y∂z=>zeros(n))
    end
    set∇̂³𝝭!.(as)
end

