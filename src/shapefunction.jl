
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

function set𝝭!(as::Vector{T}) where T<:AbstractElement
    push!(as,:𝝭)
    set𝝭!.(as)
end

function set∇𝝭!(as::Vector{T}) where T<:AbstractElement
    if T<:Element1D
        push!(as,:𝝭,:∂𝝭∂x)
    elseif T<:Element2D
        push!(as,:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    else
        push!(as,:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z)
    end
    set∇𝝭!.(as)
end
function set∇²𝝭!(as::Vector{T}) where T<:AbstractElement
    if T<:Element1D
        push!(as,:𝝭,:∂𝝭∂x,:∂²𝝭∂x²)
    elseif T<:Element2D
        push!(as,:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂y²,:∂²𝝭∂x∂y)
    else
        push!(as,:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z,:∂²𝝭∂x²,:∂²𝝭∂y²,:∂²𝝭∂z²,:∂²𝝭∂x∂y,:∂²𝝭∂x∂z,:∂²𝝭∂y∂z)
    end
    set∇²𝝭!.(as)
end
function set∇̂³𝝭!(as::Vector{T}) where T<:AbstractElement
    if T<:Element1D
        push!(as,:𝝭,:∂𝝭∂x,:∂²𝝭∂x²,:∂³𝝭∂x³)
    elseif T<:Element2D
        push!(as,:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂y²,:∂²𝝭∂x∂y,:∂³𝝭∂x³,:∂³𝝭∂y³,:∂³𝝭∂x²∂y)
    else
        push!(as,:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z,
              :∂²𝝭∂x²,:∂²𝝭∂y²,:∂²𝝭∂z²,
              :∂²𝝭∂x∂y,:∂²𝝭∂x∂z,:∂²𝝭∂y∂z,
              :∂³𝝭∂x³,:∂³𝝭∂y³,:∂³𝝭∂z³,
              :∂³𝝭∂x²∂y,:∂³𝝭∂x∂y²,
              :∂³𝝭∂x²∂z,:∂³𝝭∂x∂z²,
              :∂³𝝭∂y²∂z,:∂³𝝭∂y∂z²,:∂³𝝭∂x∂y∂z)
    end
    set∇̂³𝝭!.(as)
end

