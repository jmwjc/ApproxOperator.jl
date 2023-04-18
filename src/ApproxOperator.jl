module ApproxOperator

import Base: +, -, *, /, getindex, setindex!, getproperty, setproperty!, length, push!, fill!, issubset, intersect
import InteractiveUtils: subtypes
import Printf: @printf

abstract type AbstractElement{T} end
abstract type SpatialPartition end

include("node.jl")
include("element.jl")
# include("meshfree.jl")
include("operation.jl")

include("preprocession/integration.jl")
include("preprocession/geometry.jl")
include("preprocession/importcomsol.jl")
include("approximation/tri3.jl")
include("approximation/seg2.jl")
include("approximation/poi1.jl")
include("operation/potential.jl")
include("operation/elasticity.jl")
include("preprocession/importmsh.jl")
include("operation/Kirchhoff-Love_plate.jl")
export Operator
export set𝝭!, set∇𝝭!

end