module ApproxOperator

import Base: +, -, *, /, getindex, setindex!, getproperty, setproperty!, length, push!, fill!, issubset, intersect
import InteractiveUtils: subtypes
# import Printf: @printf

abstract type AbstractElement{T} end
abstract type SpatialPartition end

include("node.jl")
include("element.jl")
include("operation.jl")

include.(filter(contains(r".jl$"), readdir("preprocession"; join=true)))
include.(filter(contains(r".jl$"), readdir("approximation"; join=true)))
include.(filter(contains(r".jl$"), readdir("operation"; join=true)))

include("littletools.jl")

export prescribe!
export Operator
export Node, Element
export set𝝭!, set∇𝝭!, set∇²𝝭!

end