module ApproxOperator

import Base: +, -, *, /, getindex, setindex!, getproperty, setproperty!, length, push!, fill!, issubset, intersect, show, Pair
import Printf: @printf

abstract type AbstractElement end
abstract type AbstractPiecewise<:AbstractElement end
abstract type SpatialPartition end

include("node.jl")
include("element.jl")
include("operation.jl")
include("shapefunction.jl")

include("preprocession/importmsh.jl")
include("preprocession/convert.jl")

include("approximation/quad.jl")
include("approximation/quad8.jl")
include("approximation/tri3.jl")
include("approximation/tri6.jl")
include("approximation/tet4.jl")
include("approximation/hex8.jl")
include("approximation/seg2.jl")
include("approximation/seg3.jl")
include("approximation/poi1.jl")
include("approximation/seghermite.jl")
include("approximation/trihermite.jl")
include("approximation/tribell.jl")
include("approximation/meshfree.jl")
include("approximation/reproducingkernel.jl")
include("approximation/kernelfunction.jl")
include("approximation/CrouzeixRaviart.jl")
include("approximation/piecewise.jl")

include("operation/heat.jl")
include("operation/elasticity.jl")
include("operation/hyperelasticity.jl")
include("operation/curved_beam.jl")
include("operation/thin_plate.jl")
include("operation/thin_shell.jl")
include("operation/hamilton.jl")
include("operation/test.jl")
include("operation/stokes.jl")
# include("operation/phasefield.jl")
# include("operation/error_estimates.jl")


export prescribe!
export 𝑿ᵢ, 𝑿ₛ
export Element
export TRElement
export ReproducingKernel, RegularGrid
export PiecewiseParametric, PiecewisePolynomial
export set𝝭!, set∇𝝭!, set∇²𝝭!, set∇̂³𝝭!
# export RKGradientSmoothing, GRKGradientSmoothing
export getPhysicalGroups, get𝑿ᵢ, getElements, addEdgeElements, getDOFs
export getPiecewiseElements, getPiecewiseBoundaryElements
export getMacroElements, getMacroBoundaryElements, getCurvedElements, getCurvedPiecewiseElements
export Tri3toTriBell, Tri3toTriHermite, Seg2toSegHermite

end