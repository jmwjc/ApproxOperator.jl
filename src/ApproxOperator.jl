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

include("approximation/quad.jl")
include("approximation/quad8.jl")
include("approximation/tri3.jl")
include("approximation/tri6.jl")
include("approximation/tet4.jl")
include("approximation/tet10.jl")
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
include("approximation/rkwave.jl")
include("approximation/CrouzeixRaviart.jl")
include("approximation/piecewise.jl")
include("preprocession/importmsh.jl")

# `GmshImport` is a submodule; without this the `export`s below would name
# symbols that do not exist in `ApproxOperator` (Aqua flags them as undefined
# exports, and `using ApproxOperator` would hand the user nothing).
using .GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements, getCurvedElements
include("preprocession/convert.jl")

include("operation/heat.jl")
include("operation/elasticity.jl")
include("operation/hyperelasticity.jl")
include("operation/curved_beam.jl")
include("operation/beam.jl")
include("operation/thin_plate.jl")
include("operation/thin_shell.jl")
include("operation/thick_plate.jl")
include("operation/hamilton.jl")
include("operation/test.jl")
include("operation/stokes.jl")
include("operation/weighted_residual.jl")
include("operation/phasefield.jl")
include("operation/plasticity.jl")
# NOTE: error_estimates.jl has been retired. Its operators now live with the
# physics they belong to:
#   Heat.L₂_heat_flux, Heat.L₂_heat_flux_u, Heat.H₃
#   Elasticity.Hₑ_up_fem, Elasticity.Hₑ_up_mix, Elasticity.Hₑ_up_mix_Q4P1,
#   Elasticity.Hₑ_Incompressible
# Dropped as duplicates or dead code: T_error (== Heat.L₂), set∇𝑢! (no callers).


export prescribe!
export 𝑿ᵢ, 𝑿ₛ
export Element
export TRElement
export ReproducingKernel, RegularGrid
export PiecewiseParametric, PiecewisePolynomial
export set𝝭!, set∇𝝭!, set∇²𝝭!, set∇̂³𝝭!
# export RKGradientSmoothing, GRKGradientSmoothing
export getPhysicalGroups, get𝑿ᵢ, getElements, getCurvedElements, getDOFs
# `addCentroidNodes!` lives in `GmshImport` and is called as
# `GmshImport.addCentroidNodes!(nodes, dimTag)`; it is not re-exported
# because every other importer entry point is qualified too.
export Tri3toTri6!, Seg2toSeg3!, Tri3toTriBell, Tri3toTriHermite, Seg2toSegHermite, Tri3toTRTri3, Tri3toDTRTri3, Seg2toTRTri3

end
