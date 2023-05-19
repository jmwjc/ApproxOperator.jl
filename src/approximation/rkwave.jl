get𝑛𝒑(  ::ReproducingKernel{:Quadratic2D}) = 5
get𝒑(   ::ReproducingKernel{:Quadratic2D},k::NTuple{2,Float64},x::NTuple{3,Float64}) = (1.,sin(k[1]x),cos(k[1]x),sin(k[2]x),cos(k[2]x))
get∂𝒑∂x(::ReproducingKernel{:Quadratic2D},k::NTuple{2,Float64},x::NTuple{3,Float64}) = (0.,cos(k[1]x),-sin(k[1]x),0.,0.)
get∂𝒑∂y(::ReproducingKernel{:Quadratic2D},k::NTuple{2,Float64},x::NTuple{3,Float64}) = (0.,0.,0.,cos(k[1]x),-sin(k[1]x))