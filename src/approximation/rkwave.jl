get𝑛𝒑(  ::ReproducingKernel{:Quadratic2D}) = 5
get𝒑(   ::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (1.,sin(x[1]),cos(x[1]),sin(x[2]),cos(x[2]))
get∂𝒑∂x(::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (0.,cos(x[1]),-sin(x[1]),0.,0.)
get∂𝒑∂y(::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (0.,0.,0.,cos(x[1]),-sin(x[1]))