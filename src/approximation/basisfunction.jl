
get∇̃²₂𝒑(ap::ReproducingKernel,x::Any) = get∂²𝒑∂x²(ap,x), get∂²𝒑∂x∂y(ap,x), get∂²𝒑∂x∂y(ap,x), get∂²𝒑∂y²(ap,x)
get∇²𝒑(ap::ReproducingKernel,x::Any) = get𝒑(ap,x), get∂𝒑∂x(ap,x), get∂𝒑∂y(ap,x), get∂²𝒑∂x²(ap,x), get∂²𝒑∂x∂y(ap,x), get∂²𝒑∂y²(ap,x), get∂𝒑∂z(ap,x), get∂²𝒑∂x∂z(ap,x), get∂²𝒑∂y∂z(ap,x), get∂²𝒑∂z²(ap,x)
get∇∇²𝒑(ap::ReproducingKernel,x::Any) = get𝒑(ap,x), get∂³𝒑∂x³(ap,x), get∂³𝒑∂x²∂y(ap,x), get∂³𝒑∂x²∂y(ap,x), get∂³𝒑∂x∂y²(ap,x), get∂³𝒑∂x∂y²(ap,x), get∂³𝒑∂y³(ap,x)
get∇𝒑₁(ap::ReproducingKernel{:Linear1D,𝑠,𝜙,T},ξ::Any) where {𝑠,𝜙,T} = get𝒑₁(ap,ξ), get∂𝒑₁∂ξ(ap,ξ)
get∇𝒑₁(ap::ReproducingKernel{𝒑,𝑠,𝜙,:Seg2},ξ::Any) where {𝒑,𝑠,𝜙} = get𝒑₁(ap,ξ), get∂𝒑₁∂ξ(ap,ξ)
get∇𝒑₁(ap::ReproducingKernel{𝒑,𝑠,𝜙,:Tri3},ξ::Any) where {𝒑,𝑠,𝜙} = get𝒑₁(ap,ξ), get∂𝒑₁∂ξ(ap,ξ), get∂𝒑₁∂η(ap,ξ)
get∇𝒑₂(ap::ReproducingKernel{𝒑,𝑠,𝜙,:Tri3},ξ::Any) where {𝒑,𝑠,𝜙} = get𝒑₂(ap,ξ), get∂𝒑₂∂ξ(ap,ξ), get∂𝒑₂∂η(ap,ξ)
get∇²𝒑₂(ap::ReproducingKernel{𝒑,𝑠,𝜙,:Tri3},ξ::Any) where {𝒑,𝑠,𝜙} = get𝒑₂(ap,ξ), get∂𝒑₂∂ξ(ap,ξ), get∂𝒑₂∂η(ap,ξ), get∂²𝒑₂∂ξ²(ap,ξ), get∂²𝒑₂∂ξ∂η(ap,ξ), get∂²𝒑₂∂η²(ap,ξ)


