struct PiecewisePolynomial{𝑝}<:AbstractPiecewise
    𝓒::Vector{𝑿ᵢ}
    𝓖::Vector{𝑿ₛ}
end

struct PiecewiseParametric{𝑝,T}<:AbstractPiecewise
    𝓒::Vector{𝑿ᵢ}
    𝓖::Vector{𝑿ₛ}
end

get𝑛𝑝(::PiecewisePolynomial{:Constant1D}) = 1
get𝑛𝑝(::PiecewiseParametric{:Constant1D}) = 1
get𝑛𝑝(::PiecewisePolynomial{:Linear1D}) = 2
get𝑛𝑝(::PiecewiseParametric{:Linear1D}) = 2
get𝑛𝑝(::PiecewisePolynomial{:Quadratic1D}) = 3
get𝑛𝑝(::PiecewiseParametric{:Quadratic1D}) = 3
get𝑛𝑝(::PiecewisePolynomial{:Cubic1D}) = 4
get𝑛𝑝(::PiecewiseParametric{:Cubic1D}) = 4

get𝑛𝑝(::PiecewisePolynomial{:Constant2D}) = 1
get𝑛𝑝(::PiecewiseParametric{:Constant2D}) = 1
get𝑛𝑝(::PiecewisePolynomial{:Linear2D}) = 3
get𝑛𝑝(::PiecewiseParametric{:Linear2D}) = 3
get𝑛𝑝(::PiecewisePolynomial{:Quadratic2D}) = 6
get𝑛𝑝(::PiecewiseParametric{:Quadratic2D}) = 6
get𝑛𝑝(::PiecewisePolynomial{:Cubic2D}) = 10
get𝑛𝑝(::PiecewiseParametric{:Cubic2D}) = 10
get𝑛𝑝(::PiecewiseParametric{:Bubble,:Tri3}) = 1
get𝑛𝑝(::PiecewiseParametric{:Bubble,:Quad}) = 1

function set𝝭!(::PiecewisePolynomial{:Constant1D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    𝝭[1] = 1.0
end

function set𝝭!(::PiecewisePolynomial{:Linear1D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    𝝭[1] = 1.0
    𝝭[2] = 𝒙.x
end

function set∇𝝭!(::PiecewisePolynomial{:Linear1D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    𝝭[1] = 1.0
    𝝭[2] = 𝒙.x
    ∂𝝭∂x = 𝒙[:∂𝝭∂x]
    ∂𝝭∂x[1] = 0.0
    ∂𝝭∂x[2] = 1.0
end

function set𝝭!(::PiecewisePolynomial{:Quadratic1D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    x = 𝒙.x
    𝝭[1] = 1.0
    𝝭[2] = x
    𝝭[3] = x^2
end

function set∇𝝭!(ap::PiecewisePolynomial{:Quadratic1D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    x = 𝒙.x
    𝝭[1] = 1.0
    𝝭[2] = x
    𝝭[3] = x^2
    ∂𝝭∂x = 𝒙[:∂𝝭∂x]
    ∂𝝭∂x[1] = 0.0
    ∂𝝭∂x[2] = 1.0
    ∂𝝭∂x[3] = 2.0*x
end

function set𝝭!(::PiecewiseParametric{:Constant1D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    𝝭[1] = 1.0
end

function set𝝭!(::PiecewiseParametric{:Linear1D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    ξ = 0.5*(1.0+𝒙.ξ)
    𝝭[1] = 1.0
    𝝭[2] = ξ
end

function set𝝭!(::PiecewiseParametric{:Quadratic1D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    ξ = 0.5*(1.0+𝒙.ξ)
    𝝭[1] = 1.0
    𝝭[2] = ξ
    𝝭[3] = ξ^2
end

function set𝝭!(::PiecewisePolynomial{:Constant2D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    𝝭[1] = 1.0
end

function set∇𝝭!(::PiecewisePolynomial{:Linear2D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    ∂𝝭∂x = 𝒙[:∂𝝭∂x]
    ∂𝝭∂y = 𝒙[:∂𝝭∂y]
    𝝭[1] = 1.0
    𝝭[2] = 𝒙.x
    𝝭[3] = 𝒙.y
    ∂𝝭∂x[1] = 0.0
    ∂𝝭∂x[2] = 1.0
    ∂𝝭∂x[3] = 0.0
    ∂𝝭∂y[1] = 0.0
    ∂𝝭∂y[2] = 0.0
    ∂𝝭∂y[3] = 1.0
end

function set∇²𝝭!(::PiecewisePolynomial{:Linear2D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    ∂𝝭∂x = 𝒙[:∂𝝭∂x]
    ∂𝝭∂y = 𝒙[:∂𝝭∂y]
    ∂²𝝭∂x² = 𝒙[:∂²𝝭∂x²]
    ∂²𝝭∂y² = 𝒙[:∂²𝝭∂y²]
    ∂²𝝭∂x∂y = 𝒙[:∂²𝝭∂x∂y]
    𝝭[1] = 1.0
    𝝭[2] = 𝒙.x
    𝝭[3] = 𝒙.y
    ∂𝝭∂x[1] = 0.0
    ∂𝝭∂x[2] = 1.0
    ∂𝝭∂x[3] = 0.0
    ∂𝝭∂y[1] = 0.0
    ∂𝝭∂y[2] = 0.0
    ∂𝝭∂y[3] = 1.0
    ∂²𝝭∂x²[1] = 0.0
    ∂²𝝭∂x²[2] = 0.0
    ∂²𝝭∂x²[3] = 0.0 
    ∂²𝝭∂y²[1] = 0.0
    ∂²𝝭∂y²[2] = 0.0
    ∂²𝝭∂y²[3] = 0.0
    ∂²𝝭∂x∂y[1] = 0.0
    ∂²𝝭∂x∂y[2] = 0.0
    ∂²𝝭∂x∂y[3] = 0.0
end

function set𝝭!(::PiecewisePolynomial{:Linear2D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    𝝭[1] = 1.0
    𝝭[2] = 𝒙.x
    𝝭[3] = 𝒙.y
end
function set𝝭!(::PiecewisePolynomial{:Quadratic2D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    x = 𝒙.x
    y = 𝒙.y
    𝝭[1] = 1.0
    𝝭[2] = x
    𝝭[3] = y
    𝝭[4] = x^2
    𝝭[5] = x*y
    𝝭[6] = y^2
end

function set𝝭!(::PiecewiseParametric{:Constant2D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    𝝭[1] = 1.0
end

function set𝝭!(::PiecewiseParametric{:Linear2D,:Tri3},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    𝝭[1] = 1.0
    𝝭[2] = 𝒙.ξ
    𝝭[3] = 𝒙.η
end

function set𝝭!(::PiecewiseParametric{:Quadratic2D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    ξ = 𝒙.ξ
    η = 𝒙.η
    𝝭[1] = 1.0
    𝝭[2] = ξ
    𝝭[3] = η
    𝝭[4] = ξ^2
    𝝭[5] = ξ*η
    𝝭[6] = η^2
end

function set𝝭!(::PiecewiseParametric{:Bubble,:Tri3},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    ξ = 𝒙.ξ
    η = 𝒙.η
    γ = 1.0-ξ-η
    𝝭[1] = 27.0*ξ*η*γ
end

function set∇𝝭!(::PiecewiseParametric{:Bubble,:Tri3},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    ξ = 𝒙.ξ
    η = 𝒙.η
    γ = 1.0-ξ-η
    𝝭[1] = 27.0*ξ*η*γ
    ∂𝝭∂x = 𝒙[:∂𝝭∂x]
    ∂𝝭∂y = 𝒙[:∂𝝭∂y]
    ∂𝝭∂x[1] = 27.0*((γ-ξ)*𝒙.∂ξ∂x + (γ-η)*𝒙.∂η∂x)
    ∂𝝭∂y[1] = 27.0*((γ-ξ)*𝒙.∂ξ∂y + (γ-η)*𝒙.∂η∂y)
end

function set𝝭!(::PiecewiseParametric{:Bubble,:Quad},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    ξ = 𝒙.ξ
    η = 𝒙.η
    𝝭[1] = (1.0-ξ^2)*(1.0-η^2)
end

function set∇𝝭!(::PiecewiseParametric{:Bubble,:Quad},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    ξ = 𝒙.ξ
    η = 𝒙.η
    𝝭[1] = (1.0-ξ^2)*(1.0-η^2)
    ∂𝝭∂x = x[:∂𝝭∂x]
    ∂𝝭∂y = x[:∂𝝭∂y]
    ∂𝝭∂x[1] = - 2*ξ*(1.0-η^2)*𝒙.∂ξ∂x - 2*η*(1.0-ξ^2)*𝒙.∂η∂x
    ∂𝝭∂y[1] = - 2*ξ*(1.0-η^2)*𝒙.∂ξ∂y - 2*η*(1.0-ξ^2)*𝒙.∂η∂y
end
