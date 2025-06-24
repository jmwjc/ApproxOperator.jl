struct PiecewisePolynomial{𝑝}<:AbstractPiecewise
    𝓒::Vector{𝑿ᵢ}
    𝓖::Vector{𝑿ₛ}
end

struct PiecewiseParametric{𝑝,T}<:AbstractPiecewise
    𝓒::Vector{𝑿ᵢ}
    𝓖::Vector{𝑿ₛ}
end

get𝑛𝑝(::PiecewisePolynomial{:Constant}) = 1
get𝑛𝑝(::PiecewiseParametric{:Constant}) = 1
get𝑛𝑝(::PiecewisePolynomial{:Linear1D}) = 2
get𝑛𝑝(::PiecewiseParametric{:Linear1D}) = 2
get𝑛𝑝(::PiecewisePolynomial{:Quadratic1D}) = 3
get𝑛𝑝(::PiecewiseParametric{:Quadratic1D}) = 3
get𝑛𝑝(::PiecewisePolynomial{:Cubic1D}) = 4
get𝑛𝑝(::PiecewiseParametric{:Cubic1D}) = 4

get𝑛𝑝(::PiecewisePolynomial{:Linear2D}) = 3
get𝑛𝑝(::PiecewiseParametric{:Linear2D}) = 3
get𝑛𝑝(::PiecewisePolynomial{:Quadratic2D}) = 6
get𝑛𝑝(::PiecewiseParametric{:Quadratic2D}) = 6
get𝑛𝑝(::PiecewisePolynomial{:Cubic2D}) = 10
get𝑛𝑝(::PiecewiseParametric{:Cubic2D}) = 10
get𝑛𝑝(::PiecewisePolynomial{:4}) = 15
get𝑛𝑝(::PiecewisePolynomial{:5}) = 21
get𝑛𝑝(::PiecewisePolynomial{:6}) = 28
get𝑛𝑝(::PiecewiseParametric{:Bubble,:Tri3}) = 1
get𝑛𝑝(::PiecewiseParametric{:Bubble,:Quad}) = 1

function set𝝭!(::PiecewisePolynomial{:Constant},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    𝝭[1] = 1.0
end

function set∇𝝭!(::PiecewisePolynomial{:Constant},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    ∂𝝭∂x = 𝒙[:∂𝝭∂x]
    ∂𝝭∂y = 𝒙[:∂𝝭∂x]
    𝝭[1] = 1.0
    ∂𝝭∂x[1] = 0.0
    ∂𝝭∂y[1] = 0.0
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
function set∇𝝭!(::PiecewisePolynomial{:Constant2D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    ∂𝝭∂x = 𝒙[:∂𝝭∂x]
    ∂𝝭∂y = 𝒙[:∂𝝭∂y]
    𝝭[1] = 1.0
    ∂𝝭∂x[1] = 0.0
    ∂𝝭∂y = 0.0
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

function set∇𝝭!(::PiecewisePolynomial{:Quadratic2D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    x = 𝒙.x
    y = 𝒙.y
    ∂𝝭∂x = 𝒙[:∂𝝭∂x]
    ∂𝝭∂y = 𝒙[:∂𝝭∂y]
    𝝭[1] = 1.0
    𝝭[2] = x
    𝝭[3] = y
    𝝭[4] = x^2
    𝝭[5] = x*y
    𝝭[6] = y^2
    ∂𝝭∂x[1] = 0.0
    ∂𝝭∂x[2] = 1.0
    ∂𝝭∂x[3] = 0.0
    ∂𝝭∂x[4] = 2*x
    ∂𝝭∂x[5] = y
    ∂𝝭∂x[6] = 0.0
    ∂𝝭∂y[1] = 0.0
    ∂𝝭∂y[2] = 0.0
    ∂𝝭∂y[3] = 1.0
    ∂𝝭∂y[4] = 0.0
    ∂𝝭∂y[5] = x
    ∂𝝭∂y[6] = 2*y
end

function set𝝭!(::PiecewisePolynomial{:Cubic2D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    x = 𝒙.x
    y = 𝒙.y
    𝝭[1] = 1.0
    𝝭[2] = x
    𝝭[3] = y
    𝝭[4] = x^2
    𝝭[5] = x*y
    𝝭[6] = y^2
    𝝭[7] = x^3
    𝝭[8] = x^2*y
    𝝭[9] = x*y^2
    𝝭[10] = y^3
end

function set∇𝝭!(::PiecewisePolynomial{:Cubic2D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    x = 𝒙.x
    y = 𝒙.y
    ∂𝝭∂x = 𝒙[:∂𝝭∂x]
    ∂𝝭∂y = 𝒙[:∂𝝭∂y]
    𝝭[1] = 1.0
    𝝭[2] = x
    𝝭[3] = y
    𝝭[4] = x^2
    𝝭[5] = x*y
    𝝭[6] = y^2
    𝝭[7] = x^3
    𝝭[8] = x^2*y
    𝝭[9] = x*y^2
    𝝭[10] = y^3
    ∂𝝭∂x[1] = 0.0
    ∂𝝭∂x[2] = 1.0
    ∂𝝭∂x[3] = 0.0
    ∂𝝭∂x[4] = 2*x
    ∂𝝭∂x[5] = y
    ∂𝝭∂x[6] = 0.0
    ∂𝝭∂x[7] = 3*x^2
    ∂𝝭∂x[8] = 2*x*y
    ∂𝝭∂x[9] = y^2
    ∂𝝭∂x[10] = 0.0
    ∂𝝭∂y[1] = 0.0
    ∂𝝭∂y[2] = 0.0
    ∂𝝭∂y[3] = 1.0
    ∂𝝭∂y[4] = 0.0
    ∂𝝭∂y[5] = x
    ∂𝝭∂y[6] = 2*y
    ∂𝝭∂y[7] = 0.0
    ∂𝝭∂y[8] = x^2
    ∂𝝭∂y[9] = 2*x*y
    ∂𝝭∂y[10] = 3*y^2
end
function set𝝭!(::PiecewisePolynomial{:Quartic2D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    x = 𝒙.x
    y = 𝒙.y
    𝝭[1] = 1.0
    𝝭[2] = x
    𝝭[3] = y
    𝝭[4] = x^2
    𝝭[5] = x*y
    𝝭[6] = y^2
    𝝭[7] = x^3
    𝝭[8] = x^2*y
    𝝭[9] = x*y^2
    𝝭[10] = y^3
    𝝭[11] = x^4
    𝝭[12] = x^3*y
    𝝭[13] = x*y^3
    𝝭[14] = x^2*y^2
    𝝭[15] = y^4
end

function set∇𝝭!(::PiecewisePolynomial{:Quartic2D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    x = 𝒙.x
    y = 𝒙.y
    ∂𝝭∂x = 𝒙[:∂𝝭∂x]
    ∂𝝭∂y = 𝒙[:∂𝝭∂y]
    𝝭[1] = 1.0
    𝝭[2] = x
    𝝭[3] = y
    𝝭[4] = x^2
    𝝭[5] = x*y
    𝝭[6] = y^2
    𝝭[7] = x^3
    𝝭[8] = x^2*y
    𝝭[9] = x*y^2
    𝝭[10] = y^3
    𝝭[11] = x^4
    𝝭[12] = x^3*y
    𝝭[13] = x*y^3
    𝝭[14] = x^2*y^2
    𝝭[15] = y^4
    ∂𝝭∂x[1] = 0.0
    ∂𝝭∂x[2] = 1.0
    ∂𝝭∂x[3] = 0.0
    ∂𝝭∂x[4] = 2*x
    ∂𝝭∂x[5] = y
    ∂𝝭∂x[6] = 0.0
    ∂𝝭∂x[7] = 3*x^2
    ∂𝝭∂x[8] = 2*x*y
    ∂𝝭∂x[9] = y^2
    ∂𝝭∂x[10] = 0.0
    ∂𝝭∂x[11] = 4*x^3
    ∂𝝭∂x[12] = 3*x^2*y
    ∂𝝭∂x[13] = y^3
    ∂𝝭∂x[14] = 2*x*y^2
    ∂𝝭∂x[15] = 0.0
    ∂𝝭∂y[1] = 0.0
    ∂𝝭∂y[2] = 0.0
    ∂𝝭∂y[3] = 1.0
    ∂𝝭∂y[4] = 0.0
    ∂𝝭∂y[5] = x
    ∂𝝭∂y[6] = 2*y
    ∂𝝭∂y[7] = 0.0
    ∂𝝭∂y[8] = x^2
    ∂𝝭∂y[9] = 2*x*y
    ∂𝝭∂y[10] = 3*y^2
    ∂𝝭∂y[11] = 0.0
    ∂𝝭∂y[12] = x^3
    ∂𝝭∂y[13] = 3*x*y^2
    ∂𝝭∂y[14] = 2*x^2*y
    ∂𝝭∂y[15] = 4*y^3
end
function set𝝭!(::PiecewisePolynomial{:Quintic2D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    x = 𝒙.x
    y = 𝒙.y
    𝝭[1] = 1.0
    𝝭[2] = x
    𝝭[3] = y
    𝝭[4] = x^2
    𝝭[5] = x*y
    𝝭[6] = y^2
    𝝭[7] = x^3
    𝝭[8] = x^2*y
    𝝭[9] = x*y^2
    𝝭[10] = y^3
    𝝭[11] = x^4
    𝝭[12] = x^3*y
    𝝭[13] = x*y^3
    𝝭[14] = x^2*y^2
    𝝭[15] = y^4
    𝝭[16] = x^5
    𝝭[17] = x^4*y
    𝝭[18] = x^3*y^2
    𝝭[19] = x^2*y^3
    𝝭[20] = x*y^4
    𝝭[21] = y^5
end

function set∇𝝭!(::PiecewisePolynomial{:Quintic2D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    x = 𝒙.x
    y = 𝒙.y
    ∂𝝭∂x = 𝒙[:∂𝝭∂x]
    ∂𝝭∂y = 𝒙[:∂𝝭∂y]
    𝝭[1] = 1.0
    𝝭[2] = x
    𝝭[3] = y
    𝝭[4] = x^2
    𝝭[5] = x*y
    𝝭[6] = y^2
    𝝭[7] = x^3
    𝝭[8] = x^2*y
    𝝭[9] = x*y^2
    𝝭[10] = y^3
    𝝭[11] = x^4
    𝝭[12] = x^3*y
    𝝭[13] = x*y^3
    𝝭[14] = x^2*y^2
    𝝭[15] = y^4
    𝝭[16] = x^5
    𝝭[17] = x^4*y
    𝝭[18] = x^3*y^2
    𝝭[19] = x^2*y^3
    𝝭[20] = x*y^4
    𝝭[21] = y^5
    ∂𝝭∂x[1] = 0.0
    ∂𝝭∂x[2] = 1.0
    ∂𝝭∂x[3] = 0.0
    ∂𝝭∂x[4] = 2*x
    ∂𝝭∂x[5] = y
    ∂𝝭∂x[6] = 0.0
    ∂𝝭∂x[7] = 3*x^2
    ∂𝝭∂x[8] = 2*x*y
    ∂𝝭∂x[9] = y^2
    ∂𝝭∂x[10] = 0.0
    ∂𝝭∂x[11] = 4*x^3
    ∂𝝭∂x[12] = 3*x^2*y
    ∂𝝭∂x[13] = y^3
    ∂𝝭∂x[14] = 2*x*y^2
    ∂𝝭∂x[15] = 0.0
    ∂𝝭∂x[16] = 5*x^4
    ∂𝝭∂x[17] = 4*x^3*y
    ∂𝝭∂x[18] = 3*x^2*y^2
    ∂𝝭∂x[19] = 2*x*y^3
    ∂𝝭∂x[20] = y^4
    ∂𝝭∂x[21] = 0.0
    ∂𝝭∂y[1] = 0.0
    ∂𝝭∂y[2] = 0.0
    ∂𝝭∂y[3] = 1.0
    ∂𝝭∂y[4] = 0.0
    ∂𝝭∂y[5] = x
    ∂𝝭∂y[6] = 2*y
    ∂𝝭∂y[7] = 0.0
    ∂𝝭∂y[8] = x^2
    ∂𝝭∂y[9] = 2*x*y
    ∂𝝭∂y[10] = 3*y^2
    ∂𝝭∂y[11] = 0.0
    ∂𝝭∂y[12] = x^3
    ∂𝝭∂y[13] = 3*x*y^2
    ∂𝝭∂y[14] = 2*x^2*y
    ∂𝝭∂y[15] = 4*y^3
    ∂𝝭∂y[16] = 0.0
    ∂𝝭∂y[17] = x^4
    ∂𝝭∂y[18] = 2*x^3*y
    ∂𝝭∂y[19] = 3*x^2*y^2
    ∂𝝭∂y[20] = 4*x*y^3
    ∂𝝭∂y[21] = 5*y^4
end
function set𝝭!(::PiecewisePolynomial{:Sextic2D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    x = 𝒙.x
    y = 𝒙.y
    𝝭[1] = 1.0
    𝝭[2] = x
    𝝭[3] = y
    𝝭[4] = x^2
    𝝭[5] = x*y
    𝝭[6] = y^2
    𝝭[7] = x^3
    𝝭[8] = x^2*y
    𝝭[9] = x*y^2
    𝝭[10] = y^3
    𝝭[11] = x^4
    𝝭[12] = x^3*y
    𝝭[13] = x*y^3
    𝝭[14] = x^2*y^2
    𝝭[15] = y^4
    𝝭[16] = x^5
    𝝭[17] = x^4*y
    𝝭[18] = x^3*y^2
    𝝭[19] = x^2*y^3
    𝝭[20] = x*y^4
    𝝭[21] = x^5
    𝝭[22] = x^6
    𝝭[23] = x^5*y
    𝝭[24] = x^4*y^2
    𝝭[25] = x^3*y^3
    𝝭[26] = x^2*y^4
    𝝭[27] = x*y^5
    𝝭[28] = y^6

end

function set∇𝝭!(::PiecewisePolynomial{:Sextic2D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    x = 𝒙.x
    y = 𝒙.y
    ∂𝝭∂x = 𝒙[:∂𝝭∂x]
    ∂𝝭∂y = 𝒙[:∂𝝭∂y]
    𝝭[1] = 1.0
    𝝭[2] = x
    𝝭[3] = y
    𝝭[4] = x^2
    𝝭[5] = x*y
    𝝭[6] = y^2
    𝝭[7] = x^3
    𝝭[8] = x^2*y
    𝝭[9] = x*y^2
    𝝭[10] = y^3
    𝝭[11] = x^4
    𝝭[12] = x^3*y
    𝝭[13] = x*y^3
    𝝭[14] = x^2*y^2
    𝝭[15] = y^4
    𝝭[16] = x^5
    𝝭[17] = x^4*y
    𝝭[18] = x^3*y^2
    𝝭[19] = x^2*y^3
    𝝭[20] = x*y^4
    𝝭[21] = y^5
    𝝭[22] = x^6
    𝝭[23] = x^5*y
    𝝭[24] = x^4*y^2
    𝝭[25] = x^3*y^3
    𝝭[26] = x^2*y^4
    𝝭[27] = x*y^5
    𝝭[28] = y^6
    ∂𝝭∂x[1] = 0.0
    ∂𝝭∂x[2] = 1.0
    ∂𝝭∂x[3] = 0.0
    ∂𝝭∂x[4] = 2*x
    ∂𝝭∂x[5] = y
    ∂𝝭∂x[6] = 0.0
    ∂𝝭∂x[7] = 3*x^2
    ∂𝝭∂x[8] = 2*x*y
    ∂𝝭∂x[9] = y^2
    ∂𝝭∂x[10] = 0.0
    ∂𝝭∂x[11] = 4*x^3
    ∂𝝭∂x[12] = 3*x^2*y
    ∂𝝭∂x[13] = y^3
    ∂𝝭∂x[14] = 2*x*y^2
    ∂𝝭∂x[15] = 0.0
    ∂𝝭∂x[16] = 5*x^4
    ∂𝝭∂x[17] = 4*x^3*y
    ∂𝝭∂x[18] = 3*x^2*y^2
    ∂𝝭∂x[19] = 2*x*y^3
    ∂𝝭∂x[20] = y^4
    ∂𝝭∂x[21] = 0.0
    ∂𝝭∂x[22] = 6*x^5
    ∂𝝭∂x[23] = 5*x^4*y
    ∂𝝭∂x[24] = 4*x^3*y^2
    ∂𝝭∂x[25] = 3*x^2*y^3
    ∂𝝭∂x[26] = 2*x*y^4
    ∂𝝭∂x[27] = y^5
    ∂𝝭∂x[28] = 0.0
    ∂𝝭∂y[1] = 0.0
    ∂𝝭∂y[2] = 0.0
    ∂𝝭∂y[3] = 1.0
    ∂𝝭∂y[4] = 0.0
    ∂𝝭∂y[5] = x
    ∂𝝭∂y[6] = 2*y
    ∂𝝭∂y[7] = 0.0
    ∂𝝭∂y[8] = x^2
    ∂𝝭∂y[9] = 2*x*y
    ∂𝝭∂y[10] = 3*y^2
    ∂𝝭∂y[11] = 0.0
    ∂𝝭∂y[12] = x^3
    ∂𝝭∂y[13] = 3*x*y^2
    ∂𝝭∂y[14] = 2*x^2*y
    ∂𝝭∂y[15] = 4*y^3
    ∂𝝭∂y[16] = 0.0
    ∂𝝭∂y[17] = x^4
    ∂𝝭∂y[18] = 2*x^3*y
    ∂𝝭∂y[19] = 3*x^2*y^2
    ∂𝝭∂y[20] = 4*x*y^3
    ∂𝝭∂y[21] = 5*y^4
    ∂𝝭∂y[22] = 0.0
    ∂𝝭∂y[23] = x^5
    ∂𝝭∂y[24] = 2*x^4*y
    ∂𝝭∂y[25] = 3*x^3*y^2
    ∂𝝭∂y[26] = 4*x^2*y^3
    ∂𝝭∂y[27] = 5*x*y^4
    ∂𝝭∂y[28] = 6*y^5
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
    ∂𝝭∂x = 𝒙[:∂𝝭∂x]
    ∂𝝭∂y = 𝒙[:∂𝝭∂y]
    ∂𝝭∂x[1] = - 2*ξ*(1.0-η^2)*𝒙.∂ξ∂x - 2*η*(1.0-ξ^2)*𝒙.∂η∂x
    ∂𝝭∂y[1] = - 2*ξ*(1.0-η^2)*𝒙.∂ξ∂y - 2*η*(1.0-ξ^2)*𝒙.∂η∂y
end
