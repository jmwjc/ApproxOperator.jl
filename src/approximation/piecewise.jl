struct PiecewisePolynomial{𝑝,T}<:AbstractElement{T}
    𝓒::Tuple{Int,Int,Vector{Node{(:𝐼,),1}}}
    𝓖::Tuple{Int,Int,Vector{Node{(:𝐼,),1}}}
end

struct PiecewiseParametric{𝑝,T}<:AbstractElement{T}
    𝓒::Tuple{Int,Int,Vector{Node{(:𝐼,),1}}}
    𝓖::Tuple{Int,Int,Vector{Node{(:𝐼,),1}}}
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

function set𝝭!(ap::PiecewisePolynomial{:Constant1D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    𝝭[1] = 1.0
end

function set𝝭!(ap::PiecewisePolynomial{:Linear1D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    𝝭[1] = 1.0
    𝝭[2] = 𝒙.x
end

function set𝝭!(ap::PiecewisePolynomial{:Quadratic1D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    x = 𝒙.x
    𝝭[1] = 1.0
    𝝭[2] = x
    𝝭[3] = x^2
end

function set𝝭!(ap::PiecewiseParametric{:Constant1D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    𝝭[1] = 1.0
end

function set𝝭!(ap::PiecewiseParametric{:Linear1D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    𝝭[1] = 1.0
    𝝭[2] = 𝒙.ξ
end

function set𝝭!(ap::PiecewiseParametric{:Quadratic1D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    x = 𝒙.ξ
    𝝭[1] = 1.0
    𝝭[2] = ξ
    𝝭[3] = ξ^2
end

function set𝝭!(ap::PiecewisePolynomial{:Constant2D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    𝝭[1] = 1.0
end

function set𝝭!(ap::PiecewisePolynomial{:Linear2D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    𝝭[1] = 1.0
    𝝭[2] = 𝒙.x
    𝝭[3] = 𝒙.y
end

function set𝝭!(ap::PiecewisePolynomial{:Quadratic2D},𝒙::Node)
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

function set𝝭!(ap::PiecewiseParametric{:Constant2D},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    𝝭[1] = 1.0
end

function set𝝭!(ap::PiecewiseParametric{:Linear2D,:Tri3},𝒙::Node)
    𝝭 = 𝒙[:𝝭]
    𝝭[1] = 1.0
    𝝭[2] = 𝒙.ξ
    𝝭[3] = 𝒙.η
end

function set𝝭!(ap::PiecewiseParametric{:Quadratic2D},𝒙::Node)
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