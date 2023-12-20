struct PiecewisePolynomial{𝑝,T}<:AbstractElement{T}
    𝓒::Tuple{Int,Int,Vector{Node{(:𝐼,),1}}}
    𝓖::Tuple{Int,Int,Vector{Node{(:𝐼,),1}}}
end

struct PiecewiseParametric{𝑝,T}<:AbstractElement{T}
    𝓒::Tuple{Int,Int,Vector{Node{(:𝐼,),1}}}
    𝓖::Tuple{Int,Int,Vector{Node{(:𝐼,),1}}}
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