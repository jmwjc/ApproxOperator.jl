# Hermite segment (SegHermite), 4 DOFs: (w, w') at the two end nodes.
#
# Gradient form. With ξ₁ = ½(1-ξ), ξ₂ = ½(1+ξ) we have dξ₁/dx = -½/𝐽 and
# dξ₂/dx = +½/𝐽, where 𝐽 = det(∂x/∂ξ) is the arc-length scale of the element.
#
# Note: this element stores values and first derivatives as separate DOFs, so
# its shape-function vector has 4 entries while `𝓒` has 2 nodes. The tests in
# test/shape_function_consistency.jl only cover C0 elements.

function set𝝭!(::Element{:SegHermite},x::Node)
    𝝭 = x[:𝝭]
    ξ₁ = 0.5*(1-x.ξ)
    ξ₂ = 0.5*(1+x.ξ)
    𝝭[1] = 3*ξ₁^2-2*ξ₁^3
    𝝭[2] = 3*ξ₂^2-2*ξ₂^3
    𝝭[3] = ξ₁^2*ξ₂
    𝝭[4] = ξ₁*ξ₂^2
end

function set∇𝝭!(ap::Element{:SegHermite},x::Node)
    𝝭 = x[:𝝭]
    ξ₁ = 0.5*(1-x.ξ)
    ξ₂ = 0.5*(1+x.ξ)
    𝝭[1] = 3*ξ₁^2-2*ξ₁^3
    𝝭[2] = 3*ξ₂^2-2*ξ₂^3
    𝝭[3] = ξ₁^2*ξ₂
    𝝭[4] = ξ₁*ξ₂^2

    𝐽 = ap.𝐽
    ∂ξ₁∂x = -0.5/𝐽
    ∂ξ₂∂x =  0.5/𝐽
    ∂𝝭∂x = x[:∂𝝭∂x]
    ∂𝝭∂x[1] = 6*ξ₁*(1-ξ₁)*∂ξ₁∂x
    ∂𝝭∂x[2] = 6*ξ₂*(1-ξ₂)*∂ξ₂∂x
    ∂𝝭∂x[3] = 2*ξ₁*ξ₂*∂ξ₁∂x + ξ₁^2*∂ξ₂∂x
    ∂𝝭∂x[4] = 2*ξ₁*ξ₂*∂ξ₂∂x + ξ₂^2*∂ξ₁∂x
end
