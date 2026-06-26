
struct ReproducingKernel{𝑝,𝑠,𝜙}<:AbstractReproducingKernel{𝑠,𝜙}
    𝓒::Vector{𝑿ᵢ}
    𝓖::Vector{𝑿ₛ}
end

get𝑛𝒑(    ::ReproducingKernel{:Linear1D}) = 2
get𝒑(     ::ReproducingKernel{:Linear1D},x::NTuple{3,Float64}) = (1.,x[1])
get∂𝒑∂x(  ::ReproducingKernel{:Linear1D},x::NTuple{3,Float64}) = (0.,1.)
get∂²𝒑∂x²(::ReproducingKernel{:Linear1D},x::NTuple{3,Float64}) = (0.,0.)

get𝑛𝒑(    ::ReproducingKernel{:Quadratic1D}) = 3
get𝒑(     ::ReproducingKernel{:Quadratic1D},x::NTuple{3,Float64}) = (1.,x[1],x[1]^2)
get∂𝒑∂x(  ::ReproducingKernel{:Quadratic1D},x::NTuple{3,Float64}) = (0.,1.,2*x[1])
get∂²𝒑∂x²(::ReproducingKernel{:Quadratic1D},x::NTuple{3,Float64}) = (0.,0.,2.)
get∂³𝒑∂x³(::ReproducingKernel{:Quadratic1D},x::NTuple{3,Float64}) = (0.,0.,0.)

get𝑛𝒑(    ::ReproducingKernel{:Cubic1D}) = 4
get𝒑(     ::ReproducingKernel{:Cubic1D},x::NTuple{3,Float64}) = (1.,x[1],x[1]^2,x[1]^3)
get∂𝒑∂x(  ::ReproducingKernel{:Cubic1D},x::NTuple{3,Float64}) = (0.,1.,2*x[1],3*x[1]^2)
get∂²𝒑∂x²(::ReproducingKernel{:Cubic1D},x::NTuple{3,Float64}) = (0.,0.,2.,6*x[1])
get∂³𝒑∂x³(::ReproducingKernel{:Cubic1D},x::NTuple{3,Float64}) = (0.,0.,0.,6.)

get𝑛𝒑(  ::ReproducingKernel{:Linear2D}) = 3
get𝒑(   ::ReproducingKernel{:Linear2D},x::NTuple{3,Float64}) = (1.,x[1],x[2])
get∂𝒑∂x(::ReproducingKernel{:Linear2D},x::NTuple{3,Float64}) = (0.,1.,0.)
get∂𝒑∂y(::ReproducingKernel{:Linear2D},x::NTuple{3,Float64}) = (0.,0.,1.)

get𝑛𝒑(      ::ReproducingKernel{:Quadratic2D}) = 6
get𝒑(       ::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (1.,x[1],x[2],x[1]^2,x[1]*x[2],x[2]^2)
get∂𝒑∂x(    ::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (0.,1.,0.,2*x[1],x[2],0.)
get∂𝒑∂y(    ::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (0.,0.,1.,0.,x[1],2*x[2])
get∂²𝒑∂x²(  ::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (0.,0.,0.,2.,0.,0.)
get∂²𝒑∂x∂y( ::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (0.,0.,0.,0.,1.,0.)
get∂²𝒑∂y²(  ::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (0.,0.,0.,0.,0.,2.)
get∂³𝒑∂x³(  ::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (0.,0.,0.,0.,0.,0.)
get∂³𝒑∂x²∂y(::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (0.,0.,0.,0.,0.,0.)
get∂³𝒑∂x∂y²(::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (0.,0.,0.,0.,0.,0.)
get∂³𝒑∂y³(  ::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (0.,0.,0.,0.,0.,0.)

get𝑛𝒑(::ReproducingKernel{:Cubic2D}) = 10
get𝒑(::ReproducingKernel{:Cubic2D},x::NTuple{3,Float64}) =
(
    1., x[1], x[2], x[1]^2, x[1]*x[2], x[2]^2, x[1]^3, x[1]^2*x[2], x[1]*x[2]^2, x[2]^3
)
get∂𝒑∂x(::ReproducingKernel{:Cubic2D},x::NTuple{3,Float64}) =
(
    0., 1., 0., 2*x[1], x[2], 0., 3*x[1]^2, 2*x[1]*x[2], x[2]^2, 0.
)
get∂𝒑∂y(::ReproducingKernel{:Cubic2D},x::NTuple{3,Float64}) =
(
    0., 0., 1., 0., x[1], 2*x[2], 0., x[1]^2, 2*x[1]*x[2], 3*x[2]^2
)
get∂²𝒑∂x²(::ReproducingKernel{:Cubic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 2., 0., 0., 6*x[1], 2*x[2], 0., 0.
)
get∂²𝒑∂x∂y(::ReproducingKernel{:Cubic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 1., 0., 0., 2*x[1], 2*x[2], 0.
)
get∂²𝒑∂y²(::ReproducingKernel{:Cubic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 2., 0., 0., 2*x[1], 6*x[2]
)
get∂³𝒑∂x³(::ReproducingKernel{:Cubic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 0., 6., 0., 0., 0.
)
get∂³𝒑∂x²∂y(::ReproducingKernel{:Cubic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 0., 0., 2., 0., 0.
)
get∂³𝒑∂x∂y²(::ReproducingKernel{:Cubic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 0., 0., 0., 2., 0.
)
get∂³𝒑∂y³(::ReproducingKernel{:Cubic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 0., 0., 0., 0., 6.
)

get𝑛𝒑(::ReproducingKernel{:Quartic2D}) = 15
get𝒑(::ReproducingKernel{:Quartic2D},x::NTuple{3,Float64}) =
(
    1., x[1], x[2], x[1]^2, x[1]*x[2], x[2]^2, x[1]^3, x[1]^2*x[2], x[1]*x[2]^2, x[2]^3, x[1]^4, x[1]^3*x[2], x[1]^2*x[2]^2, x[1]*x[2]^3, x[2]^4
)
get∂𝒑∂x(::ReproducingKernel{:Quartic2D},x::NTuple{3,Float64}) =
(
    0., 1., 0., 2*x[1], x[2], 0., 3*x[1]^2, 2*x[1]*x[2], x[2]^2, 0., 4.0*x[1]^3, 3.0*x[1]^2*x[2], 2.0*x[1]*x[2]^2, x[2]^3, 0.0
)
get∂𝒑∂y(::ReproducingKernel{:Quartic2D},x::NTuple{3,Float64}) =
(
    0., 0., 1., 0., x[1], 2*x[2], 0., x[1]^2, 2*x[1]*x[2], 3*x[2]^2, 0.0, x[1]^3, 2.0*x[1]^2*x[2], 3.0*x[1]*x[2]^2, 4.0*x[2]^3
)
get∂²𝒑∂x²(::ReproducingKernel{:Quartic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 2., 0., 0., 6*x[1], 2*x[2], 0., 0., 12.0*x[1]^2, 6.0*x[1]*x[2], 2.0*x[2]^2, 0.0, 0.0
)
get∂²𝒑∂x∂y(::ReproducingKernel{:Quartic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 1., 0., 0., 2*x[1], 2*x[2], 0., 0.0, 3.0*x[1]^2, 4.0*x[1]*x[2], 3.0*x[2]^2, 0.0
)
get∂²𝒑∂y²(::ReproducingKernel{:Quartic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 2., 0., 0., 2*x[1], 6*x[2], 0.0, 0.0, 2.0*x[1]^2, 6.0*x[1]*x[2], 12.0*x[2]^2
)
get∂³𝒑∂x³(::ReproducingKernel{:Quartic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 0., 6., 0., 0., 0., 24.0*x[1], 6.0*x[2], 0.0, 0.0, 0.0
)
get∂³𝒑∂x²∂y(::ReproducingKernel{:Quartic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0., 6.0*x[1], 4.0*x[2], 0., 0.
)
get∂³𝒑∂x∂y²(::ReproducingKernel{:Quartic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0., 4.0*x[1], 6.0*x[2],0.
)
get∂³𝒑∂y³(::ReproducingKernel{:Quartic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 0., 0., 0., 0., 6., 0., 0., 0., 6.0*x[1], 24.0*x[2]
)


get𝑛𝒑(  ::ReproducingKernel{:Linear3D}) = 4
get𝒑(   ::ReproducingKernel{:Linear3D},x::NTuple{3,Float64}) = (1.,x[1],x[2],x[3])
get∂𝒑∂x(::ReproducingKernel{:Linear3D},x::NTuple{3,Float64}) = (0.,1.,0.,0.)
get∂𝒑∂y(::ReproducingKernel{:Linear3D},x::NTuple{3,Float64}) = (0.,0.,1.,0.)
get∂𝒑∂z(::ReproducingKernel{:Linear3D},x::NTuple{3,Float64}) = (0.,0.,0.,1.)

get𝑛𝒑(      ::ReproducingKernel{:Quadratic3D}) = 10
get𝒑(       ::ReproducingKernel{:Quadratic3D},x::NTuple{3,Float64}) = (1.,x[1],x[2],x[3],x[1]^2,x[1]*x[2],x[1]*x[3],x[2]^2,x[2]*x[3],x[3]^2)
get∂𝒑∂x(    ::ReproducingKernel{:Quadratic3D},x::NTuple{3,Float64}) = (0.,1.,0.,0.,2*x[1],x[2],x[3],0.,0.,0.)
get∂𝒑∂y(    ::ReproducingKernel{:Quadratic3D},x::NTuple{3,Float64}) = (0.,0.,1.,0.,0.,x[1],0.,2*x[2],x[3],0.)

get∂𝒑∂z(::ReproducingKernel{:Quadratic3D},x::NTuple{3,Float64}) = (0.,0.,0.,1.,0.,0.,x[1],0.,x[2],2*x[3])
get∂²𝒑∂x²(  ::ReproducingKernel{:Quadratic3D},x::NTuple{3,Float64}) = (0.,0.,0.,0.,2.,0.,0.,0.,0.,0.)
get∂²𝒑∂x∂y( ::ReproducingKernel{:Quadratic3D},x::NTuple{3,Float64}) = (0.,0.,0.,0.,0.,1.,0.,0.,0.,0.)
get∂²𝒑∂x∂z( ::ReproducingKernel{:Quadratic3D},x::NTuple{3,Float64}) = (0.,0.,0.,0.,0.,0.,1.,0.,0.,0.)
get∂²𝒑∂y²(  ::ReproducingKernel{:Quadratic3D},x::NTuple{3,Float64}) = (0.,0.,0.,0.,0.,0.,0.,2.,0.,0.)
get∂²𝒑∂y∂z( ::ReproducingKernel{:Quadratic3D},x::NTuple{3,Float64}) = (0.,0.,0.,0.,0.,0.,0.,0.,1.,0.)

get∂²𝒑∂z²(  ::ReproducingKernel{:Quadratic3D},x::NTuple{3,Float64}) = (0.,0.,0.,0.,0.,0.,0.,0.,0.,2.)

function cal𝗠!(ap::AbstractReproducingKernel,x::Node)
    𝓒 = ap.𝓒
    𝗠 = get𝗠(ap,:𝗠)
    n = get𝑛𝒑(ap)
    for xᵢ in 𝓒
        Δx = x - xᵢ
        𝒑 = get𝒑(ap,Δx)
        𝜙 = get𝜙(ap,xᵢ,Δx)
        for I in 1:n
            for J in 1:I
                𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
            end
        end
    end
    cholesky!(𝗠)
    inverse!(𝗠)
    UUᵀ!(𝗠)
    return 𝗠
end

function get𝗠(ap::AbstractReproducingKernel,s::Symbol)
    n = get𝑛𝒑(ap)
    data = getfield(ap.𝓖[1],:data)
    fill!(data[s][2],0.)
    return SymMat(n,data[s][2])
end

function set𝝭!(ap::AbstractReproducingKernel,𝒙::Node)
    𝓒 = ap.𝓒
    𝝭 = 𝒙[:𝝭]
    𝒑₀ᵀ𝗠⁻¹ = cal𝗠!(ap,𝒙)
    for (i,𝒙ᵢ) in enumerate(𝓒)
        Δ𝒙 = 𝒙 - 𝒙ᵢ
        𝒑= get𝒑(ap,Δ𝒙)
        𝜙 = get𝜙(ap,𝒙ᵢ,Δ𝒙)
        𝝭[i] = 𝒑₀ᵀ𝗠⁻¹*𝒑*𝜙
    end
end

for 𝒑 in (:(:Linear1D),:(:Quadratic1D),:(:Cubic1D),:(:Quartic1D))
    @eval begin
        get∇𝒑( ap::ReproducingKernel{$𝒑},x::Any) = get𝒑(ap,x), get∂𝒑∂x(ap,x)
        get∇²𝒑(ap::ReproducingKernel{$𝒑},x::Any) = get𝒑(ap,x), get∂𝒑∂x(ap,x), get∂²𝒑∂x²(ap,x)
        get∇³𝒑(ap::ReproducingKernel{$𝒑},x::Any) = get𝒑(ap,x), get∂𝒑∂x(ap,x), get∂²𝒑∂x²(ap,x), get∂³𝒑∂x³(ap,x)

        function cal∇𝗠!(ap::ReproducingKernel{$𝒑},x::Node)
            𝓒 = ap.𝓒
            𝗠 = get𝗠(ap,:𝗠)
            ∂𝗠∂x = get𝗠(ap,:∂𝗠∂x)
            n = get𝑛𝒑(ap)
            for xᵢ in 𝓒
                Δx = x - xᵢ
                𝒑, ∂𝒑∂x = get∇𝒑(ap,Δx)
                𝜙, ∂𝜙∂x = get∇𝜙(ap,xᵢ,Δx)
                for I in 1:n
                    for J in 1:I
                        𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
                        ∂𝗠∂x[I,J] += ∂𝜙∂x*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂x[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂x[J]
                    end
                end
            end
            cholesky!(𝗠)
            U = inverse!(𝗠)
            ∂𝗠⁻¹∂x = - UUᵀAUUᵀ!(∂𝗠∂x,U)
            𝗠⁻¹ = UUᵀ!(U)
            return 𝗠⁻¹, ∂𝗠⁻¹∂x
        end

        function set∇𝝭!(ap::ReproducingKernel{$𝒑},𝒙::Node)
            𝓒 = ap.𝓒
            𝝭 = 𝒙[:𝝭]
            ∂𝝭∂x = 𝒙[:∂𝝭∂x]
            𝒑₀ᵀ𝗠⁻¹, 𝒑₀ᵀ∂𝗠⁻¹∂x = cal∇𝗠!(ap,𝒙)
            for (i,𝒙ᵢ) in enumerate(𝓒)
                Δ𝒙 = 𝒙 - 𝒙ᵢ
                𝒑, ∂𝒑∂x = get∇𝒑(ap,Δ𝒙)
                𝜙, ∂𝜙∂x = get∇𝜙(ap,𝒙ᵢ,Δ𝒙)
                𝝭[i] = 𝒑₀ᵀ𝗠⁻¹*𝒑*𝜙
                ∂𝝭∂x[i] = 𝒑₀ᵀ∂𝗠⁻¹∂x*𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂x*𝜙 + 𝒑₀ᵀ𝗠⁻¹*𝒑*∂𝜙∂x
            end
        end
    end
end

for 𝒑 in (:(:Linear2D),:(:Quadratic2D),:(:Cubic2D))
    @eval begin
        get∇𝒑( ap::ReproducingKernel{$𝒑},x::Any) = get𝒑(ap,x),get∂𝒑∂x(ap,x),get∂𝒑∂y(ap,x)
        get∇²𝒑(ap::ReproducingKernel{$𝒑},x::Any) = get𝒑(ap,x),get∂𝒑∂x(ap,x),get∂𝒑∂y(ap,x),get∂²𝒑∂x²(ap,x),get∂²𝒑∂x∂y(ap,x),get∂²𝒑∂y²(ap,x)
        get∇³𝒑(ap::ReproducingKernel{$𝒑},x::Any) = get𝒑(ap,x),get∂𝒑∂x(ap,x),get∂𝒑∂y(ap,x),get∂²𝒑∂x²(ap,x),get∂²𝒑∂x∂y(ap,x),get∂²𝒑∂y²(ap,x),get∂³𝒑∂x³(ap,x),get∂³𝒑∂x²∂y(ap,x),get∂³𝒑∂x∂y²(ap,x),get∂³𝒑∂y³(ap,x)

        function cal∇𝗠!(ap::ReproducingKernel{$𝒑},x::Node)
            𝓒 = ap.𝓒
            𝗠 = get𝗠(ap,:𝗠)
            ∂𝗠∂x = get𝗠(ap,:∂𝗠∂x)
            ∂𝗠∂y = get𝗠(ap,:∂𝗠∂y)
            n = get𝑛𝒑(ap)
            for xᵢ in 𝓒
                Δx = x - xᵢ
                𝒑, ∂𝒑∂x, ∂𝒑∂y = get∇𝒑(ap,Δx)
                𝜙, ∂𝜙∂x, ∂𝜙∂y = get∇𝜙(ap,xᵢ,Δx)
                for I in 1:n
                    for J in 1:I
                        𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
                        ∂𝗠∂x[I,J] += ∂𝜙∂x*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂x[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂x[J]
                        ∂𝗠∂y[I,J] += ∂𝜙∂y*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂y[J]
                    end
                end
            end
            cholesky!(𝗠)
            U = inverse!(𝗠)
            ∂𝗠⁻¹∂x = - UUᵀAUUᵀ!(∂𝗠∂x,U)
            ∂𝗠⁻¹∂y = - UUᵀAUUᵀ!(∂𝗠∂y,U)
            𝗠⁻¹ = UUᵀ!(U)
            return 𝗠⁻¹, ∂𝗠⁻¹∂x, ∂𝗠⁻¹∂y
        end

        function set∇𝝭!(ap::ReproducingKernel{$𝒑},𝒙::Node)
            𝓒 = ap.𝓒
            𝝭 = 𝒙[:𝝭]
            ∂𝝭∂x = 𝒙[:∂𝝭∂x]
            ∂𝝭∂y = 𝒙[:∂𝝭∂y]
            𝒑₀ᵀ𝗠⁻¹, 𝒑₀ᵀ∂𝗠⁻¹∂x, 𝒑₀ᵀ∂𝗠⁻¹∂y = cal∇𝗠!(ap,𝒙)
            for (i,𝒙ᵢ) in enumerate(𝓒)
                Δ𝒙 = 𝒙 - 𝒙ᵢ
                𝒑, ∂𝒑∂x, ∂𝒑∂y = get∇𝒑(ap,Δ𝒙)
                𝜙, ∂𝜙∂x, ∂𝜙∂y = get∇𝜙(ap,𝒙ᵢ,Δ𝒙)
                𝝭[i] = 𝒑₀ᵀ𝗠⁻¹*𝒑*𝜙
                ∂𝝭∂x[i] = 𝒑₀ᵀ∂𝗠⁻¹∂x*𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂x*𝜙 + 𝒑₀ᵀ𝗠⁻¹*𝒑*∂𝜙∂x
                ∂𝝭∂y[i] = 𝒑₀ᵀ∂𝗠⁻¹∂y*𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂y*𝜙 + 𝒑₀ᵀ𝗠⁻¹*𝒑*∂𝜙∂y
            end
        end

        function cal∇²𝗠!(ap::ReproducingKernel{$𝒑},x::Node)
            𝓒 = ap.𝓒
            𝗠 = get𝗠(ap,:𝗠)
            ∂𝗠∂x = get𝗠(ap,:∂𝗠∂x)
            ∂𝗠∂y = get𝗠(ap,:∂𝗠∂y)
            ∂²𝗠∂x² = get𝗠(ap,:∂²𝗠∂x²)
            ∂²𝗠∂y² = get𝗠(ap,:∂²𝗠∂y²)
            ∂²𝗠∂x∂y = get𝗠(ap,:∂²𝗠∂x∂y)
            n = get𝑛𝒑(ap)
            for xᵢ in 𝓒
                Δx = x - xᵢ
                𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂²𝒑∂x², ∂²𝒑∂x∂y, ∂²𝒑∂y² = get∇²𝒑(ap,Δx)
                𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂²𝜙∂x², ∂²𝜙∂x∂y, ∂²𝜙∂y² = get∇²𝜙(ap,xᵢ,Δx)
                for I in 1:n
                    for J in 1:I
                        𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
                        ∂𝗠∂x[I,J] += ∂𝜙∂x*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂x[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂x[J]
                        ∂𝗠∂y[I,J] += ∂𝜙∂y*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂y[J]
                        ∂²𝗠∂x²[I,J] += ∂²𝜙∂x²*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂x²[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂x²[J] + 2.0*∂𝜙∂x*∂𝒑∂x[I]*𝒑[J] + 2.0*∂𝜙∂x*𝒑[I]*∂𝒑∂x[J] + 2.0*𝜙*∂𝒑∂x[I]*∂𝒑∂x[J]
                        ∂²𝗠∂y²[I,J] += ∂²𝜙∂y²*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂y²[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂y²[J] + 2.0*∂𝜙∂y*∂𝒑∂y[I]*𝒑[J] + 2.0*∂𝜙∂y*𝒑[I]*∂𝒑∂y[J] + 2.0*𝜙*∂𝒑∂y[I]*∂𝒑∂y[J]
                        ∂²𝗠∂x∂y[I,J] += ∂²𝜙∂x∂y*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂x∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂x∂y[J] + ∂𝜙∂x*∂𝒑∂y[I]*𝒑[J] + ∂𝜙∂y*∂𝒑∂x[I]*𝒑[J] + ∂𝜙∂x*𝒑[I]*∂𝒑∂y[J] + ∂𝜙∂y*𝒑[I]*∂𝒑∂x[J] + 𝜙*∂𝒑∂x[I]*∂𝒑∂y[J] + 𝜙*∂𝒑∂y[I]*∂𝒑∂x[J]
                    end
                end
            end
            cholesky!(𝗠)
            U = inverse!(𝗠)
            Uᵀ∂𝗠∂xU = UᵀAU!(∂𝗠∂x,U)
            Uᵀ∂𝗠∂yU = UᵀAU!(∂𝗠∂y,U)
            Uᵀ∂²𝗠∂x²U = UᵀAU!(∂²𝗠∂x²,U)
            Uᵀ∂²𝗠∂y²U = UᵀAU!(∂²𝗠∂y²,U)
            Uᵀ∂²𝗠∂x∂yU = UᵀAU!(∂²𝗠∂x∂y,U)
            for i in 1:n
                for j in 1:i
                    for k in 1:n
                        Uᵀ∂²𝗠∂x²U[i,j] -= 2*Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂𝗠∂xU[k,j]
                        Uᵀ∂²𝗠∂y²U[i,j] -= 2*Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂𝗠∂yU[k,j]
                        Uᵀ∂²𝗠∂x∂yU[i,j] -= Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂𝗠∂yU[k,j] + Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂𝗠∂xU[k,j]
                    end
                end
            end
    
            ∂²𝗠⁻¹∂x² = - UAUᵀ!(Uᵀ∂²𝗠∂x²U,U)
            ∂²𝗠⁻¹∂y² = - UAUᵀ!(Uᵀ∂²𝗠∂y²U,U)
            ∂²𝗠⁻¹∂x∂y = - UAUᵀ!(Uᵀ∂²𝗠∂x∂yU,U)
            ∂𝗠⁻¹∂x = - UAUᵀ!(Uᵀ∂𝗠∂xU,U)
            ∂𝗠⁻¹∂y = - UAUᵀ!(Uᵀ∂𝗠∂yU,U)
            𝗠⁻¹ = UUᵀ!(U)
            return 𝗠⁻¹, ∂𝗠⁻¹∂x, ∂𝗠⁻¹∂y, ∂²𝗠⁻¹∂x², ∂²𝗠⁻¹∂x∂y, ∂²𝗠⁻¹∂y²
        end

        function cal∇³𝗠!(ap::ReproducingKernel{$𝒑},x::Node)
            𝓒 = ap.𝓒
            𝗠 = get𝗠(ap,:𝗠)
            ∂𝗠∂x = get𝗠(ap,:∂𝗠∂x)
            ∂𝗠∂y = get𝗠(ap,:∂𝗠∂y)
            ∂²𝗠∂x² = get𝗠(ap,:∂²𝗠∂x²)
            ∂²𝗠∂x∂y = get𝗠(ap,:∂²𝗠∂x∂y)
            ∂²𝗠∂y² = get𝗠(ap,:∂²𝗠∂y²)
            ∂³𝗠∂x³ = get𝗠(ap,:∂³𝗠∂x³)
            ∂³𝗠∂x²∂y = get𝗠(ap,:∂³𝗠∂x²∂y)
            ∂³𝗠∂x∂y² = get𝗠(ap,:∂³𝗠∂x∂y²)
            ∂³𝗠∂y³ = get𝗠(ap,:∂³𝗠∂y³)
            n = get𝑛𝒑(ap)
            for xᵢ in 𝓒
                Δx = x - xᵢ
                𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂²𝒑∂x², ∂²𝒑∂x∂y, ∂²𝒑∂y², ∂³𝒑∂x³, ∂³𝒑∂x²∂y, ∂³𝒑∂x∂y², ∂³𝒑∂y³ = get∇³𝒑(ap,Δx)
                𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂²𝜙∂x², ∂²𝜙∂x∂y, ∂²𝜙∂y², ∂³𝜙∂x³, ∂³𝜙∂x²∂y, ∂³𝜙∂x∂y², ∂³𝜙∂y³ = get∇³𝜙(ap,xᵢ,Δx)
                for I in 1:n
                    for J in I:n
                        𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
                        ∂𝗠∂x[I,J] += ∂𝜙∂x*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂x[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂x[J]
                        ∂𝗠∂y[I,J] += ∂𝜙∂y*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂y[J]

                        ∂²𝗠∂x²[I,J] += ∂²𝜙∂x²*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂x²[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂x²[J] + 2.0*∂𝜙∂x*∂𝒑∂x[I]*𝒑[J] + 2.0*∂𝜙∂x*𝒑[I]*∂𝒑∂x[J] + 2.0*𝜙*∂𝒑∂x[I]*∂𝒑∂x[J]

                        ∂²𝗠∂x∂y[I,J] += ∂²𝜙∂x∂y*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂x∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂x∂y[J] + ∂𝜙∂x*∂𝒑∂y[I]*𝒑[J] + ∂𝜙∂y*∂𝒑∂x[I]*𝒑[J] + ∂𝜙∂x*𝒑[I]*∂𝒑∂y[J] + ∂𝜙∂y*𝒑[I]*∂𝒑∂x[J] + 𝜙*∂𝒑∂x[I]*∂𝒑∂y[J] + 𝜙*∂𝒑∂y[I]*∂𝒑∂x[J]

                        ∂²𝗠∂y²[I,J] += ∂²𝜙∂y²*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂y²[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂y²[J] + 2.0*∂𝜙∂y*∂𝒑∂y[I]*𝒑[J] + 2.0*∂𝜙∂y*𝒑[I]*∂𝒑∂y[J] + 2.0*𝜙*∂𝒑∂y[I]*∂𝒑∂y[J]

                        ∂³𝗠∂x³[I,J] += ∂³𝜙∂x³*𝒑[I]*𝒑[J] + 𝜙*∂³𝒑∂x³[I]*𝒑[J] + 𝜙*𝒑[I]*∂³𝒑∂x³[J] + 3.0*∂²𝜙∂x²*∂𝒑∂x[I]*𝒑[J] + 3.0*∂𝜙∂x*∂²𝒑∂x²[I]*𝒑[J] + 3.0*∂²𝜙∂x²*𝒑[I]*∂𝒑∂x[J] + 3.0*∂𝜙∂x*𝒑[I]*∂²𝒑∂x²[J] + 3.0*𝜙*∂²𝒑∂x²[I]*∂𝒑∂x[J] + 3.0*𝜙*∂𝒑∂x[I]*∂²𝒑∂x²[J] + 6.0*∂𝜙∂x*∂𝒑∂x[I]*∂𝒑∂x[J]

                        ∂³𝗠∂x²∂y[I,J] += ∂³𝜙∂x²∂y*𝒑[I]*𝒑[J] + 𝜙*∂³𝒑∂x²∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂³𝒑∂x²∂y[J] + 2.0*∂²𝜙∂x∂y*∂𝒑∂x[I]*𝒑[J] + ∂²𝜙∂x²*∂𝒑∂y[I]*𝒑[J] + 2.0*∂𝜙∂x*∂²𝒑∂x∂y[I]*𝒑[J] + ∂𝜙∂y*∂²𝒑∂x²[I]*𝒑[J] + 2.0*∂²𝜙∂x∂y*𝒑[I]*∂𝒑∂x[J] + ∂²𝜙∂x²*𝒑[I]*∂𝒑∂y[J] + 2.0*∂𝜙∂x*𝒑[I]*∂²𝒑∂x∂y[J] + ∂𝜙∂y*𝒑[I]*∂²𝒑∂x²[J] + 2.0*𝜙*∂²𝒑∂x∂y[I]*∂𝒑∂x[J] + 𝜙*∂²𝒑∂x²[I]*∂𝒑∂y[J] + 2.0*𝜙*∂𝒑∂x[I]*∂²𝒑∂x∂y[J] + 𝜙*∂𝒑∂y[I]*∂²𝒑∂x²[J] + 2.0*∂𝜙∂y*∂𝒑∂x[I]*∂𝒑∂x[J] + 2.0*∂𝜙∂x*∂𝒑∂y[I]*∂𝒑∂x[J] + 2.0*∂𝜙∂x*∂𝒑∂x[I]*∂𝒑∂y[J]

                        ∂³𝗠∂x∂y²[I,J] += ∂³𝜙∂x∂y²*𝒑[I]*𝒑[J] + 𝜙*∂³𝒑∂x∂y²[I]*𝒑[J] + 𝜙*𝒑[I]*∂³𝒑∂x∂y²[J] + 2.0*∂²𝜙∂x∂y*∂𝒑∂y[I]*𝒑[J] + ∂²𝜙∂y²*∂𝒑∂x[I]*𝒑[J] + 2.0*∂𝜙∂y*∂²𝒑∂x∂y[I]*𝒑[J] + ∂𝜙∂x*∂²𝒑∂y²[I]*𝒑[J] + 2.0*∂²𝜙∂x∂y*𝒑[I]*∂𝒑∂y[J] + ∂²𝜙∂y²*𝒑[I]*∂𝒑∂x[J] + 2.0*∂𝜙∂y*𝒑[I]*∂²𝒑∂x∂y[J] + ∂𝜙∂x*𝒑[I]*∂²𝒑∂y²[J] + 2.0*𝜙*∂²𝒑∂x∂y[I]*∂𝒑∂y[J] + 𝜙*∂²𝒑∂y²[I]*∂𝒑∂x[J] + 2.0*𝜙*∂𝒑∂y[I]*∂²𝒑∂x∂y[J] + 𝜙*∂𝒑∂x[I]*∂²𝒑∂y²[J] + 2.0*∂𝜙∂x*∂𝒑∂y[I]*∂𝒑∂y[J] + 2.0*∂𝜙∂y*∂𝒑∂x[I]*∂𝒑∂y[J] + 2.0*∂𝜙∂y*∂𝒑∂y[I]*∂𝒑∂x[J]

                        ∂³𝗠∂y³[I,J] += ∂³𝜙∂y³*𝒑[I]*𝒑[J] + 𝜙*∂³𝒑∂y³[I]*𝒑[J] + 𝜙*𝒑[I]*∂³𝒑∂y³[J] + 3.0*∂²𝜙∂y²*∂𝒑∂y[I]*𝒑[J] + 3.0*∂𝜙∂y*∂²𝒑∂y²[I]*𝒑[J] + 3.0*∂²𝜙∂y²*𝒑[I]*∂𝒑∂y[J] + 3.0*∂𝜙∂y*𝒑[I]*∂²𝒑∂y²[J] + 3.0*𝜙*∂²𝒑∂y²[I]*∂𝒑∂y[J] + 3.0*𝜙*∂𝒑∂y[I]*∂²𝒑∂y²[J] + 6.0*∂𝜙∂y*∂𝒑∂y[I]*∂𝒑∂y[J]
                    end
                end
            end
            cholesky!(𝗠)
            U = inverse!(𝗠)
            Uᵀ∂𝗠∂xU = UᵀAU!(∂𝗠∂x,U)
            Uᵀ∂𝗠∂yU = UᵀAU!(∂𝗠∂y,U)
            Uᵀ∂²𝗠∂x²U = UᵀAU!(∂²𝗠∂x²,U)
            Uᵀ∂²𝗠∂y²U = UᵀAU!(∂²𝗠∂y²,U)
            Uᵀ∂²𝗠∂x∂yU = UᵀAU!(∂²𝗠∂x∂y,U)
            Uᵀ∂³𝗠∂x³U = UᵀAU!(∂³𝗠∂x³,U)
            Uᵀ∂³𝗠∂x²∂yU = UᵀAU!(∂³𝗠∂x²∂y,U)
            Uᵀ∂³𝗠∂x∂y²U = UᵀAU!(∂³𝗠∂x∂y²,U)
            Uᵀ∂³𝗠∂y³U = UᵀAU!(∂³𝗠∂y³,U)

            for i in 1:n
                for j in 1:i
                    for k in 1:n
                        Uᵀ∂³𝗠∂x³U[i,j] -= 3*Uᵀ∂²𝗠∂x²U[i,k]*Uᵀ∂𝗠∂xU[k,j]
                        Uᵀ∂³𝗠∂x²∂yU[i,j] -= 2*Uᵀ∂²𝗠∂x∂yU[i,k]*Uᵀ∂𝗠∂xU[k,j]+Uᵀ∂²𝗠∂x²U[i,k]*Uᵀ∂𝗠∂yU[k,j]
                        Uᵀ∂³𝗠∂x∂y²U[i,j] -= 2*Uᵀ∂²𝗠∂x∂yU[i,k]*Uᵀ∂𝗠∂yU[k,j]+Uᵀ∂²𝗠∂y²U[i,k]*Uᵀ∂𝗠∂xU[k,j]
                        Uᵀ∂³𝗠∂y³U[i,j] -= 3*Uᵀ∂²𝗠∂y²U[i,k]*Uᵀ∂𝗠∂yU[k,j]
                    end
                end
                for j in 1:i
                    for k in 1:n
                        Uᵀ∂²𝗠∂x²U[i,j] -= 2*Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂𝗠∂xU[k,j]
                        Uᵀ∂²𝗠∂y²U[i,j] -= 2*Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂𝗠∂yU[k,j]
                        Uᵀ∂²𝗠∂x∂yU[i,j] -= Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂𝗠∂yU[k,j] + Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂𝗠∂xU[k,j]
                    end
                end
            end
            for i in 1:n
                for j in 1:i
                    for k in 1:n
                        Uᵀ∂³𝗠∂x³U[i,j] -= 3*Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂²𝗠∂x²U[k,j]
                        Uᵀ∂³𝗠∂x²∂yU[i,j] -= 2*Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂²𝗠∂x∂yU[k,j]+Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂²𝗠∂x²U[k,j]
                        Uᵀ∂³𝗠∂x∂y²U[i,j] -= 2*Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂²𝗠∂x∂yU[k,j]+Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂²𝗠∂y²U[k,j]
                        Uᵀ∂³𝗠∂y³U[i,j] -= 3*Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂²𝗠∂y²U[k,j]
                    end
                end
            end

            ∂³𝗠⁻¹∂x³ = - UAUᵀ!(Uᵀ∂³𝗠∂x³U,U)
            ∂³𝗠⁻¹∂x²∂y = - UAUᵀ!(Uᵀ∂³𝗠∂x²∂yU,U)
            ∂³𝗠⁻¹∂x∂y² = - UAUᵀ!(Uᵀ∂³𝗠∂x∂y²U,U)
            ∂³𝗠⁻¹∂y³ = - UAUᵀ!(Uᵀ∂³𝗠∂y³U,U)
            ∂²𝗠⁻¹∂x² = - UAUᵀ!(Uᵀ∂²𝗠∂x²U,U)
            ∂²𝗠⁻¹∂y² = - UAUᵀ!(Uᵀ∂²𝗠∂y²U,U)
            ∂²𝗠⁻¹∂x∂y = - UAUᵀ!(Uᵀ∂²𝗠∂x∂yU,U)
            ∂𝗠⁻¹∂x = - UAUᵀ!(Uᵀ∂𝗠∂xU,U)
            ∂𝗠⁻¹∂y = - UAUᵀ!(Uᵀ∂𝗠∂yU,U)
            𝗠⁻¹ = UUᵀ!(U)

            return 𝗠⁻¹, ∂𝗠⁻¹∂x, ∂𝗠⁻¹∂y, ∂²𝗠⁻¹∂x², ∂²𝗠⁻¹∂x∂y, ∂²𝗠⁻¹∂y², ∂³𝗠⁻¹∂x³, ∂³𝗠⁻¹∂x²∂y, ∂³𝗠⁻¹∂x∂y², ∂³𝗠⁻¹∂y³
        end

        function set∇²𝝭!(ap::ReproducingKernel{$𝒑},𝒙::Node)
            𝓒 = ap.𝓒
            𝝭 = 𝒙[:𝝭]
            ∂𝝭∂x = 𝒙[:∂𝝭∂x]
            ∂𝝭∂y = 𝒙[:∂𝝭∂y]
            ∂²𝝭∂x² = 𝒙[:∂²𝝭∂x²]
            ∂²𝝭∂y² = 𝒙[:∂²𝝭∂y²]
            ∂²𝝭∂x∂y = 𝒙[:∂²𝝭∂x∂y]
            𝒑₀ᵀ𝗠⁻¹, 𝒑₀ᵀ∂𝗠⁻¹∂x, 𝒑₀ᵀ∂𝗠⁻¹∂y, 𝒑₀ᵀ∂²𝗠⁻¹∂x², 𝒑₀ᵀ∂²𝗠⁻¹∂x∂y, 𝒑₀ᵀ∂²𝗠⁻¹∂y² = cal∇²𝗠!(ap,𝒙)
            for (i,𝒙ᵢ) in enumerate(𝓒)
                Δ𝒙 = 𝒙 - 𝒙ᵢ
                𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂²𝒑∂x², ∂²𝒑∂x∂y, ∂²𝒑∂y² = get∇²𝒑(ap,Δ𝒙)
                𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂²𝜙∂x², ∂²𝜙∂x∂y, ∂²𝜙∂y² = get∇²𝜙(ap,𝒙ᵢ,Δ𝒙)
                𝒑₀ᵀ𝗠⁻¹𝒑 = 𝒑₀ᵀ𝗠⁻¹*𝒑
                𝒑₀ᵀ∂𝗠⁻¹∂x𝒑 = 𝒑₀ᵀ∂𝗠⁻¹∂x*𝒑
                𝒑₀ᵀ∂𝗠⁻¹∂y𝒑 = 𝒑₀ᵀ∂𝗠⁻¹∂y*𝒑
                𝒑₀ᵀ𝗠⁻¹∂𝒑∂x = 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂x
                𝒑₀ᵀ𝗠⁻¹∂𝒑∂y = 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂y
                𝒑₀ᵀ∂²𝗠⁻¹∂x²𝒑 = 𝒑₀ᵀ∂²𝗠⁻¹∂x²*𝒑
                𝒑₀ᵀ∂²𝗠⁻¹∂y²𝒑 = 𝒑₀ᵀ∂²𝗠⁻¹∂y²*𝒑
                𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x² = 𝒑₀ᵀ𝗠⁻¹*∂²𝒑∂x²
                𝒑₀ᵀ𝗠⁻¹∂²𝒑∂y² = 𝒑₀ᵀ𝗠⁻¹*∂²𝒑∂y²
                𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂x = 𝒑₀ᵀ∂𝗠⁻¹∂x*∂𝒑∂x
                𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂y = 𝒑₀ᵀ∂𝗠⁻¹∂y*∂𝒑∂y
                𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂y = 𝒑₀ᵀ∂𝗠⁻¹∂x*∂𝒑∂y
                𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂x = 𝒑₀ᵀ∂𝗠⁻¹∂y*∂𝒑∂x
                𝒑₀ᵀ∂²𝗠⁻¹∂x∂y𝒑 = 𝒑₀ᵀ∂²𝗠⁻¹∂x∂y*𝒑
                𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x∂y = 𝒑₀ᵀ𝗠⁻¹*∂²𝒑∂x∂y

                𝝭[i] = 𝒑₀ᵀ𝗠⁻¹𝒑*𝜙
                ∂𝝭∂x[i] = 𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂𝜙∂x
                ∂𝝭∂y[i] = 𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂𝜙∂y

                ∂²𝝭∂x²[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂x²𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x²*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂x² + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂x*𝜙 + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*∂𝜙∂x + 2.0*𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*∂𝜙∂x

                ∂²𝝭∂y²[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂y²𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂y²*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂y² + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂y*𝜙 + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*∂𝜙∂y + 2.0*𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*∂𝜙∂y

                ∂²𝝭∂x∂y[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂x∂y𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x∂y*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂x∂y + 𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂y*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂x*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*∂𝜙∂y + 𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*∂𝜙∂x + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*∂𝜙∂y +𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*∂𝜙∂x
            end
        end

        function set∇³𝝭!(ap::ReproducingKernel{$𝒑},𝒙::Node)
            𝓒 = ap.𝓒
            𝝭 = 𝒙[:𝝭]
            ∂𝝭∂x = 𝒙[:∂𝝭∂x]
            ∂𝝭∂y = 𝒙[:∂𝝭∂y]
            ∂²𝝭∂x² = 𝒙[:∂²𝝭∂x²]
            ∂²𝝭∂y² = 𝒙[:∂²𝝭∂y²]
            ∂²𝝭∂x∂y = 𝒙[:∂²𝝭∂x∂y]
            ∂³𝝭∂x³ = 𝒙[:∂³𝝭∂x³]
            ∂³𝝭∂x²∂y = 𝒙[:∂³𝝭∂x²∂y]
            ∂³𝝭∂x∂y² = 𝒙[:∂³𝝭∂x∂y²]
            ∂³𝝭∂y³ = 𝒙[:∂³𝝭∂y³]
            𝒑₀ᵀ𝗠⁻¹, 𝒑₀ᵀ∂𝗠⁻¹∂x, 𝒑₀ᵀ∂𝗠⁻¹∂y, 𝒑₀ᵀ∂²𝗠⁻¹∂x², 𝒑₀ᵀ∂²𝗠⁻¹∂x∂y, 𝒑₀ᵀ∂²𝗠⁻¹∂y², 𝒑₀ᵀ∂³𝗠⁻¹∂x³, 𝒑₀ᵀ∂³𝗠⁻¹∂x²∂y, 𝒑₀ᵀ∂³𝗠⁻¹∂x∂y², 𝒑₀ᵀ∂³𝗠⁻¹∂y³ = cal∇³𝗠!(ap,𝒙)
            for (i,𝒙ᵢ) in enumerate(𝓒)
                Δ𝒙 = 𝒙 - 𝒙ᵢ
                𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂²𝒑∂x², ∂²𝒑∂x∂y, ∂²𝒑∂y², ∂³𝒑∂x³, ∂³𝒑∂x²∂y, ∂³𝒑∂x∂y², ∂³𝒑∂y³ = get∇³𝒑(ap,Δ𝒙)
                𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂²𝜙∂x², ∂²𝜙∂x∂y, ∂²𝜙∂y², ∂³𝜙∂x³, ∂³𝜙∂x²∂y, ∂³𝜙∂x∂y², ∂³𝜙∂y³ = get∇³𝜙(ap,𝒙ᵢ,Δ𝒙)
                𝒑₀ᵀ𝗠⁻¹𝒑 = 𝒑₀ᵀ𝗠⁻¹*𝒑
                𝒑₀ᵀ∂𝗠⁻¹∂x𝒑 = 𝒑₀ᵀ∂𝗠⁻¹∂x*𝒑
                𝒑₀ᵀ∂𝗠⁻¹∂y𝒑 = 𝒑₀ᵀ∂𝗠⁻¹∂y*𝒑
                𝒑₀ᵀ𝗠⁻¹∂𝒑∂x = 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂x
                𝒑₀ᵀ𝗠⁻¹∂𝒑∂y = 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂y
                𝒑₀ᵀ∂²𝗠⁻¹∂x²𝒑 = 𝒑₀ᵀ∂²𝗠⁻¹∂x²*𝒑
                𝒑₀ᵀ∂²𝗠⁻¹∂x∂y𝒑 = 𝒑₀ᵀ∂²𝗠⁻¹∂x∂y*𝒑
                𝒑₀ᵀ∂²𝗠⁻¹∂y²𝒑 = 𝒑₀ᵀ∂²𝗠⁻¹∂y²*𝒑
                𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x² = 𝒑₀ᵀ𝗠⁻¹*∂²𝒑∂x²
                𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x∂y = 𝒑₀ᵀ𝗠⁻¹*∂²𝒑∂x∂y
                𝒑₀ᵀ𝗠⁻¹∂²𝒑∂y² = 𝒑₀ᵀ𝗠⁻¹*∂²𝒑∂y²
                𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂x = 𝒑₀ᵀ∂𝗠⁻¹∂x*∂𝒑∂x
                𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂y = 𝒑₀ᵀ∂𝗠⁻¹∂x*∂𝒑∂y
                𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂x = 𝒑₀ᵀ∂𝗠⁻¹∂y*∂𝒑∂x
                𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂y = 𝒑₀ᵀ∂𝗠⁻¹∂y*∂𝒑∂y
                𝒑₀ᵀ∂³𝗠⁻¹∂x³𝒑 = 𝒑₀ᵀ∂³𝗠⁻¹∂x³*𝒑
                𝒑₀ᵀ∂³𝗠⁻¹∂x²∂y𝒑 = 𝒑₀ᵀ∂³𝗠⁻¹∂x²∂y*𝒑
                𝒑₀ᵀ∂³𝗠⁻¹∂x∂y²𝒑 = 𝒑₀ᵀ∂³𝗠⁻¹∂x∂y²*𝒑
                𝒑₀ᵀ∂³𝗠⁻¹∂y³𝒑 = 𝒑₀ᵀ∂³𝗠⁻¹∂y³*𝒑
                𝒑₀ᵀ𝗠⁻¹∂³𝒑∂x³ = 𝒑₀ᵀ𝗠⁻¹*∂³𝒑∂x³
                𝒑₀ᵀ𝗠⁻¹∂³𝒑∂x²∂y = 𝒑₀ᵀ𝗠⁻¹*∂³𝒑∂x²∂y
                𝒑₀ᵀ𝗠⁻¹∂³𝒑∂x∂y² = 𝒑₀ᵀ𝗠⁻¹*∂³𝒑∂x∂y²
                𝒑₀ᵀ𝗠⁻¹∂³𝒑∂y³ = 𝒑₀ᵀ𝗠⁻¹*∂³𝒑∂y³
                𝒑₀ᵀ∂²𝗠⁻¹∂x²∂𝒑∂x = 𝒑₀ᵀ∂²𝗠⁻¹∂x²*∂𝒑∂x
                𝒑₀ᵀ∂𝗠⁻¹∂x∂²𝒑∂x² = 𝒑₀ᵀ∂𝗠⁻¹∂x*∂²𝒑∂x²
                𝒑₀ᵀ∂²𝗠⁻¹∂x²∂𝒑∂y = 𝒑₀ᵀ∂²𝗠⁻¹∂x²*∂𝒑∂y
                𝒑₀ᵀ∂𝗠⁻¹∂y∂²𝒑∂x² = 𝒑₀ᵀ∂𝗠⁻¹∂y*∂²𝒑∂x²
                𝒑₀ᵀ∂²𝗠⁻¹∂x∂y∂𝒑∂x = 𝒑₀ᵀ∂²𝗠⁻¹∂x∂y*∂𝒑∂x
                𝒑₀ᵀ∂𝗠⁻¹∂x∂²𝒑∂x∂y = 𝒑₀ᵀ∂𝗠⁻¹∂x*∂²𝒑∂x∂y
                𝒑₀ᵀ∂²𝗠⁻¹∂y²∂𝒑∂x = 𝒑₀ᵀ∂²𝗠⁻¹∂y²*∂𝒑∂x
                𝒑₀ᵀ∂𝗠⁻¹∂x∂²𝒑∂y² = 𝒑₀ᵀ∂𝗠⁻¹∂x*∂²𝒑∂y²
                𝒑₀ᵀ∂²𝗠⁻¹∂x∂y∂𝒑∂y = 𝒑₀ᵀ∂²𝗠⁻¹∂x∂y*∂𝒑∂y
                𝒑₀ᵀ∂𝗠⁻¹∂y∂²𝒑∂x∂y = 𝒑₀ᵀ∂𝗠⁻¹∂y*∂²𝒑∂x∂y
                𝒑₀ᵀ∂²𝗠⁻¹∂y²∂𝒑∂y = 𝒑₀ᵀ∂²𝗠⁻¹∂y²*∂𝒑∂y
                𝒑₀ᵀ∂𝗠⁻¹∂y∂²𝒑∂y² = 𝒑₀ᵀ∂𝗠⁻¹∂y*∂²𝒑∂y²

                𝝭[i] = 𝒑₀ᵀ𝗠⁻¹𝒑*𝜙
                ∂𝝭∂x[i] = 𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂𝜙∂x
                ∂𝝭∂y[i] = 𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂𝜙∂y

                ∂²𝝭∂x²[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂x²𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x²*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂x² + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂x*𝜙 + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*∂𝜙∂x + 2.0*𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*∂𝜙∂x

                ∂²𝝭∂y²[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂y²𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂y²*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂y² + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂y*𝜙 + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*∂𝜙∂y + 2.0*𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*∂𝜙∂y

                ∂²𝝭∂x∂y[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂x∂y𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x∂y*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂x∂y + 𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂y*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂x*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*∂𝜙∂y + 𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*∂𝜙∂x + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*∂𝜙∂y +𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*∂𝜙∂x

                ∂³𝝭∂x³[i] = 𝒑₀ᵀ∂³𝗠⁻¹∂x³𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂³𝒑∂x³*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂³𝜙∂x³ + 3*𝒑₀ᵀ∂²𝗠⁻¹∂x²∂𝒑∂x*𝜙  + 3*𝒑₀ᵀ∂𝗠⁻¹∂x∂²𝒑∂x²*𝜙 + 3*𝒑₀ᵀ∂²𝗠⁻¹∂x²𝒑*∂𝜙∂x + 3*𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*∂²𝜙∂x² + 3*𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x²*∂𝜙∂x + 3*𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*∂²𝜙∂x² + 6*𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂x*∂𝜙∂x

                ∂³𝝭∂x²∂y[i] = 𝒑₀ᵀ∂³𝗠⁻¹∂x²∂y𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂³𝒑∂x²∂y*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂³𝜙∂x²∂y + 2*𝒑₀ᵀ∂²𝗠⁻¹∂x∂y∂𝒑∂x*𝜙 + 𝒑₀ᵀ∂²𝗠⁻¹∂x²∂𝒑∂y*𝜙 + 2*𝒑₀ᵀ∂𝗠⁻¹∂x∂²𝒑∂x∂y*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂y∂²𝒑∂x²*𝜙 + 2*𝒑₀ᵀ∂²𝗠⁻¹∂x∂y𝒑*∂𝜙∂x + 𝒑₀ᵀ∂²𝗠⁻¹∂x²𝒑*∂𝜙∂y + 2*𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*∂²𝜙∂x∂y + 𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*∂²𝜙∂x² + 2*𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x∂y*∂𝜙∂x + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x²*∂𝜙∂y + 2*𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*∂²𝜙∂x∂y + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*∂²𝜙∂x² + 2*𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂x*∂𝜙∂x + 2*𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂y*∂𝜙∂x + 2*𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂x*∂𝜙∂y

                ∂³𝝭∂x∂y²[i] = 𝒑₀ᵀ∂³𝗠⁻¹∂x∂y²𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂³𝒑∂x∂y²*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂³𝜙∂x∂y² + 2*𝒑₀ᵀ∂²𝗠⁻¹∂x∂y∂𝒑∂y*𝜙 + 𝒑₀ᵀ∂²𝗠⁻¹∂y²∂𝒑∂x*𝜙 + 2*𝒑₀ᵀ∂𝗠⁻¹∂y∂²𝒑∂x∂y*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂x∂²𝒑∂y²*𝜙 + 2*𝒑₀ᵀ∂²𝗠⁻¹∂x∂y𝒑*∂𝜙∂y + 𝒑₀ᵀ∂²𝗠⁻¹∂y²𝒑*∂𝜙∂x + 2*𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*∂²𝜙∂x∂y + 𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*∂²𝜙∂y² + 2*𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x∂y*∂𝜙∂y + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂y²*∂𝜙∂x + 2*𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*∂²𝜙∂x∂y + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*∂²𝜙∂y² + 2*𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂y*∂𝜙∂y + 2*𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂x*∂𝜙∂y + 2*𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂y*∂𝜙∂x

                ∂³𝝭∂y³[i] = 𝒑₀ᵀ∂³𝗠⁻¹∂y³𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂³𝒑∂y³*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂³𝜙∂y³ + 3*𝒑₀ᵀ∂²𝗠⁻¹∂y²∂𝒑∂y*𝜙  + 3*𝒑₀ᵀ∂𝗠⁻¹∂y∂²𝒑∂y²*𝜙 + 3*𝒑₀ᵀ∂²𝗠⁻¹∂y²𝒑*∂𝜙∂y + 3*𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*∂²𝜙∂y² + 3*𝒑₀ᵀ𝗠⁻¹∂²𝒑∂y²*∂𝜙∂y + 3*𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*∂²𝜙∂y² + 6*𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂y*∂𝜙∂y
            end
        end
    end
end

# for 𝒑 in (:(:Linear3D),:(:Quadratic3D),:(:Cubic3D))
    for 𝒑 in (:(:Linear3D),:(:Quadratic3D))
    @eval begin
        get∇𝒑(ap::ReproducingKernel{$𝒑},x::Any) = get𝒑(ap,x),get∂𝒑∂x(ap,x),get∂𝒑∂y(ap,x),get∂𝒑∂z(ap,x)

        function cal∇𝗠!(ap::ReproducingKernel{$𝒑},x::Node)
            𝓒 = ap.𝓒
            𝗠 = get𝗠(ap,:𝗠)
            ∂𝗠∂x = get𝗠(ap,:∂𝗠∂x)
            ∂𝗠∂y = get𝗠(ap,:∂𝗠∂y)
            ∂𝗠∂z = get𝗠(ap,:∂𝗠∂z)
            n = get𝑛𝒑(ap)
            for xᵢ in 𝓒
                Δx = x - xᵢ
                𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂𝒑∂z = get∇𝒑(ap,Δx)
                𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂𝜙∂z = get∇𝜙(ap,xᵢ,Δx)
                for I in 1:n
                    for J in 1:I
                        𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
                        ∂𝗠∂x[I,J] += ∂𝜙∂x*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂x[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂x[J]
                        ∂𝗠∂y[I,J] += ∂𝜙∂y*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂y[J]
                        ∂𝗠∂z[I,J] += ∂𝜙∂z*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂z[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂z[J]
                    end
                end
            end
            cholesky!(𝗠)
            U = inverse!(𝗠)
            ∂𝗠⁻¹∂x = - UUᵀAUUᵀ!(∂𝗠∂x,U)
            ∂𝗠⁻¹∂y = - UUᵀAUUᵀ!(∂𝗠∂y,U)
            ∂𝗠⁻¹∂z = - UUᵀAUUᵀ!(∂𝗠∂z,U)
            𝗠⁻¹ = UUᵀ!(U)
            return 𝗠⁻¹, ∂𝗠⁻¹∂x, ∂𝗠⁻¹∂y, ∂𝗠⁻¹∂z
        end

        function set∇𝝭!(ap::ReproducingKernel{$𝒑},𝒙::Node)
            𝓒 = ap.𝓒
            𝝭 = 𝒙[:𝝭]
            ∂𝝭∂x = 𝒙[:∂𝝭∂x]
            ∂𝝭∂y = 𝒙[:∂𝝭∂y]
            ∂𝝭∂z = 𝒙[:∂𝝭∂z]
            𝒑₀ᵀ𝗠⁻¹, 𝒑₀ᵀ∂𝗠⁻¹∂x, 𝒑₀ᵀ∂𝗠⁻¹∂y, 𝒑₀ᵀ∂𝗠⁻¹∂z= cal∇𝗠!(ap,𝒙)
            for (i,𝒙ᵢ) in enumerate(𝓒)
                Δ𝒙 = 𝒙 - 𝒙ᵢ
                𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂𝒑∂z = get∇𝒑(ap,Δ𝒙)
                𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂𝜙∂z = get∇𝜙(ap,𝒙ᵢ,Δ𝒙)
                𝝭[i] = 𝒑₀ᵀ𝗠⁻¹*𝒑*𝜙
                ∂𝝭∂x[i] = 𝒑₀ᵀ∂𝗠⁻¹∂x*𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂x*𝜙 + 𝒑₀ᵀ𝗠⁻¹*𝒑*∂𝜙∂x
                ∂𝝭∂y[i] = 𝒑₀ᵀ∂𝗠⁻¹∂y*𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂y*𝜙 + 𝒑₀ᵀ𝗠⁻¹*𝒑*∂𝜙∂y
                ∂𝝭∂z[i] = 𝒑₀ᵀ∂𝗠⁻¹∂z*𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂z*𝜙 + 𝒑₀ᵀ𝗠⁻¹*𝒑*∂𝜙∂z
            end
        end

        function cal∇²𝗠!(ap::ReproducingKernel{$𝒑},x::Node)
            𝓒 = ap.𝓒
            𝗠 = get𝗠(ap,:𝗠)
            ∂𝗠∂x = get𝗠(ap,:∂𝗠∂x)
            ∂𝗠∂y = get𝗠(ap,:∂𝗠∂y)
            ∂𝗠∂z = get𝗠(ap,:∂𝗠∂z)
            ∂²𝗠∂x² = get𝗠(ap,:∂²𝗠∂x²)
            ∂²𝗠∂y² = get𝗠(ap,:∂²𝗠∂y²)
            ∂²𝗠∂z² = get𝗠(ap,:∂²𝗠∂z²)
            ∂²𝗠∂x∂y = get𝗠(ap,:∂²𝗠∂x∂y)
            ∂²𝗠∂x∂z = get𝗠(ap,:∂²𝗠∂x∂z)
            ∂²𝗠∂y∂z = get𝗠(ap,:∂²𝗠∂y∂z)
            n = get𝑛𝒑(ap)
            for xᵢ in 𝓒
                Δx = x - xᵢ
                𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂²𝒑∂x², ∂²𝒑∂x∂y, ∂²𝒑∂y², ∂𝒑∂z, ∂²𝒑∂x∂z, ∂²𝒑∂y∂z, ∂²𝒑∂z² = get∇²𝒑(ap,Δx)
                𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂²𝜙∂x², ∂²𝜙∂x∂y, ∂²𝜙∂y², ∂𝜙∂z, ∂²𝜙∂x∂z, ∂²𝜙∂y∂z, ∂²𝜙∂z² = get∇²𝜙(ap,xᵢ,Δx)
                for I in 1:n
                    for J in 1:I
                        𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
                        ∂𝗠∂x[I,J] += ∂𝜙∂x*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂x[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂x[J]
                        ∂𝗠∂y[I,J] += ∂𝜙∂y*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂y[J]
                        ∂𝗠∂z[I,J] += ∂𝜙∂z*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂z[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂z[J]

                        ∂²𝗠∂x²[I,J] += ∂²𝜙∂x²*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂x²[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂x²[J] + 2.0*∂𝜙∂x*∂𝒑∂x[I]*𝒑[J] + 2.0*∂𝜙∂x*𝒑[I]*∂𝒑∂x[J] + 2.0*𝜙*∂𝒑∂x[I]*∂𝒑∂x[J]

                        ∂²𝗠∂y²[I,J] += ∂²𝜙∂y²*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂y²[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂y²[J] + 2.0*∂𝜙∂y*∂𝒑∂y[I]*𝒑[J] + 2.0*∂𝜙∂y*𝒑[I]*∂𝒑∂y[J] + 2.0*𝜙*∂𝒑∂y[I]*∂𝒑∂y[J]

                        ∂²𝗠∂z²[I,J] += ∂²𝜙∂z²*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂z²[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂z²[J] + 2.0*∂𝜙∂z*∂𝒑∂z[I]*𝒑[J] + 2.0*∂𝜙∂z*𝒑[I]*∂𝒑∂z[J] + 2.0*𝜙*∂𝒑∂z[I]*∂𝒑∂z[J]

                        ∂²𝗠∂x∂y[I,J] += ∂²𝜙∂x∂y*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂x∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂x∂y[J] + ∂𝜙∂x*∂𝒑∂y[I]*𝒑[J] + ∂𝜙∂y*∂𝒑∂x[I]*𝒑[J] + ∂𝜙∂x*𝒑[I]*∂𝒑∂y[J] + ∂𝜙∂y*𝒑[I]*∂𝒑∂x[J] + 𝜙*∂𝒑∂x[I]*∂𝒑∂y[J] + 𝜙*∂𝒑∂y[I]*∂𝒑∂x[J]

                        ∂²𝗠∂x∂z[I,J] += ∂²𝜙∂x∂z*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂x∂z[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂x∂z[J] + ∂𝜙∂x*∂𝒑∂z[I]*𝒑[J] + ∂𝜙∂z*∂𝒑∂x[I]*𝒑[J] + ∂𝜙∂x*𝒑[I]*∂𝒑∂z[J] + ∂𝜙∂z*𝒑[I]*∂𝒑∂x[J] + 𝜙*∂𝒑∂x[I]*∂𝒑∂z[J] + 𝜙*∂𝒑∂z[I]*∂𝒑∂x[J]

                        ∂²𝗠∂y∂z[I,J] += ∂²𝜙∂y∂z*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂y∂z[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂y∂z[J] + ∂𝜙∂y*∂𝒑∂z[I]*𝒑[J] + ∂𝜙∂z*∂𝒑∂y[I]*𝒑[J] + ∂𝜙∂y*𝒑[I]*∂𝒑∂z[J] + ∂𝜙∂z*𝒑[I]*∂𝒑∂y[J] + 𝜙*∂𝒑∂y[I]*∂𝒑∂z[J] + 𝜙*∂𝒑∂z[I]*∂𝒑∂y[J]
                    end
                end
            end
            cholesky!(𝗠)
            U = inverse!(𝗠)
            Uᵀ∂𝗠∂xU = UᵀAU!(∂𝗠∂x,U)
            Uᵀ∂𝗠∂yU = UᵀAU!(∂𝗠∂y,U)
            Uᵀ∂𝗠∂zU = UᵀAU!(∂𝗠∂z,U)
            Uᵀ∂²𝗠∂x²U = UᵀAU!(∂²𝗠∂x²,U)
            Uᵀ∂²𝗠∂y²U = UᵀAU!(∂²𝗠∂y²,U)
            Uᵀ∂²𝗠∂z²U = UᵀAU!(∂²𝗠∂z²,U)
            Uᵀ∂²𝗠∂x∂yU = UᵀAU!(∂²𝗠∂x∂y,U)
            Uᵀ∂²𝗠∂x∂zU = UᵀAU!(∂²𝗠∂x∂z,U)
            Uᵀ∂²𝗠∂y∂zU = UᵀAU!(∂²𝗠∂y∂z,U)
            for i in 1:n
                for j in 1:i
                    for k in 1:n
                        Uᵀ∂²𝗠∂x²U[i,j] -= 2*Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂𝗠∂xU[k,j]
                        Uᵀ∂²𝗠∂y²U[i,j] -= 2*Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂𝗠∂yU[k,j]
                        Uᵀ∂²𝗠∂z²U[i,j] -= 2*Uᵀ∂𝗠∂zU[i,k]*Uᵀ∂𝗠∂zU[k,j]
                        Uᵀ∂²𝗠∂x∂yU[i,j] -= Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂𝗠∂yU[k,j] + Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂𝗠∂xU[k,j]
                        Uᵀ∂²𝗠∂x∂zU[i,j] -= Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂𝗠∂zU[k,j] + Uᵀ∂𝗠∂zU[i,k]*Uᵀ∂𝗠∂xU[k,j]
                        Uᵀ∂²𝗠∂y∂zU[i,j] -= Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂𝗠∂zU[k,j] + Uᵀ∂𝗠∂zU[i,k]*Uᵀ∂𝗠∂yU[k,j]
                    end
                end
            end

            ∂²𝗠⁻¹∂x² = - UAUᵀ!(Uᵀ∂²𝗠∂x²U,U)
            ∂²𝗠⁻¹∂y² = - UAUᵀ!(Uᵀ∂²𝗠∂y²U,U)
            ∂²𝗠⁻¹∂z² = - UAUᵀ!(Uᵀ∂²𝗠∂z²U,U)
            ∂²𝗠⁻¹∂x∂y = - UAUᵀ!(Uᵀ∂²𝗠∂x∂yU,U)
            ∂²𝗠⁻¹∂x∂z = - UAUᵀ!(Uᵀ∂²𝗠∂x∂zU,U)
            ∂²𝗠⁻¹∂y∂z = - UAUᵀ!(Uᵀ∂²𝗠∂y∂zU,U)
            ∂𝗠⁻¹∂x = - UAUᵀ!(Uᵀ∂𝗠∂xU,U)
            ∂𝗠⁻¹∂y = - UAUᵀ!(Uᵀ∂𝗠∂yU,U)
            ∂𝗠⁻¹∂z = - UAUᵀ!(Uᵀ∂𝗠∂zU,U)
            𝗠⁻¹ = UUᵀ!(U)
            return 𝗠⁻¹, ∂𝗠⁻¹∂x, ∂𝗠⁻¹∂y, ∂²𝗠⁻¹∂x², ∂²𝗠⁻¹∂x∂y, ∂²𝗠⁻¹∂y², ∂𝗠⁻¹∂z, ∂²𝗠⁻¹∂x∂z, ∂²𝗠⁻¹∂y∂z, ∂²𝗠⁻¹∂z²
        end

        function set∇²𝝭!(ap::ReproducingKernel{$𝒑},𝒙::Node)
            𝓒 = ap.𝓒
            𝝭 = 𝒙[:𝝭]
            ∂𝝭∂x = 𝒙[:∂𝝭∂x]
            ∂𝝭∂y = 𝒙[:∂𝝭∂y]
            ∂𝝭∂z = 𝒙[:∂𝝭∂z]
            ∂²𝝭∂x² = 𝒙[:∂²𝝭∂x²]
            ∂²𝝭∂y² = 𝒙[:∂²𝝭∂y²]
            ∂²𝝭∂z² = 𝒙[:∂²𝝭∂z²]
            ∂²𝝭∂x∂y = 𝒙[:∂²𝝭∂x∂y]
            ∂²𝝭∂x∂z = 𝒙[:∂²𝝭∂x∂z]
            ∂²𝝭∂y∂z = 𝒙[:∂²𝝭∂y∂z]
            𝒑₀ᵀ𝗠⁻¹, 𝒑₀ᵀ∂𝗠⁻¹∂x, 𝒑₀ᵀ∂𝗠⁻¹∂y, 𝒑₀ᵀ∂²𝗠⁻¹∂x², 𝒑₀ᵀ∂²𝗠⁻¹∂x∂y, 𝒑₀ᵀ∂²𝗠⁻¹∂y², 𝒑₀ᵀ∂𝗠⁻¹∂z, 𝒑₀ᵀ∂²𝗠⁻¹∂x∂z, 𝒑₀ᵀ∂²𝗠⁻¹∂y∂z, 𝒑₀ᵀ∂²𝗠⁻¹∂z² = cal∇²𝗠!(ap,𝒙)
            for (i,𝒙ᵢ) in enumerate(𝓒)
                Δ𝒙 = 𝒙 - 𝒙ᵢ
                𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂²𝒑∂x², ∂²𝒑∂x∂y, ∂²𝒑∂y², ∂𝒑∂z, ∂²𝒑∂x∂z, ∂²𝒑∂y∂z, ∂²𝒑∂z² = get∇²𝒑(ap,Δ𝒙)
                𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂²𝜙∂x², ∂²𝜙∂x∂y, ∂²𝜙∂y², ∂𝜙∂z, ∂²𝜙∂x∂z, ∂²𝜙∂y∂z, ∂²𝜙∂z² = get∇²𝜙(ap,𝒙ᵢ,Δ𝒙)
                𝒑₀ᵀ𝗠⁻¹𝒑 = 𝒑₀ᵀ𝗠⁻¹*𝒑
                𝒑₀ᵀ∂𝗠⁻¹∂x𝒑 = 𝒑₀ᵀ∂𝗠⁻¹∂x*𝒑
                𝒑₀ᵀ∂𝗠⁻¹∂y𝒑 = 𝒑₀ᵀ∂𝗠⁻¹∂y*𝒑
                𝒑₀ᵀ∂𝗠⁻¹∂z𝒑 = 𝒑₀ᵀ∂𝗠⁻¹∂z*𝒑
                𝒑₀ᵀ𝗠⁻¹∂𝒑∂x = 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂x
                𝒑₀ᵀ𝗠⁻¹∂𝒑∂y = 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂y
                𝒑₀ᵀ𝗠⁻¹∂𝒑∂z = 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂z
                𝒑₀ᵀ∂²𝗠⁻¹∂x²𝒑 = 𝒑₀ᵀ∂²𝗠⁻¹∂x²*𝒑
                𝒑₀ᵀ∂²𝗠⁻¹∂y²𝒑 = 𝒑₀ᵀ∂²𝗠⁻¹∂y²*𝒑
                𝒑₀ᵀ∂²𝗠⁻¹∂z²𝒑 = 𝒑₀ᵀ∂²𝗠⁻¹∂z²*𝒑
                𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x² = 𝒑₀ᵀ𝗠⁻¹*∂²𝒑∂x²
                𝒑₀ᵀ𝗠⁻¹∂²𝒑∂y² = 𝒑₀ᵀ𝗠⁻¹*∂²𝒑∂y²
                𝒑₀ᵀ𝗠⁻¹∂²𝒑∂z² = 𝒑₀ᵀ𝗠⁻¹*∂²𝒑∂z²
                𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂x = 𝒑₀ᵀ∂𝗠⁻¹∂x*∂𝒑∂x
                𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂y = 𝒑₀ᵀ∂𝗠⁻¹∂y*∂𝒑∂y
                𝒑₀ᵀ∂𝗠⁻¹∂z∂𝒑∂z = 𝒑₀ᵀ∂𝗠⁻¹∂z*∂𝒑∂z
                𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂y = 𝒑₀ᵀ∂𝗠⁻¹∂x*∂𝒑∂y
                𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂x = 𝒑₀ᵀ∂𝗠⁻¹∂y*∂𝒑∂x
                𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂z = 𝒑₀ᵀ∂𝗠⁻¹∂x*∂𝒑∂z
                𝒑₀ᵀ∂𝗠⁻¹∂z∂𝒑∂x = 𝒑₀ᵀ∂𝗠⁻¹∂z*∂𝒑∂x
                𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂z = 𝒑₀ᵀ∂𝗠⁻¹∂y*∂𝒑∂z
                𝒑₀ᵀ∂𝗠⁻¹∂z∂𝒑∂y = 𝒑₀ᵀ∂𝗠⁻¹∂z*∂𝒑∂y
                𝒑₀ᵀ∂²𝗠⁻¹∂x∂y𝒑 = 𝒑₀ᵀ∂²𝗠⁻¹∂x∂y*𝒑
                𝒑₀ᵀ∂²𝗠⁻¹∂x∂z𝒑 = 𝒑₀ᵀ∂²𝗠⁻¹∂x∂z*𝒑
                𝒑₀ᵀ∂²𝗠⁻¹∂y∂z𝒑 = 𝒑₀ᵀ∂²𝗠⁻¹∂y∂z*𝒑
                𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x∂y = 𝒑₀ᵀ𝗠⁻¹*∂²𝒑∂x∂y
                𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x∂z = 𝒑₀ᵀ𝗠⁻¹*∂²𝒑∂x∂z
                𝒑₀ᵀ𝗠⁻¹∂²𝒑∂y∂z = 𝒑₀ᵀ𝗠⁻¹*∂²𝒑∂y∂z

                𝝭[i] = 𝒑₀ᵀ𝗠⁻¹𝒑*𝜙
                ∂𝝭∂x[i] = 𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂𝜙∂x
                ∂𝝭∂y[i] = 𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂𝜙∂y
                ∂𝝭∂z[i] = 𝒑₀ᵀ∂𝗠⁻¹∂z𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂z*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂𝜙∂z

                ∂²𝝭∂x²[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂x²𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x²*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂x² + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂x*𝜙 + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*∂𝜙∂x + 2.0*𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*∂𝜙∂x

                ∂²𝝭∂y²[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂y²𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂y²*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂y² + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂y*𝜙 + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*∂𝜙∂y + 2.0*𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*∂𝜙∂y

                ∂²𝝭∂z²[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂z²𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂z²*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂z² + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂z∂𝒑∂z*𝜙 + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂z𝒑*∂𝜙∂z + 2.0*𝒑₀ᵀ𝗠⁻¹∂𝒑∂z*∂𝜙∂z

                ∂²𝝭∂x∂y[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂x∂y𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x∂y*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂x∂y + 𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂y*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂x*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*∂𝜙∂y + 𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*∂𝜙∂x + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*∂𝜙∂y +𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*∂𝜙∂x

                ∂²𝝭∂x∂z[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂x∂z𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x∂z*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂x∂z + 𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂z*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂z∂𝒑∂x*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*∂𝜙∂z + 𝒑₀ᵀ∂𝗠⁻¹∂z𝒑*∂𝜙∂x + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*∂𝜙∂z +𝒑₀ᵀ𝗠⁻¹∂𝒑∂z*∂𝜙∂x

                ∂²𝝭∂y∂z[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂y∂z𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂y∂z*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂y∂z + 𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂z*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂z∂𝒑∂y*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*∂𝜙∂z + 𝒑₀ᵀ∂𝗠⁻¹∂z𝒑*∂𝜙∂y + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*∂𝜙∂z +𝒑₀ᵀ𝗠⁻¹∂𝒑∂z*∂𝜙∂y
            end
        end
    end
end

function set∇̂³𝝭!(ap::ReproducingKernel,𝒙::Node)
    𝓒 = ap.𝓒
    𝝭 = 𝒙[:𝝭]
    ∂𝝭∂x = 𝒙[:∂𝝭∂x]
    ∂𝝭∂y = 𝒙[:∂𝝭∂y]
    ∂²𝝭∂x² = 𝒙[:∂²𝝭∂x²]
    ∂²𝝭∂y² = 𝒙[:∂²𝝭∂y²]
    ∂²𝝭∂x∂y = 𝒙[:∂²𝝭∂x∂y]
    ∂³𝝭∂x³ = 𝒙[:∂³𝝭∂x³]
    ∂³𝝭∂x²∂y = 𝒙[:∂³𝝭∂x²∂y]
    ∂³𝝭∂x∂y² = 𝒙[:∂³𝝭∂x∂y²]
    ∂³𝝭∂y³ = 𝒙[:∂³𝝭∂y³]

    n = get𝑛𝒑(ap)
    𝒑₀ᵀ𝗠⁻¹ = cal𝗠!(ap,𝒙)
    for (i,𝒙ᵢ) in enumerate(𝓒)
        Δ𝒙 = 𝒙 - 𝒙ᵢ
        𝒑= get𝒑(ap,Δ𝒙)
        𝜙 = get𝜙(ap,𝒙ᵢ,Δ𝒙)
        ∂𝝭∂x_ = 0.0
        ∂𝝭∂y_ = 0.0
        ∂²𝝭∂x²_ = 0.0
        ∂²𝝭∂x∂y_ = 0.0
        ∂²𝝭∂y²_ = 0.0
        ∂³𝝭∂x³_ = 0.0
        ∂³𝝭∂x²∂y_ = 0.0
        ∂³𝝭∂x∂y²_ = 0.0
        ∂³𝝭∂y³_ = 0.0
        for j in 1:n
            ∂𝝭∂x_ -= 𝒑₀ᵀ𝗠⁻¹[2,j]*𝒑[j]*𝜙
            ∂𝝭∂y_ -= 𝒑₀ᵀ𝗠⁻¹[3,j]*𝒑[j]*𝜙
            ∂²𝝭∂x²_ += 2*𝒑₀ᵀ𝗠⁻¹[4,j]*𝒑[j]*𝜙
            ∂²𝝭∂x∂y_ += 𝒑₀ᵀ𝗠⁻¹[5,j]*𝒑[j]*𝜙
            ∂²𝝭∂y²_ += 2*𝒑₀ᵀ𝗠⁻¹[6,j]*𝒑[j]*𝜙
            ∂³𝝭∂x³_ -= 6*𝒑₀ᵀ𝗠⁻¹[7,j]*𝒑[j]*𝜙
            ∂³𝝭∂x²∂y_ -= 2*𝒑₀ᵀ𝗠⁻¹[8,j]*𝒑[j]*𝜙
            ∂³𝝭∂x∂y²_ -= 2*𝒑₀ᵀ𝗠⁻¹[9,j]*𝒑[j]*𝜙
            ∂³𝝭∂y³_ -= 6*𝒑₀ᵀ𝗠⁻¹[10,j]*𝒑[j]*𝜙
        end
        𝝭[i] = 𝒑₀ᵀ𝗠⁻¹*𝒑*𝜙
        ∂𝝭∂x[i] = ∂𝝭∂x_
        ∂𝝭∂y[i] = ∂𝝭∂y_
        ∂²𝝭∂x²[i] = ∂²𝝭∂x²_
        ∂²𝝭∂x∂y[i] = ∂²𝝭∂x∂y_
        ∂²𝝭∂y²[i] = ∂²𝝭∂y²_
        ∂³𝝭∂x³[i] = ∂³𝝭∂x³_
        ∂³𝝭∂x²∂y[i] = ∂³𝝭∂x²∂y_
        ∂³𝝭∂x∂y²[i] = ∂³𝝭∂x∂y²_
        ∂³𝝭∂y³[i] = ∂³𝝭∂y³_
    end
end