
get𝑛𝒑(  ::ReproducingKernel{:Wave2D}) = 5
get𝒑(   ::ReproducingKernel{:Wave2D},x::NTuple{3,Float64}) = (1.,sin(x[1]),cos(x[1]),sin(x[2]),cos(x[2]))
# get𝒑(   ::ReproducingKernel{:Wave2D},x::NTuple{3,Float64}) = (1.,x[1],x[2],x[1]^2,x[2]^2)
get∂𝒑∂x(::ReproducingKernel{:Wave2D},x::NTuple{3,Float64}) = (0.,cos(x[1]),-sin(x[1]),0.,0.)
get∂𝒑∂y(::ReproducingKernel{:Wave2D},x::NTuple{3,Float64}) = (0.,0.,0.,cos(x[2]),-sin(x[2]))

function set𝝭!(ap::ReproducingKernel{:Wave2D},𝒙::Node)
    𝓒 = ap.𝓒
    𝝭 = 𝒙[:𝝭]
    𝗠⁻¹ = cal𝗠!(ap,𝒙)
    𝒑₀ = get𝒑(ap,(0.0,0.0,0.0))
    n = get𝑛𝒑(ap)
    for (i,𝒙ᵢ) in enumerate(𝓒)
        Δ𝒙 = 𝒙 - 𝒙ᵢ
        𝒑 = get𝒑(ap,Δ𝒙)
        𝜙 = get𝜙(ap,𝒙ᵢ,Δ𝒙)
        for j in 1:n
            for k in 1:n
                𝝭[i] += 𝒑₀[j]*𝗠⁻¹[j,k]*𝒑[k]*𝜙
            end
        end
    end
end
get∇𝒑( ap::ReproducingKernel{:Wave2D},x::Any) = get𝒑(ap,x),get∂𝒑∂x(ap,x),get∂𝒑∂y(ap,x)
function cal∇𝗠!(ap::ReproducingKernel{:Wave2D},x::Node)
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

function set∇𝝭!(ap::ReproducingKernel{:Wave2D},𝒙::Node)
    𝓒 = ap.𝓒
    𝝭 = 𝒙[:𝝭]
    ∂𝝭∂x = 𝒙[:∂𝝭∂x]
    ∂𝝭∂y = 𝒙[:∂𝝭∂y]
    𝗠⁻¹, ∂𝗠⁻¹∂x, ∂𝗠⁻¹∂y = cal∇𝗠!(ap,𝒙)
    𝒑₀ = get𝒑(ap,(0.0,0.0,0.0))
    n = get𝑛𝒑(ap)
    for (i,𝒙ᵢ) in enumerate(𝓒)
        Δ𝒙 = 𝒙 - 𝒙ᵢ
        𝒑, ∂𝒑∂x, ∂𝒑∂y = get∇𝒑(ap,Δ𝒙)
        𝜙, ∂𝜙∂x, ∂𝜙∂y = get∇𝜙(ap,𝒙ᵢ,Δ𝒙)
        for j in 1:n
            for k in 1:n
                𝝭[i] += 𝒑₀[j]*𝗠⁻¹[j,k]*𝒑[k]*𝜙
                ∂𝝭∂x[i] += 𝒑₀[j]*∂𝗠⁻¹∂x[j,k]*𝒑[k]*𝜙 + 𝒑₀[j]𝗠⁻¹[j,k]*∂𝒑∂x[k]*𝜙 + 𝒑₀[j]𝗠⁻¹[j,k]*𝒑[k]*∂𝜙∂x
                ∂𝝭∂y[i] += 𝒑₀[j]*∂𝗠⁻¹∂y[j,k]*𝒑[k]*𝜙 + 𝒑₀[j]𝗠⁻¹[j,k]*∂𝒑∂y[k]*𝜙 + 𝒑₀[j]𝗠⁻¹[j,k]*𝒑[k]*∂𝜙∂y
            end
        end
    end
end