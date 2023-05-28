
get𝑛𝒑(  ::ReproducingKernel{:Wave2D}) = 5
get𝒑(   ::ReproducingKernel{:Wave2D},x::NTuple{3,Float64}) = (1.,sin(x[1]),cos(x[1]),sin(x[2]),cos(x[2]))
get∂𝒑∂x(::ReproducingKernel{:Wave2D},x::NTuple{3,Float64}) = (0.,cos(x[1]),-sin(x[1]),0.,0.)
get∂𝒑∂y(::ReproducingKernel{:Wave2D},x::NTuple{3,Float64}) = (0.,0.,0.,cos(x[2]),-sin(x[2]))

function set𝝭!(ap::ReproducingKernel{Wave2D},𝒙::Node)
    𝓒 = ap.𝓒
    𝝭 = 𝒙[:𝝭]
    𝗠⁻¹ = cal𝗠!(ap,𝒙)
    𝒑₀ = get𝒑(ap,(0.0,0.0,0.0))
    n = get𝑛𝒑(ap)
    for (i,𝒙ᵢ) in enumerate(𝓒)
        Δ𝒙 = 𝒙 - 𝒙ᵢ
        𝒑= get𝒑(ap,Δ𝒙)
        𝜙 = get𝜙(ap,𝒙ᵢ,Δ𝒙)
        for j in 1:n
            for k in 1:n
                𝝭[i] = 𝒑₀[j]*𝗠⁻¹[j,k]*𝒑[k]*𝜙
            end
        end
    end
end

function set∇𝝭!(ap::ReproducingKernel{Wave2D},𝒙::Node)
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
        𝜙, ∂𝜙∂x, ∂𝜙∂x = get∇𝜙(ap,𝒙ᵢ,Δ𝒙)
        for j in 1:n
            for k in 1:n
                𝝭[i] = 𝒑₀[j]*𝗠⁻¹[j,k]*𝒑[k]*𝜙
                ∂𝝭∂x[i] = 𝒑₀[j]*∂𝗠⁻¹∂x[j,k]*𝒑[k]*𝜙 + 𝒑₀[j]𝗠⁻¹[j,k]*∂𝒑∂x[k]*𝜙 + 𝒑₀[j]𝗠⁻¹[j,k]*𝒑[k]*∂𝜙∂x
                ∂𝝭∂y[i] = 𝒑₀[j]*∂𝗠⁻¹∂y[j,k]*𝒑[k]*𝜙 + 𝒑₀[j]𝗠⁻¹[j,k]*∂𝒑∂y[k]*𝜙 + 𝒑₀[j]𝗠⁻¹[j,k]*𝒑[k]*∂𝜙∂y
            end
        end
    end
end