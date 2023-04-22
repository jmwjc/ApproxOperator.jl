
abstract type AbstractReproducingKernel{𝑠,𝜙,T}<:AbstractElement{T} end



for set𝝭 in (:set𝝭!,:set∇𝝭!,:set∇²𝝭!,:set∇³𝝭!)
    @eval begin
        function $set𝝭(aps::Vector{T}) where T<:AbstractReproducingKernel
            for ap in aps
                𝓖 = ap.𝓖
                for 𝒙 in 𝓖
                    $set𝝭(ap,𝒙)
                end
            end
        end
    end
end

for set𝝭 in (:set∇̃𝝭!,:set∇̃²𝝭!,:set∇∇̃²𝝭!)
    @eval begin
        function $set𝝭(gps::Vector{T},aps::Vector{S}) where {T<:ReproducingKernel,S<:ReproducingKernel}
            if length(gps) ≠ length(aps)
                error("Miss match element numbers")
            else
                for i in 1:length(gps)
                    $set𝝭(gps[i],aps[i])
                end
            end
        end
    end
end

function set∇̃𝝭!(cps::Vector{T},gps::Vector{T},aps::Vector{T}) where T<:ReproducingKernel
    if length(gps) ≠ length(aps)
        error("Miss match element numbers")
    else
        for i in 1:length(cps)
            set∇̃𝝭!(cps[i],gps[i],aps[i])
        end
    end
end

for set𝝭 in (:set∇̄𝝭!,:set∇̃₁𝝭!)
    @eval begin
        function $set𝝭(aps::Vector{T}) where T<:ReproducingKernel
            for ap in aps
                $set𝝭(ap)
            end
        end
    end
end

function set∇̄²𝝭!(aps::Vector{T};Γᵍ::Vector{T}=T[],Γᶿ::Vector{T}=T[],Γᴾ::Vector{T}=T[]) where T<:ReproducingKernel
    for ap in aps
        set∇̄²𝝭!(ap,Γᵍ=Γᵍ,Γᶿ=Γᶿ,Γᴾ=Γᴾ)
    end
end

function set∇∇̄²𝝭!(aps::Vector{T};Γᵍ::Vector{T}=T[],Γᶿ::Vector{T}=T[],Γᴾ::Vector{T}=T[]) where T<:ReproducingKernel
    for i in 1:length(aps)
        isempty(Γᵍ) ? a = nothing : a = Γᵍ[i]
        isempty(Γᶿ) ? b = nothing : b = Γᶿ[i]
        set∇∇̄²𝝭!(aps[i],Γᵍ=a,Γᶿ=b,Γᴾ=Γᴾ)
    end
end

"""
SymMat
"""
struct SymMat
    n::Int
    m::Vector{Float64}
end
SymMat(n::Int) = SymMat(n,zeros(Int(n*(n+1)/2)))

@inline function getindex(A::SymMat,i::Int,j::Int)
    i > j ? A.m[Int(j+i*(i-1)/2)] : A.m[Int(i+j*(j-1)/2)]
end

@inline function setindex!(A::SymMat,val::Float64,i::Int,j::Int)
    i > j ? A.m[Int(j+i*(i-1)/2)] = val : A.m[Int(i+j*(j-1)/2)] = val
end
@inline function setindex!(A::SymMat,val::Float64,i::Int)
    A.m[i] = val
end

@inline *(A::SymMat,v::NTuple{N,Float64}) where N = sum(A[1,i]*v[i] for i in 1:N)
@inline function *(v::NTuple{N,Float64},A::SymMat) where N
    for j in 1:N
        A[1,j] = sum(v[i]*A[i,j] for i in 1:N)
    end
    return A
end
@inline function -(A::SymMat)
    A.m .= .-A.m
    return A
end
@inline fill!(A::SymMat,val::Float64) = fill!(A.m,val)
function inverse!(A::SymMat)
    n = A.n
    for i in 1:n
        A[i,i] = 1.0/A[i,i]
        for j in i+1:n
            A[i,j] = - sum(A[i,k]*A[k,j] for k in i:j-1)/A[j,j]
        end
    end
    return A
end

function UUᵀ!(A::SymMat)
    n = A.n
    for i in 1:n
        for j in 1:i
            A[i,j] = sum(A[i,k]*A[k,j] for k in i:n)
        end
    end
    return A
end

function UᵀAU!(A::SymMat,U::SymMat)
    n = A.n
    for i in n:-1:1
        for j in n:-1:i
            A[i,j] = sum(U[k,i]*A[k,l]*U[l,j] for k in 1:i for l in 1:j)
        end
    end
    return A
end

function UAUᵀ!(A::SymMat,U::SymMat)
    n = A.n
    for i in 1:n
        for j in i:n
            A[i,j] = sum(U[i,k]*A[k,l]*U[j,l] for k in i:n for l in j:n)
        end
    end
    return A
end

function UUᵀAUUᵀ!(A::SymMat,U::SymMat)
    UᵀAU!(A,U)
    UAUᵀ!(A,U)
    return A
end

function cholesky!(A::SymMat)
    n = A.n
    for i in 1:n
        for k in 1:i-1
            A[i,i] -= A[k,i]^2
        end
        A[i,i] = A[i,i]^0.5
        for j in i+1:n
            for k in 1:i-1
                A[i,j] -= A[k,i]A[k,j]
            end
            A[i,j] = A[i,j]/A[i,i]
        end
    end
    return A
end