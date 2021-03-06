"""
Operator
"""
struct Operator{T,D}
    type::Val{T}
    data::Dict{Symbol,D}
end
Operator(t::Symbol) = Operator(Val(t),Dict{Symbol,Any}())
function Operator(t::Symbol,d::Pair{Symbol,D}...) where D<:Any
    return Operator(Val(t),Dict(d))
end


## General Functions
push!(op::Operator,d::Pair{Symbol,D}...) where D<:Any = push!(op.data,d...)
@inline getproperty(op::Operator,f::Symbol) = hasfield(Operator,f) ? getfield(op,f) : getfield(op,:data)[f]
@inline function setproperty!(op::Operator,f::Symbol,x)
    getfield(op,:data)[f] = x
end

@inline function (op::Operator)(aps::Vector{T},gps::Vector{S},k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    for i in 1:length(aps)
        @inbounds op(aps[i],gps[i],k,f)
    end
end

@inline function (op::Operator)(aps::Vector{T},k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    for ap in aps
        op(ap,k,f)
    end
end
@inline function (op::Operator)(aps::Vector{T},k::AbstractMatrix{Float64}) where T<:AbstractElement
    for ap in aps
        op(ap,k)
    end
end
@inline function (op::Operator)(aps::Vector{T},f::AbstractVector{Float64}) where T<:AbstractElement
    for ap in aps
        op(ap,f)
    end
end

@inline function (op::Operator)(aps::Vector{T},s::Symbol) where T<:AbstractElement
    for ap in aps
        op(ap,s)
    end
end
@inline function (op::Operator)(aps::Vector{T}) where T<:AbstractElement
    for ap in aps
        op(ap)
    end
end

function prescribe!(ap::T,sf::Pair{Symbol,F}) where {T<:AbstractElement,F<:Function}
    š = ap.š
    s,f = sf
    for Ī¾ in š
        š = (Ī¾.x,Ī¾.y,Ī¾.z)
        if applicable(f,š...)
            v = f(š...)
        elseif applicable(f,š...,Ī¾.nā)
            v = f(š...,Ī¾.nā)
        elseif applicable(f,š...,Ī¾.nā,Ī¾.nā)
            v = f(š...,Ī¾.nā,Ī¾.nā)
        elseif applicable(f,š...,Ī¾.nā,Ī¾.nā,Ī¾.nā)
            v = f(š...,Ī¾.nā,Ī¾.nā,Ī¾.nā)
        end
        setproperty!(Ī¾,s,v)
    end
end

function prescribe!(aps::Vector{T},sf::Pair{Symbol,F}) where {T<:AbstractElement,F<:Function}
    s,f = sf
    n = length(getfield(aps[1].š[1],:data)[:x][2])
    haskey(getfield(aps[1].š[1],:data),s) ? nothing : push!(getfield(aps[1].š[1],:data),s=>(2,zeros(n)))
    for ap in aps
        prescribe!(ap,sf)
    end
end

"""
# Potential Problem
"""
function (op::Operator{:šš£})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    for Ī¾ in š
        š¤ = Ī¾.š¤
        N = Ī¾.š­
        u = Ī¾.u
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            f[I] += N[i]*u*š¤
        end
    end
end

function (op::Operator{:ā«āvāuvbdĪ©})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    for Ī¾ in š
        N = Ī¾[:š­]
        Bā = Ī¾[:āš­āx]
        Bā = Ī¾[:āš­āy]
        Bā = Ī¾[:āš­āz]
        š¤ = Ī¾.š¤
        b = Ī¾.b
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            for (j,xā±¼) in enumerate(š)
                J = xā±¼.š¼
                k[I,J] += op.k*(Bā[i]*Bā[j] + Bā[i]*Bā[j] + Bā[i]*Bā[j])*š¤
            end
            f[I] += N[i]*b*š¤
        end
    end
end

function (op::Operator{:ā«āvāudĪ©})(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    for Ī¾ in š
        Bā = Ī¾[:āš­āx]
        Bā = Ī¾[:āš­āy]
        Bā = Ī¾[:āš­āz]
        š¤ = Ī¾.š¤
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            for (j,xā±¼) in enumerate(š)
                J = xā±¼.š¼
                k[I,J] += op.k*(Bā[i]*Bā[j] + Bā[i]*Bā[j] + Bā[i]*Bā[j])*š¤
            end
        end
    end
end

function (op::Operator{:ā«vbdĪ©})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    for Ī¾ in š
        š¤ = Ī¾.š¤
        N = Ī¾[:š­]
        b = Ī¾.b
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            f[I] += N[i]*b*š¤
        end
    end
end

function (op::Operator{:ā«vtdĪ})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    for Ī¾ in š
        š¤ = Ī¾.š¤
        N = Ī¾[:š­]
        t = Ī¾.t
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            f[I] += N[i]*t*š¤
        end
    end
end

function (op::Operator{:ā«vgdĪ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    for Ī¾ in š
        š¤ = Ī¾.š¤
        N = Ī¾[:š­]
        g = Ī¾.g
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            for (j,xā±¼) in enumerate(š)
                J = xā±¼.š¼
                k[I,J] += op.Ī±*N[i]*N[j]*š¤
            end
            f[I] += op.Ī±*N[i]*g*š¤
        end
    end
end

function (op::Operator{:ā«Ī»gdĪ})(ap1::T,ap2::S,g::AbstractMatrix{Float64},q::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    for j in 1:length(ap1.š)
        Ī¾ā = ap1.š[j]
        Ī¾ā = ap2.š[j]
        š¤ = Ī¾ā.š¤
        N = Ī¾ā[:š­]
        NĢ = Ī¾ā[:š­]
        gĢ = Ī¾ā.g
        for (k,xā) in enumerate(ap2.š)
            K = xā.š¼
            for (i,xįµ¢) in enumerate(ap1.š)
                I = xįµ¢.š¼
                g[I,K] -= N[i]*NĢ[k]*š¤
            end
            q[K] -= NĢ[k]*gĢ*š¤
        end
    end
end

function (op::Operator{:ā«āšvgdĪ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    for Ī¾ in ap.š
        N = Ī¾[:š­]
        Bā = Ī¾[:āš­āx]
        Bā = Ī¾[:āš­āy]
        Bā = Ī¾[:āš­āz]
        š¤ = Ī¾.š¤
        nā = Ī¾.nā
        nā = Ī¾.nā
        nā = Ī¾.nā
        g = Ī¾.g
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            for (j,xā±¼) in enumerate(š)
                J = xā±¼.š¼
                k[I,J] -= op.k*((Bā[i]*nā+Bā[i]*nā+Bā[i]*nā)*N[j] + N[i]*(Bā[j]*nā+Bā[j]*nā+Bā[j]*nā))*š¤
            end
            f[I] -= op.k*(Bā[i]*nā+Bā[i]*nā+Bā[i]*nā)*g*š¤
        end
    end
end

function (op::Operator{:ā«āĢšvgdĪ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    for Ī¾ in š
        š¤ = Ī¾.š¤
        N = Ī¾[:š­]
        Bā = Ī¾[:āš­āx_]
        Bā = Ī¾[:āš­āy_]
        Bā = Ī¾[:āš­āz_]
        nā = Ī¾.nā
        nā = Ī¾.nā
        nā = Ī¾.nā
        g = Ī¾.g
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            for (j,xā±¼) in enumerate(š)
                J = xā±¼.š¼
                k[I,J] += (Bā[i]*nā+Bā[i]*nā+Bā[i]*nā)*N[j]*š¤
            end
            f[I] += (Bā[i]*nā+Bā[i]*nā+Bā[i]*nā)*g*š¤
        end
    end
end

function (op::Operator{:g})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64};dof::Symbol=:d) where T<:AbstractElement{:Poi1}
    x = ap.š[1]
    j = x.š¼
    g = getproperty(x,dof)
    for i in 1:length(f)
        f[i] -= k[i,j]*g
    end
    k[j,:] .= 0.
    k[:,j] .= 0.
    k[j,j] = 1.
    f[j] = g
end

"""
Plane Strain
"""
function (op::Operator{:ā«ā«Īµįµ¢ā±¼Ļįµ¢ā±¼vįµ¢bįµ¢dxdy})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    for Ī¾ in š
        N = Ī¾[:š­]
        Bā = Ī¾[:āš­āx]
        Bā = Ī¾[:āš­āy]
        š¤ = Ī¾.š¤
        E = op.E
        Ī½ = op.Ī½
        Cįµ¢įµ¢įµ¢įµ¢ = E/(1-Ī½^2)
        Cįµ¢įµ¢ā±¼ā±¼ = E*Ī½/(1-Ī½^2)
        Cįµ¢ā±¼įµ¢ā±¼ = E/2/(1+Ī½)
        bā = Ī¾.bā
        bā = Ī¾.bā
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            for (j,xā±¼) in enumerate(š)
                J = xā±¼.š¼
                k[2*I-1,2*J-1] += (Cįµ¢įµ¢įµ¢įµ¢*Bā[i]*Bā[j] + Cįµ¢ā±¼įµ¢ā±¼*Bā[i]*Bā[j])*š¤
                k[2*I-1,2*J]   += (Cįµ¢įµ¢ā±¼ā±¼*Bā[i]*Bā[j] + Cįµ¢ā±¼įµ¢ā±¼*Bā[i]*Bā[j])*š¤
                k[2*I,2*J-1]   += (Cįµ¢įµ¢ā±¼ā±¼*Bā[i]*Bā[j] + Cįµ¢ā±¼įµ¢ā±¼*Bā[i]*Bā[j])*š¤
                k[2*I,2*J]     += (Cįµ¢įµ¢įµ¢įµ¢*Bā[i]*Bā[j] + Cįµ¢ā±¼įµ¢ā±¼*Bā[i]*Bā[j])*š¤
            end
            f[2*I-1] += N[i]*bā*š¤
            f[2*I]   += N[i]*bā*š¤
        end
    end
end

function (op::Operator{:ā«ā«Īµįµ¢ā±¼Ļįµ¢ā±¼dxdy})(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    for Ī¾ in š
        Bā = Ī¾[:āš­āx]
        Bā = Ī¾[:āš­āy]
        š¤ = Ī¾.š¤
        E = op.E
        Ī½ = op.Ī½
        Cįµ¢įµ¢įµ¢įµ¢ = E/(1-Ī½^2)
        Cįµ¢įµ¢ā±¼ā±¼ = E*Ī½/(1-Ī½^2)
        Cįµ¢ā±¼įµ¢ā±¼ = E/2/(1+Ī½)
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            for (j,xā±¼) in enumerate(š)
                J = xā±¼.š¼
                k[2*I-1,2*J-1] += (Cįµ¢įµ¢įµ¢įµ¢*Bā[i]*Bā[j] + Cįµ¢ā±¼įµ¢ā±¼*Bā[i]*Bā[j])*š¤
                k[2*I-1,2*J]   += (Cįµ¢įµ¢ā±¼ā±¼*Bā[i]*Bā[j] + Cįµ¢ā±¼įµ¢ā±¼*Bā[i]*Bā[j])*š¤
                k[2*I,2*J-1]   += (Cįµ¢įµ¢ā±¼ā±¼*Bā[i]*Bā[j] + Cįµ¢ā±¼įµ¢ā±¼*Bā[i]*Bā[j])*š¤
                k[2*I,2*J]     += (Cįµ¢įµ¢įµ¢įµ¢*Bā[i]*Bā[j] + Cįµ¢ā±¼įµ¢ā±¼*Bā[i]*Bā[j])*š¤
            end
        end
    end
end

function (op::Operator{:ā«ā«Īµįµįµ¢ā±¼Ļįµįµ¢ā±¼dxdy})(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    for Ī¾ in š
        Bā = Ī¾[:āš­āx]
        Bā = Ī¾[:āš­āy]
        š¤ = Ī¾.š¤
        E = op.E
        Ī½ = op.Ī½
        Cįµ = E/(1-2*Ī½)
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            for (j,xā±¼) in enumerate(š)
                J = xā±¼.š¼
                k[2*I-1,2*J-1] += Cįµ/3*Bā[i]*Bā[j]*š¤
                k[2*I-1,2*J]   += Cįµ/3*Bā[i]*Bā[j]*š¤
                k[2*I,2*J-1]   += Cįµ/3*Bā[i]*Bā[j]*š¤
                k[2*I,2*J]     += Cįµ/3*Bā[i]*Bā[j]*š¤
            end
        end
    end
end

function (op::Operator{:ā«ā«Īµįµįµ¢ā±¼Ļįµįµ¢ā±¼dxdy})(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    for Ī¾ in š
        Bā = Ī¾[:āš­āx]
        Bā = Ī¾[:āš­āy]
        š¤ = Ī¾.š¤
        E = op.E
        Ī½ = op.Ī½
        Cįµ = E/(1+Ī½)
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            for (j,xā±¼) in enumerate(š)
                J = xā±¼.š¼
                k[2*I-1,2*J-1] += Cįµ*( 2/3*Bā[i]*Bā[j]+1/2*Bā[i]*Bā[j])*š¤
                k[2*I-1,2*J]   += Cįµ*(-1/3*Bā[i]*Bā[j]+1/2*Bā[i]*Bā[j])*š¤
                k[2*I,2*J-1]   += Cįµ*(-1/3*Bā[i]*Bā[j]+1/2*Bā[i]*Bā[j])*š¤
                k[2*I,2*J]     += Cįµ*( 2/3*Bā[i]*Bā[j]+1/2*Bā[i]*Bā[j])*š¤
            end
        end
    end
end

function (op::Operator{:ā«ā«vįµ¢bįµ¢dxdy})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    for Ī¾ in š
        N = Ī¾[:š­]
        š¤ = Ī¾.š¤
        bā = Ī¾.bā
        bā = Ī¾.bā
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            f[2*I-1] += N[i]*bā*š¤
            f[2*I]   += N[i]*bā*š¤
        end
    end
end

function (op::Operator{:ā«vįµ¢tįµ¢ds})(ap::T,f::Vector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    for Ī¾ in š
        N = Ī¾[:š­]
        š¤ = Ī¾.š¤
        tā = Ī¾.tā
        tā = Ī¾.tā
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            f[2*I-1] += N[i]*tā*š¤
            f[2*I]   += N[i]*tā*š¤
        end
    end
end

function (op::Operator{:ā«Ī»įµ¢gįµ¢ds})(ap1::T,ap2::S,g::AbstractMatrix{Float64},q::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    for j in 1:length(ap1.š)
        Ī¾ā = ap1.š[j]
        Ī¾ā = ap2.š[j]
        š¤ = Ī¾ā.š¤
        N = Ī¾ā[:š­]
        NĢ = Ī¾ā[:š­]
        gā = Ī¾.gā
        gā = Ī¾.gā
        nāā = Ī¾.nāā
        nāā = Ī¾.nāā
        nāā = Ī¾.nāā
        for (k,xā) in enumerate(ap2.š)
            K = xā.š¼
            for (i,xįµ¢) in enumerate(ap1.š)
                I = xįµ¢.š¼
                NĢāNįµ¢ = NĢ[k]*N[i]
                g[2*K-1,2*I-1] -= nāā*NĢāNįµ¢*š¤
                g[2*K-1,2*I]   -= nāā*NĢāNįµ¢*š¤
                g[2*K,2*I-1]   -= nāā*NĢāNįµ¢*š¤
                g[2*K,2*I]     -= nāā*NĢāNįµ¢*š¤
            end
            q[2*K-1] -= NĢ[k]*(nāā*gā+nāā*gā)*š¤
            q[2*K]   -= NĢ[k]*(nāā*gā+nāā*gā)*š¤
        end
    end
end

function (op::Operator{:ā«Ļįµ¢ā±¼nā±¼gįµ¢ds})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    for Ī¾ in š
        š¤ = Ī¾.š¤
        N = Ī¾[:š­]
        Bā = Ī¾[:āš­āx]
        Bā = Ī¾[:āš­āy]
        E = op.E
        Ī½ = op.Ī½
        Cįµ¢įµ¢įµ¢įµ¢ = E/(1-Ī½^2)
        Cįµ¢įµ¢ā±¼ā±¼ = E*Ī½/(1-Ī½^2)
        Cįµ¢ā±¼įµ¢ā±¼ = E/2/(1+Ī½)
        nāā = Ī¾.nāā
        nāā = Ī¾.nāā
        nāā = Ī¾.nāā
        gā = Ī¾.gā
        gā = Ī¾.gā
        nā = Ī¾.nā
        nā = Ī¾.nā
        Cāāā = Cįµ¢įµ¢įµ¢įµ¢*nā*nāā+Cįµ¢įµ¢ā±¼ā±¼*nā*nāā
        Cāāā = Cįµ¢įµ¢įµ¢įµ¢*nā*nāā+Cįµ¢įµ¢ā±¼ā±¼*nā*nāā
        Cāāā = Cįµ¢įµ¢ā±¼ā±¼*nā*nāā+Cįµ¢įµ¢įµ¢įµ¢*nā*nāā
        Cāāā = Cįµ¢įµ¢ā±¼ā±¼*nā*nāā+Cįµ¢įµ¢įµ¢įµ¢*nā*nāā
        Cāāā = Cįµ¢ā±¼įµ¢ā±¼*(nā*nāā+nā*nāā)
        Cāāā = Cįµ¢ā±¼įµ¢ā±¼*(nā*nāā+nā*nāā)
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            for (j,xā±¼) in enumerate(š)
                J = xā±¼.š¼
                k[2*I-1,2*J-1] -= (Cāāā*(N[i]*Bā[j]+Bā[i]*N[j]) + Cāāā*(N[i]*Bā[j]+Bā[i]*N[j]))*š¤
                k[2*I-1,2*J]   -= (Cāāā*N[i]*Bā[j] + Cāāā*Bā[i]*N[j] + Cāāā*N[i]*Bā[j] + Cāāā*Bā[i]*N[j])*š¤
                k[2*I,2*J-1]   -= (Cāāā*N[i]*Bā[j] + Cāāā*Bā[i]*N[j] + Cāāā*N[i]*Bā[j] + Cāāā*Bā[i]*N[j])*š¤
                k[2*I,2*J]     -= (Cāāā*(N[i]*Bā[j]+Bā[i]*N[j]) + Cāāā*(N[i]*Bā[j]+Bā[i]*N[j]))*š¤
            end
            f[2*I-1] -= ((Cāāā*Bā[i]+Cāāā*Bā[i])*gā + (Cāāā*Bā[i]+Cāāā*Bā[i])*gā)*š¤
            f[2*I]   -= ((Cāāā*Bā[i]+Cāāā*Bā[i])*gā + (Cāāā*Bā[i]+Cāāā*Bā[i])*gā)*š¤
        end
    end
end

function (op::Operator{:ā«Ļįµ¢ā±¼nā±¼gįµ¢vįµ¢gįµ¢ds})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    for Ī¾ in š
        š¤ = Ī¾.š¤
        N = Ī¾[:š­]
        Bā = Ī¾[:āš­āx]
        Bā = Ī¾[:āš­āy]
        E = op.E
        Ī½ = op.Ī½
        Ī± = op.Ī±
        Cįµ¢įµ¢įµ¢įµ¢ = E/(1-Ī½^2)
        Cįµ¢įµ¢ā±¼ā±¼ = E*Ī½/(1-Ī½^2)
        Cįµ¢ā±¼įµ¢ā±¼ = E/2/(1+Ī½)
        nāā = Ī¾.nāā
        nāā = Ī¾.nāā
        nāā = Ī¾.nāā
        gā = Ī¾.gā
        gā = Ī¾.gā
        nā = Ī¾.nā
        nā = Ī¾.nā
        Cāāā = Cįµ¢įµ¢įµ¢įµ¢*nā*nāā+Cįµ¢įµ¢ā±¼ā±¼*nā*nāā
        Cāāā = Cįµ¢įµ¢įµ¢įµ¢*nā*nāā+Cįµ¢įµ¢ā±¼ā±¼*nā*nāā
        Cāāā = Cįµ¢įµ¢ā±¼ā±¼*nā*nāā+Cįµ¢įµ¢įµ¢įµ¢*nā*nāā
        Cāāā = Cįµ¢įµ¢ā±¼ā±¼*nā*nāā+Cįµ¢įµ¢įµ¢įµ¢*nā*nāā
        Cāāā = Cįµ¢ā±¼įµ¢ā±¼*(nā*nāā+nā*nāā)
        Cāāā = Cįµ¢ā±¼įµ¢ā±¼*(nā*nāā+nā*nāā)
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            for (j,xā±¼) in enumerate(š)
                J = xā±¼.š¼
                k[2*I-1,2*J-1] -= (Cāāā*(N[i]*Bā[j]+Bā[i]*N[j]) + Cāāā*(N[i]*Bā[j]+Bā[i]*N[j]) - Ī±*N[i]*nāā*N[j])*š¤
                k[2*I-1,2*J]   -= (Cāāā*N[i]*Bā[j] + Cāāā*Bā[i]*N[j] + Cāāā*N[i]*Bā[j] + Cāāā*Bā[i]*N[j] -Ī±*N[i]*nāā*N[j])*š¤
                k[2*I,2*J-1]   -= (Cāāā*N[i]*Bā[j] + Cāāā*Bā[i]*N[j] + Cāāā*N[i]*Bā[j] + Cāāā*Bā[i]*N[j] -Ī±*N[i]*nāā*N[j])*š¤
                k[2*I,2*J]     -= (Cāāā*(N[i]*Bā[j]+Bā[i]*N[j]) + Cāāā*(N[i]*Bā[j]+Bā[i]*N[j]) -Ī±*N[i]*nāā*N[j])*š¤
            end
            f[2*I-1] -= ((Cāāā*Bā[i]+Cāāā*Bā[i]-Ī±*N[i]*nāā)*gā + (Cāāā*Bā[i]+Cāāā*Bā[i]-Ī±*N[i]*nāā)*gā)*š¤
            f[2*I]   -= ((Cāāā*Bā[i]+Cāāā*Bā[i]-Ī±*N[i]*nāā)*gā + (Cāāā*Bā[i]+Cāāā*Bā[i]-Ī±*N[i]*nāā)*gā)*š¤
        end
    end
end

function (op::Operator{:ā«ĻĢįµ¢ā±¼nā±¼gįµ¢ds})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    for Ī¾ in š
        š¤ = Ī¾.š¤
        N = Ī¾[:š­]
        Bā = Ī¾[:āš­āx_]
        Bā = Ī¾[:āš­āy_]
        E = op.E
        Ī½ = op.Ī½
        Cįµ¢įµ¢įµ¢įµ¢ = E/(1-Ī½^2)
        Cįµ¢įµ¢ā±¼ā±¼ = E*Ī½/(1-Ī½^2)
        Cįµ¢ā±¼įµ¢ā±¼ = E/2/(1+Ī½)
        nāā = Ī¾.nāā
        nāā = Ī¾.nāā
        nāā = Ī¾.nāā
        gā = Ī¾.gā
        gā = Ī¾.gā
        nā = Ī¾.nā
        nā = Ī¾.nā
        Cāāā = Cįµ¢įµ¢įµ¢įµ¢*nā*nāā+Cįµ¢įµ¢ā±¼ā±¼*nā*nāā
        Cāāā = Cįµ¢įµ¢įµ¢įµ¢*nā*nāā+Cįµ¢įµ¢ā±¼ā±¼*nā*nāā
        Cāāā = Cįµ¢įµ¢ā±¼ā±¼*nā*nāā+Cįµ¢įµ¢įµ¢įµ¢*nā*nāā
        Cāāā = Cįµ¢įµ¢ā±¼ā±¼*nā*nāā+Cįµ¢įµ¢įµ¢įµ¢*nā*nāā
        Cāāā = Cįµ¢ā±¼įµ¢ā±¼*(nā*nāā+nā*nāā)
        Cāāā = Cįµ¢ā±¼įµ¢ā±¼*(nā*nāā+nā*nāā)
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            for (j,xā±¼) in enumerate(š)
                J = xā±¼.š¼
                k[2*I-1,2*J-1] += (Cāāā*Bā[i]*N[j] + Cāāā*Bā[i]*N[j])*š¤
                k[2*I-1,2*J]   += (Cāāā*Bā[i]*N[j] + Cāāā*Bā[i]*N[j])*š¤
                k[2*I,2*J-1]   += (Cāāā*Bā[i]*N[j] + Cāāā*Bā[i]*N[j])*š¤
                k[2*I,2*J]     += (Cāāā*Bā[i]*N[j] + Cāāā*Bā[i]*N[j])*š¤
            end
            f[2*I-1] += ((Cāāā*Bā[i] + Cāāā*Bā[i])*gā + (Cāāā*Bā[i] + Cāāā*Bā[i])*gā)*š¤
            f[2*I]   += ((Cāāā*Bā[i] + Cāāā*Bā[i])*gā + (Cāāā*Bā[i] + Cāāā*Bā[i])*gā)*š¤
        end
    end
end

function (op::Operator{:ā«ĻĢĢįµ¢ā±¼nā±¼gįµ¢ds})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    E = op.E
    Ī½ = op.Ī½
    for Ī¾ in š
        š¤ = Ī¾.š¤
        N = Ī¾[:š­]
        Bā = Ī¾[:āš­āx]
        Bā = Ī¾[:āš­āy]
        BĢā = Ī¾[:āš­āx_]
        BĢā = Ī¾[:āš­āy_]
        Cįµ¢įµ¢įµ¢įµ¢ = E/(1-Ī½^2)
        Cįµ¢įµ¢ā±¼ā±¼ = E*Ī½/(1-Ī½^2)
        Cįµ¢ā±¼įµ¢ā±¼ = E/2/(1+Ī½)
        nāā = Ī¾.nāā
        nāā = Ī¾.nāā
        nāā = Ī¾.nāā
        gā = Ī¾.gā
        gā = Ī¾.gā
        nā = Ī¾.nā
        nā = Ī¾.nā
        Cāāā = Cįµ¢įµ¢įµ¢įµ¢*nā*nāā+Cįµ¢įµ¢ā±¼ā±¼*nā*nāā
        Cāāā = Cįµ¢įµ¢įµ¢įµ¢*nā*nāā+Cįµ¢įµ¢ā±¼ā±¼*nā*nāā
        Cāāā = Cįµ¢įµ¢ā±¼ā±¼*nā*nāā+Cįµ¢įµ¢įµ¢įµ¢*nā*nāā
        Cāāā = Cįµ¢įµ¢ā±¼ā±¼*nā*nāā+Cįµ¢įµ¢įµ¢įµ¢*nā*nāā
        Cāāā = Cįµ¢ā±¼įµ¢ā±¼*(nā*nāā+nā*nāā)
        Cāāā = Cįµ¢ā±¼įµ¢ā±¼*(nā*nāā+nā*nāā)
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            ĪBā = Bā[i]-BĢā[i]
            ĪBā = Bā[i]-BĢā[i]
            for (j,xā±¼) in enumerate(š)
                J = xā±¼.š¼
                k[2*I-1,2*J-1] -= (Cāāā*(N[i]*Bā[j]+ĪBā*N[j]) + Cāāā*(N[i]*Bā[j]+ĪBā*N[j]))*š¤
                k[2*I-1,2*J]   -= (Cāāā*N[i]*Bā[j] + Cāāā*ĪBā*N[j] + Cāāā*N[i]*Bā[j] + Cāāā*ĪBā*N[j])*š¤
                k[2*I,2*J-1]   -= (Cāāā*N[i]*Bā[j] + Cāāā*ĪBā*N[j] + Cāāā*N[i]*Bā[j] + Cāāā*ĪBā*N[j])*š¤
                k[2*I,2*J]     -= (Cāāā*(N[i]*Bā[j]+ĪBā*N[j]) + Cāāā*(N[i]*Bā[j]+ĪBā*N[j]))*š¤
            end
            f[2*I-1] -= ((Cāāā*ĪBā+Cāāā*ĪBā)*gā + (Cāāā*ĪBā+Cāāā*ĪBā)*gā)*š¤
            f[2*I]   -= ((Cāāā*ĪBā+Cāāā*ĪBā)*gā + (Cāāā*ĪBā+Cāāā*ĪBā)*gā)*š¤
        end
    end
end

function (op::Operator{:ā«vįµ¢gįµ¢ds})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    for Ī¾ in š
        š¤ = Ī¾.š¤
        N = Ī¾[:š­]
        nāā = Ī¾.nāā
        nāā = Ī¾.nāā
        nāā = Ī¾.nāā
        gā = Ī¾.gā
        gā = Ī¾.gā
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            for (j,xā±¼) in enumerate(š)
                J = xā±¼.š¼
                k[2*I-1,2*J-1] += op.Ī±*N[i]*nāā*N[j]*š¤
                k[2*I,2*J-1]   += op.Ī±*N[i]*nāā*N[j]*š¤
                k[2*I-1,2*J]   += op.Ī±*N[i]*nāā*N[j]*š¤
                k[2*I,2*J]     += op.Ī±*N[i]*nāā*N[j]*š¤
            end
            f[2*I-1] += op.Ī±*N[i]*(nāā*gā+nāā*gā)*š¤
            f[2*I]   += op.Ī±*N[i]*(nāā*gā+nāā*gā)*š¤
        end
    end
end


"""
Kirchhoff-Love plate
"""
function (op::Operator{:ā«Īŗįµ¢ā±¼Mįµ¢ā±¼dĪ©})(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    D = op.D
    Ī½ = op.Ī½
    for Ī¾ in š
        Bāā = Ī¾[:āĀ²š­āxĀ²]
        Bāā = Ī¾[:āĀ²š­āxāy]
        Bāā = Ī¾[:āĀ²š­āyĀ²]
        š¤ = Ī¾.š¤
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            for (j,xā±¼) in enumerate(š)
                J = xā±¼.š¼
                k[I,J] += D*(Bāā[i]*Bāā[j] + Ī½*(Bāā[i]*Bāā[j] + Bāā[i]*Bāā[j]) + Bāā[i]*Bāā[j] + 2*(1-Ī½)*Bāā[i]*Bāā[j])*š¤
            end
        end
    end
end

function (op::Operator{:ā«wqdĪ©})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    for Ī¾ in š
        š¤ = Ī¾.š¤
        N = Ī¾[:š­]
        q = Ī¾.q
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            f[I] += N[i]*q*š¤
        end
    end
end

function (op::Operator{:ā«wVdĪ})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    for Ī¾ in š
        š¤ = Ī¾.š¤
        N = Ī¾[:š­]
        V = Ī¾.V
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            f[I] += N[i]*V*š¤
        end
    end
end

function (op::Operator{:ā«ĪøāMāādĪ})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    for Ī¾ in š
        š¤ = Ī¾.š¤
        Bā = Ī¾[:āš­āx]
        Bā = Ī¾[:āš­āy]
        nā = Ī¾.nā
        nā = Ī¾.nā
        M = Ī¾.M
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            f[I] -= (Bā[i]*nā+Bā[i]*nā)*M*š¤
        end
    end
end

function (op::Operator{:ā«āwMdĪ})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    nā,nā = getš(ap)
    for Ī¾ in š
        š¤ = Ī¾.š¤
        Bā = Ī¾[:āš­āx]
        Bā = Ī¾[:āš­āy]
        Mā = Ī¾.Mā
        Mā = Ī¾.Mā
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            f[I] -= (Bā[i]*Mā+Bā[i]*Mā)*š¤
        end
    end
end

function (op::Operator{:ā«āšvĪødĪ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    for Ī¾ in š
        š¤ = Ī¾.š¤
        nā = Ī¾.nā
        nā = Ī¾.nā
        Bā = Ī¾[:āš­āx]
        Bā = Ī¾[:āš­āy]
        Īø = Ī¾.Īø
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            Īøįµ¢ = Bā[i]*nā+Bā[i]*nā
            for (j,xā±¼) in enumerate(š)
                J = xā±¼.š¼
                Īøā±¼ = Bā[j]*nā+Bā[j]*nā
                k[I,J] += op.Ī±*Īøįµ¢*Īøā±¼*š¤
            end
            f[I] += op.Ī±*Īøįµ¢*Īø*š¤
        end
    end
end

function (op::Operator{:ā«MāāĪødĪ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    Ī± = op.Ī±
    D = op.D
    Ī½ = op.Ī½
    for Ī¾ in š
        š¤ = Ī¾.š¤
        nā = Ī¾.nā
        nā = Ī¾.nā
        Dāā = -D*(nā^2+Ī½*nā^2)
        Dāā = -2*D*nā*nā*(1-Ī½)
        Dāā = -D*(Ī½*nā^2+nā^2)
        Bā = Ī¾[:āš­āx]
        Bā = Ī¾[:āš­āy]
        Bāā = Ī¾[:āĀ²š­āxĀ²]
        Bāā = Ī¾[:āĀ²š­āxāy]
        Bāā = Ī¾[:āĀ²š­āyĀ²]
        Īø = Ī¾.Īø
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            Īøįµ¢ = Bā[i]*nā + Bā[i]*nā
            Mįµ¢ = Dāā*Bāā[i] + Dāā*Bāā[i] + Dāā*Bāā[i]
            for (j,xā±¼) in enumerate(š)
                J = xā±¼.š¼
                Īøā±¼ = Bā[j]*nā + Bā[j]*nā
                Mā±¼ = Dāā*Bāā[j] + Dāā*Bāā[j] + Dāā*Bāā[j]
                k[I,J] += (Mįµ¢*Īøā±¼+Īøįµ¢*Mā±¼+Ī±*Īøįµ¢*Īøā±¼)*š¤
            end
            f[I] += (Mįµ¢+Ī±*Īøįµ¢)*Īø*š¤
        end
    end
end

function (op::Operator{:ā«VgdĪ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    Ī± = op.Ī±
    D = op.D
    Ī½ = op.Ī½
    for Ī¾ in š
        š¤ = Ī¾.š¤
        nā = Ī¾.nā
        nā = Ī¾.nā
        sā = -Ī¾.nā
        sā = Ī¾.nā
        Dāāā = -D*(nā + nā*sā*sā + Ī½*nā*sā*sā)
        Dāāā = -D*(nā + nā*sā*sā + 2*nā*sā*sā + (nā*sā*sā - nā*sā*sā - nā*sā*sā)*Ī½)
        Dāāā = -D*(nā + nā*sā*sā + 2*nā*sā*sā + (nā*sā*sā - nā*sā*sā - nā*sā*sā)*Ī½)
        Dāāā = -D*(nā + nā*sā*sā + Ī½*nā*sā*sā)
        N = Ī¾[:š­]
        Bāāā = Ī¾[:āĀ³š­āxĀ³]
        Bāāā = Ī¾[:āĀ³š­āxĀ²āy]
        Bāāā = Ī¾[:āĀ³š­āxāyĀ²]
        Bāāā = Ī¾[:āĀ³š­āyĀ³]
        g = Ī¾.g
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            Vįµ¢ = Dāāā*Bāāā[i] + Dāāā*Bāāā[i] + Dāāā*Bāāā[i] + Dāāā*Bāāā[i]
            for (j,xā±¼) in enumerate(š)
                J = xā±¼.š¼
                Vā±¼ = Dāāā*Bāāā[j] + Dāāā*Bāāā[j] + Dāāā*Bāāā[j] + Dāāā*Bāāā[j]
                k[I,J] += (-Vįµ¢*N[j]-N[i]*Vā±¼+Ī±*N[i]*N[j])*š¤
            end
            f[I] += (-Vįµ¢+Ī±*N[i])*g*š¤
        end
    end
end

function (op::Operator{:ā«MĢāāĪødĪ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    nā = š[1].nā
    nā = š[1].nā
    sā = š[1].sā
    sā = š[1].sā
    D = op.D
    Ī½ = op.Ī½
    Dāā = -D*(nā^2+Ī½*nā^2)
    Dāā = -2*D*nā*nā*(1-Ī½)
    Dāā = -D*(Ī½*nā^2+nā^2)
    for Ī¾ in š
        š¤ = Ī¾.š¤
        Bā = Ī¾[:āš­āx]
        Bā = Ī¾[:āš­āy]
        Bāā = Ī¾[:āĀ²š­āxĀ²]
        Bāā = Ī¾[:āĀ²š­āxāy]
        Bāā = Ī¾[:āĀ²š­āyĀ²]
        BĢāā = Ī¾[:āĀ²š­āxĀ²_]
        BĢāā = Ī¾[:āĀ²š­āxāy_]
        BĢāā = Ī¾[:āĀ²š­āyĀ²_]
        Īø = Ī¾.Īø
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            Īøįµ¢ = Bā[i]*nā + Bā[i]*nā
            Mįµ¢ = Dāā*Bāā[i] + Dāā*Bāā[i] + Dāā*Bāā[i]
            MĢįµ¢ = Dāā*BĢāā[i] + Dāā*BĢāā[i] + Dāā*BĢāā[i]
            for (j,xā±¼) in enumerate(š)
                J = xā±¼.š¼
                Īøā±¼ = Bā[j]*nā + Bā[j]*nā
                Mā±¼ = Dāā*Bāā[j] + Dāā*Bāā[j] + Dāā*Bāā[j]
                k[I,J] += (Mįµ¢*Īøā±¼+Īøįµ¢*Mā±¼-MĢįµ¢*Īøā±¼)*š¤
            end
            f[I] += (Mįµ¢-MĢįµ¢)*Īø*š¤
        end
    end
end

function (op::Operator{:ā«VĢgdĪ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    nā = š[1].nā
    nā = š[1].nā
    sā = š[1].sā
    sā = š[1].sā
    D = op.D
    Ī½ = op.Ī½
    Dāāā = -D*(nā + nā*sā*sā + Ī½*nā*sā*sā)
    Dāāā = -D*(nā*sā*sā + (nā*sā*sā + nā)*Ī½)
    Dāāā = -D*(nā + nā*sā*sā + nā*sā*sā + (-nā - nā*sā*sā - nā*sā*sā)*Ī½)
    Dāāā = -D*(nā + nā*sā*sā + nā*sā*sā + (-nā - nā*sā*sā - nā*sā*sā)*Ī½)
    Dāāā = -D*(nā*sā*sā + (nā*sā*sā + nā)*Ī½)
    Dāāā = -D*(nā + nā*sā*sā + Ī½*nā*sā*sā)
    for Ī¾ in š
        š¤ = Ī¾.š¤
        N = Ī¾[:š­]
        Bāāā = Ī¾[:āāĀ²š­āxĀ²āx]
        Bāāā = Ī¾[:āāĀ²š­āxĀ²āy]
        Bāāā = Ī¾[:āāĀ²š­āxāyāx]
        Bāāā = Ī¾[:āāĀ²š­āxāyāy]
        Bāāā = Ī¾[:āāĀ²š­āyĀ²āx]
        Bāāā = Ī¾[:āāĀ²š­āyĀ²āy]
        BĢāāā = Ī¾[:āāĀ²š­āxĀ²āx_]
        BĢāāā = Ī¾[:āāĀ²š­āxĀ²āy_]
        BĢāāā = Ī¾[:āāĀ²š­āxāyāx_]
        BĢāāā = Ī¾[:āāĀ²š­āxāyāy_]
        BĢāāā = Ī¾[:āāĀ²š­āyĀ²āx_]
        BĢāāā = Ī¾[:āāĀ²š­āyĀ²āy_]
        g = Ī¾.g
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            Vįµ¢ = Dāāā*Bāāā[i] + Dāāā*Bāāā[i] + Dāāā*Bāāā[i] + Dāāā*Bāāā[i] + Dāāā*Bāāā[i] + Dāāā*Bāāā[i]
            VĢįµ¢ = Dāāā*BĢāāā[i] + Dāāā*BĢāāā[i] + Dāāā*BĢāāā[i] + Dāāā*BĢāāā[i] + Dāāā*BĢāāā[i] + Dāāā*BĢāāā[i]
            for (j,xā±¼) in enumerate(š)
                J = xā±¼.š¼
                Vā±¼ = Dāāā*Bāāā[j] + Dāāā*Bāāā[j] + Dāāā*Bāāā[j] + Dāāā*Bāāā[j] + Dāāā*Bāāā[j] + Dāāā*Bāāā[j]
                k[I,J] -= (Vįµ¢*N[j]+N[i]*Vā±¼-VĢįµ¢*N[j])*š¤
            end
            f[I] -= (Vįµ¢-VĢįµ¢)*g*š¤
        end
    end
end

function (op::Operator{:wĪMāā})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    for Ī¾ in š
        N = Ī¾[:š­]
        ĪM = Ī¾.ĪM
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            f[I] -= N[i]*ĪM
        end
    end
end

function (op::Operator{:ĪMāāg})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    D = op.D
    Ī½ = op.Ī½
    Ī± = op.Ī±
    for Ī¾ in š
        Īnāsā = Ī¾.Īnāsā
        Īnāsānāsā = Ī¾.Īnāsānāsā
        Īnāsā = Ī¾.Īnāsā
        Dāā = - D*(Īnāsā+Īnāsā*Ī½)
        Dāā = - D*(1-Ī½)*Īnāsānāsā
        Dāā = - D*(Īnāsā*Ī½+Īnāsā)
        N = Ī¾[:š­]
        Bāā = Ī¾[:āĀ²š­āxĀ²]
        Bāā = Ī¾[:āĀ²š­āxāy]
        Bāā = Ī¾[:āĀ²š­āyĀ²]
        g = Ī¾.g
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            ĪMāāįµ¢ = Dāā*Bāā[i] + Dāā*Bāā[i] + Dāā*Bāā[i]
            for (j,xā±¼) in enumerate(š)
                J = xā±¼.š¼
                ĪMāāā±¼ = Dāā*Bāā[j] + Dāā*Bāā[j] + Dāā*Bāā[j]
                k[I,J] += ĪMāāįµ¢*N[j] + N[i]*ĪMāāā±¼ + Ī±*N[i]*N[j]
            end
            f[I] += (ĪMāāįµ¢ + Ī±*N[i])*g
        end
    end
end

function (op::Operator{:ĪMĢāāg})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    š = ap.š; š = ap.š
    D = op.D
    Ī½ = op.Ī½
    for Ī¾ in š
        Īnāsā = Ī¾.Īnāsā
        Īnāsānāsā = Ī¾.Īnāsānāsā
        Īnāsā = Ī¾.Īnāsā
        Dāā = - D*(Īnāsā+Īnāsā*Ī½)
        Dāā = - D*(1-Ī½)*Īnāsānāsā
        Dāā = - D*(Īnāsā*Ī½+Īnāsā)
        N = Ī¾[:š­]
        Bāā = Ī¾[:āĀ²š­āxĀ²]
        Bāā = Ī¾[:āĀ²š­āxāy]
        Bāā = Ī¾[:āĀ²š­āyĀ²]
        BĢāā = Ī¾[:āĀ²š­āxĀ²_]
        BĢāā = Ī¾[:āĀ²š­āxāy_]
        BĢāā = Ī¾[:āĀ²š­āyĀ²_]
        g = Ī¾.g
        for (i,xįµ¢) in enumerate(š)
            I = xįµ¢.š¼
            ĪMāāįµ¢ = Dāā*Bāā[i] + Dāā*Bāā[i] + Dāā*Bāā[i]
            ĪMĢāāįµ¢ = Dāā*BĢāā[i] + Dāā*BĢāā[i] + Dāā*BĢāā[i]
            for (j,xā±¼) in enumerate(š)
                J = xā±¼.š¼
                ĪMāāā±¼ = Dāā*Bāā[j] + Dāā*Bāā[j] + Dāā*Bāā[j]
                k[I,J] += ĪMāāįµ¢*N[j] + N[i]*ĪMāāā±¼ + ĪMĢāāįµ¢*N[j]
            end
            f[I] += (ĪMāāįµ¢ + ĪMĢāāįµ¢)*g
        end
    end
end

"""
Error Estimates
"""
function (op::Operator{:Lā})(ap::T) where T<:AbstractElement
    ĪuĀ²= 0
    uĢĀ² = 0
    for Ī¾ in ap.š
        š¤ = Ī¾.š¤
        N = Ī¾[:š­]
        uĢįµ¢ = Ī¾.u
        uįµ¢ = 0
        for (i,xįµ¢) in enumerate(š)
            uįµ¢ += N[i]*xįµ¢.d
        end
        ĪuĀ² += (uįµ¢ - uĢįµ¢)^2*š¤
        uĢĀ²  += uĢįµ¢^2*š¤
    end
    return ĪuĀ², uĢĀ²
end

function (op::Operator{:Lā})(aps::Vector{T}) where T<:AbstractElement
    LāNorm_ĪuĀ²= 0
    LāNorm_uĢĀ² = 0
    for ap in aps
        ĪuĀ², uĢĀ² = op(ap)
        LāNorm_ĪuĀ² += ĪuĀ²
        LāNorm_uĢĀ²  += uĢĀ²
    end
    return (LāNorm_ĪuĀ²/LāNorm_uĢĀ²)^0.5
end

function (op::Operator{:Hā})(ap::T) where T<:AbstractElement
    ĪāuĀ²= 0
    āuĢĀ² = 0
    ĪuĀ²= 0
    uĢĀ² = 0
    for Ī¾ in ap.š
        š¤ = Ī¾.š¤
        N = Ī¾[:š­]
        Bā = Ī¾[:āš­āx]
        Bā = Ī¾[:āš­āy]
        Bā = Ī¾[:āš­āz]
        uĢįµ¢ = Ī¾.u
        āuĢįµ¢āx = Ī¾.āuāx
        āuĢįµ¢āy = Ī¾.āuāy
        āuĢįµ¢āz = Ī¾.āuāz
        uįµ¢ = 0.
        āuįµ¢āx = 0.
        āuįµ¢āy = 0.
        āuįµ¢āz = 0.
        for (i,xįµ¢) in enumerate(ap.š)
            uįµ¢ += N[i]*xįµ¢.d
            āuįµ¢āx += Bā[i]*xįµ¢.d
            āuįµ¢āy += Bā[i]*xįµ¢.d
            āuįµ¢āz += Bā[i]*xįµ¢.d
        end
        ĪāuĀ² += ((āuįµ¢āx - āuĢįµ¢āx)^2 + (āuįµ¢āy - āuĢįµ¢āy)^2 + (āuįµ¢āz - āuĢįµ¢āz)^2)*š¤
        āuĢĀ² += (āuĢįµ¢āx^2 + āuĢįµ¢āy^2 + āuĢįµ¢āz^2)*š¤
        ĪuĀ² += (uįµ¢ - uĢįµ¢)^2*š¤
        uĢĀ² += uĢįµ¢^2*š¤
    end
    return ĪāuĀ², āuĢĀ², ĪuĀ², uĢĀ²
end

function (op::Operator{:Hā})(aps::Vector{T}) where T<:AbstractElement
    HāNorm_ĪuĀ²= 0
    HāNorm_uĢĀ² = 0
    LāNorm_ĪuĀ²= 0
    LāNorm_uĢĀ² = 0
    for ap in aps
        ĪāuĀ², āuĢĀ², ĪuĀ², uĢĀ² = op(ap)
        HāNorm_ĪuĀ² += ĪuĀ² + ĪāuĀ²
        HāNorm_uĢĀ²  += uĢĀ² + āuĢĀ²
        LāNorm_ĪuĀ² += ĪuĀ²
        LāNorm_uĢĀ²  += uĢĀ²
    end
    return (HāNorm_ĪuĀ²/HāNorm_uĢĀ²)^0.5, (LāNorm_ĪuĀ²/LāNorm_uĢĀ²)^0.5
end

function (op::Operator{:Hā_PlaneStress})(ap::T) where T<:AbstractElement
    ĪWĀ²= 0
    WĢĀ² = 0
    ĪuĀ²= 0
    uĢĀ² = 0
    E = op.E
    Ī½ = op.Ī½
    Cįµ¢įµ¢įµ¢įµ¢ = E/(1-Ī½^2)
    Cįµ¢įµ¢ā±¼ā±¼ = E*Ī½/(1-Ī½^2)
    Cįµ¢ā±¼įµ¢ā±¼ = E/2/(1+Ī½)
    for Ī¾ in ap.š
        š¤ = Ī¾.š¤
        N = Ī¾[:š­]
        Bā = Ī¾[:āš­āx]
        Bā = Ī¾[:āš­āy]
        uĢā = Ī¾.u
        uĢā = Ī¾.v
        āuĢāāx = Ī¾.āuāx
        āuĢāāy = Ī¾.āuāy
        āuĢāāx = Ī¾.āvāx
        āuĢāāy = Ī¾.āvāy
        ĪµĢāā = āuĢāāx
        ĪµĢāā = āuĢāāy
        ĪµĢāā = āuĢāāy + āuĢāāx
        ĻĢāā = Cįµ¢įµ¢įµ¢įµ¢*ĪµĢāā + Cįµ¢įµ¢ā±¼ā±¼*ĪµĢāā
        ĻĢāā = Cįµ¢įµ¢įµ¢įµ¢*ĪµĢāā + Cįµ¢įµ¢ā±¼ā±¼*ĪµĢāā
        ĻĢāā = Cįµ¢ā±¼įµ¢ā±¼*ĪµĢāā
        uā = 0.
        uā = 0.
        Īµāā = 0.
        Īµāā = 0.
        Īµāā = 0.
        for (i,xįµ¢) in enumerate(ap.š)
            uā += N[i]*xįµ¢.dā
            uā += N[i]*xįµ¢.dā
            Īµāā += Bā[i]*xįµ¢.dā
            Īµāā += Bā[i]*xįµ¢.dā
            Īµāā += Bā[i]*xįµ¢.dā + Bā[i]*xįµ¢.dā
        end
        Ļāā = Cįµ¢įµ¢įµ¢įµ¢*Īµāā + Cįµ¢įµ¢ā±¼ā±¼*Īµāā
        Ļāā = Cįµ¢įµ¢įµ¢įµ¢*Īµāā + Cįµ¢įµ¢ā±¼ā±¼*Īµāā
        Ļāā = Cįµ¢ā±¼įµ¢ā±¼*Īµāā
        ĪWĀ² += 0.5*((Ļāā-ĻĢāā)*(Īµāā-ĪµĢāā) + (Ļāā-ĻĢāā)*(Īµāā-ĪµĢāā) + (Ļāā-ĻĢāā)*(Īµāā-ĪµĢāā))*š¤
        WĢĀ² += 0.5*(Ļāā*Īµāā + Ļāā*Īµāā + Ļāā*Īµāā)*š¤
        ĪuĀ² += ((uā - uĢā)^2 + (uā - uĢā)^2)*š¤
        uĢĀ² += (uĢā^2 + uĢā^2)*š¤
    end
    return ĪWĀ², WĢĀ², ĪuĀ², uĢĀ²
end

function (op::Operator{:Hā_PlaneStress})(aps::Vector{T}) where T<:AbstractElement
    HāNorm_ĪWĀ²= 0.0
    HāNorm_WĢĀ² = 0.0
    LāNorm_ĪuĀ²= 0.0
    LāNorm_uĢĀ² = 0.0
    for ap in aps
        ĪWĀ², WĢĀ², ĪuĀ², uĢĀ² = op(ap)
        HāNorm_ĪWĀ² += ĪWĀ²
        HāNorm_WĢĀ²  += WĢĀ²
        LāNorm_ĪuĀ² += ĪuĀ²
        LāNorm_uĢĀ²  += uĢĀ²
    end
    return (HāNorm_ĪWĀ²/HāNorm_WĢĀ²)^0.5, (LāNorm_ĪuĀ²/LāNorm_uĢĀ²)^0.5
end

function setāš¢!(ap::T) where T<:AbstractElement
    for Ī¾ in ap.š
        N = Ī¾[:š­]
        Bā = Ī¾[:āš­āx_]
        Bā = Ī¾[:āš­āy_]
        Bā = Ī¾[:āš­āz_]
        š = (Ī¾.x,Ī¾.y,Ī¾.z)
        u = 0.
        āuāx = 0.
        āuāy = 0.
        āuāz = 0.
        for (i,xįµ¢) in enumerate(ap.š)
            u += N[i]*x.d
            āuāx += Bā[i]*x.d
            āuāy += Bā[i]*x.d
            āuāz += Bā[i]*x.d
        end
        Ī¾.x = š[1]
        Ī¾.y = š[2]
        Ī¾.z = š[3]
        Ī¾.u = u
        Ī¾.āuāx = āuāx
        Ī¾.āuāy = āuāy
        Ī¾.āuāz = āuāz
    end
end

function (op::Operator{:Hā})(aps::Vector{T}) where T<:AbstractElement
    HāNorm_ĪuĀ²= 0
    HāNorm_uĢĀ² = 0
    HāNorm_ĪuĀ²= 0
    HāNorm_uĢĀ² = 0
    HāNorm_ĪuĀ²= 0
    HāNorm_uĢĀ² = 0
    LāNorm_ĪuĀ²= 0
    LāNorm_uĢĀ² = 0
    for ap in aps
        ĪāĀ³uĀ², āĀ³uĢĀ²,ĪāĀ²uĀ², āĀ²uĢĀ²,ĪāuĀ², āuĢĀ², ĪuĀ², uĢĀ² = op(ap)
        HāNorm_ĪuĀ² += ĪuĀ² + ĪāuĀ² + ĪāĀ²uĀ² + ĪāĀ³uĀ²
        HāNorm_uĢĀ²  += uĢĀ² + āuĢĀ² + āĀ²uĢĀ² + āĀ³uĢĀ²
        HāNorm_ĪuĀ² += ĪuĀ² + ĪāuĀ² + ĪāĀ²uĀ²
        HāNorm_uĢĀ²  += uĢĀ² + āuĢĀ² + āĀ²uĢĀ²
        HāNorm_ĪuĀ² += ĪuĀ² + ĪāuĀ²
        HāNorm_uĢĀ²  += uĢĀ² + āuĢĀ²
        LāNorm_ĪuĀ² += ĪuĀ²
        LāNorm_uĢĀ²  += uĢĀ²
    end
    return (HāNorm_ĪuĀ²/HāNorm_uĢĀ²)^0.5, (HāNorm_ĪuĀ²/HāNorm_uĢĀ²)^0.5, (HāNorm_ĪuĀ²/HāNorm_uĢĀ²)^0.5, (LāNorm_ĪuĀ²/LāNorm_uĢĀ²)^0.5
end
function (op::Operator{:Hā})(ap::T) where T<:AbstractElement
    ĪāĀ³uĀ²= 0
    āĀ³uĢĀ² = 0
    ĪāĀ²uĀ²= 0
    āĀ²uĢĀ² = 0
    ĪāuĀ²= 0
    āuĢĀ² = 0
    ĪuĀ²= 0
    uĢĀ² = 0
    for Ī¾ in ap.š
        š¤ = getš¤(ap,Ī¾)
        N = Ī¾[:š­]
        Bā = Ī¾[:āš­āx]
        Bā = Ī¾[:āš­āy]
        Bāā = Ī¾[:āĀ²š­āxĀ²]
        Bāā = Ī¾[:āĀ²š­āxāy]
        Bāā = Ī¾[:āĀ²š­āyĀ²]
        Bāāā = Ī¾[:āĀ³š­āxĀ³]
        Bāāā = Ī¾[:āĀ³š­āxĀ²āy]
        Bāāā = Ī¾[:āĀ³š­āxāyĀ²]
        Bāāā = Ī¾[:āĀ³š­āyĀ³]
        uĢįµ¢ = Ī¾.u
        āuĢįµ¢āx = Ī¾.āuāx
        āuĢįµ¢āy = Ī¾.āuāy
        āĀ²uĢįµ¢āxĀ² = Ī¾.āĀ²uāxĀ²
        āĀ²uĢįµ¢āxāy = Ī¾.āĀ²uāxāy
        āĀ²uĢįµ¢āyĀ² = Ī¾.āĀ²uāyĀ²
        āĀ³uĢįµ¢āxĀ³ = Ī¾.āĀ³uāxĀ³
        āĀ³uĢįµ¢āxĀ²āy = Ī¾.āĀ³uāxĀ²āy
        āĀ³uĢįµ¢āxāyĀ² = Ī¾.āĀ³uāxāyĀ²
        āĀ³uĢįµ¢āyĀ³ = Ī¾.āĀ³uāyĀ³
        uįµ¢ = 0.
        āuįµ¢āx = 0.
        āuįµ¢āy = 0.
        āĀ²uįµ¢āxĀ² = 0.
        āĀ²uįµ¢āxāy = 0.
        āĀ²uįµ¢āyĀ² = 0.
        āĀ³uįµ¢āxĀ³ = 0.
        āĀ³uįµ¢āxĀ²āy = 0.
        āĀ³uįµ¢āxāyĀ² = 0.
        āĀ³uįµ¢āyĀ³ = 0.
        for i in 1:length(ap.š)
            xįµ¢ = ap.š[i]
            I = xįµ¢.id
            uįµ¢ += N[i]*xįµ¢.d
            āuįµ¢āx += Bā[i]*xįµ¢.d
            āuįµ¢āy += Bā[i]*xįµ¢.d
            āĀ²uįµ¢āxĀ² += Bāā[i]*xįµ¢.d
            āĀ²uįµ¢āxāy += Bāā[i]*xįµ¢.d
            āĀ²uįµ¢āyĀ² += Bāā[i]*xįµ¢.d
            āĀ³uįµ¢āxĀ³ += Bāāā[i]*xįµ¢.d
            āĀ³uįµ¢āxĀ²āy += Bāāā[i]*xįµ¢.d
            āĀ³uįµ¢āxāyĀ² += Bāāā[i]*xįµ¢.d
            āĀ³uįµ¢āyĀ³ += Bāāā[i]*xįµ¢.d
        end
        ĪāĀ³uĀ² += ((āĀ³uįµ¢āxĀ³ - āĀ³uĢįµ¢āxĀ³)^2 + (āĀ³uįµ¢āxĀ²āy - āĀ³uĢįµ¢āxĀ²āy)^2 + (āĀ³uįµ¢āxāyĀ² - āĀ³uĢįµ¢āxāyĀ²)^2 + (āĀ³uįµ¢āyĀ³ - āĀ³uĢįµ¢āyĀ³)^2)*š¤
        āĀ³uĢĀ² += (āĀ³uĢįµ¢āxĀ³^2 + āĀ³uĢįµ¢āxĀ²āy^2  + āĀ³uĢįµ¢āxāyĀ²^2+ āĀ³uĢįµ¢āyĀ³^2)*š¤
        ĪāĀ²uĀ² += ((āĀ²uįµ¢āxĀ² - āĀ²uĢįµ¢āxĀ²)^2 + (āĀ²uįµ¢āxāy - āĀ²uĢįµ¢āxāy)^2 + (āĀ²uįµ¢āyĀ² - āĀ²uĢįµ¢āyĀ²)^2)*š¤
        āĀ²uĢĀ² += (āĀ²uĢįµ¢āxĀ²^2 + āĀ²uĢįµ¢āxāy^2 + āĀ²uĢįµ¢āyĀ²^2)*š¤
        ĪāuĀ² += ((āuįµ¢āx - āuĢįµ¢āx)^2 + (āuįµ¢āy - āuĢįµ¢āy)^2)*š¤
        āuĢĀ² += (āuĢįµ¢āx^2 + āuĢįµ¢āy^2)*š¤
        ĪuĀ² += (uįµ¢ - uĢįµ¢)^2*š¤
        uĢĀ² += uĢįµ¢^2*š¤
    end
    return ĪāĀ³uĀ², āĀ³uĢĀ², ĪāĀ²uĀ², āĀ²uĢĀ², ĪāuĀ², āuĢĀ², ĪuĀ², uĢĀ²
end

function setāš¢!(aps::Vector{T}) where T<:AbstractElement
    for ap in aps
        setāš¢!(ap)
    end
end
