"""
Element{T}<:AbstractElement{T}
"""
struct Element{T} <: AbstractElement{T}
    𝓒::Tuple{Int,Int,Vector{Node{(:𝐼,),1}}}
    𝓖::Tuple{Int,Int,Vector{Node{(:𝑔, :𝐺, :𝐶, :𝑠),4}}}
end

function Base.getproperty(a::T,s::Symbol) where T<:AbstractElement
    if s∈(:𝓒,:𝓖)
        𝓐 =  getfield(a,s)
        return (𝓐[3][𝓐[1]+i] for i in 1:𝓐[2])
    else
        𝓖 = getfield(a,:𝓖)
        ξ = 𝓖[3][𝓖[1]+1]
        return getproperty(ξ,s)
    end
end

function Base.setproperty!(ap::T,s::Symbol,v::Float64) where T<:AbstractElement
    𝓖 = getfield(ap,:𝓖)
    ξ = 𝓖[3][𝓖[1]+1]
    setproperty!(ξ,s,v)
end

for set𝝭 in (:set𝝭!,:set∇𝝭!,:set∇²𝝭!)
    @eval begin
        function $set𝝭(a::T) where T<:AbstractElement
            𝓖 = a.𝓖
            for x in 𝓖
                $set𝝭(a,x)
            end
        end
        $set𝝭(as::Vector{T}) where T<:AbstractElement = $set𝝭.(as)
    end
end

function Base.show(io::IO,::MIME"text/plain",a::T) where T<:AbstractElement
    𝓒 = a.𝓒
    𝓖 = a.𝓖
    println(T)
    println("𝓒:")
    for (i,p) in enumerate(𝓒)
        print("$i. ")
        printinfo(p)
    end
    println("𝓖:")
    for (j,p) in enumerate(𝓖)
        print("$j. ")
        shapes = printinfo(p)
        @printf "         "
        for shape in shapes
            @printf "%13s" string(shape)
        end
        @printf "\n"
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            @printf "  (𝐼 = %i)" I
            for shape in shapes
                @printf " %e" p[shape][i]
            end
            @printf "\n"
        end
    end
end