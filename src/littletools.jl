
const listπ­ = (:π­,)
const listβπ­ = (:π­,:βπ­βx,:βπ­βy,:βπ­βz)
const listββπ­ = (:π­,:βπ­βx,:βπ­βy)
const listβΒ²π­ = (:π­,:βπ­βx,:βπ­βy,:βΒ²π­βxΒ²,:βΒ²π­βxβy,:βΒ²π­βyΒ²,:βπ­βz,:βΒ²π­βxβz,:βΒ²π­βyβz,:βΒ²π­βzΒ²)
const listβΒ²βπ­ = (:π­,:βπ­βx,:βπ­βy,:βΒ²π­βxΒ²,:βΒ²π­βxβy,:βΒ²π­βyΒ²)
const listβΒ³π­ = (:π­,:βπ­βx,:βπ­βy,:βΒ²π­βxΒ²,:βΒ²π­βxβy,:βΒ²π­βyΒ²,:βΒ³π­βxΒ³,:βΒ³π­βxΒ²βy,:βΒ³π­βxβyΒ²,:βΒ³π­βyΒ³)
for (π­,π,list) in ((:checkπ­,:getπ,:listπ­),
                   (:checkβπ­,:getβπ,:listβπ­),
                   (:checkββπ­,:getββπ,:listββπ­),
                   (:checkβΒ²π­,:getβΒ²π,:listβΒ²π­),
                   (:checkβΒ²βπ­,:getβΒ²βπ,:listβΒ²βπ­),
                   (:checkβΒ³π­,:getβΒ³π,:listβΒ³π­))
    @eval begin
        function $π­(a::T,f::Matrix{Float64},π::Matrix{Float64},πΚ°::Matrix{Float64}) where T<:AbstractElement
            n = getππ(a)
            for ΞΎ in a.π
                π = getπ(a,ΞΎ)
                π€ = ΞΎ.π€
                πs = $π(a,π)
                for i in 1:n
                    for (j,π_) in enumerate(πs)
                        π[i,j] = π_[i]
                    end
                end
                fill!(πΚ°,0.0)
                for (k,πα΅’) in enumerate(a.π)
                    πα΅’ = getπ(a,(πα΅’.x,πα΅’.y,πα΅’.z))
                    for i in 1:n
                        for (j,s) in enumerate($list)
                            πΚ°[i,j] += ΞΎ[s][k]*πα΅’[i]
                        end
                    end
                end
                f .+= (π .- πΚ°).^2 .* π€
            end
        end

        function $π­(as::Vector{T}) where T<:ReproducingKernel
            nα΅ = getππ(as[1])
            n = length($list)
            f = zeros(nα΅,n)
            π = zeros(nα΅,n)
            πΚ° = zeros(nα΅,n)
            for a in as
                $π­(a,f,π,πΚ°)
            end
            return f.^0.5
        end
    end
end
