# finite element analysis for 1D bar problem
# tuthor: @wujc
# problem: EA*d²u/dx² = x,   x∈(0,1)
#          u(0) = 0.
#          EAdu/dx(1) = 1.

using ApproxOperator, Printf

# length of bar
Lb = 1.
# material coefficients
EA = 1.

# num of nodes
nₚ = 11

# num of cells
nₑ = nₚ - 1

# nodes
x = zeros(nₚ)
for i in 1:nₑ
    x[i+1] = i*Lb/nₑ
end
nodes = ApproxOperator.Node(:x=>x,:y=>zeros(nₚ),:z=>zeros(nₚ))

# elements
elements = Dict{String,Any}()
elements["Ω"] = [ApproxOperator.Element{:Seg2}([nodes[i],nodes[i+1]]) for i in 1:nₑ]
elements["Γᵍ"] = [ApproxOperator.Element{:Poi1}([nodes[1]])]
elements["Γᵗ"] = [ApproxOperator.Element{:Poi1}([nodes[nₚ]])]

# set ingeration points
set𝓖!(elements["Ω"],:SegGI2)
set𝓖!(elements["Γᵗ"],:PoiGI1)
set𝓖!(elements["Γᵍ"],:PoiGI1)

# set shape functions
set_memory_𝝭!(elements["Ω"],:𝝭,:∂𝝭∂x)
set_memory_𝝭!(elements["Γᵗ"],:𝝭)
set_memory_𝝭!(elements["Γᵍ"],:𝝭)
set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵗ"])
set𝝭!(elements["Γᵍ"])

# prescribe
prescribe!(elements["Ω"],:σₙ=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:αₙ=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:εᵖₙ=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:Δεₙ=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:ε=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:g=>(x,y,z)->0.0)

# set operator
ops = [
    Operator{:∫vₓσdx}(:E=>100.0,:K=>100.0,:σy=>1.0,:tol=>1e14),
    Operator{:∫vtdΓ}(),
    Operator{:∫vgdΓ}(:α=>1e15)
]

# assembly
k = zeros(nₚ,nₚ)
fint = zeros(nₚ)
fext = zeros(nₚ)
d = zeros(nₚ)
Δd = zeros(nₚ)
push!(nodes,:d=>d)
push!(nodes,:Δd=>Δd)

total_steps = 100
F = 2.0
for n in 1:total_steps
    fill!(k,0.0)
    fill!(fint,0.0)
    fill!(fext,0.0)

    prescribe!(elements["Γᵗ"],:t=>(x,y,z)->F*n/total_steps)
    ops[2](elements["Γᵗ"],fext)
    ops[1](elements["Ω"],k,fint)
    ops[3](elements["Γᵍ"],k,fext)

    Δd .= k\(fext-fint)
    d .+= Δd
    # update Δεₙ₊₁
    for ap in elements["Ω"]
        𝓒 = ap.𝓒;𝓖 = ap.𝓖
        for ξ in 𝓖
            Δε = 0.0
            ε = 0.0
            B = ξ[:∂𝝭∂x]
            for (i,xᵢ) in enumerate(𝓒)
                Δε += B[i]*xᵢ.Δd
                ε += B[i]*xᵢ.d
            end
            ξ.Δεₙ = Δε
            ξ.ε = ε
        end
    end

    a = elements["Ω"][5]
    ξ = a.𝓖[1]
    @printf  "Load step=%i,f=%e, σₙ=%e, εₙ₊₁=%e \n" n F*n/total_steps ξ.σₙ ξ.ε
end
