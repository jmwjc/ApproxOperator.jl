# finite element analysis for 1D bar problem
# tuthor: @wujc
# problem: EA*d²u/dx² = x,   x∈(0,1)
#          u(0) = 0.
#          EAdu/dx(1) = 1.

using ApproxOperator

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
set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵗ"])

# prescribe
prescribe!(elements["Ω"],:b=>(x,y,z)->x)
prescribe!(elements["Γᵗ"],:t=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:g=>(x,y,z)->0.0)

# set operator
ops = [
    Operator{:∫vₓuₓdx}(:EA=>1.0),
    Operator{:∫vbdΩ}(),
    Operator{:∫vtdΓ}(),
    Operator{:g}()
]

# assembly
k = zeros(nₚ,nₚ)
f = zeros(nₚ)
d = zeros(nₚ)
push!(nodes,:d=>d)
ops[1](elements["Ω"],k)
ops[2](elements["Ω"],f)
ops[3](elements["Γᵗ"],f)
ops[4](elements["Γᵍ"],k,f)

# solve
d = k\f