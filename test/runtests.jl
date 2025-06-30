using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy

using Test
import Gmsh: gmsh

@testset "ApproxOperator.jl" begin
    𝑢(x,y,z) = x + y
    
    gmsh.initialize()
    gmsh.open("patchtest.msh")
    
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    
    nₚ = length(nodes)
    k = zeros(nₚ,nₚ)
    f = zeros(nₚ)
    
    elements = getElements(nodes, entities["Ω"])
    prescribe!(elements, :k=>1.0)
    set∇𝝭!(elements)
    𝑎 = ∫∫∇v∇udxdy=>elements
    𝑎(k)
    
    elements = getElements(nodes, entities["Γᵍ"])
    prescribe!(elements, :g=>𝑢)
    prescribe!(elements, :α=>1e7)
    set𝝭!(elements)
    𝑓 = ∫vgdΓ=>elements
    𝑓(k,f)
    
    d = k\f
    push!(nodes, :d=>d)
    
    elements = getElements(nodes, entities["Ω"], 10)
    prescribe!(elements, :u=>𝑢)
    set∇𝝭!(elements)
    L₂error = L₂(elements)
    println("L₂ error: ", L₂error)
    @test L₂error ≈ 0.0 atol=1e-6
    gmsh.finalize()
end
