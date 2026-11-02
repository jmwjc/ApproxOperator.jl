# Patch test: the canonical smoke test for the whole assembly pipeline.
#
# This file doubles as the minimal end-to-end example of the package:
#   1. read a mesh from Gmsh,
#   2. build integration elements,
#   3. prescribe material data and evaluate shape functions,
#   4. assemble a stiffness matrix and a load vector through operator pairs,
#   5. solve and measure the L2 error against the exact solution.
#
# The declared function u = x + y is reproduced exactly by the (multi-)linear
# approximation space, so the discrete solution must be exact up to round-off.

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy

using Test
import Gmsh: gmsh

function run_patchtest(io::IO)
    𝑢(x,y,z) = x + y

    gmsh.open(joinpath(@__DIR__,"patchtest.msh"))

    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()

    nₚ = length(nodes)
    k = zeros(nₚ,nₚ)
    f = zeros(nₚ)

    # Bilinear form a(v,u) = ∫ k ∇v·∇u dΩ on the domain.
    # The conductivity is a material constant and is passed as a keyword.
    elements = getElements(nodes, entities["Ω"])
    set∇𝝭!(elements)
    𝑎 = ∫∫∇v∇udxdy=>elements
    𝑎(k; kᶜ=1.0)

    # Penalised essential boundary condition on Γᵍ.
    # The boundary value is prescribed data; the penalty coefficient is a
    # material constant and goes in as a keyword.
    elements = getElements(nodes, entities["Γᵍ"])
    prescribe!(elements, :g=>𝑢)
    set𝝭!(elements)
    𝑓 = ∫vgdΓ=>elements
    𝑓(k,f; α=1e7)

    d = k\f
    push!(nodes, :d=>d)

    # L2 error, sampled on a finer quadrature rule.
    elements = getElements(nodes, entities["Ω"], 10)
    prescribe!(elements, :u=>𝑢)
    set∇𝝭!(elements)
    L₂error = L₂(elements)
    println(io, "L2 error: ", L₂error)
    return L₂error
end

@testset "patch test (heat, Tri3, penalty BC)" begin
    gmsh.initialize()
    try
        L₂error = run_patchtest(stdout)
        # A linear solution is in the approximation space, so this is a
        # consistency check rather than an accuracy measurement.
        @test L₂error < 1.0e-6
    finally
        gmsh.finalize()
    end
end
