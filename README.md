# ApproxOperator.jl

Weak-form assembly for finite-element and meshfree approximations in Julia.

The package turns a mesh into integration elements and assembles weak forms into
global matrices and vectors. Operators are plain Julia functions, so the weak
form you read in a paper is the weak form you write in code:

```julia
(∫∫εᵢⱼσᵢⱼdxdy => elements)(k; E = 210e9, ν = 0.3)
```

- [Installation](#installation)
- [A complete example](#a-complete-example)
- [The data model](#the-data-model)
- [Writing a weak form](#writing-a-weak-form)
- [Adding an element type](#adding-an-element-type)
- [What is implemented](#what-is-implemented)
- [Testing](#testing)

---

## Installation

Requires Julia 1.9 or newer and a working Gmsh installation (supplied by the
`Gmsh.jl` artifact, no system install needed).

```julia
using Pkg
Pkg.develop(path = "path/to/ApproxOperator.jl")
```

```julia
using ApproxOperator          # getPhysicalGroups, get𝑿ᵢ, getElements are exported
import Gmsh: gmsh
```

The importer entry points are re-exported, so a qualified
`import ApproxOperator.GmshImport: getElements` is not needed.

See [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) for the environment the published
results were produced with, and [`RELEASE_NOTES.md`](RELEASE_NOTES.md) for what
changed in v0.2.0 — including the breaking changes to the operator signatures.

---

## A complete example

A patch test: `-∇²u = 0` on the unit square with a penalty-enforced Dirichlet
condition. The exact solution `u = x + y` lies in the approximation space, so the
discrete solution is exact up to round-off.

```julia
using ApproxOperator
import ApproxOperator.Heat: ∫vgdΓ, ∫∫∇v∇udxdy, L₂
using Gmsh: gmsh

𝑢(x, y, z) = x + y

gmsh.initialize()
gmsh.open("patchtest.msh")          # any .msh with physical groups Ω and Γᵍ

entities = getPhysicalGroups()
nodes    = get𝑿ᵢ()

k = zeros(length(nodes), length(nodes))
f = zeros(length(nodes))

# ∫ k ∇v·∇u dΩ  on the domain.  The conductivity is a keyword argument.
elements = getElements(nodes, entities["Ω"])
set∇𝝭!(elements)
(∫∫∇v∇udxdy => elements)(k; kᶜ = 1.0)

# α ∫ (u - g) v dΓ  on the boundary.  The boundary value `g` is prescribed data.
elements = getElements(nodes, entities["Γᵍ"])
prescribe!(elements, :g => 𝑢)
set𝝭!(elements)
(∫vgdΓ => elements)(k, f; α = 1e7)

d = k \ f
push!(nodes, :d => d)

# L2 error, evaluated on a finer quadrature rule.
elements = getElements(nodes, entities["Ω"], 10)
prescribe!(elements, :u => 𝑢)
set∇𝝭!(elements)
println("L2 error: ", L₂(elements))     # 6.6e-8

gmsh.finalize()
```

The three-step rhythm — **build elements → set the shape functions → assemble** —
is the same for every problem in the package.

### A vector problem: 2-D elasticity

The same rhythm with a vector-valued unknown, two degrees of freedom per node,
and material constants as keyword arguments. This example has been run; it
prints the deflection under the load.

```julia
using ApproxOperator
import ApproxOperator.Elasticity: ∫∫εᵢⱼσᵢⱼdxdy
using Gmsh: gmsh

E, ν = 210.0e9, 0.3

gmsh.initialize()
gmsh.open("square.msh")             # physical groups: Ω, Γᵍ

entities = getPhysicalGroups()
nodes    = get𝑿ᵢ()

n = length(nodes)
K = zeros(2n, 2n)
F = zeros(2n)

# ∫ ε(v):C:ε(u) dΩ  — the element stiffness.  `E` and `ν` are keywords.
elements = getElements(nodes, entities["Ω"], 2)
set∇𝝭!(elements)
(∫∫εᵢⱼσᵢⱼdxdy => elements)(K; E = E, ν = ν)

# A downward body force, assembled with the same loop over elements.
# `prescribe!` supplies the load per quadrature point, the operator reads it.
import ApproxOperator.Elasticity: ∫∫vᵢbᵢdxdy
prescribe!(elements, :b₁ => (x, y) -> 0.0, :b₂ => (x, y) -> -1.0e6)
(∫∫vᵢbᵢdxdy => elements)(F)

# Clamped bottom edge: eliminate the constrained DOFs.
fixed = Int[]
for xᵢ in nodes
    if xᵢ.y ≈ 0.0 atol = 1.0e-12
        push!(fixed, 2xᵢ.𝐼 - 1, 2xᵢ.𝐼)     # both components
    end
end
free = setdiff(eachindex(F), fixed)

U = zeros(length(F))
U[free] = K[free, free] \ F[free]

push!(nodes, :d₁ => U[1:2:end], :d₂ => U[2:2:end])
println("max |u₂| = ", maximum(abs, U[2:2:end]))   # 2.36e-6

gmsh.finalize()
```

Notes on this example:

- Both components of the load are prescribed separately (`:b₁`, `:b₂`); the
  operator reads one scalar per component and per quadrature point.
- The stiffness operator and the load operator use the **same** `elements`, so
  `set∇𝝭!` is done once.
- Traction boundary conditions (`∫σᵢⱼnⱼgᵢds`) additionally need the **surface
  normal derivative of the shape functions** (`ξ.n₁₁`, `ξ.n₁₂`, `ξ.n₂₂`), which
  are produced when the boundary elements are built with `normal = true` and the
  displacement elements are available; see the operator docstring.

### Which shape-function setter do I need?

Every operator docstring states this. In short:

| setter | fills | needed by |
|---|---|---|
| `set𝝭!` | `ξ[:𝝭]` | mass matrices, load vectors, penalty terms, mixed coupling |
| `set∇𝝭!` | first derivatives | second-order PDEs: elasticity, heat, Stokes, beams, Mindlin plates |
| `set∇²𝝭!` | second derivatives | fourth-order PDEs: Kirchhoff plates and shells |
| `set∇̂³𝝭!` | third derivatives | reproducing-kernel operators |

See [`docs/src/operators.md`](docs/src/operators.md) for the per-operator list.

---

## The data model

There are two node types and one shared dictionary.

```
Element{:Tri3}
 ├── 𝓒 :: Vector{𝑿ᵢ}     support points — the nodes carrying the shape functions
 └── 𝓖 :: Vector{𝑿ₛ}     quadrature points — where the integrand is evaluated

Node
 ├── index :: (𝑔, 𝐺, 𝐶, 𝑠)   or (𝐼,) for a support point
 └── data  :: Dict{Symbol, Tuple{Int, Vector{Float64}}}
                                                slot   values
```

Every node of an element points at the **same** `data` dictionary. A field read
`ξ.name` resolves to

```
data[name][2][  node.index[ data[name][1] ] ]
```

so the slot number selects which component of the index tuple picks the value:

| slot | for `𝑿ᵢ` (support) | for `𝑿ₛ` (quadrature) |
|---|---|---|
| 1 | indexed by `𝐼` — per node | indexed by `𝑔` — per quadrature point *within* an element |
| 2 | — | indexed by `𝐺` — per quadrature point, globally |
| 3 | — | indexed by `𝐶` — per element |
| 4 | — | indexed by `𝑠` — offset into the per-element DOF list |

This is why `ξ.x` is the coordinate of a quadrature point while `𝒙ᵢ.x` is the
coordinate of a support node, and why `ξ[:𝝭]` is a *view* (`RV`) into the shared
vector rather than a per-element array.

`prescribe!(elements, :E => value)` writes into that dictionary, which is how
boundary data and exact solutions are supplied. **Material constants are not
prescribed** — they are keyword arguments of the operator.

---

## Writing a weak form

An operator is a function that accumulates into `k` (and optionally `f`):

```julia
"""
    ∫∫εᵢⱼσᵢⱼdxdy(ap, k; E, ν)

Two-dimensional linear-elastic stiffness, `∫ ε(u):σ(v) dΩ`.

Plane stress:

    C₁₁₁₁ = E/(1-ν²),   C₁₁₂₂ = Eν/(1-ν²),   C₁₂₁₂ = E/(2(1+ν))

# Keyword arguments
- `E` : Young's modulus
- `ν` : Poisson ratio

# Element requirements
`set∇𝝭!` — reads `ξ[:∂𝝭∂x]`, `ξ[:∂𝝭∂y]`, `ξ.𝑤`.
"""
function ∫∫εᵢⱼσᵢⱼdxdy(ap::T, k::AbstractMatrix{Float64}; E::Real, ν::Real) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤  = ξ.𝑤
        Cᵢᵢᵢᵢ = E/(1-ν^2)
        Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2I-1,2J-1] += (Cᵢᵢᵢᵢ*B₁[i]*B₁[j] + Cᵢⱼᵢⱼ*B₂[i]*B₂[j])*𝑤
                k[2I-1,2J]   += (Cᵢᵢⱼⱼ*B₁[i]*B₂[j] + Cᵢⱼᵢⱼ*B₂[i]*B₁[j])*𝑤
                k[2I,2J-1]   += (Cᵢᵢⱼⱼ*B₂[i]*B₁[j] + Cᵢⱼᵢⱼ*B₁[i]*B₂[j])*𝑤
                k[2I,2J]     += (Cᵢᵢᵢᵢ*B₂[i]*B₂[j] + Cᵢⱼᵢⱼ*B₁[i]*B₁[j])*𝑤
            end
        end
    end
end
```

Rules that the whole package follows:

1. **Material constants are keyword arguments**, typed `::Real`. A missing one
   raises `UndefKeywordError` at the call site.
2. **Boundary values and exact fields are read from `ξ`** and supplied with
   `prescribe!`.
3. **`k` is the stiffness matrix, `f`/`fint` the load vector** — always
   positional, never keyword.
4. If the name `k` is also needed for a material constant, expose it as `kᶜ`
   (conductance); this is the convention in `heat.jl`.
5. Document the quadrature-point fields the operator reads, so a caller knows
   whether `set𝝭!`, `set∇𝝭!` or `set∇²𝝭!` is required.

Several operators can be applied in one call:

```julia
import ApproxOperator.KirchhoffPlate: ∫κᵢⱼMᵢⱼdΩ, ∫VgdΓ

[
    ∫κᵢⱼMᵢⱼdΩ => elements_plate,
    ∫VgdΓ     => elements_edge,
](k; D = D, ν = 0.3)
```

---

## Adding an element type

1. **Shape functions.** Add `set𝝭!`, `set∇𝝭!` (and `set∇²𝝭!` if the element is
   `C¹`) in `src/approximation/`. Follow an existing file of the same family —
   `tri6.jl` for a 2-D simplex, `hex8.jl` for an isoparametric brick.

2. **Register the type** in the `Element2D` / `Element3D` union in
   `src/shapefunction.jl` if operators need to dispatch on it.

3. **Teach the importer the Gmsh tag** — add the entry to `GMSH_ELEMENT_TYPES`
   in `src/preprocession/elements.jl`. The tag is the Gmsh/VTK element
   identifier.

4. **Write a consistency test.** This is not optional; it is the only thing that
   catches a permuted node ordering or a Jacobian used the wrong way round.
   `test/shape_function_consistency.jl` checks, for every element family:

   - partition of unity, `Σᵢ Nᵢ(𝒙) = 1`
   - nodal delta, `Nᵢ(𝒙ⱼ) = δᵢⱼ`
   - reproduction of a linear field and its gradient

   ```julia
   @testset "MyElement" begin
       aps, dg, n, ng = prepare(:Quad8, REF_QUAD8, [(0.0,0.0,0.0)]; gradient=true)
       @test partition_of_unity_error(aps) < 1.0e-14
       @test nodal_delta_error(aps, dg, REF_QUAD8; gradient=true) < 1.0e-14
       Δu, Δ∇u = linear_field_error(aps)
       @test Δu  < 1.0e-14
       @test Δ∇u < 1.0e-13
   end
   ```

5. **Document the `:𝐽` convention** the element relies on. `:𝐽` is
   `det(∂𝒙/∂𝝃)`, so a 1-D element divides its gradient by it, a 2-D simplex
   also divides (the difference product already carries the reference-simplex
   determinant ½), and an isoparametric element never touches it — it inverts the
   full Jacobian instead.

---

## What is implemented

| module | contents |
|---|---|
| `Elasticity` | 3-D and 2-D elasticity, plane stress/strain, mixed `u/p`, compliance form, energy norms |
| `Heat` | conduction, convection–diffusion, Nitsche and penalty boundaries, L₂/H¹/H³ norms, heat-flux norms |
| `KirchhoffPlate` | bending stiffness, corner forces, Kirchhoff shear, penalty boundaries |
| `MindlinPlate` | Reissner–Mindlin, MITC4 and DSG locking-free variants, mixed `w/φ/Q` |
| `KirchhoffShell` | membrane and bending, Nitsche and Hellinger–Reissner boundaries |
| `Beam` | 3-D beam, rotation matrices, block assembly |
| `Hamilton` | space–time rod, least-squares stabilisation, truncation error |
| `Hyperelasticity` | Saint-Venant–Kirchhoff and Neo-Hookean, force and consistent tangent |
| `Plasticity` | J2, Tresca, Mohr–Coulomb, Coulomb contact, phase-field coupling |
| `Phasefield` | 1-D and 2-D AT2 damage, history update, crack normal |
| `Stokes` | viscous term |
| `WeightedResidual` | space–time weighted residual |

**264 operators** with docstrings; see [`docs/src/operators.md`](docs/src/operators.md)
for the generated index. Element families: `Seg2/3`, `Tri3/6`, `Quad4/8`,
`Tet4/10`, `Hex8`, `TriBell`, `TriHermite`, `SegHermite`, `CrouzeixRaviart`,
plus reproducing-kernel approximations (`Linear1D/2D/3D`, `Quadratic*`,
`Cubic*`, `Quartic2D`, `Wave2D`).

Reproducing-kernel operators additionally accept a spatial partition, e.g.
`RegularGrid`, so each quadrature point gets the cloud of nodes inside its
support:

```julia
sp = RegularGrid(nodes.x, nodes.y, nodes.z; n = 2)
elements = getElements(nodes, entities["Ω"], RK, 2, sp)
```

---

## Testing

```julia
using Pkg; Pkg.test()
```

The suite is organised from "does the pipeline run" to "is each element
mathematically correct":

| file | what it pins |
|---|---|
| `test/patchtest.jl` | the end-to-end assembly pipeline |
| `test/importer_characterisation.jl` | the importer's output — element count, data keys, weight sums, centroid nodes, curved geometry |
| `test/shape_function_consistency.jl` | partition of unity, nodal delta, linear-field reproduction for every element family |
| `test/operators_material_parameter.jl` | material constants arrive via keyword arguments; the keyword path equals the `prescribe!` path |
| `test/restored_modules.jl` | the plasticity and phase-field operators are reachable and callable |

---

## Citing

If you use this package in academic work, please cite the papers listed in
[`INTRODUCTION.md`](INTRODUCTION.md).

## Licence

See [`LICENSE`](LICENSE).
