# ApproxOperator.jl

Weak-form assembly for finite-element and meshfree approximations in Julia.

The package turns a mesh into integration elements and assembles weak forms into
global matrices and vectors. Operators are plain Julia functions, so the weak
form you read in a paper is the weak form you write in code:

```julia
(∫∫εᵢⱼσᵢⱼdxdy => elements)(k; E = 210e9, ν = 0.3)
```

The user guide lives in the repository
[`README.md`](https://github.com/jmwjc/ApproxOperator.jl#readme). This site
contains the generated [operator index](@ref) — every weak-form operator with the
shape-function setter it requires and the material constants it takes.

## Quick reference

The three-step rhythm behind every problem:

```julia
elements = getElements(nodes, entities["Ω"], 2)   # 1. build elements
set∇𝝭!(elements)                                  # 2. shape functions
(∫∫∇v∇udxdy => elements)(k; kᶜ = 1.0)             # 3. assemble
```

| setter | needed by |
|---|---|
| `set𝝭!` | mass matrices, loads, penalty terms, mixed coupling |
| `set∇𝝭!` | second-order PDEs: elasticity, heat, Stokes, beams, Mindlin plates |
| `set∇²𝝭!` | fourth-order PDEs: Kirchhoff plates and shells |
| `set∇̂³𝝭!` | reproducing-kernel operators |

## Modules

`Elasticity`, `Heat`, `KirchhoffPlate`, `MindlinPlate`, `KirchhoffShell`, `Beam`,
`Hamilton`, `Hyperelasticity`, `Plasticity`, `Phasefield`, `Stokes`,
`WeightedResidual`, `Test`.
