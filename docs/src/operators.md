# Operator index

Every weak-form operator in `src/operation/`, generated from the source
by `scripts/build_operator_index.py`. Do not edit by hand.

**264 operators** across 14 modules.

## How to read this table

| column | meaning |
|---|---|
| *operator* | the Julia function to put on the left of `form => elements` |
| *requires* | which shape-function setter must run before assembly |
| *material constants* | keyword arguments the operator takes |
| *description* | the first line of the operator's docstring |

```julia
set∇𝝭!(elements)                       # whichever the table says
(∫∫εᵢⱼσᵢⱼdxdy => elements)(k; E = 210e9, ν = 0.3)
```

## Contents

- [3-D beam](#3-d-beam) — 15 operators
- [Curved_beam](#curved_beam) — 6 operators
- [Elasticity](#elasticity) — 60 operators
- [Hamiltonian / space–time rod](#hamiltonian-/-spacetime-rod) — 18 operators
- [Heat conduction and convection](#heat-conduction-and-convection) — 27 operators
- [Hyperelasticity](#hyperelasticity) — 16 operators
- [Kirchhoff plate](#kirchhoff-plate) — 14 operators
- [Kirchhoff shell](#kirchhoff-shell) — 28 operators
- [Mindlin (thick) plate](#mindlin-thick-plate) — 57 operators
- [Phase-field fracture](#phase-field-fracture) — 8 operators
- [Plasticity and contact](#plasticity-and-contact) — 6 operators
- [Stokes flow](#stokes-flow) — 1 operators
- [Consistency checks](#consistency-checks) — 4 operators
- [Space–time weighted residual](#spacetime-weighted-residual) — 4 operators

## 3-D beam

`ApproxOperator.Beam`

| operator | requires | material constants | description |
|---|---|---|---|
| `L₂d` | `set𝝭!` | — | Relative L₂ error of the corresponding field against an exact solution, `‖e‖/‖u‖`. The exact field is supplied with `prescribe!`. *(2 overloads)* |
| `∫uᵢgᵢdΓ` | `set𝝭!` | `α` | Penalty enforcement of a displacement boundary condition on a beam, `α ∫ (u-g)·v dΓ`. |
| `∫εᵘNᵘdΩ` | — | — | Element-set overload of `∫εᵘNᵘdΩ`: the same integrand, summed over a vector of elements. See the element-level method for the formula. |
| `∫εᵘNᵠdΩ` | — | — | Element-set overload of `∫εᵘNᵠdΩ`: the same integrand, summed over a vector of elements. See the element-level method for the formula. |
| `∫εᵠNᵘdΩ` | — | — | Element-set overload of `∫εᵠNᵘdΩ`: the same integrand, summed over a vector of elements. See the element-level method for the formula. |
| `∫εᵠNᵠdΩ` | — | — | Element-set overload of `∫εᵠNᵠdΩ`: the same integrand, summed over a vector of elements. See the element-level method for the formula. |
| `∫φᵢgᵢdΓ` | `set𝝭!` | `α` | Penalty enforcement of a rotation boundary condition on a beam, `α ∫ (φ-φ̄)·θ dΓ`. |
| `𝐊ᵘᵘ!` | `set∇𝝭!` | `E` | 3-D beam stiffness, displacement–displacement block, `Kᵘᵘ`. |
| `𝐊ᵘᵠ!` | `set∇𝝭!`, `set𝝭!` | `E` | 3-D beam stiffness, displacement–rotation coupling block, `Kᵘᵠ` (row `u`, column `φ`). |
| `𝐊ᵠᵘ!` | `set∇𝝭!`, `set𝝭!` | `E` | 3-D beam stiffness, rotation–displacement coupling block, `Kᵠᵘ` (row `φ`, column `u`). This is the transpose of [`𝐊ᵘᵠ!`](@ref) for a symmetric element. |
| `𝐊ᵠᵠ!` | `set∇𝝭!`, `set𝝭!` | `E` | 3-D beam stiffness, rotation–rotation block, `Kᵠᵠ`. |
| `𝐐` | — | — | Transformation matrix between the local and global beam DOFs. |
| `𝐑` | — | — | 3×3 rotation matrix from the element axis to the global axes, built from the first two support nodes. |
| `𝐑ᵀ𝐊𝐑` | — | — | Rotate a 3×3 element stiffness from local to global axes, `𝐑ᵀ 𝐊 𝐑`. |

## Curved_beam

`ApproxOperator.Curved_beam`

| operator | requires | material constants | description |
|---|---|---|---|
| `Bᵢvᵢ` | `set∇²𝝭!`, `set∇𝝭!`, `set𝝭!` | — | Curved-beam strain matrix applied to the displacement, assembling `∫ Bᵀ v ds` into the right-hand side. |
| `∫vᵢθᵢds` | `set𝝭!` | `α` | Curved-beam penalty/constraint term over the reference surface, |
| `∫δMEI⁻¹Mds` | `set∇²𝝭!`, `set∇𝝭!`, `set𝝭!` | — | Element-set overload of `∫δMEI⁻¹Mds`: the same integrand, summed over a vector of elements. See the element-level method for the formula. |
| `∫δNEA⁻¹Nds` | `set∇𝝭!`, `set𝝭!` | — | Element-set overload of `∫δNEA⁻¹Nds`: the same integrand, summed over a vector of elements. See the element-level method for the formula. |
| `∫κεds` | `set∇²𝝭!`, `set∇𝝭!` | — | Element-set overload of `∫κεds`: the same integrand, summed over a vector of elements. See the element-level method for the formula. |
| `∫κεγds` | `set∇𝝭!`, `set𝝭!` | `EA`, `EI` | Curved-beam stiffness from membrane strain `ε`, curvature `κ` and shear `γ`, |

## Elasticity

`ApproxOperator.Elasticity`

| operator | requires | material constants | description |
|---|---|---|---|
| `Hₑ` | `set∇𝝭!`, `set𝝭!` | `E`, `ν` | Relative energy-norm error against an exact displacement field. *(2 overloads)* |
| `Hₑ_Incompressible` | `set∇𝝭!`, `set𝝭!` | `E`, `ν` | Relative errors for a purely deviatoric (incompressible-limit) material. *(2 overloads)* |
| `Hₑ_PlaneStrain_Deviatoric` | `set∇𝝭!` | `E`, `ν` | Relative error of the **deviatoric** part of the energy norm. *(2 overloads)* |
| `Hₑ_PlaneStress` | `set∇𝝭!`, `set𝝭!` | `E`, `ν` | Relative energy-norm error for 2-D **plane stress** problems. *(2 overloads)* |
| `Hₑ_up_fem` | `set∇𝝭!`, `set𝝭!` | `E`, `ν` | Relative errors of a displacement–pressure formulation with the **FEM** (continuous) pressure approximation. *(2 overloads)* |
| `Hₑ_up_mix` | `set∇𝝭!`, `set𝝭!` | `E`, `ν` | Relative errors of a displacement–pressure formulation with an independent **mixed** pressure approximation. *(2 overloads)* |
| `Hₑ_up_mix_Q4P1` | `set∇𝝭!`, `set𝝭!` | `E`, `ν` | Relative errors for the Q4/P1 pair (bilinear displacement, elementwise-constant pressure). *(2 overloads)* |
| `L₂` | `set𝝭!` | — | Relative L₂ error of the corresponding field against an exact solution, `‖e‖/‖u‖`. The exact field is supplied with `prescribe!`. *(2 overloads)* |
| `L₂𝑝` | `set𝝭!` | — | Relative L₂ error of the corresponding field against an exact solution, `‖e‖/‖u‖`. The exact field is supplied with `prescribe!`. *(2 overloads)* |
| `getσₙ` | — | — | Principal stresses and their unit directions of a 2-D stress state. *(2 overloads)* |
| `g₂` | — | — | Element-set overload of `g₂`: the same integrand, summed over a vector of elements. See the element-level method for the formula. |
| `∫pnᵢgᵢds` | `set𝝭!` | — | Boundary pressure load `∫ p (n·g) ds`. |
| `∫pnᵢuᵢds` | `set𝝭!` | — | Boundary pressure term `∫ p (n·u) ds`. |
| `∫p∇udΩ` | `set∇𝝭!`, `set𝝭!` | — | Pressure–displacement coupling in 3-D, `∫ p (∇·u) dΩ`. |
| `∫qpdΩ` | `set𝝭!` | `E`, `ν` | Divergence–pressure coupling block, `∫ q (∇·u) dΩ` (or its transpose). |
| `∫sᵢⱼnⱼgᵢds` | `set𝝭!` | — | Deviatoric traction load `∫ (s·n)·g ds`. |
| `∫sᵢⱼnⱼuᵢds` | `set𝝭!` | — | Deviatoric traction coupling `∫ (s·n)·u ds`. |
| `∫vᵢbᵢdΩ` | `set𝝭!` | — | Body force load in 3-D, `f[I] += ∫ b·v dΩ`. |
| `∫vᵢgᵢds` | `set𝝭!` | `α` | Penalty enforcement of an essential boundary condition, `α ∫ (u-g)·v ds`. *(2 overloads)* |
| `∫vᵢgᵢdΓ` | `set𝝭!` | — | Element-set overload of `∫vᵢgᵢdΓ`: the same integrand, summed over a vector of elements. See the element-level method for the formula. |
| `∫vᵢtᵢds` | `set𝝭!` | — | Surface traction load in 2-D, `f[I] += ∫ t·v ds`. |
| `∫vᵢtᵢdΓ` | `set𝝭!` | — | Surface traction load in 3-D, `f[I] += ∫ t·v dΓ`. |
| `∫εᵈᵢⱼσᵈᵢⱼdΩ` | `set∇𝝭!` | `E`, `ν` | Three-dimensional **deviatoric** contribution to the stiffness, `∫ εᵈ(u):σᵈ(v) dΩ`. |
| `∫εᵢⱼσᵢⱼdΩ` | `set∇𝝭!` | `E`, `ν` | Three-dimensional linear-elastic stiffness, `∫ ε(u):σ(v) dΩ`. |
| `∫λᵢgᵢds` | `set𝝭!` | — | Lagrange-multiplier boundary term `∫ λ·g ds` for a mixed essential condition. |
| `∫σᵢⱼnⱼgᵢds` | `set∇𝝭!`, `set𝝭!` | `E`, `ν` | Traction boundary term, `∫ (σ(u)·n)·v ds`, for a prescribed traction `g`. *(2 overloads)* |
| `∫σᵢⱼnⱼgᵢvᵢgᵢds` | `set∇𝝭!`, `set𝝭!` | — | Element-set overload of `∫σᵢⱼnⱼgᵢvᵢgᵢds`: the same integrand, summed over a vector of elements. See the element-level method for the formula. |
| `∫σᵢⱼnⱼuᵢds` | `set𝝭!` | — | Traction–displacement coupling on the boundary, `∫ (σ(u)·n)·v ds`. |
| `∫∫p∇udxdy` | `set∇𝝭!`, `set𝝭!` | — | Pressure–displacement coupling in 2-D, `∫ p (∇·u) dΩ`. |
| `∫∫sᵢⱼsᵢⱼdxdy` | `set𝝭!` | `E`, `ν` | Deviation norm of the stress, `∫ sᵢⱼ sᵢⱼ dΩ` — a post-processing operator, not a stiffness block. |
| `∫∫sᵢⱼεᵢⱼdxdy` | `set∇𝝭!`, `set𝝭!` | — | Deviatoric stress–strain coupling `∫ s : ε dΩ`. |
| `∫∫vᵢbᵢdxdy` | `set𝝭!` | — | Body force load in 2-D, `f[I] += ∫ b·v dΩ`. `b` is supplied with `prescribe!(el, :bᵢ => f)`. |
| `∫∫εᵈᵢⱼσᵈᵢⱼdxdy` | `set∇𝝭!` | `E`, `ν` | Deviatoric contribution to the 2-D stiffness, `∫ εᵈ(u):σᵈ(v) dΩ`. *(2 overloads)* |
| `∫∫εᵛᵢⱼσᵛᵢⱼdxdy` | `set∇𝝭!` | `E`, `ν` | Volumetric (dilatational) contribution to the 2-D stiffness, `∫ εᵛ(u):σᵛ(v) dΩ`. |
| `∫∫εᵢⱼσᵢⱼdxdy` | `set∇𝝭!` | `E`, `ν` | Two-dimensional linear-elastic stiffness, `∫ ε(u):σ(v) dΩ`. *(2 overloads)* |
| `∫∫εᵢⱼσᵢⱼdxdy_PlaneStrian` | `set∇𝝭!` | — | Element-set overload of `∫∫εᵢⱼσᵢⱼdxdy_PlaneStrian`: the same integrand, summed over a vector of elements. See the element-level method for the formula. |
| `∫∫ρvᵢuᵢdxdy` | `set𝝭!` | `ρ` | Consistent mass matrix, `∫ ρ vᵢ uᵢ dΩ`. |
| `∫∫σᵢⱼσₖₗdxdy` | `set𝝭!` | `E`, `ν` | Compliance form, `∫ σ(u):C⁻¹:σ(v) dΩ`. |
| `∫∫σᵢⱼσₖₗdxdy_PlaneStrian` | `set𝝭!` | — | Element-set overload of `∫∫σᵢⱼσₖₗdxdy_PlaneStrian`: the same integrand, summed over a vector of elements. See the element-level method for the formula. |
| `∫∫τ∇q∇pdxdy` | `set∇𝝭!` | — | Stabilisation `τ ∫ ∇q·∇p dΩ` with the local time/streamline scale `τ`. |
| `∫∫τ∇sᵢⱼ∇pdxdy` | `set∇𝝭!` | — | Stress–pressure stabilisation coupling `τ ∫ ∇s:∇p dΩ`. |
| `∫∫τ∇sᵢⱼ∇sᵢₖdxdy` | `set∇𝝭!` | — | Stabilisation of the deviatoric stress `τ ∫ ∇s:∇s dΩ`. |
| `∫∫τ∇σᵢⱼ∇σᵢₖdxdy` | `set∇𝝭!` | — | Stabilisation of the full stress `τ ∫ ∇σ:∇σ dΩ`. |
| `∫∫∇puᵢdxdy` | `set∇𝝭!`, `set𝝭!` | — | Transpose pressure coupling `∫ ∇p·u dΩ`. |
| `∫∫∇sᵢⱼuᵢdxdy` | `set∇𝝭!`, `set𝝭!` | — | Divergence of the deviatoric stress `∫ (∇·s)·u dΩ`. |
| `∫∫∇σᵢⱼuᵢdxdy` | `set∇𝝭!`, `set𝝭!` | — | Divergence of the stress against the displacement, `∫ (∇·σ)·v dΩ`. |

## Hamiltonian / space–time rod

`ApproxOperator.Hamilton`

| operator | requires | material constants | description |
|---|---|---|---|
| `Hₑ` | `set∇𝝭!` | `EA`, `c²`, `ρA` | Relative energy-norm error of the 2-D space–time rod. *(3 overloads)* |
| `stabilization_bar_LSG` | `set∇²𝝭!` | `ρA`, `α`, `EA` | Least-squares Galerkin (LSG) stabilisation for the space–time rod. |
| `stabilization_bar_LSG_Γ` | `set∇𝝭!` | — | LSG stabilisation applied on the boundary of the space–time rod (the `_Γ` counterpart of `stabilization_bar_LSG`). |
| `truncation_error` | `set∇𝝭!` | `c` | Truncation-error measures of the space–time discretisation. *(2 overloads)* |
| `∫ppdΩ` | `set𝝭!` | `ρA` | Mass matrix of the pressure variable, `∫ ρA p p dΩ`. |
| `∫pudΩ` | `set∇𝝭!`, `set𝝭!` | — | Element-set overload of `∫pudΩ`: the same integrand, summed over a vector of elements. See the element-level method for the formula. |
| `∫qkpdΩ` | `set𝝭!` | `EA` | Axial stiffness term of a rod, `∫ EA q p dΩ`. |
| `∫qmpdΩ` | `set𝝭!` | `ρA` | Mass-like term of a rod, `∫ ρA q p dΩ`. |
| `∫q∇𝑛pds` | `set∇𝝭!`, `set𝝭!` | `ρA`, `EA` | Space–time boundary term `∫ q ∂ₙp ds` with the wave coefficients. |
| `∫uudΩ` | `set∇𝝭!` | `EA` | Stiffness matrix of the displacement variable, `∫ EA u u dΩ`. |
| `∫∫q̇mṗqkpdxdt` | `set∇𝝭!`, `set𝝭!` | `m`, `kᶜ` | Space–time operator `∫∫ (m q̇ ṗ - k q p) dx dt`. |
| `∫∫αqṗdxdt` | `set∇𝝭!`, `set𝝭!` | `α` | Space–time stabilisation term `α ∫∫ q ṗ dx dt`. |
| `∫∫∇q∇pdxdt` | `set∇𝝭!` | `EA`, `ρA` | Space–time wave operator `∫∫ (ρA q,ₜ p,ₜ - EA q,ₓ p,ₓ) dx dt`. *(2 overloads)* |
| `∫∫∇q∇pdΩdt` | `set∇𝝭!` | `c²` | Space–time operator `∫∫ (q,ₜ p,ₜ - c² ∇q·∇p) dΩ dt`. |

## Heat conduction and convection

`ApproxOperator.Heat`

| operator | requires | material constants | description |
|---|---|---|---|
| `H₁` | `set∇𝝭!`, `set𝝭!` | — | Relative H¹ error of the corresponding field against an exact solution, `‖e‖/‖u‖`. The exact field is supplied with `prescribe!`. *(2 overloads)* |
| `L₂` | `set𝝭!` | — | Relative L₂ error of the corresponding field against an exact solution, `‖e‖/‖u‖`. The exact field is supplied with `prescribe!`. *(2 overloads)* |
| `L₂𝒑` | `set𝝭!` | — | Relative L₂ error of the corresponding field against an exact solution, `‖e‖/‖u‖`. The exact field is supplied with `prescribe!`. *(2 overloads)* |
| `g` | — | — | Assemble the prescribed boundary values into the right-hand side, `f[I] += Nᵢ g w`. The nodal field read is selected by `dof` (default `:d`). |
| `∫pᵢnᵢgds` | `set𝝭!` | — | Boundary flux load `∫ (p·n) g ds`. |
| `∫pᵢnᵢuds` | `set𝝭!` | — | Boundary flux coupling `∫ (p·n) u ds`. |
| `∫uds` | — | — | Element-set overload of `∫uds`: the same integrand, summed over a vector of elements. See the element-level method for the formula. |
| `∫vbdΩ` | `set𝝭!` | — | Volumetric heat source load, `f[I] += ∫ b Nᵢ dΩ`. |
| `∫vgdΓ` | `set𝝭!` | `α` | Penalty enforcement of a Dirichlet condition, `α ∫ (u-g) v dΓ`. |
| `∫vtdΓ` | `set𝝭!` | — | Boundary heat flux load, `f[I] += ∫ t Nᵢ dΓ`. |
| `∫vᵢnᵢgds` | `set∇𝝭!` | `kᶜ`, `f` | Boundary flux load, `k ∫ (∂ₙu) g ds` evaluated into `f`. |
| `∫vᵢnᵢuds` | `set∇𝝭!`, `set𝝭!` | `kᶜ` | Two-set boundary flux coupling, `k ∫ (∂ₙu) v ds`, with test functions on `a₁` and trial functions on `a₂`. |
| `∫λgdΓ` | `set𝝭!` | — | Lagrange-multiplier boundary term `∫ λ g dΓ`. The single-element overload takes `λ` from the same set as the boundary data. *(2 overloads)* |
| `∫λudΓ` | `set𝝭!` | — | Lagrange-multiplier coupling `∫ λ u dΓ`. |
| `∫∇v∇udx` | `set∇𝝭!` | `kᶜ` | One-dimensional conduction stiffness, `∫ k ∇v ∇u dx`. |
| `∫∇v∇udΩ` | `set∇𝝭!` | `kᶜ` | Three-dimensional conduction stiffness, `∫ k ∇v·∇u dΩ`. |
| `∫∇𝑛vgds` | `set∇𝝭!`, `set𝝭!` | `α`, `kᶜ` | Nitsche-type weak boundary condition for convection–diffusion, `∫ (k ∂ₙu - α(u-g)) v ds`. |
| `∫∇𝑛vgdΓ` | `set∇𝝭!`, `set𝝭!` | — | Boundary normal-derivative load `∫ (∂ₙu) g dΓ`. |
| `∫∫qᵢpᵢdxdy` | `set𝝭!` | — | Mass matrix of the heat-flux variable, `∫ q·p dΩ`. |
| `∫∫∇v∇udxdy` | `set∇𝝭!` | `kᶜ` | Two-dimensional conduction stiffness, `∫ k ∇v·∇u dΩ`. |
| `∫∫∇𝒑bdxdy` | `set∇𝝭!` | — | Source coupling for the mixed heat formulation, `∫ (∇·p) b dΩ`. |
| `∫∫∇𝒑udxdy` | `set∇𝝭!`, `set𝝭!` | — | Transpose flux coupling `∫ (∇·p) u dΩ`. |
| `∫∫𝒑∇udxdy` | `set∇𝝭!`, `set𝝭!` | — | Flux–gradient coupling `∫ p·∇u dΩ`. |

## Hyperelasticity

`ApproxOperator.Hyperelasticity`

| operator | requires | material constants | description |
|---|---|---|---|
| `Δ∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean` | `set∇𝝭!` | `E`, `ν` | Deviatoric part of the tangent of a 2-D Neo-Hookean material. |
| `Δ∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean2` | `set∇𝝭!` | `E`, `ν` | Deviatoric part of the tangent of a 2-D Neo-Hookean material (the `_2` variant). |
| `Δ∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean` | `set∇𝝭!` | `E`, `ν` | Volumetric part of the tangent of a 2-D Neo-Hookean material. |
| `Δ∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean2` | `set∇𝝭!` | `E`, `ν` | Volumetric part of the tangent of a 2-D Neo-Hookean material (the `_2` variant). |
| `Δ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean` | `set∇𝝭!` | `E`, `ν` | Full part of the tangent of a 2-D Neo-Hookean material. |
| `Δ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean2` | `set∇𝝭!` | `E`, `ν` | Full part of the tangent of a 2-D Neo-Hookean material (the `_2` variant). |
| `Δ∫∫EᵢⱼSᵢⱼdxdy_SaintVenantKirchhoff` | `set∇𝝭!` | `E`, `ν` | Incremental tangent stiffness of a 2-D Saint-Venant–Kirchhoff material. |
| `∫ESdx` | `set∇𝝭!` | `E` | 1-D hyperelastic bar: tangent stiffness and internal force from the Saint-Venant–Kirchhoff stored energy. |
| `∫vᵢuᵢds` | `set𝝭!` | — | Element-set overload of `∫vᵢuᵢds`: the same integrand, summed over a vector of elements. See the element-level method for the formula. |
| `∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean` | `set∇𝝭!` | `E`, `ν` | Deviatoric part of the internal force of a 2-D Neo-Hookean material. |
| `∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean2` | `set∇𝝭!` | `E`, `ν` | Deviatoric part of the internal force of a 2-D Neo-Hookean material (the `_2` variant). |
| `∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean` | `set∇𝝭!` | `E`, `ν` | Volumetric part of the internal force of a 2-D Neo-Hookean material. |
| `∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean2` | `set∇𝝭!` | `E`, `ν` | Volumetric part of the internal force of a 2-D Neo-Hookean material (the `_2` variant). |
| `∫∫EᵢⱼSᵢⱼdxdy_NeoHookean` | `set∇𝝭!` | `E`, `ν` | Full part of the internal force of a 2-D Neo-Hookean material. |
| `∫∫EᵢⱼSᵢⱼdxdy_NeoHookean2` | `set∇𝝭!` | `E`, `ν` | Full part of the internal force of a 2-D Neo-Hookean material (the `_2` variant). |
| `∫∫EᵢⱼSᵢⱼdxdy_SaintVenantKirchhoff` | `set∇𝝭!` | `E`, `ν` | Incremental internal force of a 2-D Saint-Venant–Kirchhoff material. |

## Kirchhoff plate

`ApproxOperator.KirchhoffPlate`

| operator | requires | material constants | description |
|---|---|---|---|
| `wΔMₙₛ` | `set𝝭!` | — | Corner-force load of the Kirchhoff plate, `∫ w ΔMₙₛ`. |
| `ΔM̃ₙₛg` | `set∇²𝝭!`, `set𝝭!` | — | Element-set overload of `ΔM̃ₙₛg`: the same integrand, summed over a vector of elements. See the element-level method for the formula. |
| `ΔMₙₛg` | `set∇²𝝭!`, `set𝝭!` | `D`, `ν` | Corner-force terms of the Kirchhoff plate, `∫ ΔMₙₛ ḡ`. |
| `∫M̃ₙₙθdΓ` | `set∇²𝝭!`, `set∇𝝭!` | — | Element-set overload of `∫M̃ₙₙθdΓ`: the same integrand, summed over a vector of elements. See the element-level method for the formula. |
| `∫MₙₙθdΓ` | `set∇²𝝭!`, `set∇𝝭!` | `D`, `ν` | Kirchhoff plate boundary term `∫ (Mₙₙ(w) + α θ) θ̄ dΓ` — normal bending moment plus penalty, as used in the weak (Nitsche/penalty) free-edge condition. |
| `∫VgdΓ` | `set∇³𝝭!`, `set𝝭!` | `D`, `ν` | Kirchhoff plate Kirchhoff-shear boundary term, `∫ (Vₙ(w) + α w) ḡ dΓ`. |
| `∫ṼgdΓ` | `set𝝭!` | — | Element-set overload of `∫ṼgdΓ`: the same integrand, summed over a vector of elements. See the element-level method for the formula. |
| `∫wVdΓ` | `set𝝭!` | — | Kirchhoff shear boundary load `∫ V Nᵢ dΓ`. |
| `∫wqdΩ` | `set𝝭!` | — | Transverse load vector, `f[I] += ∫ q Nᵢ dΩ`. |
| `∫θₙMₙₙdΓ` | `set∇𝝭!` | — | Rotation–moment boundary coupling `∫ θₙ Mₙₙ dΓ`. |
| `∫κ̃ᵢⱼM̃ᵢⱼdΩ` | `set∇²𝝭!` | `D`, `ν` | Kirchhoff plate bending stiffness with the **symmetric** mixed-derivative treatment. |
| `∫κᵢⱼMᵢⱼdΩ` | `set∇²𝝭!` | `D`, `ν` | Kirchhoff plate bending stiffness, `∫ κᵢⱼ(w) Mᵢⱼ(v) dΩ`. |
| `∫∇wMdΓ` | `set∇𝝭!` | — | Deflection-gradient boundary term `∫ (∂ₙw) M dΓ`. |
| `∫∇𝑛vθdΓ` | `set∇𝝭!` | — | Boundary penalty on the normal derivative, `α ∫ (∂ₙw) θ dΓ`. |

## Kirchhoff shell

`ApproxOperator.KirchhoffShell`

| operator | requires | material constants | description |
|---|---|---|---|
| `Hₑ_ThinShell` | `set∇²𝝭!`, `set∇𝝭!`, `set𝝭!` | — | Relative energy-norm error of a Kirchhoff shell solution, `‖e‖_E / ‖u‖_E` over membrane plus bending energy. *(2 overloads)* |
| `ScordelisLoRoof_𝐴` | `set𝝭!` | — | Geometry helper for the Scordelis–Lo roof benchmark: the physical coordinates of a quadrature point on the shell surface. |
| `SphericalShell_𝐴` | `set𝝭!` | — | Geometry helper for the pinched spherical shell benchmark: physical coordinates on the sphere. |
| `SphericalShell_𝐵` | `set𝝭!` | — | Geometry helper for the pinched spherical shell benchmark: the second coordinate of the paramaterisation. |
| `ΔMₙₛ𝒂₃𝒈_HR` | `set𝝭!` | — | Overload or Nitsche/Hellinger–Reissner sibling of the same shell term; see the documented operators with a similar name. |
| `ΔMₙₛ𝒂₃𝒈_Nitsche` | `set∇²𝝭!`, `set∇𝝭!`, `set𝝭!` | — | Overload or Nitsche/Hellinger–Reissner sibling of the same shell term; see the documented operators with a similar name. |
| `ΔMₙₛ𝒂₃𝒗` | `set𝝭!` | — | Overload or Nitsche/Hellinger–Reissner sibling of the same shell term; see the documented operators with a similar name. |
| `∫MₙₙθₙdΓ_HR` | `set∇𝝭!`, `set𝝭!` | — | Overload or Nitsche/Hellinger–Reissner sibling of the same shell term; see the documented operators with a similar name. |
| `∫MₙₙθₙdΓ_Nitsche` | `set∇²𝝭!`, `set∇𝝭!`, `set𝝭!` | — | Overload or Nitsche/Hellinger–Reissner sibling of the same shell term; see the documented operators with a similar name. |
| `∫δMκdΩ` | `set𝝭!` | — | Bending coupling `∫ δM : κ dΩ`. |
| `∫δNεdΩ` | `set𝝭!` | — | Membrane coupling `∫ δN : ε dΩ`, the stress–strain dual of [`∫δεCεdΩ`](@ref). |
| `∫δεCεdΩ` | `set𝝭!` | — | Membrane energy in strain form, `∫ δε : C ε dΩ`. |
| `∫δθθdΓ` | `set∇𝝭!` | — | Boundary term coupling the rotation trace, `∫ δθ θ dΓ`. |
| `∫δκCκdΩ` | `set𝝭!` | — | Bending energy in curvature form, `∫ δκ : C κ dΩ`. |
| `∫εᵢⱼNᵢⱼdΩ` | `set∇𝝭!` | — | Membrane (stretching) part of the Kirchhoff shell, `∫ εᵢⱼ(u) Nᵢⱼ(v) dΩ` with `N = Dᵐ C ε`, `Dᵐ = Eh/(1-ν²)`. |
| `∫εᵢⱼNᵢⱼκᵢⱼMᵢⱼdΩ` | `set∇²𝝭!`, `set∇𝝭!` | — | Kirchhoff shell stiffness, membrane plus bending, `∫ (εᵢⱼNᵢⱼ + κᵢⱼMᵢⱼ) dΩ`. |
| `∫κᵢⱼMᵢⱼdΩ` | `set∇²𝝭!`, `set∇𝝭!` | — | Bending part of the Kirchhoff shell, `∫ κᵢⱼ(θ) Mᵢⱼ(θ̄) dΩ` with `M = Dᵇ C κ`, `Dᵇ = Eh³/(12(1-ν²))`. |
| `∫∇𝑴∇𝒂₃𝒗dΩ` | `set∇²𝝭!`, `set∇𝝭!`, `set𝝭!` | — | Bending-moment/deflection coupling in gradient form. |
| `∫∇𝑴𝒏𝒂₃𝒈dΓ_HR` | `set∇𝝭!`, `set𝝭!` | — | Overload or Nitsche/Hellinger–Reissner sibling of the same shell term; see the documented operators with a similar name. |
| `∫∇𝑴𝒏𝒂₃𝒈dΓ_Nitsche` | `set∇²𝝭!`, `set∇³𝝭!`, `set∇𝝭!`, `set𝝭!` | — | Overload or Nitsche/Hellinger–Reissner sibling of the same shell term; see the documented operators with a similar name. |
| `∫∇𝑴𝒏𝒂₃𝒗dΓ` | `set∇𝝭!`, `set𝝭!` | — | Overload or Nitsche/Hellinger–Reissner sibling of the same shell term; see the documented operators with a similar name. |
| `∫∇𝑵𝒗dΩ` | `set∇𝝭!`, `set𝝭!` | — | Divergence coupling `∫ ∇·𝑵 v dΩ` for the shell. |
| `∫𝑴ₙₙ𝜽ₙdΓ` | `set∇𝝭!`, `set𝝭!` | — | Overload or Nitsche/Hellinger–Reissner sibling of the same shell term; see the documented operators with a similar name. |
| `∫𝒂ᵢvᵢgᵢdΓ` | `set𝝭!` | — | Overload or Nitsche/Hellinger–Reissner sibling of the same shell term; see the documented operators with a similar name. |
| `∫𝒏𝑵𝒈dΓ_HR` | `set𝝭!` | — | Overload or Nitsche/Hellinger–Reissner sibling of the same shell term; see the documented operators with a similar name. |
| `∫𝒏𝑵𝒈dΓ_Nitsche` | `set∇𝝭!`, `set𝝭!` | — | Nitsche weak boundary condition on the shell, `∫ 𝒏·𝑵 𝑔 dΓ`. |
| `∫𝒏𝑵𝒗dΓ` | `set𝝭!` | — | Overload or Nitsche/Hellinger–Reissner sibling of the same shell term; see the documented operators with a similar name. |

## Mindlin (thick) plate

`ApproxOperator.MindlinPlate`

| operator | requires | material constants | description |
|---|---|---|---|
| `H₁` | — | — | Relative H¹ error of the mixed `(w, φ)` solution. |
| `H₁w` | `set∇𝝭!`, `set𝝭!` | — | Relative H¹ error of the deflection (L₂ plus the gradient part). |
| `H₁φ` | `set∇𝝭!`, `set𝝭!` | — | Relative H¹ error of the rotation field, combining the L₂ and gradient parts. |
| `Hₑ` | `set∇𝝭!`, `set𝝭!` | `E`, `h`, `ν` | Relative energy-norm error of a Mindlin plate solution against an exact field: *(4 overloads)* |
| `L₂Q` | `set𝝭!` | — | Relative L₂ error of the shear force, `‖Q-Q̄‖/‖Q̄‖`. *(2 overloads)* |
| `L₂w` | `set𝝭!` | — | Relative L₂ error of the deflection. *(2 overloads)* |
| `L₂γ` | `set∇𝝭!`, `set𝝭!` | — | Relative L₂ error of the transverse shear strain, `‖γ-γ̄‖/‖γ̄‖`. *(4 overloads)* |
| `L₂φ` | `set𝝭!` | — | Relative L₂ error of the rotation field. *(2 overloads)* |
| `∫MMdΩ` | `set𝝭!` | `E`, `ν`, `h` | Bending-moment mass term `∫ M·M dΩ`, the dual counterpart of [`∫QQdΩ`](@ref) in mixed plate formulations. |
| `∫MφdΓ` | `set𝝭!` | — | Bending moment on the boundary, `∫ (M·n) φ dΓ`. *(2 overloads)* |
| `∫M∇φdΩ` | `set∇𝝭!`, `set𝝭!` | — | Mixed bending moment coupling, `∫ M:∇φ dΩ`. |
| `∫QQdΩ` | `set𝝭!` | `E`, `ν`, `h` | Shear-force mass term `∫ Q·Q dΩ`, used in mixed and stabilised Mindlin formulations. |
| `∫QwdΓ` | `set𝝭!` | — | Mixed shear force on the boundary, `∫ (Q·n) w dΓ`; the overloads differ in whether the boundary value and the two element sets are supplied. *(4 overloads)* |
| `∫QφdΩ` | `set𝝭!` | — | Mixed shear–rotation coupling, `∫ Q·φ dΩ`. *(2 overloads)* |
| `∫Q∇wdΩ` | `set∇𝝭!`, `set𝝭!` | — | Mixed shear force coupling, `∫ Q·∇w dΩ`. The two-element-set overload takes the shear and deflection approximations separately. *(2 overloads)* |
| `∫wVdΓ` | `set𝝭!` | — | Kirchhoff shear boundary load, `f[I] += ∫ V Nᵢ dΓ`. |
| `∫wqdΩ` | `set𝝭!` | — | Transverse load vector, `f[I] += ∫ q Nᵢ dΩ`. The load `q` is prescribed with `prescribe!(el, :q => f)`. |
| `∫wwdΩ` | `set𝝭!` | — | Bending stiffness of the Mindlin plate in deflection form, `∫ κ(δw)·κ(w) dΩ` with `κ = ∇²w`. |
| `∫αwwdΓ` | `set𝝭!` | `α` | Penalty term on the deflection, `α ∫ w w̄ dΓ`, for a weak simply-supported or clamped edge. *(2 overloads)* |
| `∫αφφdΓ` | `set𝝭!` | `α` | Penalty term on the rotation, `α ∫ φ φ̄ dΓ`, for a weakly clamped edge. *(2 overloads)* |
| `∫κMγQdΩ` | `set∇𝝭!`, `set𝝭!` | `E`, `ν`, `h` | Reissner–Mindlin plate bilinear form, `∫ (κ(w,φ)·M(κ) + γ(w,φ)·Q(γ)) dΩ`. |
| `∫κκdΩ` | `set∇𝝭!` | `E`, `ν`, `h` | Bending-only block `∫ κ(δφ)·κ(φ) dΩ` of the Mindlin plate, carrying the bending stiffness. |
| `∫ρwwdΩ` | `set𝝭!` | `ρ`, `h` | Deflection mass matrix of the Mindlin plate, `∫ ρ w w / h² dΩ`. |
| `∫ρφφdΩ` | `set𝝭!` | `ρ` | Rotary-inertia mass matrix, `∫ ρ φ·φ dΩ`. |
| `∫φMdΓ` | `set𝝭!` | — | Boundary moment load, `f[I] += ∫ M Nᵢ dΓ`. |
| `∫φmdΩ` | `set𝝭!` | — | Distributed couple load, `f[I] += ∫ m Nᵢ dΩ`. |
| `∫φwdΩ` | `set∇𝝭!`, `set𝝭!` | `E`, `h`, `ν` | Rotation–deflection coupling `∫ φ·∇w dΩ` of the Mindlin plate. *(2 overloads)* |
| `∫φwdΩ_DSG` | `set∇𝝭!` | `E`, `ν`, `h` | Rotation–deflection coupling of the **DSG** triangle. |
| `∫φwdΩ_MITC` | — | `E`, `ν`, `h` | Rotation–deflection shear coupling of the **MITC4** element (mid-side sampling). |
| `∫φφdΩ` | `set𝝭!` | `E`, `ν`, `h` | Rotation–rotation shear block `∫ φ·φ dΩ` of the Mindlin plate. |
| `∫φφdΩ_DSG` | — | `E`, `ν`, `h` | Rotation–rotation block of the **DSG** (Discrete Shear Gap) triangle. The shear gap is integrated exactly on the element boundary, which avoids shear locking on `Tri3`. |
| `∫φφdΩ_MITC` | — | `E`, `ν`, `h` | Rotation–rotation shear block of the **MITC4** element (mid-side sampling). |
| `∫∇MφdΩ` | `set∇𝝭!`, `set𝝭!` | — | Mixed bending moment coupling in gradient form, `∫ ∇·M · φ dΩ`. |
| `∫∇QwdΩ` | `set∇𝝭!`, `set𝝭!` | — | Mixed shear force coupling in gradient form, `∫ ∇Q · w dΩ`. The transpose of [`∫Q∇wdΩ`](@ref) for the rotated variable ordering. *(2 overloads)* |
| `∫∇wσ∇wdΩ` | `set∇𝝭!` | `h` | Deflection-gradient block of a `u/p` or `u/σ` mixed formulation, `∫ ∇w·σ·∇w dΩ`, where `σ` comes from its own element set. |
| `∫∇w∇wdΩ` | `set∇𝝭!` | `E`, `ν`, `h` | Shear block `∫ ∇w·∇w dΩ` of the Mindlin plate. Together with [`∫φφdΩ`](@ref) and [`∫φwdΩ`](@ref) it forms the shear contribution. |
| `∫∇w∇wdΩ_MITC` | — | `E`, `ν`, `h` | Shear block of the **MITC4** element: the transverse shear is sampled at the mid-side points instead of the Gauss points, which removes shear locking. |
| `∫∇φσ∇φdΩ` | `set∇𝝭!` | `h` | Rotation-gradient block of a mixed formulation, `∫ ∇φ·σ·∇φ dΩ`. |

## Phase-field fracture

`ApproxOperator.Phasefield`

| operator | requires | material constants | description |
|---|---|---|---|
| `CRACK_NORMAL` | — | — | Compute the crack normal at each quadrature point from the damage gradient, `n = -∇v/\|∇v\|`, and write it back into the quadrature point. |
| `UPDATE_PFM_1D` | `set∇𝝭!`, `set𝝭!` | — | Advance the 1-D phase-field history variable: `ℋ ← max(ℋ, (ε-1)²)` at every quadrature point. Stateful — call once per load step before assembling. |
| `UPDATE_PFM_2D` | `set∇𝝭!`, `set𝝭!` | — | Advance the 2-D phase-field history variable from the strain energy, `ℋ ← max(ℋ, ψ⁺(ε))`. Stateful — call once per load step. |
| `∫v²uₓuₓdx` | `set∇𝝭!`, `set𝝭!` | — | Element-set overload of `∫v²uₓuₓdx`: the same integrand, summed over a vector of elements. See the element-level method for the formula. |
| `∫vₓvₓvvdx` | `set∇𝝭!`, `set𝝭!` | — | Element-set overload of `∫vₓvₓvvdx`: the same integrand, summed over a vector of elements. See the element-level method for the formula. *(2 overloads)* |
| `∫vₓvₓvvdx_hard_device` | `set∇𝝭!`, `set𝝭!` | `l` | 1-D phase-field damage equation, AT2 model. |
| `∫∫∇v∇vvvdxdy` | `set∇𝝭!`, `set𝝭!` | `l` | 2-D phase-field damage equation, AT2 model. |

## Plasticity and contact

`ApproxOperator.Plasticity`

| operator | requires | material constants | description |
|---|---|---|---|
| `∫vᵢσdΩ_frictional_contact` | `set∇𝝭!`, `set𝝭!` | `E`, `ν` | Coulomb frictional contact with a phase-field damage coupling. |
| `∫vᵢσdΩ_frictional_contact_2` | `set∇𝝭!`, `set𝝭!` | — | frictional-contact2 |
| `∫vᵢσdΩ_mc_phasefield` | `set∇𝝭!`, `set𝝭!` | `E`, `ν` | Mohr–Coulomb plasticity coupled to a phase field. The friction angle is fixed at `𝜙 = π/3` inside the operator, and the damage enters the tangent through the `(1-v²+η)` degradation factor, with `η` a small residual-stiffness parameter. |
| `∫vᵢσdΩ_mohr_coulomb` | `set∇𝝭!` | `λ`, `μ` | 2-D Mohr–Coulomb elastoplastic tangent and internal force. |
| `∫vᵢσdΩ_tresca` | `set∇𝝭!` | `λ`, `μ`, `c`, `c` | Tresca elastoplastic tangent and internal force. |
| `∫vₓσdx` | `set∇𝝭!` | `E`, `K`, `σy` | 1-D elastoplastic tangent stiffness and internal force, J2 with linear isotropic hardening. The quadrature-point state `ξ.σₙ`, `ξ.εᵖₙ`, `ξ.αₙ` is updated in place, so the operator is **stateful**. |

## Stokes flow

`ApproxOperator.Stokes`

| operator | requires | material constants | description |
|---|---|---|---|
| `∫∫μ∇u∇vdxdy` | `set∇𝝭!` | `μ` | Viscous term of the Stokes problem, `2μ ∫ ∇u:∇v dΩ`. |

## Consistency checks

`ApproxOperator.Test`

| operator | requires | material constants | description |
|---|---|---|---|
| `cc∇𝝭` | `set∇𝝭!` | — | Consistency check of the shape-function gradients: relative L₂ difference between the reproduced `∇u` and the prescribed `ξ.∂u∂x`. *(2 overloads)* |
| `cc𝝭` | `set𝝭!` | — | Consistency check of the shape functions: relative L₂ difference between the reproduced field `Σᵢ Nᵢ dᵢ` and the prescribed `ξ.u`. *(2 overloads)* |

## Space–time weighted residual

`ApproxOperator.WeightedResidual`

| operator | requires | material constants | description |
|---|---|---|---|
| `∫NṄdxdt` | `set∇𝝭!`, `set𝝭!` | — | Space–time weighted-residual operator `∫∫ Nᵢ ∂Nⱼ/∂t dx dt` with unit coefficient. |
| `∫c²B₁B₁dxdt` | `set∇𝝭!` | `c` | Space–time weighted-residual stiffness, `c² ∫∫ (∂ₓNᵢ)(∂ₓNⱼ) dx dt`. |
| `∫kNNdxdt` | `set𝝭!` | `kᶜ` | Space–time weighted-residual mass-like term `-k ∫∫ Nᵢ Nⱼ dx dt`. |
| `∫kNṄdxdt` | `set∇𝝭!`, `set𝝭!` | `kᶜ` | Space–time weighted-residual operator `k ∫∫ Nᵢ ∂Nⱼ/∂t dx dt`. |

