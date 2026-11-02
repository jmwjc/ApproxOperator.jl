# Optimisation Progress — ApproxOperator.jl

**Priority order（优先级）:** 可读性 > 使用便利性 > 计算效率
**Last updated:** S0 – S6 complete — `242/242` tests green

---

## 1. Status overview（总览）

| Phase | Scope | Status | Tests |
|---|---|---|---|
| **S0** | Test infrastructure | ✅ Done | — |
| **S1** | Mathematical correctness | ✅ Done | 99 pass |
| **S1.5** | Retire `error_estimates.jl`, redistribute its operators | ✅ Done | — |
| **S2** | Material parameters → explicit keyword arguments + docstrings | ✅ Done | 103 pass |
| **S3** | `importmsh.jl` rewrite (`ElementSpec`) | ✅ Done | 125 pass |
| **S3.7** | `addCentroidNodes!` + `get𝑿ᵢ(centroid = …)` option | ✅ Done | 193 pass |
| **S4** | Docstrings + `docs/operators.md` | ✅ Done | 199 pass |
| **S5** | `README` + `docs/` + dead code revival | ✅ Done | 227 pass |
| **S6** | `Project.toml` compat, `Aqua`, ambiguity and export fixes | ✅ Done | 242 pass |

---

## 2. S0 — Test infrastructure ✅

| Item | Change | File |
|---|---|---|
| Dead tests removed | `Timoshenko` module does not exist anywhere in `src/`; both files referenced `ApproxOperator.Timoshenko` and could never run | `test/timoshenko_static_ss.jl` ❌ deleted<br>`test/timoshenko_geometric_stiffness.jl` ❌ deleted |
| Patch test extracted | Became a reusable end-to-end example instead of being buried inside a `@testset` | `test/patchtest.jl` ✨ new |
| Mesh path fixed | `gmsh.open("patchtest.msh")` → `joinpath(@__DIR__, "patchtest.msh")`; no longer CWD-dependent | `test/patchtest.jl` |
| All tests wired up | `runtests.jl` now `include`s every test file instead of containing them inline | `test/runtests.jl` |

**Rationale:** a test file that is not `include`d is a test that silently rots. Two of the three original test files were in that state.

---

## 3. S1 — Mathematical correctness ✅

### 3.1 The `op` global — 138 dangling references

Seven operator modules read material constants from a global `op` that was **deleted** from the package (it used to be a `struct Operator`, removed in an earlier commit). Every one of those operators raised `UndefVarError: op not defined` on first call.

| File | `op.` references at HEAD | Now | Bound via |
|---|---:|---:|---|
| `thin_plate.jl` | 20 | **0** | `ξ.D`, `ξ.ν`, `ξ.α`, `𝓖[1].D`, `𝓖[1].ν` |
| `thin_shell.jl` | 33 | **0** | `ξ.E`, `ξ.ν`, `ξ.h` |
| `hyperelasticity.jl` | 30 | **0** | `ξ.E`, `ξ.ν`, `ap.E`, `ap.ν` |
| `curved_beam.jl` | 11 | **0** | `ξ.EI`, `ξ.EA`, `ξ.kGA`, `ξ.R`, `ξ.α`, `ξ.a₁`, `ξ.a₃` |
| `error_estimates.jl` | 9 | **0** | `apu.E`, `apu.ν`, `ap.E`, `ap.ν` |
| `phasefield.jl` | 11 | **0** | `ξ.k`, `ξ.l`, `ξ.η` |
| `plasticity.jl` | 24 | **0** | `ξ.E`, `ξ.ν`, `ξ.λ`, `ξ.μ`, `ξ.𝜙`, `ξ.σy`, `ap.E`, `ap.ν` |
| **Total** | **138** | **0** | |

> Note: `plasticity.jl` is **24**, not 23. Line 143 held `op.𝜙` where `𝜙` is
> U+1D719 (Mathematical Italic Small Phi), which is outside the Greek block and
> was missed by an early regex. Verified by byte inspection.

**Convention after the change:** identical to `elasticity.jl`, `stokes.jl` and `thick_plate.jl` — the constant is read from the quadrature point and supplied by the caller:

```julia
elements = getElements(nodes, entities["Ω"], 2)
prescribe!(elements, :D => D, :ν => ν)
(∫κᵢⱼMᵢⱼdΩ => elements)(k)
```

**Do not revive `struct Operator{T}`.** It would create a second, parallel API
next to the existing `Pair{Function,Vector}` assembler and make the package
*less* readable, which violates the stated priority.

### 3.2 Bugs fixed

| # | Bug | File | Verification |
|---|---|---|---|
| 1 | `Seg3` gradient **multiplied** by `J` instead of dividing | `src/approximation/seg3.jl` | see §3.3 |
| 2 | `PiecewisePolynomial{:Constant}`: `∂𝝭∂y = 𝒙[:∂𝝭∂x]` (copy-paste) | `src/approximation/piecewise.jl` | see §3.4 |
| 3 | `SegHermite` `set∇𝝭!` referenced an undefined `ap` | `src/approximation/seghermite.jl` | parses + loads |
| 4 | `Tri6` stale `# FIXME: The shape function order is wrong!` | `src/approximation/tri6.jl` | see §3.5 |
| 5 | Statements naming `E`/`ν`/`h` **before** any binding, exposed once `op.` was removed (46 sites) | `hyperelasticity.jl`, `error_estimates.jl`, `plasticity.jl`, `thin_plate.jl`, `thin_shell.jl` | guarded by a new static test |

> **Bug 5 is the important one.** The first `op.` removal pass deleted the alias
> lines but did not re-bind the names they fed. A naive "count `op.` = 0" check
> reported success while `thin_shell.jl` still had 33 unbound identifiers.
> Detected by scanning each function body for pre-loop uses of `E`/`ν`/`h`
> with no binding, then fixed and locked with a regression test.

### 3.3 `Seg3` — `∂N/∂x = ∂N/∂ξ / J`

Segment of length `L = 3` mapped to `ξ ∈ [-1,1]`, so `det J = L/2 = 1.5`.
For `u(x) = x` the nodal values are the node coordinates, so
`Σᵢ (∂Nᵢ/∂x)·uᵢ` must equal `1` everywhere.

```
 xi      du/dx    expected
 -0.75   1.0      1.0
 -0.50   1.0      1.0
  0.00   1.0      1.0
  0.50   1.0      1.0
  0.75   1.0      1.0
max |du/dx − 1| = 0.0        PASS
```

Multiplying by `J` (the old code) yields `L² = 9`.

**Reverse-verified:** restoring `*𝐽` makes the consistency test fail with
`1.1555 < 1e-12`, so the test genuinely guards this bug.

### 3.4 `PiecewisePolynomial{:Constant}` — `∂𝝭∂y` independence

```julia
-    ∂𝝭∂y = 𝒙[:∂𝝭∂x]      # HEAD: aliases the ∂𝝭∂x array
+    ∂𝝭∂y = 𝒙[:∂𝝭∂y]      # fixed
```

Runtime proof — poison `∂𝝭∂x` with 42.0 and read `∂𝝭∂y`:

```
after writing ∂𝝭∂x = 42, ∂𝝭∂y reads 0.0     PASS (independent arrays)
```

### 3.5 `Tri6` node ordering — FIXME was stale

The reference triangle is `(0,0), (1,0), (0,1)` with midside nodes on
`(1,2), (2,3), (3,1)` — exactly what `Tri3toTri6!` produces. Evaluating all
six shape functions at each of the six reference nodes gives the identity
matrix:

```
pt 1  N = [1, 0, 0, 0, 0, 0]
pt 2  N = [0, 1, 0, 0, 0, 0]
pt 3  N = [0, 0, 1, 0, 0, 0]
pt 4  N = [0, 0, 0, 1, 0, 0]
pt 5  N = [0, 0, 0, 0, 1, 0]
pt 6  N = [0, 0, 0, 0, 0, 1]
```

### 3.6 `:𝐽` semantics documented

`:𝐽` is `det(∂𝒙/∂𝝃)`. A 50-line block at the top of
`src/preprocession/importmsh.jl` now records:

* the data-dictionary slot layout (`slot 1` = per-node geometry for `𝑿ᵢ`,
  `slot 2` = per-quadrature-point, `slot 4` = DOF offset via `ξ.𝑠`);
* why 1-D code **divides** by `𝐽`;
* why 2-D/3-D simplex code also divides (the difference product already
  carries the reference-simplex determinant 1/2 or 1/6);
* why isoparametric elements (`Quad4/8`, `Hex8`) never touch `𝐽`.

---

## 3.7 S1.5 — retire `error_estimates.jl` ✅

`src/operation/error_estimates.jl` was never `include`d (commented out in
`ApproxOperator.jl`) and used the long-dead `Operator{:tag}` dispatch style.
Its 17 functions were redistributed to the physics modules they belong to.

### Migrated (7 operators, 14 methods)

| Operator | Destination | Notes |
|---|---|---|
| `L₂_heat_flux` (el + vec) | `Heat` | heat-flux norm vs. exact flux `ξ.p₁/ξ.p₂` |
| `L₂_heat_flux_u` (el + vec) | `Heat` | flux recovered as `-∇u` |
| `H₃` (el + vec) | `Heat` | third-order Sobolev norm; **fixed 2 bugs** — see below |
| `Hₑ_up_fem` (el + vec) | `Elasticity` | |
| `Hₑ_up_mix` (el + vec) | `Elasticity` | |
| `Hₑ_up_mix_Q4P1` (el + vec) | `Elasticity` | |
| `Hₑ_Incompressible` (el + vec) | `Elasticity` | **kept** — see the near-miss below |

All were converted from `(op::Operator{:tag})(args)` to plain
`name(args)` with an explicit overload for the vector form, so callers no
longer need an `Operator` object:

```julia
# before
L₂Norm = (Operator(:L₂_heat_flux))(elements)
# after
L₂Norm = Heat.L₂_heat_flux(elements)
```

### Dropped (3 items)

| Item | Reason |
|---|---|
| `T_error` (el + vec) | Bit-for-bit the same computation as `Heat.L₂`, which compares `ξ.T`/`xᵢ.T`. Straight duplicate. |
| `set∇𝑢!` | No callers anywhere in `src/` or `test/`; also read `ξ[:∂𝝭∂x_]` with a trailing underscore, a field nothing produces. |
| *(nothing else)* | — |

### Two bugs fixed during the move

`H₃` referenced two things that do not exist in the package:

* `𝑤 = get𝑤(ap,ξ)` — `get𝑤` is defined nowhere → replaced with `𝑤 = ξ.𝑤`;
* `I = xᵢ.id` — `Node` has no `id` field (it is `𝐼`); the variable was never
  used → line deleted.

### Near-miss worth recording

`Hₑ_Incompressible` was initially classified as a duplicate of
`Elasticity.Hₑ_PlaneStrain_Deviatoric` and nearly deleted. They use different
shear conventions:

```
Incompressible:           ε₁₂ = 0.5*(u,y + v,x),  σ₁₂ = Cᵈ·ε₁₂
PlaneStrain_Deviatoric:   ε₁₂ =      u,y + v,x,    σ₁₂ = Cᵈ·ε₁₂/2
```

so `σ₁₂` agrees but the product `σ₁₂·ε₁₂` differs by a **factor of 2**. The
`σ₁₁`/`σ₂₂` terms do *not* scale, so the 2 does **not** cancel in
`√(ΔW²/W̄²)`:

```
Hₑ (Incompressible)          = 1.61995
Hₑ (PlaneStrain_Deviatoric)  = 1.38250
difference                   = 0.23745
```

`Hₑ_Incompressible` is therefore a distinct estimator and was migrated, not
dropped. Also fixed while moving: its `E = op.E` / `ν = op.ν` bindings became
`𝓖[1].E` / `𝓖[1].ν` (the last surviving `op.` reference, caught by the static
guard).

---

## 3.8 S2 — material constants as keyword arguments 🟡

Design: `S2_DESIGN.md`. Revised before execution — the proposed `resolve`
fallback was **dropped**; keywords are mandatory.

### Decision: `k` renamed to `kᶜ` (option A)

Julia rejects a signature where one name is both positional and keyword
(`ERROR: syntax: function argument name not unique: "k"`). Every heat /
weighted-residual / hamilton operator takes the stiffness matrix as `k`, and
the thermal conductivity was also `k`. Per review, option **A** was chosen: the
material constant is exposed as **`kᶜ`** (conductance), matching the local name
the code already used (`kᶜ = ξ.k`). The other names keep their spelling; a
docstring states the meaning where `α` / `k` / `c` differ between modules.

### Result

**123 of 263 operators** now take material constants as keyword arguments.

| File | operators | with kwargs |
|---|---:|---:|
| `elasticity.jl` | 60 | 25 |
| `thick_plate.jl` | 57 | 22 |
| `thin_shell.jl` | 28 | 13 |
| `heat.jl` | 27 | 7 |
| `hamilton.jl` | 18 | 15 |
| `hyperelasticity.jl` | 16 | 15 |
| `beam.jl` | 15 | 6 |
| `thin_plate.jl` | 14 | 6 |
| `phasefield.jl` | 7 | 3 |
| `curved_beam.jl` | 6 | 2 |
| `plasticity.jl` | 6 | 5 |
| `weighted_residual.jl` | 4 | 3 |
| `stokes.jl` | 1 | 1 |
| `test.jl` | 4 | 0 |
| **Total** | **263** | **123** |

The remainder are operators that need no material constant, plus the
element-only helpers.

### Shape of the change

```julia
# before
function ∫∫εᵢⱼσᵢⱼdxdy(ap::T, k::AbstractMatrix{Float64}) where T<:AbstractElement
    for ξ in 𝓖
        E = ξ.E
        ν = ξ.ν
        Cᵢᵢᵢᵢ = E/(1-ν^2)

# after
function ∫∫εᵢⱼσᵢⱼdxdy(ap::T, k::AbstractMatrix{Float64}; E, ν) where T<:AbstractElement
    for ξ in 𝓖
        Cᵢᵢᵢᵢ = E/(1-ν^2)

# call site
(∫∫εᵢⱼσᵢⱼdxdy => elements)(k; E = 210e9, ν = 0.3)
```

`prescribe!` keeps its role for boundary data and loads (`g`, `t`, `b`, `q`)
and for solution / exact fields (`u`, `d`, `T`); it is no longer used for
material constants.

### Verification

* **`103/103` tests pass.** The patch test's L2 error is unchanged at
  `6.594194188924542e-8`, so the migration did not alter the numerics.
* **keyword path == `prescribe!` path** — both assembly routes give the same
  matrix, and doubling `μ` doubles `k`.
* **missing keyword** — `@test_throws UndefKeywordError`, i.e. the failure is
  now at the call site and names the argument, instead of
  `KeyError: key :E not found` from a `Dict` lookup inside `Node.getproperty`.
* static scan: **0** operators still read a material constant from `ξ` / `ap`
  without a binding, and **0** dangling references.

### Two defects found and fixed during the conversion

1. `thin_shell.jl` `Hₑ_ThinShell` used `E`/`ν`/`h` before any binding
   (surfaced by the dedent, caught by the static guard).
2. The static guard itself had to learn that a material constant bound as a
   **keyword in the signature** is as valid as one bound inside the loop.

### Tooling

* `scripts/list_material_ops.py` — lists operators and the constants they read.
* `scripts/convert_to_keywords.py` — performs the conversion. Handles the
  `kᶜ` rename, `ξ`/`ξ₁`/`ξᵤ` quadrature-node spellings, an existing `;` in the
  signature, and processes blocks back-to-front so line numbers stay valid.

### S2.7 — docstrings ✅

All **263** operators in `src/operation/` now carry a docstring. Each one states

* the mathematical form of the weak term;
* which quadrature-point fields it reads, so it is discoverable whether
  `set𝝭!`, `set∇𝝭!` or `set∇²𝝭!` must be called first;
* each keyword argument and its meaning;
* where the same name means different things in different modules
  (`α` penalty vs. stiffness scaling; `k`/`kᶜ`; `c`/`c²`).

| File | documented |
|---|---:|
| `elasticity.jl` | 60/60 |
| `thick_plate.jl` | 57/57 |
| `thin_shell.jl` | 28/28 |
| `heat.jl` | 27/27 |
| `hamilton.jl` | 18/18 |
| `hyperelasticity.jl` | 16/16 |
| `beam.jl` | 15/15 |
| `thin_plate.jl` | 14/14 |
| `phasefield.jl` | 7/7 |
| `curved_beam.jl` | 6/6 |
| `plasticity.jl` | 6/6 |
| `weighted_residual.jl` | 4/4 |
| `test.jl` | 4/4 |
| `stokes.jl` | 1/1 |
| **Total** | **263/263** |

Example of the resulting entry:

```julia
"""
    ∫κᵢⱼMᵢⱼdΩ(ap, k; D, ν)

Kirchhoff plate bending stiffness, `∫ κᵢⱼ(w) Mᵢⱼ(v) dΩ`.

    a(w,v) = ∫_Ω D[∇²w ∇²v + ν(∇²w v,ₓᵧ + w,ₓᵧ ∇²v) + 2(1-ν) w,ₓᵧ v,ₓᵧ] dΩ

# Keyword arguments
- `D` : bending stiffness `E h³ / (12(1-ν²))`
- `ν` : Poisson ratio

# Element requirements
`set∇²𝝭!` — reads `ξ[:∂²𝝭∂x²]`, `ξ[:∂²𝝭∂x∂y]`, `ξ[:∂²𝝭∂y²]`, `ξ.𝑤`.
"""
```

### S2.8 — every material constant is typed ✅

All **234** keyword parameters across the 123 converted operators now carry an
explicit type:

```julia
# before
function ∫κᵢⱼMᵢⱼdΩ(ap::T, k::AbstractMatrix{Float64}; ν, D)
# after
function ∫κᵢⱼMᵢⱼdΩ(ap::T, k::AbstractMatrix{Float64}; ν::Real, D::Real)
```

**Why `Real` and not `Float64`:** `thick_plate.jl` accumulates some error norms
in `BigFloat`, and callers legitimately pass integer literals such as
`prescribe!(el, :ν => 0)`. Verified behaviour:

| argument | `::Real` | `::Float64` |
|---|---|---|
| `2.5` | accepted | accepted |
| `3` (Int) | accepted | **TypeError** |
| `big"2.5"` | accepted | **TypeError** |

So `Real` keeps every current call site working while still rejecting nonsense
(`"x"` → `TypeError`) at the signature instead of failing deep in a multiply.

### A second defect: `k` declared after `;` made the operators uncallable

Several modules — 78 occurrences across `thin_shell` (46), `hyperelasticity`
(10), `curved_beam` (10), `phasefield` (6), `plasticity` (5), `heat` (1) — wrote

```julia
function ∫εᵢⱼNᵢⱼκᵢⱼMᵢⱼdΩ(ap::T;k::AbstractMatrix{Float64})
```

with `k` **after the semicolon**, making the stiffness matrix a *required
keyword* rather than a positional argument. The `Pair` assembler calls
`form(elm, k)`, so none of those operators could be used through the documented
`(form => elements)(k)` path — a `MethodError` with no obvious cause.

Fixed by moving `k`, `f`, `fint`, `g` and `𝐑` back into the positional list, so
all 263 operators now share one calling convention. While doing so a
follow-up script over-matched and pulled `kᶜ` (the conductance, a genuine
material constant) into the positional list as well; that was caught by the
test suite and reverted.

### A defect found while writing the docstrings

The docstring generator matched the operator name inside an *existing*
docstring body instead of at `^function`, which nested 7 docstrings inside
others. Julia then failed with `cannot document the following expression`. The
file was reverted and re-applied with the anchor fixed; a check now scans for
orphan `"""` blocks (0 found).

---

## 3.9 S3 — rewrite `importmsh.jl` ✅

Design: `S3_DESIGN.md`. Decisions taken: **D1c** (delete the three piecewise
APIs), **D2b** (split into three files), **D3** (characterisation tests first).

### Where the code lives now

| file | lines | contents |
|---|---:|---|
| `preprocession/importmsh.jl` | 245 | `getPhysicalGroups`, `get𝑿ᵢ`, and the legacy signatures as thin forwarders |
| `preprocession/elements.jl` | 285 | `ElementSpec`, `resolve_*`, `entity_list`, `generate_*`, `log_summary`, `getElements(nodes, spec)` |
| `preprocession/coordinates.jl` | 350 | `initial_data`, `fill_coordinates!`, `fill_edge_coordinates!`, `fill_curved_coordinates!`, `fill_jacobian!`, `fill_normal!` |

`importmsh.jl` went from **1049 → 245 lines**. The 16 `quote` fragments and the
`@eval begin … end` block that spliced them into 17 function bodies are gone;
there is now **one** implementation with four dispatch points.

### S3.0 — characterisation tests

`test/importer_characterisation.jl`, 22 assertions, written **before** the
rewrite and green against the old code. They pin:

| observable | why |
|---|---|
| element count | catches a lost or duplicated entity |
| sorted `data` keys | catches a missing or extra field |
| `Σ ξ.𝑤` | catches a changed quadrature rule or Jacobian convention |
| quadrature centroid | catches swapped coordinates |
| `𝓖[1].𝑔` = 1, `𝓖[end].𝑔` = ng | catches the `𝑔` counter leaking across elements |
| `:ξ`, `:η` length = ng; `:x` length = ng·ne | catches local coordinates no longer being written |
| one shared `data` object | catches per-entity dictionaries, which break `prescribe!` |

**These tests earned their keep.** During the port they caught two regressions
that the patch test alone did not:

1. `𝑔` was advanced across elements instead of restarting at 1 for each one, so
   `ξ.ξ` indexed past the end of a length-1 array;
2. `:ξ`/`:η` stopped being written at all, because the responsibility moved out
   of the old `prequote` fragment without being picked up by
   `fill_coordinates!`.

Both were fixed before the suite was declared green.

### S3.1 / S3.5 — deletions

| removed | quantity |
|---|---|
| commented-out `getElements` (old QuoteNode version) | 38 lines |
| commented-out `getMacroCurvedElements` / `getMacroBoundaryElements` ×4 | 70 lines |
| `generateForMarco` fragment (0 users) | 22 lines |
| `getPiecewiseElements` ×2 | per D1c |
| `getPiecewiseBoundaryElements` ×2 | per D1c |
| `getCurvedPiecewiseElements` ×2 | per D1c |

The three piecewise APIs were deleted because they call `get𝑛𝑝(::Element)`, a
method that has never existed in this package (`get𝑛𝑝` is only defined for
`PiecewisePolynomial` and friends in `approximation/piecewise.jl`), so they could
not run; they also had zero callers. Their exports were removed from
`ApproxOperator.jl` as well.

### S3.3 — `getCurvedElements` and its unfixable signature

During the port the missing `fill_jacobian!` call was added to the curved path.
The API still cannot run, and **the legacy code could not either**:

```julia
getCurvedElements(nodes, dimTag, cs::Function, order)
```

requires `cs::Function`, but the implementation reads `cs.𝐽`, `cs.a₁₁`, … as
*fields*. A `Function` cannot carry those fields, so no argument can satisfy both
requirements. Verified: the old code raises `UndefVarError`, the new code raises
`FieldError`. `test/importer_characterisation.jl` records this with
`@test_throws` rather than pretending the API works, and the docstring says so.
Fixing it needs a decision about what `cs` should be (a struct with callable
metrics, or a `NamedTuple` of functions with no `Function` restriction) — that
is a separate change, not part of a behaviour-preserving refactor.

### S3.6 — docstrings

Every public entry point and helper carries an English docstring:
`ElementSpec`, `getElements(nodes, spec)`, `initial_data`, `fill_coordinates!`,
`fill_edge_coordinates!`, `fill_curved_coordinates!`, `fill_jacobian!`,
`fill_normal!`, `resolve_element`, `resolve_quadrature`, `entity_list`,
`generate_fem!`, `generate_neighbour!`, `log_summary`, plus all nine legacy
forwarders.

`fill_jacobian!` also documents a **preserved quirk**: the original appended the
`z`-direction inverse-Jacobian entries in an `elseif dim >= 3` branch that can
never run, because `dim >= 2` is tested first. Those arrays are never filled for
3-D elements. The rewrite keeps the behaviour so results do not change.

### S3.7 — centroid nodes

The feature is exposed in two equivalent ways:

```julia
# one call
nodes = get𝑿ᵢ(centroid = entities["Ω"])

# or two steps, when the nodes already exist
nodes = get𝑿ᵢ()
GmshImport.addCentroidNodes!(nodes, entities["Ω"])

# either way:
elements = getElements(nodes, entities["Ω"], 2)
```

`get𝑿ᵢ()` with no argument returns only the mesh nodes, so every existing call
site is unaffected — verified by the fact that all six call sites in `src/` and
`test/` use the no-argument form and the suite stays green.

Appends one node at the centroid of every element of a Gmsh entity and returns
the extended node list. The existing nodes are untouched and keep their `𝐼`; each
appended node takes the next free `𝐼` and the **same `data` dictionary**, so
`prescribe!` and `getElements` keep working on the enlarged list.

Two design points worth recording:

* **Corner mean, not all-node mean.** `corner_count(elementType)` gives the
  number of nodes that span the element geometry, and the centroid averages only
  those. Averaging every node of a `Tri6` or `Quad9` would pull the centroid
  towards the mid-side nodes, which is not the geometric centre. Verified on an
  order-2 quad: the centroid lands at exactly `(0.5, 0.5)`, not at the mean of
  the nine nodes.
* **Unknown element types** make `corner_count` return `-1`, in which case the
  function falls back to averaging all nodes rather than failing.

The optional `elementType` keyword restricts the operation to one Gmsh element
type; `:auto` (the default) processes every type Gmsh reports for the entity.

Verified on the patch-test mesh: 165 nodes → 453 nodes (`+288` triangle
centroids), and `getElements` still builds 288 elements from the extended list.
`get𝑿ᵢ(centroid = …)` and `get𝑿ᵢ()` + `addCentroidNodes!` produce indistinguishable
results (same count, same `𝐼`, same coordinates to `1e-12`), which is asserted in
`test/importer_characterisation.jl`.

### Test count

| | before S3 | after S3 | after S3.7 |
|---|---:|---:|---:|
| `Pkg.test()` | 103 | 125 | **193** |

The S3 increase is the 22 characterisation assertions; S3.7 added 68 more for
`addCentroidNodes!`, `corner_count` and the `get𝑿ᵢ(centroid = …)` keyword. `Pkg.test()` also works from a clean
environment now: `Project.toml` gained `LinearAlgebra` in the test target (the
material-parameter tests use `norm`) and the manifest was re-resolved.

---

## 4. Test suite---

## 4. Test suite

`99` assertions, all passing.

| File | Contents |
|---|---|
| `test/patchtest.jl` | End-to-end heat problem; `L2 error ≈ 6.59e-8` |
| `test/shape_function_consistency.jl` | 8 testsets: Seg2/Seg3, Tri3/Tri6, Quad4, Hex8/Tet4, Quad8, Tet10, TriBell, TriHermite |
| `test/operators_material_parameter.jl` | `KirchhoffPlate.∫κᵢⱼMᵢⱼdΩ`, `Heat.∫∫∇v∇udxdy`, plus 2 static guards |

Both static guards proved their worth during S1.5: the `op.` guard caught a
surviving `op.E` in the migrated `Hₑ_Incompressible`.

### Construction routes

* **Gmsh-backed** — Seg2, Seg3, Tri3, Tri6, Quad4, Hex8, Tet4 read through
  `GmshImport.getElements`, exercising the real importer layout.
* **Explicit** — Quad8, Tet10, TriBell, TriHermite built from explicit node
  coordinates via `build_element`, with **no `.msh` file, no `setOrder`, no
  Gmsh**. This removes the dependency on a mesher's node-ordering choices.

### What every element is checked for

1. **Partition of unity** — `Σᵢ Nᵢ(𝒙) = 1`
2. **Nodal delta** — `Nᵢ(𝒙ⱼ) = δᵢⱼ`
3. **Linear-field reproduction** — `u = c₁ + c₂x + c₃y + c₄z` and its gradient,
   which simultaneously validates the chain rule `∂N/∂x = ∂N/∂ξ · ∂ξ/∂x`

Tolerances: `≤ 1e-14` for the directly built elements, `≤ 1e-11` where
isoparametric Jacobians accumulate round-off.

### Static guards (catch silent regressions)

* no `op.` field access remains in `src/operation/` (dispatch arguments
  `op::Operator` excepted);
* no function uses `E`/`ν`/`h` before or inside its quadrature loop without
  binding them.

---

## 5. Remaining work

### S2 — material parameters → explicit keyword arguments

`prescribe!(elements, :D => D, :ν => ν)` works and is uniform, but the operator
signature still says nothing about what the caller must supply. Proposed:

```julia
∫κᵢⱼMᵢⱼdΩ(ap, k; D, ν)
```

with the `Pair` assembler forwarding `kw...`. Keeps the current path as a
compatibility fallback.

### S3 — `importmsh.jl` rewrite

999 lines → target ~350. Currently 13 near-identical `@eval`'d `getElements`
variants built from `QuoteNode` fragments; parameter combinations already
differ (some omit `cal_jacobe` or `normal`). Replace with one `ElementSpec`
struct and a single dispatch function.

### S4 — operator docstrings and the operator index

**What is done**

| item | state |
|---|---|
| every operator has an English docstring | ✅ 263/263 |
| every docstring states the quadrature-point fields it reads | ✅ |
| every docstring lists its keyword arguments | ✅ |
| every docstring names the meaning where a symbol differs between modules (`α`, `k`/`kᶜ`, `c`/`c²`) | ✅ |
| `docs/operators.md` generated from the source | ✅ 250 operators, 12 modules |
| generator script, so the table cannot drift | ✅ `scripts/build_operator_index.py` |
| mathematical content in the docstring | 260/263 (98%) |
| a **standalone** formula block | 67/263 (25%) |

**Accepted: inline formulas stay inline**

Most docstrings state the weak form *inline* in the summary, e.g.

```
Two-dimensional conduction stiffness, `∫ k ∇v·∇u dΩ`.
```

rather than as a display block:

```
Two-dimensional conduction stiffness.

    a(u,v) = ∫_Ω k ∇u·∇v dΩ
```

**This was reviewed and accepted as-is.** Both forms carry the same mathematical
content; the inline form is unambiguous for a single-term operator, and 67
operators — those whose content genuinely needs unfolding (energy norms,
phase-field updates, shell and plate functionals, elastoplastic yield surfaces) —
already use the display form. Reformatting the remaining ~190 would be pure
churn: the integrand differs for every operator, so no script can derive it, and
a scripted attempt produced duplicate and misfiled formulas that had to be
reverted.

The gap this leaves is *typographic*, not semantic: `docs/operators.md` carries a
description for all 250 operators, and every docstring states the weak form, the
required setter and the keyword arguments.

**Scripts added**

| script | purpose |
|---|---|
| `scripts/build_operator_index.py` | regenerates `docs/operators.md` from the docstrings |
| `scripts/add_math_forms.py` | the formula insertions that were applied (44 of them), kept for provenance |

### S5 — README, docs and dead code ✅

**README**

`README.md` went from one line to a full user guide: installation, a complete
runnable patch test, the data model (the two node types and the shared `data`
dictionary with its four slots), how to write a weak form (with the six rules the
package follows), how to add an element type (including the obligatory
consistency test), the module inventory, and how the test suite is organised.

The main example in the README was executed to confirm it prints
`L2 error: 6.594194188924542e-8` as documented.

**`docs/`**

| file | purpose |
|---|---|
| `docs/make.jl` | Documenter.jl build; `warnonly` set so the 264 operator docstrings do not fail on cross-reference warnings |
| `docs/src/index.md` | landing page: quick reference and the setter table |
| `docs/src/operators.md` | generated operator index (264 entries, 14 modules) |

`scripts/build_operator_index.py` now writes to `docs/src/operators.md`, so the
index and the Documenter site are the same file.

**Keyword `getElements`** — not done, deliberately. The positional forms
(`getElements(nodes, dimTag, 2)`, `getElements(nodes, dimTag, Element{:Tri6}, 3)`)
are already unambiguous and are what every script uses; adding a parallel
keyword API would create two ways to do one thing. The `ElementSpec` constructor
*is* the keyword API for anyone who needs it:

```julia
getElements(nodes, ElementSpec(dimTag; element = Element{:Tri6}, normal = true))
```

**Dead code — resolved rather than deleted**

| file | before | after |
|---|---|---|
| `approximation/grkgsi.jl` | not included; depended on `AbstractReproducingKernel`, **a type that has never existed in the package** | **deleted** |
| `approximation/rkwave.jl` | not included; extends `ReproducingKernel{:Wave2D}`, which does exist | **included**, tested |
| `operation/phasefield.jl` | not included; dispatched on the deleted `Operator{:tag}` | **converted and included** |
| `operation/plasticity.jl` | same | **converted and included** |

The two operator files were converted from `(op::Operator{:tag})(args)` to plain
functions, as in S1.5.

**Nineteen operators had never been executed, and it showed.** Bringing them
into the build exposed a run of genuine defects:

| defect | file |
|---|---|
| `ε`, `v` read before assignment | `phasefield.jl` `UPDATE_PFM_1D` |
| `η` (residual stiffness) read but never declared | `phasefield.jl`, 4 operators |
| `kc` (fracture toughness) read but never declared | `phasefield.jl`, 3 operators |
| `λ`, `μ`, `c` (Tresca) not in the signature | `plasticity.jl` |
| `c` (Mohr–Coulomb cohesion) not in the signature | `plasticity.jl` |
| `K`, `σy` (hardening, yield) not in the signature | `plasticity.jl` |
| `tol` (return-mapping tolerance) not in the signature | `plasticity.jl`, 3 operators |
| `getσₙ` called unqualified across module boundary | `plasticity.jl` |
| three title-only docstrings (`"""Tresca"""`) orphaned the real docstring | `plasticity.jl` |

None of these could be seen while the files were excluded from the build.

**Test count**

| | before S5 | after S5 |
|---|---:|---:|
| `Pkg.test()` | 199 | **227** |

The 28 new assertions are `test/restored_modules.jl`: the 13 restored operators
are defined, the `Wave2D` basis is trigonometric at `x = 0`, the plasticity
operators assemble a finite matrix on a real quad mesh, and the phase-field
update writes a non-negative history variable.

---

## 3.11 S6 — packaging and hygiene ✅

Design: `S6_DESIGN.md`. Decisions taken: **D1** `julia = "1.9"`, **D2**
`REPRODUCIBILITY.md` instead of committing the manifest, **D3** re-export the
importer entry points, **D4** Aqua.jl (your choice; the design recommended the
hand-written checks, Aqua renders most of them unnecessary).

### `Project.toml`

```toml
[deps]                # before: Gmsh, Printf, WriteVTK
Gmsh   = "..."        # after:  Gmsh, Printf
Printf = "..."

[compat]              # before: WriteVTK only, julia = "1"
Aqua = "0.8"; Gmsh = "0.3"; LinearAlgebra = "1"
Printf = "1"; Test = "1"; julia = "1.9"
```

`WriteVTK` was a direct dependency that `src/` never imported — dropped. Every
remaining dep, including the test-only ones, now has a compat bound.

### Aqua.jl — six checks, three real defects

Added to `[extras]` and the test target; `test/package_hygiene.jl` runs it.
Three of the six checks failed on the tree as it stood:

| check | before | defect |
|---|---|---|
| `test_undefined_exports` | **5 failures** | `getPhysicalGroups`, `get𝑿ᵢ`, `getElements`, `getCurvedElements` were exported but only defined in the `GmshImport` submodule, so `using ApproxOperator` handed the user nothing; **`addEdgeElements` was exported and defined nowhere at all** |
| `test_ambiguities` | **2 failures** | `push!(aps::Vector{<:AbstractElement}, ::Symbol...)` and `push!(aps, ::Pair...)` were both variadic, so the empty call `push!(aps)` could not be resolved |
| `test_deps_compat` | **1 failure** | `[compat]` covered only `WriteVTK`, the unused package |
| `test_stale_deps` | pass | — |
| `test_piracies` | pass | the `Base` overloads all act on types the package owns |
| `test_persistent_tasks` | pass | — |

The fixes:

```julia
# ApproxOperator.jl — make the exports name something
using .GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements, getCurvedElements
export getPhysicalGroups, get𝑿ᵢ, getElements, getCurvedElements, getDOFs
# `addEdgeElements` removed: it does not exist anywhere in the package
```

```julia
# element.jl — require at least one argument, which removes the empty-call ambiguity
Base.push!(aps::Vector{T}, s::Symbol, ss::Symbol...; index=:𝑠)
Base.push!(aps::Vector{T}, sv::Pair{Symbol,Vector{Float64}},
           svs::Pair{Symbol,Vector{Float64}}...; index::Symbol=:𝑠)
```

A behavioural test also asserts that the five importer names resolve after
`using ApproxOperator`, and that `addEdgeElements` does **not** — so the export
list cannot silently drift back.

### `REPRODUCIBILITY.md`

`Manifest.toml` stays untracked (the usual convention: a committed manifest makes
`Pkg.develop`/`Pkg.add` unpredictable downstream). Instead `REPRODUCIBILITY.md`
records the Julia version, the `gmsh_jll` tree hash, and — importantly — the
**numerical values that must reproduce**, each with the test file that asserts
it:

| quantity | value |
|---|---|
| patch-test L2 error | `6.594194188924542e-8` |
| 2×2 Quad4 mesh, Σ ξ.𝑤 | `1.0` |
| 1-D boundary mesh, Σ ξ.𝑤 | `4.0` |
| flat curved map / doubled map | `1.0` / `2.0` |
| patch-test element count | `288` |

This matters because the integration weights come from
`gmsh.model.mesh.getIntegrationPoints`, so a different `gmsh_jll` can change the
last digits — and the tests will say so.

### Test count

| | before S6 | after S6 |
|---|---:|---:|
| `Pkg.test()` | 227 | **242** |

The 15 new assertions are `test/package_hygiene.jl`: the six Aqua checks plus the
five reachability assertions and four negative/positive export checks.



```
 M src/approximation/piecewise.jl
 M src/approximation/seg3.jl
 M src/approximation/seghermite.jl
 M src/approximation/tri6.jl
 M src/operation/curved_beam.jl
 M src/operation/hyperelasticity.jl
 M src/operation/phasefield.jl
 M src/operation/plasticity.jl
 M src/operation/thin_plate.jl
 M src/operation/thin_shell.jl
 M src/preprocession/importmsh.jl
 M src/ApproxOperator.jl
 M test/runtests.jl
 D src/operation/error_estimates.jl
 D test/timoshenko_geometric_stiffness.jl
 D test/timoshenko_static_ss.jl
?? test/patchtest.jl
?? test/shape_function_consistency.jl
?? test/operators_material_parameter.jl
?? scripts/fix_op_globals.py
```

`src/operation/thin_shell.jl` shows a large raw diff because indentation was
normalised; `git diff -w` reduces it to ~90 lines of real change.

All comments and docstrings are in **English**.
