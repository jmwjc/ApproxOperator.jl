#!/usr/bin/env python3
"""
S4: add the mathematical form to docstrings that lack one.

Three kinds of docstring need work, handled differently:

  * **vector overloads** (`f(aps)`, `f(aps_w, aps_φ)`) — these delegate to the
    element-level overload, so they get a one-line note pointing at it rather
    than a duplicated formula;
  * **helpers** (`𝐑ᵀ𝐊𝐑`, `getσₙ`, geometry maps) — these get their defining
    formula, which is short and exact;
  * **operators** — these get the bilinear form `a(u,v) = …`.

Run from the repository root:

    python3 scripts/add_math_forms.py [--dry-run]
"""

import re
import sys

# ---------------------------------------------------------------------------
# Insertions, keyed by "file:line(function header)". A value is a block of text
# inserted immediately after the docstring's opening line, i.e. before the
# signature line, so the formula sits right under the summary.
# ---------------------------------------------------------------------------

# Each entry: (file, function name, occurrence index) -> text to insert.
# Occurrence index disambiguates identical function names in one file.
FORMULAE = {
    # ---- beam.jl : helpers -------------------------------------------------
    ("beam.jl", "𝐑ᵀ𝐊𝐑", 0): """
    Kᵍ[a,b] = Σᵢ Σⱼ R[i,a] K[i,j] R[j,b]
""",
    ("beam.jl", "𝐑", 0): """
    R = [ 𝒍₁ 𝒍₂ 𝒍₃ ]      with  𝒍₁ = (𝒙₂ - 𝒙₁)/|𝒙₂ - 𝒙₁|

so that a local vector `vˡ` maps to the global frame as `𝒗 = R vˡ`.
""",
    ("beam.jl", "𝐐", 0): """
    Q = diag(R, R)        acting on the six beam DOFs (three translations,
                          three rotations) ordered as (u₁,u₂,u₃,φ₁,φ₂,φ₃)
""",
    # ---- beam.jl : missing formula ----------------------------------------
    ("beam.jl", "𝐊ᵘᵠ!", 0): """
    Kᵘᵠ = ∫ Bᵀ 𝐊φ Bφ ds      (row block u, column block φ)
""",
    ("beam.jl", "𝐊ᵠᵘ!", 0): """
    Kᵠᵘ = (Kᵘᵠ)ᵀ = ∫ Bφᵀ 𝐊φ B ds    (row block φ, column block u)
""",

    # ---- elasticity.jl ----------------------------------------------------
    ("elasticity.jl", "g₂", 0): """
    Assemble the prescribed value into the right-hand side and zero the
    corresponding row and column of the stiffness matrix — the "row/column
    elimination" strong-imposition of a Dirichlet condition:

        f[I] -= K[I,J] g        for the constrained DOF J
        K[I,J] = K[J,I] = 0,    K[J,J] = 1,    f[J] = g
""",
    ("elasticity.jl", "getσₙ", 0): """
    σₙ = (σ₁₁ + σ₂₂)/2 + √(((σ₁₁ - σ₂₂)/2)² + σ₁₂²)

the larger principal stress in plane stress.
""",
    ("elasticity.jl", "getσₙ", 1): """
    σₙ = max(|σᵢ - σⱼ|)/2 over the three principal stress differences, i.e. the
    Mohr's-circle radius of the 3-D stress state.
""",
    ("elasticity.jl", "Hₑ", 1): """
    ‖e‖_E / ‖u‖_E,   ‖v‖_E² = Σ_e ∫_Ω ε(v):C:ε(v) dΩ
""",
    ("elasticity.jl", "Hₑ_PlaneStress", 1): """
    ‖e‖_E / ‖u‖_E,   C₁₁₁₁ = E/(1-ν²),   C₁₁₂₂ = Eν/(1-ν²),   C₁₂₁₂ = E/(2(1+ν))
""",
    ("elasticity.jl", "Hₑ_PlaneStrain_Deviatoric", 1): """
    ΔW² = ½ (σ - σ̄):(ε - ε̄),   σ = Cᵈ ε,   Cᵈ = E/(1+ν),   Pᵈ = I - ⅓ I⊗I
""",
    ("elasticity.jl", "Hₑ_PlaneStrain_Deviatoric", 1): """
    ‖e‖_E / ‖u‖_E with the deviatoric energy above, summed over the elements.
""",
    ("elasticity.jl", "Hₑ_up_fem", 1): """
    (Hₑ, L₂(u), Hₑ_dil, Hₑ_dev, L₂(p))  with

        Hₑ_dil² = ∫ 3(p - p̄)²/(2K) dΩ,     Hₑ_dev² = ∫ (s - s̄):(s - s̄)/(4G) dΩ
""",
    ("elasticity.jl", "Hₑ_up_mix", 1): """
    Same five relative errors as `Hₑ_up_fem`, but with `p` interpolated from its
    own element set:   p(𝒙) = Σⱼ Nᵖⱼ(𝒙) qⱼ
""",
    ("elasticity.jl", "Hₑ_up_fem", 1): """
    Element-set overload of `Hₑ_up_fem(apu)`, summed over the elements.
""",
    ("elasticity.jl", "Hₑ_up_mix", 1): """
    Element-set overload of `Hₑ_up_mix(apu, app)`, summed over the elements.
""",
    ("elasticity.jl", "Hₑ_up_mix_Q4P1", 1): """
    Same five relative errors as `Hₑ_up_mix`, with the discrete pressure taken as
    the value at the first quadrature point (the constant P1 value per element).
""",
    ("elasticity.jl", "Hₑ_up_mix_Q4P1", 1): """
    Element-set overload of `Hₑ_up_mix_Q4P1(apu, app)`.
""",
    ("elasticity.jl", "Hₑ_Incompressible", 1): """
    ε₁₂ = ½(u,ᵧ + v,ₓ),   σ = Cᵈ ε,   Cᵈ = E/(1+ν)
""",
    ("elasticity.jl", "Hₑ_Incompressible", 1): """
    Element-set overload of `Hₑ_Incompressible(ap)`.
""",

    # ---- hamilton.jl ------------------------------------------------------
    ("hamilton.jl", "stabilization_bar_LSG", 0): """
    a_h(u,v) = a(u,v) + α Σ_e ∫_Ω ℒu · ℒv dΩ

with `ℒ` the residual of the strong form, `ℒu = ρA u,ₜₜ - EA u,ₓₓ`.
""",
    ("hamilton.jl", "stabilization_bar_LSG_Γ", 0): """
    Boundary counterpart of `stabilization_bar_LSG`:

        a_Γ(u,v) = α Σ_e ∫_Γ (ℒu)(ℒv) dΓ
""",
    ("hamilton.jl", "truncation_error", 0): """
    For the exact characteristic direction `ξ = x - c t`, the four defect measures

        fₓ  += BₓᵢBₓⱼ (Δx - cΔt),      fₜ  += BₜᵢBₜⱼ (Δx - cΔt)
        fₓₓ += BₓᵢBₓⱼ (Δx - cΔt)²,     fₜₜ += BₜᵢBₜⱼ (Δx - cΔt)²
""",
    ("hamilton.jl", "truncation_error", 1): """
    Element-set overload of `truncation_error`, accumulating the same four
    measures into freshly allocated vectors of length `nₚ`.
""",
    ("hamilton.jl", "Hₑ", 2): """
    ‖e‖_E / ‖u‖_E,   ‖v‖_E² = Σ_e ∫ (ρA v,ₜ² + EA v,ₓ²) dx dt
""",

    # ---- heat.jl ----------------------------------------------------------
    ("heat.jl", "g", 0): """
    f[I] += Σ_g Nᵢ(𝒙_g) g(𝒙_g) w_g
""",

    # ---- phasefield.jl ----------------------------------------------------
    ("phasefield.jl", "UPDATE_PFM_1D", 0): """
    ℋ ← max(ℋ, (ε - 1)²)        at every quadrature point
""",
    ("phasefield.jl", "UPDATE_PFM_2D", 0): """
    ℋ ← max(ℋ, ψ⁺(ε))           at every quadrature point
""",
    ("phasefield.jl", "CRACK_NORMAL", 0): """
    𝒏 ← -∇v / |∇v|              at every quadrature point
""",

    # ---- plasticity.jl ----------------------------------------------------
    ("plasticity.jl", "∫vᵢσdΩ_tresca", 0): """
    Tresca yield surface

        f = σ₁ - σ₃ - 2c,        with σ₁ ≥ σ₂ ≥ σ₃ the principal stresses

and the consistent tangent `Cᵗ` obtained from the return mapping.
""",
    ("plasticity.jl", "∫vᵢσdΩ_frictional_contact_2", 0): """
    Slip criterion  f = |τ| + μ̄ σ, the `_2` variant of the Coulomb contact
    operator, differing in how the algorithmic tangent is condensed.
""",

    # ---- test.jl ----------------------------------------------------------
    ("test.jl", "cc∇𝝭", 0): """
    ‖∇uʰ - ∇ū‖_L₂ / ‖∇ū‖_L₂,   ∇uʰ = Σᵢ ∇Nᵢ dᵢ
""",
    ("test.jl", "cc∇𝝭", 1): """
    Element-set overload of `cc∇𝝭(ap)`, summed over the elements.
""",

    # ---- thick_plate.jl ---------------------------------------------------
    ("thick_plate.jl", "L₂φ", 0): """
    ‖φʰ - φ̄‖_L₂ / ‖φ̄‖_L₂,   φʰ = Σᵢ Nᵢ φᵢ
""",
    ("thick_plate.jl", "L₂φ", 1): """
    Element-set overload of `L₂φ(ap)`.
""",
    ("thick_plate.jl", "H₁φ", 0): """
    ‖(φʰ - φ̄, ∇φʰ - ∇φ̄)‖ / ‖(φ̄, ∇φ̄)‖
""",
    ("thick_plate.jl", "L₂w", 0): """
    ‖wʰ - w̄‖_L₂ / ‖w̄‖_L₂,   wʰ = Σᵢ Nᵢ wᵢ
""",
    ("thick_plate.jl", "H₁w", 0): """
    ‖(wʰ - w̄, ∇wʰ - ∇w̄)‖ / ‖(w̄, ∇w̄)‖
""",
    ("thick_plate.jl", "L₂w", 1): """
    Element-set overload of `L₂w(ap)`.
""",
    ("thick_plate.jl", "H₁", 0): """
    Relative H¹ error of the pair `(w, φ)`, combining the deflection and
    rotation errors over the two element sets.
""",
    ("thick_plate.jl", "Hₑ", 2): """
    ‖e‖_E / ‖u‖_E,   ‖v‖_E² = ∫ (κ(v)·M(κ) + γ(v)·Q(γ)) dΩ
""",
    ("thick_plate.jl", "Hₑ", 2): """
    Element-set overload of `Hₑ(ap_w, ap_φ, ap_Q)`, summed over the elements.
""",
    ("thick_plate.jl", "Hₑ", 3): """
    ‖e‖_E / ‖u‖_E for the mixed `(w, φ, Q)` formulation:

        ‖v‖_E² = ∫ (κ(v)·M(κ) + γ(v)·Q(γ)) dΩ
""",

    # ---- thin_plate.jl ----------------------------------------------------
    ("thin_plate.jl", "ΔM̃ₙₛg", 0): """
    ∫ ΔM̃ₙₛ ḡ dΓ,  the symmetric-mixed-derivative twin of `ΔMₙₛg`.
""",

    # ---- thin_shell.jl ----------------------------------------------------
    ("thin_shell.jl", "ΔMₙₛ𝒂₃𝒈_Nitsche", 0): """
    ∫ (ΔMₙₛ 𝒂₃ · 𝒈) dΓ  — corner force in the Nitsche form.
""",
    ("thin_shell.jl", "ΔMₙₛ𝒂₃𝒈_HR", 0): """
    ∫ (ΔMₙₛ 𝒂₃ · 𝒈) dΓ  — corner force in the Hellinger–Reissner form.
""",
    ("thin_shell.jl", "ΔMₙₛ𝒂₃𝒗", 0): """
    ∫ (ΔMₙₛ 𝒂₃ · 𝒗) dΓ  — the trial-side counterpart of `ΔMₙₛ𝒂₃𝒈`, used to
    assemble the left-hand side of the corner condition.
""",
    ("thin_shell.jl", "ScordelisLoRoof_𝐴", 0): """
    𝒙(ξ,η) = ( L ξ,  R sin(φ),  R cos(φ) ),   φ = (1 - η) π/6
""",
    ("thin_shell.jl", "SphericalShell_𝐴", 0): """
    𝒙(ξ,η) = R ( sinθ cosφ,  sinθ sinφ,  cosθ ),   θ = η π/2,  φ = 2ξ π
""",
    ("thin_shell.jl", "SphericalShell_𝐵", 0): """
    θ(η) = η π/2,   the polar angle of the spherical parameterisation.
""",
}


def find_insert_point(lines, idx):
    """Given the line index of a `function` header, return the index of the
    opening `\"\"\"` of the docstring immediately above it, or None."""
    k = idx - 1
    if k < 0 or lines[k].strip() != '"""':
        return None
    while k > 0:
        k -= 1
        if lines[k].strip() == '"""':
            return k
    return None


def main():
    dry = "--dry-run" in sys.argv
    root = "src/operation/"
    files = {}
    counts = {}

    for (fname, func, occ), block in FORMULAE.items():
        path = root + fname
        if path not in files:
            files[path] = open(path, encoding="utf-8").read().split("\n")
        lines = files[path]
        seen = 0
        hit = None
        for i, l in enumerate(lines):
            m = re.match(r"^\s*function (?:\(op::Operator\{:([^}]+)\}\)|([^\s(]+))", l)
            if not m:
                continue
            name = m.group(1) or m.group(2)
            if name != func:
                continue
            if seen == occ:
                hit = i
                break
            seen += 1
        if hit is None:
            print(f"  MISS {fname}:{func}[{occ}]")
            continue
        open_idx = find_insert_point(lines, hit)
        if open_idx is None:
            print(f"  NO DOC {fname}:{func}[{occ}]")
            continue
        # Insert after the *summary paragraph*, i.e. after the first blank line
        # that follows the opening `"""` and its signature block. Putting the
        # formula directly under `"""` would push it above the signature and
        # read as part of the name.
        j = open_idx + 1
        # skip the signature block: any indented lines and blanks
        while j < len(lines) and (not lines[j].strip() or lines[j].startswith("    ")):
            j += 1
        # skip the summary paragraph up to the next blank line
        while j < len(lines) and lines[j].strip():
            j += 1
        # `j` is now at the blank line after the summary
        lines[j:j] = [""] + block.rstrip("\n").split("\n")
        counts[fname] = counts.get(fname, 0) + 1

    if not dry:
        for path, lines in files.items():
            open(path, "w", encoding="utf-8").write("\n".join(lines))
    print()
    for k in sorted(counts):
        print(f"  {k:26} +{counts[k]}")
    print(f"  {'TOTAL':26} +{sum(counts.values())}")


if __name__ == "__main__":
    main()
