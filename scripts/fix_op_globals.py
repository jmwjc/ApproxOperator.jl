#!/usr/bin/env python3
"""
S1: remove the undefined global `op` from src/operation/*.jl.

Background
----------
Seven operator modules (thin_plate, thin_shell, hyperelasticity, curved_beam,
error_estimates, phasefield, plasticity) read their material/solver constants
from a global variable `op` that no longer exists in the package: it used to be
an `Operator` struct that was deleted at some point, leaving 138 dangling
`op.<field>` references. Every one of those operators raised

    UndefVarError: op not defined

as soon as it was called.

After this script the constants are read from the quadrature-point node `ξ`,
exactly like elasticity.jl (`ξ.E`, `ξ.ν`), stokes.jl (`ξ.μ`) and
thick_plate.jl (`ap.E`) already do. The caller supplies the values through the
existing mechanism:

    prescribe!(elements, :D => D, :ν => 0.3)
    (∫κᵢⱼMᵢⱼdΩ => elements)(k)

Transformation
--------------
For each function:

    function f(ap::T, ...)
        𝓒 = ap.𝓒; 𝓖 = ap.𝓖
        D = op.D                    # <- removed
        ν = op.ν                    # <- removed
        for ξ in 𝓖
            ...  D    ...           # <- becomes ξ.D
            ...  ν    ...           # <- becomes ξ.ν

i.e. the pre-loop aliases are dropped and each use site inside the quadrature
loop is rewritten to read the value from `ξ`. Identifier uses are matched with
word boundaries so that `α` inside `αβ` or `D₁₁` are left untouched.
"""

import re
import sys
from pathlib import Path

TARGETS = [
    "thin_plate.jl",
    "thin_shell.jl",
    "hyperelasticity.jl",
    "curved_beam.jl",
    "error_estimates.jl",
    "phasefield.jl",
    "plasticity.jl",
]

ROOT = Path(__file__).resolve().parent.parent

# `name = op.field`  (the whole right-hand side must be `op.field`)
PRELOOP = re.compile(
    r"^(?P<indent>[ \t]+)(?P<name>[A-Za-zΑ-Ωα-ω_][A-Za-z0-9Α-Ωα-ω_₀-₉⁰-⁹ᵢⱼₖₗᵘᵛᵠ̄ᵉ¹²³]*)\s*=\s*op\.(?P<field>\S+?)\s*$"
)

# Any surviving `op.<field>`
OPREF = re.compile(r"\bop\.(?P<field>[A-Za-zΑ-Ωα-ω_][A-Za-z0-9Α-Ωα-ω_₀-₉⁰-⁹ᵢⱼₖₗᵘᵛᵠ̄ᵉ¹²³]*)")


def transform(text: str) -> tuple[str, list[tuple[str, str]]]:
    """Return (new_text, [(alias, field), ...]) with aliases lifted into ξ."""
    lines = text.split("\n")
    out: list[str] = []
    lifted: list[tuple[str, str]] = []
    pending: list[tuple[str, str]] = []   # aliases collected for the next loop

    in_loop = False           # inside `for ξ in 𝓖` / `for (ii,ξ) in enumerate(𝓖)`
    loop_depth = 0            # block depth relative to the loop header
    seen_loop = False         # aliases declared before a loop are collectable

    for line in lines:
        stripped = line.strip()

        # --- new function: aliases are function-local ------------------------
        if re.match(r"^(function\b|\s*function\b)", line):
            lifted = []
            pending = []

        # --- entering a quadrature loop -------------------------------------
        if re.match(r"^for\s+(\([^)]*\)|[A-Za-zΑ-Ωα-ω_][\w]*)\s+in\s+(enumerate\()?𝓖\b", stripped):
            in_loop = True
            seen_loop = True
            loop_depth = 0
            lifted = pending
            pending = []
            out.append(line)
            # Nothing is injected here: the use sites below are rewritten to
            # `ξ.<field>` directly, which needs no local alias at all.
            continue

        if info := re.match(r"^(\s*)for\b", line):
            if in_loop:
                loop_depth += 1
            out.append(line)
            continue
        if stripped.startswith("if ") or stripped.startswith("elseif ") or stripped.startswith("while "):
            if in_loop:
                loop_depth += 1
            out.append(line)
            continue
        if stripped == "else" or stripped.startswith("else "):
            out.append(line)
            continue
        if stripped == "end":
            if in_loop:
                if loop_depth == 0:
                    in_loop = False
                else:
                    loop_depth -= 1
            out.append(line)
            continue

        # --- pre-loop alias declarations ------------------------------------
        if not in_loop:
            m = PRELOOP.match(line)
            if m:
                pending.append((m.group("name"), m.group("field")))
                continue      # drop the alias, uses are rewritten to ξ.<field>

        out.append(line)

    # Remove the insertion markers and rewrite identifier uses inside loops.
    final: list[str] = []
    in_loop = False
    loop_depth = 0
    for line in out:
        stripped = line.strip()
        if re.match(r"^for\s+(\([^)]*\)|[A-Za-zΑ-Ωα-ω_][\w]*)\s+in\s+(enumerate\()?𝓖\b", stripped):
            in_loop = True
            loop_depth = 0
            final.append(line)
            continue
        if in_loop:
            if re.match(r"^\s*for\b", line) or stripped.startswith("if ") or stripped.startswith("while "):
                loop_depth += 1
            elif stripped == "end":
                if loop_depth == 0:
                    in_loop = False
                else:
                    loop_depth -= 1
            else:
                for name, field in lifted:
                    # Only rewrite uses of the lifted alias, never `op.` itself
                    # and never a substring of a longer identifier.
                    line = re.sub(rf"(?<![\w.]){re.escape(name)}(?![\w])", f"ξ.{field}", line)
        final.append(line)

    return "\n".join(final), lifted


def main() -> int:
    problems: list[str] = []
    for name in TARGETS:
        path = ROOT / "src" / "operation" / name
        if not path.exists():
            print(f"skip (missing): {path}")
            continue
        original = path.read_text()
        updated, lifted = transform(original)
        path.write_text(updated)
        fields = sorted({f for _, f in lifted})
        print(f"{name}: {len(lifted)} pre-loop aliases lifted -> ξ.{' , ξ.'.join(fields) if fields else '?'}")
        for leftover in OPREF.finditer(updated):
            problems.append(f"{name}: {leftover.group(0)}")

    print()
    if problems:
        print("REMAINING `op.` references that need manual review:")
        for p in problems:
            print("   ", p)
        return 1
    print("no `op.` reference left in src/operation/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
