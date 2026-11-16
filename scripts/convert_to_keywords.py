#!/usr/bin/env python3
"""
S2 conversion: turn material constants read from the quadrature point into
keyword arguments.

For each operator that reads `ξ.<name>` where `<name>` is a material constant,
the script

  * appends `; <kw>` (or `; <kw>...`) to the function signature, using the
    keyword spelling from RENAME for the symbols that collide with a positional
    argument;
  * deletes the `name = ξ.name` alias line inside the quadrature loop, because
    the name is now a parameter;
  * rewrites inline `ξ.name` uses to the bare `name`.

`k` is special: every operator takes the stiffness matrix as a positional
argument called `k`, and Julia rejects a signature where the same name is both
positional and keyword ("function argument name not unique"). The thermal
conductivity / stiffness `ξ.k` is therefore exposed as `kᶜ` (conductance),
matching the local name the code already used.
"""

import re
import sys

# Material constants that can become keywords under their own name.
PLAIN = ["E", "ν", "h", "D", "ρ", "ρA", "m", "μ", "λ", "σy", "K", "tol",
         "α", "c", "c²", "EA", "EI", "kGA", "R", "l"]

# Symbols that must be renamed because the operator already has a positional
# argument of that name.
RENAME = {"k": "kᶜ"}

ALL = PLAIN + list(RENAME)

# A quadrature-node identifier: ξ, ξ₁, ξ₂, ξᵤ, ξₚ, ξₛ, ...
QNODE = r"ξ[\w₀-₉ᵘᵖˢʷᵠ]*"

# Where each symbol is read from on the quadrature point.
SOURCE = {s: s for s in PLAIN}
SOURCE["k"] = "k"


OPEN = re.compile(r"^(function|for|while|if|try|let|begin|struct|quote|macro|module)\b")
NEUTRAL = re.compile(r"^(else|elseif|catch|finally)\b")
CLOSE = re.compile(r"^end\b")


def find_operators(text):
    """Yield (start, end) line spans of top-level `function` blocks.

    Depth is tracked with an explicit counter; `else`/`elseif`/`catch`/`finally`
    do not change it, and a line that opens and closes in the same statement
    (e.g. `for i in 1:n; x end`) is still counted correctly because only the
    presence of a leading keyword matters.
    """
    lines = text.split("\n")
    i = 0
    while i < len(lines):
        if not lines[i].startswith("function "):
            i += 1
            continue
        start = i
        depth = 0
        j = i
        closed = False
        while j < len(lines):
            code = lines[j].strip().split("#")[0].strip()
            if OPEN.match(code):
                depth += 1
            elif NEUTRAL.match(code):
                pass
            elif CLOSE.match(code):
                depth -= 1
                if depth == 0:
                    yield start, j
                    closed = True
                    break
            j += 1
        i = (j + 1) if closed else (i + 1)


def strip_comments(s):
    return "\n".join(re.sub(r"#.*$", "", ln) for ln in s.split("\n"))


def convert(path, dry_run=False):
    text = open(path, encoding="utf-8").read()
    lines = text.split("\n")
    spans = list(find_operators(text))

    changed = []

    # Process from the last block to the first: replacing a later block shifts
    # only the line numbers after it, which we have already handled.
    for start, end in reversed(spans):
        header = lines[start]
        pre_where = header.split("where")[0]

        # Already converted?  A `;` followed by bare names (no `::`) means yes;
        # a `;` followed by `name::Type` is just the positional separator.
        if ";" in pre_where:
            after = pre_where.split(";", 1)[1]
            if "::" not in after:
                continue

        body = "\n".join(lines[start + 1:end + 1])
        code = strip_comments(body)

        used = []
        for name in ALL:
            pat = QNODE + r"\." + re.escape(name) + r"(?![\w])"
            alias = (r"^\s*[A-Za-zΑ-Ωα-ω_]\w*\s*=\s*" + QNODE + r"\."
                     + re.escape(name) + r"\s*$")
            if re.search(pat, code, re.M) or re.search(alias, code, re.M):
                used.append(name)
        if not used:
            continue

        kws = [RENAME.get(n, n) for n in used]

        if dry_run:
            changed.append((header.strip()[:60], kws))
            continue

        # --- signature: insert before the closing paren of the argument list
        close = header.rfind(")")
        if ";" in pre_where:
            # `f(a::T; k::AbstractMatrix)` — the single `;` starts the keyword
            # list, so the new bare names join that same list. Julia allows only
            # one `;` per argument list.
            new_header = (header[:close].rstrip() + ", " + ", ".join(kws) + ")"
                          + header[close + 1:]) if False else \
                         (header[:close].rstrip() + ", " + ", ".join(kws) + ")"
                          + header[close + 1:])
        else:
            # Add the first semicolon. `pre_where` ends with `)` and may be
            # followed by " where ...", so rebuild from the original header.
            new_header = (header[:close].rstrip() + "; " + ", ".join(kws) + ")"
                          + header[close + 1:])

        # --- body: drop the alias line, rewrite inline reads to the keyword
        new_body = "\n".join(lines[start + 1:end + 1])
        for name in used:
            kw = RENAME.get(name, name)
            src = SOURCE[name]
            new_body = re.sub(
                r"^[ \t]*[A-Za-zΑ-Ωα-ω_][\w₀-₉]*[ \t]*=[ \t]*"
                + QNODE + r"\." + re.escape(src) + r"[ \t]*\n",
                "", new_body, flags=re.M)
            new_body = re.sub(
                r"(?<![\w.])" + QNODE + r"\." + re.escape(src) + r"(?![\w])",
                kw, new_body)

        lines[start:end + 1] = [new_header] + new_body.split("\n")
        changed.append((header.strip()[:60], kws))

    if not dry_run:
        open(path, "w", encoding="utf-8").write("\n".join(lines))
    return list(reversed(changed))


if __name__ == "__main__":
    dry = "--dry-run" in sys.argv
    for path in [a for a in sys.argv[1:] if not a.startswith("--")]:
        ch = convert(path, dry_run=dry)
        print(f"{path}: {len(ch)} operators")
        for h, kws in ch:
            print(f"    {h:62} ; {kws}")
