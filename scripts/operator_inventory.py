#!/usr/bin/env python3
"""
S2.7 helper: for each operator, report the quadrature-point fields it reads and
the keyword arguments it declares.

The "requires" list is what tells a reader whether `set𝝭!`, `set∇𝝭!` or
`set∇²𝝭!` must be called before assembling, which is currently discoverable
only by reading the body.
"""

import re
import sys
import glob

# Quadrature-point field names used by `set𝝭!` and friends.
SETTER = {
    ":𝝭": "set𝝭!",
    ":∂𝝭∂x": "set∇𝝭!",
    ":∂𝝭∂y": "set∇𝝭!",
    ":∂𝝭∂z": "set∇𝝭!",
    ":∂²𝝭∂x²": "set∇²𝝭!",
    ":∂²𝝭∂x∂y": "set∇²𝝭!",
    ":∂²𝝭∂y²": "set∇²𝝭!",
    ":∂³𝝭∂x³": "set∇³𝝭!",
    ":∂³𝝭∂x²∂y": "set∇³𝝭!",
    ":∂³𝝭∂x∂y²": "set∇³𝝭!",
    ":∂³𝝭∂y³": "set∇³𝝭!",
}


def find_operators(text):
    lines = text.split("\n")
    i = 0
    while i < len(lines):
        if not lines[i].startswith("function "):
            i += 1
            continue
        start = i
        depth = 0
        j = i
        ok = False
        while j < len(lines):
            code = lines[j].strip().split("#")[0]
            if re.match(r"^(function|for|while|if|try|let|begin)\b", code):
                depth += 1
            elif re.match(r"^(else|elseif)\b", code):
                pass
            elif re.match(r"^end\b", code):
                depth -= 1
                if depth == 0:
                    yield start, j
                    ok = True
                    break
            j += 1
        i = j + 1 if ok else i + 1


def keywords_of(header):
    close = header.rfind(")")
    if close < 0 or ";" not in header[:close]:
        return []
    after = header[:close].split(";", 1)[1]
    out = []
    for p in after.split(","):
        p = p.strip()
        if re.fullmatch(r"\w+", p):
            out.append(p)
    return out


def fields_of(body):
    code = "\n".join(re.sub(r"#.*$", "", ln) for ln in body.split("\n"))
    found = [f for f in SETTER if f + "]" in code or f + "]" in code]
    # also catch ξ[:field] spelled with any quadrature-node name
    for m in re.finditer(r"ξ[\w₀-₉ᵘᵖˢʷᵠ]*\[(:[^\]]+)\]", code):
        f = ":" + m.group(1)[1:] if not m.group(1).startswith(":") else m.group(1)
        if f not in found:
            found.append(f)
    return sorted(set(found))


def main(paths):
    for path in paths:
        lines = open(path, encoding="utf-8").read().split("\n")
        print("=" * 78)
        print(path)
        print("=" * 78)
        for start, end in find_operators("\n".join(lines)):
            header = lines[start]
            kws = keywords_of(header)
            if not kws:
                continue
            body = "\n".join(lines[start:end + 1])
            fields = fields_of(body)
            setters = sorted({SETTER[f] for f in fields if f in SETTER})
            print(f"\n  {header[:96]}")
            print(f"    kwargs  : {kws}")
            print(f"    fields  : {fields}")
            print(f"    setters : {setters}")


if __name__ == "__main__":
    args = sys.argv[1:] or sorted(glob.glob("src/operation/*.jl"))
    main(args)
