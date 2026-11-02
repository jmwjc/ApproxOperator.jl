#!/usr/bin/env python3
"""
S4: generate `docs/operators.md` — a single table of every weak-form operator,
its module, the quadrature-point fields it requires, its material constants and
a one-line description.

The table is produced from the source, not hand-maintained, so it cannot drift
from the code. Run it after changing any operator docstring:

    python3 scripts/build_operator_index.py
"""

import glob
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Quadrature-point field -> the setter that must be called before assembly.
SETTER = {
    ":𝝭": "set𝝭!",
    ":∂𝝭∂x": "set∇𝝭!", ":∂𝝭∂y": "set∇𝝭!", ":∂𝝭∂z": "set∇𝝭!",
    ":∂²𝝭∂x²": "set∇²𝝭!", ":∂²𝝭∂x∂y": "set∇²𝝭!", ":∂²𝝭∂y²": "set∇²𝝭!",
    ":∂³𝝭∂x³": "set∇³𝝭!", ":∂³𝝭∂x²∂y": "set∇³𝝭!",
    ":∂³𝝭∂x∂y²": "set∇³𝝭!", ":∂³𝝭∂y³": "set∇³𝝭!",
}

MODULE_TITLE = {
    "Elasticity": "Elasticity",
    "Heat": "Heat conduction and convection",
    "KirchhoffPlate": "Kirchhoff plate",
    "MindlinPlate": "Mindlin (thick) plate",
    "KirchhoffShell": "Kirchhoff shell",
    "Beam": "3-D beam",
    "Hamilton": "Hamiltonian / space–time rod",
    "Hyperelasticity": "Hyperelasticity",
    "Plasticity": "Plasticity and contact",
    "Phasefield": "Phase-field fracture",
    "Stokes": "Stokes flow",
    "WeightedResidual": "Space–time weighted residual",
    "Test": "Consistency checks",
}


def split_functions(text):
    """Yield (docstring, signature, body) for each top-level function."""
    lines = text.split("\n")
    i = 0
    pending_doc = None
    while i < len(lines):
        if lines[i].strip() == '"""':
            j = i + 1
            buf = []
            while j < len(lines) and lines[j].strip() != '"""':
                buf.append(lines[j])
                j += 1
            pending_doc = "\n".join(buf)
            i = j + 1
            continue
        stripped = lines[i].strip()
        if not stripped.startswith("function"):
            # a blank line between a docstring and a function resets the pairing
            if stripped and stripped != '"""':
                pass
            i += 1
            continue
        doc = pending_doc
        pending_doc = None
        # gather the body
        start = i
        depth = 0
        j = i
        closed = False
        while j < len(lines):
            code = lines[j].strip().split("#")[0]
            if re.match(r"^(function|for|while|if|try|let|begin|struct|quote)\b", code):
                depth += 1
            elif re.match(r"^(else|elseif)\b", code):
                pass
            elif re.match(r"^end\b", code):
                depth -= 1
                if depth == 0:
                    closed = True
                    break
            j += 1
        body = "\n".join(lines[start:j + 1])
        yield doc, lines[start].strip(), body
        i = (j + 1) if closed else (i + 1)


def parse_doc(doc):
    """Extract the one-line summary and the keyword list from a docstring.

    The convention used throughout `src/operation/` is:

        <signature lines, indented by 4>
        <blank>
        <summary paragraph, not indented>
        <blank>
        # Keyword arguments
        - `name` : meaning
    """
    if not doc:
        return "", []
    raw = doc.split("\n")

    def indented(l):
        return l.startswith("    ") or l.startswith("\t")

    i = 0
    # 1) skip the leading signature block (indented lines + blanks)
    while i < len(raw) and (not raw[i].strip() or indented(raw[i])):
        i += 1
    # 2) collect the summary paragraph (unindented, until a blank line)
    parts = []
    while i < len(raw) and raw[i].strip() and not indented(raw[i]):
        parts.append(raw[i].strip())
        i += 1
    summary = " ".join(parts)

    # 3) keyword arguments: "- `name` : meaning" under "# Keyword arguments"
    kws = []
    in_kw = False
    for l in raw:
        t = l.strip()
        if t.lower().startswith("# keyword arguments"):
            in_kw = True
            continue
        if t.startswith("# ") and not t.lower().startswith("# keyword"):
            in_kw = False
        if in_kw:
            m = re.match(r"^-\s*`([^`]+)`\s*:?\s*(.*)$", t)
            if m:
                kws.append(m.group(1))
    return summary, kws


def signature_name(sig):
    m = re.search(r"function (?:\(op::Operator\{:([^}]+)\}\)|([^\s(]+))", sig)
    if not m:
        return None
    return m.group(1) or m.group(2)


def fields_of(body):
    code = "\n".join(re.sub(r"#.*$", "", l) for l in body.split("\n"))
    found = set()
    for m in re.finditer(r"ξ[\w₀-₉ᵘᵖˢʷᵠ]*\[(:[^\]]+)\]", code):
        found.add(m.group(1))
    return sorted(found)


def main():
    rows = []
    for path in sorted(glob.glob(os.path.join(ROOT, "src", "operation", "*.jl"))):
        module = None
        text = open(path, encoding="utf-8").read()
        m = re.search(r"^module\s+(\w+)", text, re.M)
        if m:
            module = m.group(1)
        else:
            module = os.path.basename(path)[:-3].capitalize()
        for doc, sig, body in split_functions(text):
            if re.search(r"\bop::Operator\b", sig):
                continue                # legacy dispatch form, not documented
            name = signature_name(sig)
            if not name:
                continue
            summary, kws = parse_doc(doc)
            fields = fields_of(body)
            setters = sorted({SETTER[f] for f in fields if f in SETTER})
            rows.append(dict(module=module, name=name, summary=summary,
                             kws=kws, setters=setters, documented=bool(doc)))

    by_module = {}
    for r in rows:
        by_module.setdefault(r["module"], []).append(r)

    out = []
    out.append("# Operator index\n")
    out.append("Every weak-form operator in `src/operation/`, generated from the source\n"
               "by `scripts/build_operator_index.py`. Do not edit by hand.\n")
    out.append(f"**{len(rows)} operators** across {len(by_module)} modules.\n")
    out.append("## How to read this table\n")
    out.append("| column | meaning |\n|---|---|")
    out.append("| *operator* | the Julia function to put on the left of `form => elements` |")
    out.append("| *requires* | which shape-function setter must run before assembly |")
    out.append("| *material constants* | keyword arguments the operator takes |")
    out.append("| *description* | the first line of the operator's docstring |\n")
    out.append("```julia\n"
               "set∇𝝭!(elements)                       # whichever the table says\n"
               "(∫∫εᵢⱼσᵢⱼdxdy => elements)(k; E = 210e9, ν = 0.3)\n"
               "```\n")
    out.append("## Contents\n")
    for module in sorted(by_module):
        title = MODULE_TITLE.get(module, module)
        anchor = title.lower().replace(" ", "-").replace("–", "").replace("(", "").replace(")", "")
        out.append(f"- [{title}](#{anchor}) — {len(by_module[module])} operators")
    out.append("")

    for module in sorted(by_module):
        title = MODULE_TITLE.get(module, module)
        out.append(f"## {title}\n")
        out.append(f"`ApproxOperator.{module}`\n")
        out.append("| operator | requires | material constants | description |")
        out.append("|---|---|---|---|")
        merged = {}
        n_over = {}
        for r in by_module[module]:
            key = r["name"]
            n_over[key] = n_over.get(key, 0) + 1
            m = merged.get(key)
            if m is None:
                merged[key] = dict(r)
                continue
            m["setters"] = sorted(set(m["setters"]) | set(r["setters"]))
            m["kws"]     = sorted(set(m["kws"]) | set(r["kws"]))
            # prefer the documented overload's summary
            if not m["documented"] and r["documented"]:
                m["summary"], m["documented"] = r["summary"], True
        for key in sorted(merged):
            r = merged[key]
            req = ", ".join(f"`{s}`" for s in r["setters"]) if r["setters"] else "—"
            kws = ", ".join(f"`{k}`" for k in r["kws"]) if r["kws"] else "—"
            desc = r["summary"].replace("|", "\\|")
            if n_over[key] > 1:
                desc += f" *({n_over[key]} overloads)*"
            out.append(f"| `{r['name']}` | {req} | {kws} | {desc} |")
        out.append("")

    dest = os.path.join(ROOT, "docs", "src", "operators.md")
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    open(dest, "w", encoding="utf-8").write("\n".join(out) + "\n")
    print(f"wrote {dest}  ({len(rows)} operators, {len(by_module)} modules)")


if __name__ == "__main__":
    main()
