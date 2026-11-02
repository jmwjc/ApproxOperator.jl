#!/usr/bin/env python3
"""
S3 analysis: compare the generated `getElements`-family variants.

Each variant is written as

    @eval function NAME(...)
        $prequote
        ...
        <some of the quote fragments>
        ...
    end

so the interesting part is *which fragments each variant includes*, in which
order, with which arguments. This script prints that matrix.
"""

import re
import sys

FRAGMENTS = [
    "prequote", "preForEdge", "typeForFEM", "integrationByGmsh",
    "integrationByManual", "coordinates", "coordinatesForEdges",
    "curvilinearCoordinates", "cal_jacobe", "cal_normal",
    "generateForFEM", "generateForNeighbor", "generateForMarco",
    "generateForPiecewise", "generateForPiecewiseBoundary", "generateSummary",
]


def variants(path):
    lines = open(path, encoding="utf-8").read().split("\n")
    i = 0
    while i < len(lines):
        if not lines[i].startswith("function get"):
            i += 1
            continue
        sig = lines[i]
        start = i
        depth = 0
        j = i
        while j < len(lines):
            code = lines[j].strip().split("#")[0]
            if re.match(r"^(function|for|while|if|try|let|begin)\b", code):
                depth += 1
            elif re.match(r"^(else|elseif)\b", code):
                pass
            elif re.match(r"^end\b", code):
                depth -= 1
                if depth == 0:
                    break
            j += 1
        body = "\n".join(lines[start:j + 1])
        used = [f for f in FRAGMENTS if "$" + f in body]
        yield sig, used
        i = j + 1


def main(path):
    rows = list(variants(path))
    name_w = max(len(s.split("(")[0].replace("function ", "")) for s, _ in rows)
    print(f"{'variant':<{name_w}}  fragments")
    print("-" * (name_w + 60))
    for sig, used in rows:
        name = sig.split("(")[0].replace("function ", "")
        extra = sig.split("(")[1].split(")")[0]
        print(f"{name:<{name_w}}  {', '.join(used)}")
        print(f"{'':<{name_w}}  args: {extra}")
    print()
    # Which fragments are used by nobody / by everybody?
    counts = {f: sum(1 for _, u in rows if f in u) for f in FRAGMENTS}
    print("fragment usage:")
    for f, c in counts.items():
        mark = "" if 0 < c < len(rows) else ("  <-- always" if c == len(rows) else "  <-- unused")
        print(f"  {f:<32} {c}/{len(rows)}{mark}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "src/preprocession/importmsh.jl")
