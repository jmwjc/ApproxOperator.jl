#!/usr/bin/env python3
"""List every operator in a src/operation/*.jl file and the material constants
it reads. Used to drive the S2 conversion (material constants -> keywords)."""
import re
import sys

MAT = ["E", "ν", "h", "D", "ρ", "ρA", "m", "μ", "λ", "σy", "K", "tol",
       "α", "k", "c", "c²", "EA", "EI", "kGA", "R", "l"]

MAT_RE = {m: re.compile(r"(?<![\w.])" + re.escape(m) + r"(?![\w])") for m in MAT}


def strip_comments(text):
    """Drop `#` comments so that names mentioned only in prose are not counted."""
    return "\n".join(re.sub(r"#.*$", "", ln) for ln in text.split("\n"))


def operators(path):
    """Yield (header_line, body_text) for every top-level `function`."""
    lines = open(path, encoding="utf-8").read().split("\n")
    cur, body = None, []
    for ln in lines:
        if ln.startswith("function "):
            if cur is not None:
                yield cur, "\n".join(body)
            cur, body = ln, []
        elif cur is not None:
            body.append(ln)
    if cur is not None:
        yield cur, "\n".join(body)


def main(paths):
    for path in paths:
        print("=" * 78)
        print(path)
        print("=" * 78)
        total = 0
        for hdr, body in operators(path):
            code = strip_comments(body)
            # A name bound by a `for` loop is a loop variable, not material data.
            loopvars = set(re.findall(r"for\s*\(([^)]*)\)", code))
            loopvars = {v.strip() for group in loopvars for v in group.split(",")}
            loopvars |= set(re.findall(r"for\s+(\w+)\s+in", code))
            used = [m for m in MAT if MAT_RE[m].search(code) and m not in loopvars]
            if not used:
                continue
            total += 1
            print(f"  {hdr[:66]:68} -> {used}")
        print(f"  ({total} operators need material constants)\n")


if __name__ == "__main__":
    main(sys.argv[1:])
