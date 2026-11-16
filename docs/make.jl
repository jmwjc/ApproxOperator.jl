using Documenter
using ApproxOperator

makedocs(
    sitename = "ApproxOperator.jl",
    modules  = [ApproxOperator],
    format   = Documenter.HTML(prettyurls = false),
    pages = [
        "Home"             => "index.md",
        "Operator index"   => "operators.md",
    ],
    # The package has 263 operator docstrings over 14 modules; a warning about a
    # missing cross-reference is not actionable, so only errors are fatal.
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(repo = "github.com/jmwjc/ApproxOperator.jl", push_preview = true)
