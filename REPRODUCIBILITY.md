# Reproducibility

The published results of this package depend on the **mesh generator** and on
the **quadrature rules** Gmsh returns for a given element type. Both come from
`gmsh_jll`, so a different `gmsh_jll` version can produce different integration
weights for the same `.msh` file and the same `integrationOrder`.

`Manifest.toml` is **not committed** (it is in `.gitignore`), following the usual
package convention, because a committed manifest makes `Pkg.develop` / `Pkg.add`
behave unpredictably for anyone using this package as a dependency. This file
records the environment instead.

## Environment used for the results in this repository

| item | value |
|---|---|
| Julia | 1.12.6 |
| `julia` compat declared | `1.9` |
| Gmsh.jl | see `[compat]` `Gmsh = "0.3"` |
| `gmsh_jll` | tree hash in the local `Manifest.toml` |
| platform | macOS arm64 (Apple silicon) |

Reproduce the exact environment with:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.status()
```

and compare the `Pkg.status()` output against the versions above. If
`gmsh_jll` differs, expect the last digits of the numerical regression values
below to differ too.

## Numerical values that must reproduce

These are asserted by the test suite, so a drift is caught automatically.

| quantity | value | asserted in |
|---|---|---|
| patch-test L2 error | `6.594194188924542e-8` | `test/patchtest.jl` |
| 2×2 Quad4 mesh, Σ ξ.𝑤 | `1.0 ± 1e-12` | `test/importer_characterisation.jl` |
| 1-D boundary mesh, Σ ξ.𝑤 | `4.0 ± 1e-12` | `test/importer_characterisation.jl` |
| flat curved-surface map, Σ ξ.𝑤 | `1.0 ± 1e-12` | `test/importer_characterisation.jl` |
| same map scaled by 2, Σ ξ.𝑤 | `2.0 ± 1e-12` | `test/importer_characterisation.jl` |
| element count, patch-test mesh | `288` triangles | `test/patchtest.jl` |

## What is *not* pinned

- the exact Gmsh meshing algorithm for meshes generated inside the tests — the
  patch-test mesh is a committed file, so it is fixed, but meshes built by
  `gmsh.model.mesh.generate` with default algorithms may differ between
  `gmsh_jll` releases;
- BLAS/LAPACK results in `k \ f` beyond the `atol` used in the tests.

If you need bit-level reproducibility for a publication, pin the environment
yourself:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
cp("Manifest.toml", "Manifest-for-paper.toml")
```

and commit `Manifest-for-paper.toml` alongside the paper's results, or install
the recorded `gmsh_jll` version explicitly.
