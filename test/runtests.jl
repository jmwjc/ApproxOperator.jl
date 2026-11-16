# Test suite for ApproxOperator.jl.
#
# Run with `julia --project=. -e 'using Pkg; Pkg.test()'` or, inside the test
# directory, `julia --project=.. runtests.jl`.
#
# The suite is intentionally organised from "does the pipeline work at all"
# (patch test) to "is each shape function mathematically consistent"
# (consistency tests). Every test file must be listed here: a test that is not
# included is a test that will silently rot.

using Test

@testset "ApproxOperator.jl" begin
    include("patchtest.jl")
    include("importer_characterisation.jl")
    include("shape_function_consistency.jl")
    include("operators_material_parameter.jl")
    include("restored_modules.jl")
    include("package_hygiene.jl")
end
