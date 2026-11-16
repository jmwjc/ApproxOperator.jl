# Package hygiene, checked by Aqua.jl.
#
# These are not numerical tests. They guard the things that decide whether
# someone else can install and trust the package:
#
#   * no `export` of a name that does not exist in the module that exports it;
#   * no ambiguous method dispatch;
#   * no stale dependency (declared but never imported);
#   * `[compat]` covers every dependency and the Julia version;
#   * no type piracy (adding methods to functions or types the package does not own).

using Aqua
using Test
using ApproxOperator

@testset "package hygiene (Aqua)" begin

    @testset "no undefined exports" begin
        # Four importer entry points used to be exported but only defined in the
        # `GmshImport` submodule, so `using ApproxOperator` handed the user
        # nothing; `addEdgeElements` was exported and defined nowhere at all.
        Aqua.test_undefined_exports(ApproxOperator)
    end

    @testset "no method ambiguities" begin
        # `push!(aps::Vector{<:AbstractElement}, ::Symbol...)` and
        # `push!(aps, ::Pair...)` were both variadic, so the empty call
        # `push!(aps)` could not be resolved. Both now require one argument.
        Aqua.test_ambiguities(ApproxOperator)
    end

    @testset "no stale dependencies" begin
        # `WriteVTK` was a direct dependency that `src/` never imported.
        Aqua.test_stale_deps(ApproxOperator)
    end

    @testset "compatibility bounds are declared" begin
        Aqua.test_deps_compat(ApproxOperator)
    end

    @testset "no type piracy" begin
        Aqua.test_piracies(ApproxOperator)
    end

    @testset "no persistent tasks" begin
        Aqua.test_persistent_tasks(ApproxOperator)
    end

    @testset "the importer entry points really are reachable" begin
        # The behavioural counterpart of `test_undefined_exports`: after
        # `using ApproxOperator` these names must resolve in `Main`.
        @test isdefined(ApproxOperator, :getPhysicalGroups)
        @test isdefined(ApproxOperator, :get𝑿ᵢ)
        @test isdefined(ApproxOperator, :getElements)
        @test isdefined(ApproxOperator, :getCurvedElements)
        @test isdefined(ApproxOperator, :getDOFs)
        # Removed: it was exported but never defined.
        @test !isdefined(ApproxOperator, :addEdgeElements)
    end
end
