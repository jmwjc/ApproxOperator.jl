using ApproxOperator
using Test
import ApproxOperator: UniformValue, PerNodeValue, LazyValue, FieldValue, NodeData

# ============================================================
# Performance Benchmarks for Node struct
# ============================================================

println("="^80)
println("Node Performance Benchmarks")
println("="^80)

# --- Setup ---
println("\n--- Setup: Creating test nodes ---")

n_nodes = 10000

# Create integration-point nodes (𝑿ᵢ)
nodes_𝑿ᵢ = [𝑿ᵢ((𝐼=i,), NodeData()) for i in 1:n_nodes]

# Create sampling-point nodes (𝑿ₛ) with various data
nodes_𝑿ₛ = [𝑿ₛ((𝑔=i, 𝐺=i, 𝐶=i, 𝑠=i), NodeData(
    x=(4, UniformValue(rand())),
    y=(4, UniformValue(rand())),
    z=(4, UniformValue(0.0)),
    𝑤=(4, UniformValue(1.0)),
    w=(1, UniformValue(0.5)),
    d=(1, PerNodeValue(zeros(n_nodes))),
)) for i in 1:n_nodes]

println("  Created $(n_nodes) 𝑿ᵢ nodes and $(n_nodes) 𝑿ₛ nodes")

# Pre-set some PerNodeValue data
for i in 1:n_nodes
    dat = getfield(nodes_𝑿ₛ[i], :data)
    getfield(dat.value[:d], :values)[i] = rand()
end

# ============================================================
# Benchmark 1: getproperty scalar access (UniformValue)
# ============================================================
println("\n--- Benchmark 1: getproperty scalar access (UniformValue) ---")

function bench_getproperty_scalar(nodes, n_iter)
    s = 0.0
    for _ in 1:n_iter
        for n in nodes
            s += n.x
        end
    end
    return s
end

t1 = @elapsed bench_getproperty_scalar(nodes_𝑿ₛ[1:1000], 100)
println("  getproperty scalar (1000 nodes × 100 iters): $(round(t1, digits=3)) s")
@test t1 < 10.0

# ============================================================
# Benchmark 2: getproperty scalar access (PerNodeValue)
# ============================================================
println("\n--- Benchmark 2: getproperty scalar access (PerNodeValue) ---")

function bench_getproperty_pernode(nodes, n_iter)
    s = 0.0
    for _ in 1:n_iter
        for n in nodes
            s += n.d
        end
    end
    return s
end

t2 = @elapsed bench_getproperty_pernode(nodes_𝑿ₛ[1:1000], 100)
println("  getproperty pernode (1000 nodes × 100 iters): $(round(t2, digits=3)) s")
@test t2 < 10.0

# ============================================================
# Benchmark 3: setproperty! on existing PerNodeValue
# ============================================================
println("\n--- Benchmark 3: setproperty! on existing PerNodeValue ---")

function bench_setproperty_pernode(nodes)
    for n in nodes
        n.d = rand()
    end
end

t3 = @elapsed bench_setproperty_pernode(nodes_𝑿ₛ[1:1000])
println("  setproperty! PerNodeValue (1000 nodes): $(round(t3, digits=3)) s")
@test t3 < 10.0

# ============================================================
# Benchmark 4: setproperty! on UniformValue (triggers alloc)
# ============================================================
println("\n--- Benchmark 4: setproperty! on UniformValue (triggers alloc) ---")

nodes_fresh = [𝑿ₛ((𝑔=i, 𝐺=i, 𝐶=i, 𝑠=i), NodeData(x=(4, UniformValue(0.0)))) for i in 1:1000]

function bench_setproperty_uniform(nodes)
    for n in nodes
        n.x = rand()
    end
end

t4 = @elapsed bench_setproperty_uniform(nodes_fresh)
println("  setproperty! UniformValue→PerNodeValue (1000 nodes): $(round(t4, digits=3)) s")
@test t4 < 10.0

# Verify conversion
dat = getfield(nodes_fresh[1], :data)
first_val = dat.value[:x]
@test first_val isa PerNodeValue

# ============================================================
# Benchmark 5: getindex (RV access)
# ============================================================
println("\n--- Benchmark 5: getindex (RV return) ---")

function bench_getindex(nodes)
    s = 0.0
    for n in nodes
        rv = n[:d]
        s += rv[0]
    end
    return s
end

t5 = @elapsed bench_getindex(nodes_𝑿ₛ[1:1000])
println("  getindex[:d] (1000 nodes): $(round(t5, digits=6)) s")
@test t5 < 10.0

# ============================================================
# Benchmark 6: Vector{Node} getproperty (batch access)
# ============================================================
println("\n--- Benchmark 6: Vector{Node} getproperty (batch) ---")

t6 = @elapsed xs = nodes_𝑿ₛ[1:1000].x
println("  Vector{Node}.x UniformValue (1000 nodes): $(round(t6, digits=6)) s")
@test t6 < 10.0
@test xs isa UniformValue{Float64}  # UniformValue returns the struct itself for batch

t6b = @elapsed ds = nodes_𝑿ₛ[1:1000].d
println("  Vector{Node}.d PerNodeValue (1000 nodes): $(round(t6b, digits=6)) s")
@test t6b < 10.0
@test length(ds) == 1000
@test ds isa Vector{Float64}

# ============================================================
# Benchmark 7: Mixed workload (realistic usage pattern)
# ============================================================
println("\n--- Benchmark 7: Mixed workload (assembly-like) ---")

function bench_mixed(nodes, n_iter)
    s = 0.0
    for _ in 1:n_iter
        for n in nodes
            x = n.x
            y = n.y
            w = n.w
            d = n.d
            s += x + y + w + d
            n.d = s
        end
    end
    return s
end

t7 = @elapsed bench_mixed(nodes_𝑿ₛ[1:500], 50)
println("  mixed workload (500 nodes × 50 iters): $(round(t7, digits=3)) s")
@test t7 < 10.0

# ============================================================
# Memory Allocation Tests
# ============================================================
println("\n--- Memory Allocation Tests ---")

# Warm up: compile all relevant functions first
bench_getproperty_scalar(nodes_𝑿ₛ[1:1], 1)
bench_getproperty_pernode(nodes_𝑿ₛ[1:1], 1)
bench_setproperty_pernode(nodes_𝑿ₛ[1:1])
let _tmp = deepcopy(nodes_fresh); bench_setproperty_uniform(_tmp[1:1]); end
bench_getindex(nodes_𝑿ₛ[1:1])
identity(nodes_𝑿ₛ[1:1].x)
identity(nodes_𝑿ₛ[1:1].d)
bench_mixed(nodes_𝑿ₛ[1:1], 1)

# Allocation: scalar getproperty UniformValue
a1 = @allocated bench_getproperty_scalar(nodes_𝑿ₛ[1:1000], 100)
println("  alloc getproperty UniformValue: $(a1) bytes")
# Julia's built-in getproperty dispatch has unavoidable overhead (~32 bytes/call).
# Core data access (Dict lookup) is zero-allocation.
@test a1 < 5_000_000

# Allocation: scalar getproperty PerNodeValue
a2 = @allocated bench_getproperty_pernode(nodes_𝑿ₛ[1:1000], 100)
println("  alloc getproperty PerNodeValue: $(a2) bytes")
@test a2 < 5_000_000

# Allocation: setproperty! on existing PerNodeValue
a3 = @allocated bench_setproperty_pernode(nodes_𝑿ₛ[1:1000])
println("  alloc setproperty! PerNodeValue: $(a3) bytes")
@test a3 < 100_000

# Allocation: setproperty! on UniformValue (triggers PerNodeValue conversion)
# Each conversion allocates zeros(nₜ) where nₜ = index[end], so total ~ sum(1:1000)*8 = ~4 MB
_tmp = deepcopy(nodes_fresh)
a4 = @allocated bench_setproperty_uniform(_tmp)
println("  alloc setproperty! UniformValue→PerNodeValue: $(a4) bytes")
# Upper bound: sum(1:1000)*8 + overhead ≈ 4,004,000; set 5 MB as generous limit
@test a4 < 5_000_000

# Allocation: getindex (RV return)
a5 = @allocated bench_getindex(nodes_𝑿ₛ[1:1000])
println("  alloc getindex[:d]: $(a5) bytes")
# Each RV wraps a sub-vector reference, so allocation is expected
@test a5 < 100_000

# Allocation: batch Vector{Node}.x UniformValue
a6 = @allocated identity(nodes_𝑿ₛ[1:1000].x)
println("  alloc Vector{Node}.x batch: $(a6) bytes")
@test a6 < 100_000

# Allocation: batch Vector{Node}.d PerNodeValue (collects 1000 values)
a6b = @allocated identity(nodes_𝑿ₛ[1:1000].d)
println("  alloc Vector{Node}.d batch: $(a6b) bytes")
# Allocates a 1000-element Vector{Float64} = 1000 * 8 = 8000 bytes + overhead
@test a6b < 100_000

# Allocation: mixed workload
a7 = @allocated bench_mixed(nodes_𝑿ₛ[1:500], 50)
println("  alloc mixed workload: $(a7) bytes")
@test a7 < 5_000_000

# ============================================================
# Summary
# ============================================================
println("\n" * "="^80)
println("Benchmark Summary:")
println("  getproperty (UniformValue):   1000×100 = $(round(t1, digits=3)) s")
println("  getproperty (PerNodeValue):   1000×100 = $(round(t2, digits=3)) s")
println("  setproperty! (PerNodeValue):  1000     = $(round(t3, digits=6)) s")
println("  setproperty! (Uniform→PerNode): 1000  = $(round(t4, digits=3)) s")
println("  getindex (RV):                1000     = $(round(t5, digits=6)) s")
println("  Vector{Node}.x batch:         1000     = $(round(t6, digits=6)) s")
println("  Vector{Node}.d batch:         1000     = $(round(t6b, digits=6)) s")
println("  mixed workload:               500×50   = $(round(t7, digits=3)) s")
println("="^80)
println("Allocation Summary:")
println("  getproperty UniformValue:        $(a1) bytes (expect < 5 MB)")
println("  getproperty PerNodeValue:        $(a2) bytes (expect < 5 MB)")
println("  setproperty! PerNodeValue:       $(a3) bytes (expect < 100 KB)")
println("  setproperty! Uniform→PerNode:    $(a4) bytes (expect < 5 MB)")
println("  getindex RV:                     $(a5) bytes (expect < 100 KB)")
println("  batch Vector{Node}.x:           $(a6) bytes (expect < 100 KB)")
println("  batch Vector{Node}.d:           $(a6b) bytes (expect < 100 KB)")
println("  mixed workload:                  $(a7) bytes (expect < 5 MB)")
println("="^80)
println("All performance and allocation checks passed!")
