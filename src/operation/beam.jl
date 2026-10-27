module Beam

using ..ApproxOperator: AbstractElement
const 𝐈₃ = [1.0 0.0 0.0;
            0.0 1.0 0.0;
            0.0 0.0 1.0]

# ===================== 算子命名约定 =====================
# 与 thin_shell.jl 的 ∫εᵢⱼNᵢⱼdΩ 风格一致，用广义应变/广义内力命名：
#
#   ∫εᵘNᵘdΩ  → Kᵘᵘ：位移-位移块
#   ∫εᵘNᵠdΩ  → Kᵘᵠ：位移-转角块（行 u，列 φ）
#   ∫εᵠNᵘdΩ  → Kᵠᵘ：转角-位移块（行 φ，列 u）
#   ∫εᵠNᵠdΩ  → Kᵠᵠ：转角-转角块
#
# 上标 u、φ 表示该块对应的自由度块，不是某个具体的应变分量。
# 每个算子只装配一个 3nₚ×3nₚ 的子块；整体刚度矩阵由使用者在算例中拼装：
#
#   K = [Kᵘᵘ Kᵘᵠ;
#        Kᵠᵘ Kᵠᵠ]
#
# 局部 3×3 块 → 整体 3×3 块：Kᵍ = Rᵀ Kˡ R
function 𝐑ᵀ𝐊𝐑(𝐊::AbstractMatrix{Float64},𝐑::AbstractMatrix{Float64})
    𝐊ᵍ = zeros(3,3)
    for i in 1:3, j in 1:3
        for a in 1:3, b in 1:3
            𝐊ᵍ[i,j] += 𝐑[a,i]*𝐊[a,b]*𝐑[b,j]
        end
    end
    return 𝐊ᵍ
end

function 𝐑(ap::T) where T<:AbstractElement
    𝓒 = ap.𝓒
    x₁ = 𝓒[1]
    x₂ = 𝓒[2]
    L = sqrt((x₂.x-x₁.x)^2 + (x₂.y-x₁.y)^2 + (x₂.z-x₁.z)^2)
    l₁ = (x₂.x-x₁.x)/L
    l₂ = (x₂.y-x₁.y)/L
    l₃ = (x₂.z-x₁.z)/L
    s = sqrt(l₁^2+l₂^2)
    if s > 1.0e-12
        m₁ = -l₁*l₃/s
        m₂ = -l₂*l₃/s
        m₃ = s
        n₁ = l₂/s
        n₂ = -l₁/s
        n₃ = 0.0
    else
        # eₓ 平行于整体 Z 轴，改用整体 X 轴作为参考向量
        m₁ = 1.0
        m₂ = 0.0
        m₃ = 0.0
        n₁ = 0.0
        n₂ = l₃
        n₃ = 0.0
    end
    return [l₁ l₂ l₃;
            m₁ m₂ m₃;
            n₁ n₂ n₃]
end

function 𝐐(ap::T) where T<:AbstractElement
    𝐑ₑ = 𝐑(ap)
    𝐐ᵢ = zeros(6,6)
    for i in 1:3, a in 1:3
        𝐐ᵢ[i,a] = 𝐑ₑ[a,i]
        𝐐ᵢ[i+3,a+3] = 𝐑ₑ[a,i]
    end
    return 𝐐ᵢ
end

# ===================== Kᵘᵘ =====================
function ∫εᵘNᵘdΩ(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝐊ᵘᵘ!(ap,k,𝐑(ap))
end

function 𝐊ᵘᵘ!(ap::T,k::AbstractMatrix{Float64},𝐑::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B = ξ[:∂𝝭∂x]
        𝑤 = ξ.𝑤
        E = ξ.E
        G = ξ.G
        A = ξ.A
        k₂ = ξ.k₂
        k₃ = ξ.k₃
        # Kᵘᵘ = BB·𝐊，𝐊 = diag(EA, k₂GA, k₃GA)
        𝐊 = 𝐑ᵀ𝐊𝐑([E*A      0.0      0.0;
                  0.0    k₂*G*A     0.0;
                  0.0      0.0    k₃*G*A],𝐑)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                𝐵𝐵 = B[i]*B[j]*𝑤
                for a in 1:3, b in 1:3
                    k[3*I-3+a,3*J-3+b] += 𝐊[a,b]*𝐵𝐵
                end
            end
        end
    end
end

# ===================== Kᵘᵠ =====================
function ∫εᵘNᵠdΩ(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝐊ᵘᵠ!(ap,k,𝐑(ap))
end

function 𝐊ᵘᵠ!(ap::T,k::AbstractMatrix{Float64},𝐑::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B = ξ[:∂𝝭∂x]
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        E = ξ.E
        G = ξ.G
        A = ξ.A
        Iy = ξ.Iy
        Iz = ξ.Iz
        k₂ = ξ.k₂
        k₃ = ξ.k₃
        # Kᵘᵠ = BB·𝐊¹ + BN·𝐊²
        𝐊¹ = 𝐑ᵀ𝐊𝐑([0.0    E*Iy   -E*Iz;
                  -G*Iy     0.0     0.0;
                   G*Iz     0.0     0.0],𝐑)
        𝐊² = 𝐑ᵀ𝐊𝐑([   0.0      0.0       0.0;
                       0.0      0.0   -k₂*G*A;
                       0.0   k₃*G*A       0.0],𝐑)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                𝐵𝐵 = B[i]*B[j]*𝑤
                𝐵𝑁 = B[i]*N[j]*𝑤
                for a in 1:3, b in 1:3
                    k[3*I-3+a,3*J-3+b] += 𝐵𝐵*𝐊¹[a,b] + 𝐵𝑁*𝐊²[a,b]
                end
            end
        end
    end
end

# ===================== Kᵠᵘ = (Kᵘᵠ)ᵀ =====================
function ∫εᵠNᵘdΩ(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝐊ᵠᵘ!(ap,k,𝐑(ap))
end

function 𝐊ᵠᵘ!(ap::T,k::AbstractMatrix{Float64},𝐑::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B = ξ[:∂𝝭∂x]
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        E = ξ.E
        G = ξ.G
        A = ξ.A
        Iy = ξ.Iy
        Iz = ξ.Iz
        k₂ = ξ.k₂
        k₃ = ξ.k₃
        # Kᵠᵘ = BB·𝐊¹ᵀ + NB·𝐊²ᵀ
        𝐊¹ = 𝐑ᵀ𝐊𝐑([0.0    E*Iy   -E*Iz;
                  -G*Iy     0.0     0.0;
                   G*Iz     0.0     0.0],𝐑)
        𝐊² = 𝐑ᵀ𝐊𝐑([   0.0      0.0       0.0;
                       0.0      0.0   -k₂*G*A;
                       0.0   k₃*G*A       0.0],𝐑)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                𝐵𝐵 = B[i]*B[j]*𝑤
                𝑁𝐵 = N[i]*B[j]*𝑤
                for a in 1:3, b in 1:3
                    k[3*I-3+a,3*J-3+b] += 𝐵𝐵*𝐊¹[b,a] + 𝑁𝐵*𝐊²[b,a]
                end
            end
        end
    end
end

# ===================== Kᵠᵠ =====================
function ∫εᵠNᵠdΩ(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝐊ᵠᵠ!(ap,k,𝐑(ap))
end

function 𝐊ᵠᵠ!(ap::T,k::AbstractMatrix{Float64},𝐑::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        B = ξ[:∂𝝭∂x]
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        E = ξ.E
        G = ξ.G
        A = ξ.A
        Iy = ξ.Iy
        Iz = ξ.Iz
        Iyy = ξ.Iyy
        Izz = ξ.Izz
        Iyz = ξ.Iyz
        k₂ = ξ.k₂
        k₃ = ξ.k₃
        # Kᵠᵠ = BB·𝐊¹ + BN·𝐊² + NB·𝐊²ᵀ + NN·𝐊³
        #   𝐊¹：扭转 G(Iyy+Izz) 与弯曲 EIyy、EIzz、EIyz
        #   𝐊²：扭转-弯曲翘曲耦合 GIz、GIy
        #   𝐊³：剪切柔度 k₃GA(φ₂)、k₂GA(φ₃)
        𝐊¹ = 𝐑ᵀ𝐊𝐑([G*(Iyy+Izz)     0.0       0.0;
                        0.0       E*Iyy    -E*Iyz;
                        0.0      -E*Iyz     E*Izz],𝐑)
        𝐊² = 𝐑ᵀ𝐊𝐑([0.0   G*Iz   G*Iy;
                    0.0    0.0    0.0;
                    0.0    0.0    0.0],𝐑)
        𝐊³ = 𝐑ᵀ𝐊𝐑([   0.0       0.0       0.0;
                       0.0    k₃*G*A      0.0;
                       0.0       0.0    k₂*G*A],𝐑)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                𝐵𝐵 = B[i]*B[j]*𝑤
                𝐵𝑁 = B[i]*N[j]*𝑤
                𝑁𝐵 = N[i]*B[j]*𝑤
                𝑁𝑁 = N[i]*N[j]*𝑤
                for a in 1:3, b in 1:3
                    k[3*I-3+a,3*J-3+b] += 𝐵𝐵*𝐊¹[a,b] + 𝐵𝑁*𝐊²[a,b] +
                                          𝑁𝐵*𝐊²[b,a] + 𝑁𝑁*𝐊³[a,b]
                end
            end
        end
    end
end

# 罚函数施加本质边界条件（位移部分）：∫α uᵢ(gᵢ−uᵢ)dΓ
# 只写入位移块 Kᵘᵘ（3nₚ×3nₚ）与位移荷载 fᵘ（3nₚ）
function ∫uᵢgᵢdΓ(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        α = ξ.α
        g₁ = ξ.g₁; g₂ = ξ.g₂; g₃ = ξ.g₃
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += α*N[i]*N[j]*𝑤
                k[3*I-1,3*J-1] += α*N[i]*N[j]*𝑤
                k[3*I,3*J]     += α*N[i]*N[j]*𝑤
            end
            f[3*I-2] += α*N[i]*g₁*𝑤
            f[3*I-1] += α*N[i]*g₂*𝑤
            f[3*I]   += α*N[i]*g₃*𝑤
        end
    end
end

# 罚函数施加本质边界条件（转角部分）：∫α φᵢ(gᵢ−φᵢ)dΓ
# 只写入转角块 Kᵠᵠ（3nₚ×3nₚ）与转角荷载 fᵠ（3nₚ）
function ∫φᵢgᵢdΓ(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        α = ξ.α
        g₄ = ξ.g₄; g₅ = ξ.g₅; g₆ = ξ.g₆
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += α*N[i]*N[j]*𝑤
                k[3*I-1,3*J-1] += α*N[i]*N[j]*𝑤
                k[3*I,3*J]     += α*N[i]*N[j]*𝑤
            end
            f[3*I-2] += α*N[i]*g₄*𝑤
            f[3*I-1] += α*N[i]*g₅*𝑤
            f[3*I]   += α*N[i]*g₆*𝑤
        end
    end
end

# ===================== 域内 L₂ 误差（6 个自由度分量）=====================
# 用法（与 sample.jl 的 L₂w、L₂φ 相同的套路）：
#   1. 用高阶积分生成误差单元：
#        elements = getElements(nodes, entities["Ω"], 10)
#   2. 把精确场 prescribe! 到积分点上（函数或常数均可）：
#        prescribe!(elements, :u₁=>u₁, :u₂=>u₂, :u₃=>u₃,
#                             :φ₁=>φ₁, :φ₂=>φ₂, :φ₃=>φ₃)
#   3. 计算形函数：
#        set𝝭!(elements)
#   4. 把有限元解写回节点（d₁..d₃ = u₁..u₃，d₄..d₆ = φ₁..φ₃）：
#        push!(nodes, :d₁=>[...], ..., :d₆=>[...])
#   5. 调用：
#        Δ², d̄² = L₂d(elements)
#      Δ²[c] = ∫(dᶜ_FE − dᶜ_exact)² dΩ，d̄²[c] = ∫(dᶜ_exact)² dΩ；
#      各分量相对 L₂ 误差为 sqrt.(Δ² ./ d̄²)。
function L₂d(ap::T) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    Δ² = zeros(6)
    d̄² = zeros(6)
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        d̄ = [ξ.u₁, ξ.u₂, ξ.u₃, ξ.φ₁, ξ.φ₂, ξ.φ₃]   # 精确解（积分点）
        d  = zeros(6)                              # 有限元解 Σ N_I d_I
        for (i,xᵢ) in enumerate(𝓒)
            d[1] += N[i]*xᵢ.d₁
            d[2] += N[i]*xᵢ.d₂
            d[3] += N[i]*xᵢ.d₃
            d[4] += N[i]*xᵢ.d₄
            d[5] += N[i]*xᵢ.d₅
            d[6] += N[i]*xᵢ.d₆
        end
        Δ² .+= (d .- d̄).^2 .* 𝑤
        d̄² .+= d̄.^2 .* 𝑤
    end
    return Δ², d̄²
end

function L₂d(aps::Vector{T}) where T<:AbstractElement
    Δ² = zeros(6)
    d̄² = zeros(6)
    for ap in aps
        Δ²ₑ, d̄²ₑ = L₂d(ap)
        Δ² .+= Δ²ₑ
        d̄² .+= d̄²ₑ
    end
    return Δ², d̄²
end

end
