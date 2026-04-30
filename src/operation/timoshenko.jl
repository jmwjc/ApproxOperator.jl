module Timoshenko

using ..ApproxOperator: AbstractElement

function ∫κEIκds(ap::T, k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        B = ξ[:∂𝝭∂x]
        𝑤 = ξ.𝑤
        EI = ξ.EI
        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j, xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2 * I, 2 * J] += EI * B[i] * B[j] * 𝑤
            end
        end
    end
end

# Combined form: 對應剪切能 kAG (w_x - φ)^2 展開後的 4 個 block
function ∫γkGAγds(ap::T, k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        B = ξ[:∂𝝭∂x]
        𝑤 = ξ.𝑤
        kGA = ξ.kGA
        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j, xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2 * I - 1, 2 * J - 1] += kGA * B[i] * B[j] * 𝑤
                k[2 * I - 1, 2 * J] -= kGA * B[i] * N[j] * 𝑤
                k[2 * I, 2 * J - 1] -= kGA * N[i] * B[j] * 𝑤
                k[2 * I, 2 * J] += kGA * N[i] * N[j] * 𝑤
            end
        end
    end
end

# Combined form: 對應外載 work 項 ∫ v q dx (只作用在 w 自由度)
function ∫vqds(ap::T, f::AbstractVector) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        q = ξ.q
        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2 * I - 1] += N[i] * q * 𝑤
        end
    end
end

# Block form: 對應 K_φφ 中的彎曲部分 ∫ EI N_I,x N_J,x dx
function ∫κκdΩ(ap::T, k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        B = ξ[:∂𝝭∂x]
        𝑤 = ξ.𝑤
        EI = ξ.EI
        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j, xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I, J] += EI * B[i] * B[j] * 𝑤
            end
        end
    end
end

# Block form: 對應 K_φφ 中的剪切部分 ∫ kAG N_I N_J dx
function ∫φφdΩ(ap::T, k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        kGA = ξ.kGA
        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j, xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I, J] += kGA * N[i] * N[j] * 𝑤
            end
        end
    end
end

# Block form: 對應 K_ww = ∫ kAG N_I,x N_J,x dx
function ∫wwdΩ(ap::T, k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        B = ξ[:∂𝝭∂x]
        𝑤 = ξ.𝑤
        kGA = ξ.kGA
        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j, xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I, J] += kGA * B[i] * B[j] * 𝑤
            end
        end
    end
end

# Geometric stiffness (Engesser), 1D beam specialization of
# K^g_{KL} = P_{αβ} ∫_Ω N_{K,α} N_{L,β} dΩ.
# For a beam, α = β = x and P := P_xx, so
# K^g_{IJ} = P ∫_Ω N_{I,x} N_{J,x} dΩ.
# This routine assembles only the kernel ∫_Ω N_{I,x} N_{J,x} dΩ;
# the axial-force factor P must be applied outside this function.
function ∫wwGdΩ(ap::T, k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        B = ξ[:∂𝝭∂x]
        𝑤 = ξ.𝑤
        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j, xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I, J] += B[i] * B[j] * 𝑤
            end
        end
    end
end

# Geometric stiffness for 2D Mindlin plate buckling.
# Assembles σᵣₑf,αβ w_,α δw_,β only in the w-w block of [w, φ₁, φ₂].
function ∫wwGdΩ2D(ap::T, k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        σ₁₁ = ξ.σ₁₁
        σ₂₂ = ξ.σ₂₂
        σ₁₂ = ξ.σ₁₂
        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j, xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3 * I - 2, 3 * J - 2] += (
                    σ₁₁ * B₁[i] * B₁[j] +
                    σ₂₂ * B₂[i] * B₂[j] +
                    σ₁₂ * (B₁[i] * B₂[j] + B₂[i] * B₁[j])
                ) * 𝑤
            end
        end
    end
end

# Block form: 對應 K_φw = -∫ kAG N_I N_J,x dx
function ∫φwdΩ(ap::T, k::AbstractMatrix) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        B = ξ[:∂𝝭∂x]
        𝑤 = ξ.𝑤
        kGA = ξ.kGA
        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j, xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I, J] -= kGA * N[i] * B[j] * 𝑤
            end
        end
    end
end

# Block RHS (實作展開): 由 w 方程弱式的外力項離散得到，PDF 未以此函數名單列
function ∫wqdΩ(ap::T, f::AbstractVector) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        q = ξ.q
        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] += N[i] * q * 𝑤
        end
    end
end

# Block RHS (實作展開): 轉角方程的外加源項，屬數值模型擴充，PDF 主線未單列
function ∫φmdΩ(ap::T, f::AbstractVector) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        m = ξ.m
        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] += N[i] * m * 𝑤
        end
    end
end

# Penalty BC on w (數值施加邊界): α∫(w-g)v dΓ，屬實作技巧，PDF 未顯式給出
function ∫αwwdΓ(ap::T, k::AbstractMatrix, f::AbstractVector) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    α = ap.α
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        g = ξ.g
        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j, xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I, J] += α * N[i] * N[j] * 𝑤
            end
            f[I] += α * N[i] * g * 𝑤
        end
    end
end

# Penalty BC on φ (數值施加邊界): α∫(φ-g1)ψ dΓ，屬實作技巧，PDF 未顯式給出
function ∫αφφdΓ(ap::T, k::AbstractMatrix, f::AbstractVector) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    α = ap.α
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        g₁ = ξ.g₁
        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j, xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I, J] += α * N[i] * N[j] * 𝑤
            end
            f[I] += α * N[i] * g₁ * 𝑤
        end
    end
end

# Natural BC on w equation (邊界功項): ∫ v V dΓ，為弱式邊界項的離散寫法
function ∫wVdΓ(ap::T, f::AbstractVector) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        V = ξ.V
        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] += N[i] * V * 𝑤
        end
    end
end

# Natural BC on φ equation (邊界功項): -∫ ψ M dΓ，為弱式邊界項的離散寫法
function ∫φMdΓ(ap::T, f::AbstractVector) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        M = ξ.M
        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] -= N[i] * M * 𝑤
        end
    end
end

# L2 for displacement field u (對應報告表10.1的場誤差口徑)
function L₂(ap::T) where T<:AbstractElement
    Δu² = 0.0
    ū² = 0.0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        ū = ξ.u
        u = 0.0
        for (i, xᵢ) in enumerate(ap.𝓒)
            u += N[i] * xᵢ.d
        end
        Δu² += (u - ū)^2 * 𝑤
        ū² += ū^2 * 𝑤
    end
    return Δu², ū²
end

function L₂(aps::Vector{T}) where T<:AbstractElement
    L₂norm_Δu² = 0.0
    L₂norm_ū² = 0.0
    for ap in aps
        Δu², ū² = L₂(ap)
        L₂norm_Δu² += Δu²
        L₂norm_ū² += ū²
    end
    return (L₂norm_Δu² / L₂norm_ū²)^0.5
end

# L2 for rotation field φ (同樣採場積分正規化口徑)
function L₂φ(ap::T) where T<:AbstractElement
    Δφ² = 0.0
    φ̄² = 0.0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        φ̄ = ξ.φ
        φ = 0.0
        for (i, xᵢ) in enumerate(ap.𝓒)
            φ += N[i] * xᵢ.φ
        end
        Δφ² += (φ - φ̄)^2 * 𝑤
        φ̄² += φ̄^2 * 𝑤
    end
    return Δφ², φ̄²
end

function L₂φ(aps::Vector{T}) where T<:AbstractElement
    L₂norm_Δφ² = 0.0
    L₂norm_φ̄² = 0.0
    for ap in aps
        Δφ², φ̄² = L₂φ(ap)
        L₂norm_Δφ² += Δφ²
        L₂norm_φ̄² += φ̄²
    end
    return (L₂norm_Δφ² / L₂norm_φ̄²)^0.5
end

# SS exact displacement: w = w_b + w_s (報告 2.2.4(b), SS)
function w_exact_ss(x::Float64, E::Float64, I::Float64, κ::Float64, G::Float64, A::Float64, L::Float64, q::Float64)
    w_b = q * x * (L^3 - 2.0 * L * x^2 + x^3) / (24.0 * E * I)
    w_s = q * x * (L - x) / (2.0 * κ * G * A)
    return w_b + w_s
end

# CF exact displacement: w = w_b + w_s (報告 2.2.4(b), CF)
function w_exact_cf(x::Float64, E::Float64, I::Float64, κ::Float64, G::Float64, A::Float64, L::Float64, q::Float64)
    w_b = q * x^2 * (6.0 * L^2 - 4.0 * L * x + x^2) / (24.0 * E * I)
    w_s = q * x * (2.0 * L - x) / (2.0 * κ * G * A)
    return w_b + w_s
end

# SS exact rotation φ(x): 對應 SS 位移解的轉角閉式
function φ_exact_ss(x::Float64, E::Float64, I::Float64, L::Float64, q::Float64)
    return q * (L^3 - 6.0 * L * x^2 + 4.0 * x^3) / (24.0 * E * I)
end

# CF exact rotation φ(x): 對應 CF 位移解的轉角閉式
function φ_exact_cf(x::Float64, E::Float64, I::Float64, L::Float64, q::Float64)
    return q * x * (3.0 * L^2 - 3.0 * L * x + x^2) / (6.0 * E * I)
end

end
