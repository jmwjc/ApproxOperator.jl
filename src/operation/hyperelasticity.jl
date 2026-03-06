module Hyperelasticity
    
using ..ApproxOperator: AbstractElement
using LinearAlgebra,Printf
function ∫ESdx(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖 
    Eᵉ = op.E
    for ξ in 𝓖
        B = ξ[:∂𝝭∂x]
        𝑤 = ξ.𝑤
        F = 1.0
        for (i,xᵢ) in enumerate(𝓒)
            F += B[i]*xᵢ.d
        end
        E = 0.5*(F^2-1.0)
        S = Eᵉ*E
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += B[i]*(S+F^2*Eᵉ)*B[j]*𝑤
            end
            f[I] += B[i]*Eᵉ*F*E*𝑤
        end
    end
end


function Δ∫∫δSᵢⱼXᵢⱼdxdy_HR_SaintVenantKirchhoff(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
   
    for ξ in 𝓖
       E = ξ.E
       ν = ξ.ν
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        C⁻¹ᵢᵢᵢᵢ = 1/E
        C⁻¹ᵢᵢⱼⱼ = -ν/E
        C⁻¹ᵢⱼᵢⱼ = 2*(1+ν)/E
      
        
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += N[i]*C⁻¹ᵢᵢᵢᵢ*N[j]*𝑤
                k[3*I-2,3*J-1] += N[i]*C⁻¹ᵢᵢⱼⱼ*N[j]*𝑤
                k[3*I-1,3*J-2] += N[i]*C⁻¹ᵢᵢⱼⱼ*N[j]*𝑤
                k[3*I-1,3*J-1] += N[i]*C⁻¹ᵢᵢᵢᵢ*N[j]*𝑤
                k[3*I,3*J]     += N[i]*C⁻¹ᵢⱼᵢⱼ*N[j]*𝑤

            end
        end
    end
end

function Δ∭EᵢⱼSᵢⱼdxdydz_HR_SaintVenantKirchhoff(ap::T, k::AbstractMatrix{Float64}) where {T<:AbstractElement}
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖

    for ξ in 𝓖
        E = ξ.E
        ν = ξ.ν
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤   

        
        C⁻¹ᵢᵢᵢᵢ = 1/E
        C⁻¹ᵢᵢⱼⱼ = -ν/E
        C⁻¹ᵢⱼᵢⱼ = 2*(1+ν)/E

        for (a, xᵃ) in enumerate(𝓒)
            I = xᵃ.𝐼
            for (b, xᵇ) in enumerate(𝓒)
                J = xᵇ.𝐼

                NiNjw = N[a] * N[b] * 𝑤

                # --- 正应力块 (σxx, σyy, σzz) ---
                # σxx - σxx
                k[6I-5, 6J-5] += N[i]*N[j]*w * C⁻¹ᵢᵢᵢᵢ
                # σxx - σyy
                k[6I-5, 6J-4] += N[i]*N[j]*w * C⁻¹ᵢᵢⱼⱼ
                # σxx - σzz
                k[6I-5, 6J-3] += N[i]*N[j]*w * C⁻¹ᵢᵢⱼⱼ

                # σyy - σxx
                k[6I-4, 6J-5] += N[i]*N[j]*w * C⁻¹ᵢᵢⱼⱼ
                # σyy - σyy
                k[6I-4, 6J-4] += N[i]*N[j]*w * C⁻¹ᵢᵢᵢᵢ
                # σyy - σzz
                k[6I-4, 6J-3] += N[i]*N[j]*w * C⁻¹ᵢᵢⱼⱼ

                # σzz - σxx
                k[6I-3, 6J-5] += N[i]*N[j]*w * C⁻¹ᵢᵢⱼⱼ
                # σzz - σyy
                k[6I-3, 6J-4] += N[i]*N[j]*w * C⁻¹ᵢᵢⱼⱼ
                # σzz - σzz
                k[6I-3, 6J-3] += N[i]*N[j]*w * C⁻¹ᵢᵢᵢᵢ

                # --- 剪应力块 (σxy, σyz, σzx) ---
                # σxy - σxy
                k[6I-2, 6J-2] += N[i]*N[j]*w * C⁻¹ᵢⱼᵢⱼ
                # σyz - σyz
                k[6I-1, 6J-1] += N[i]*N[j]*w * C⁻¹ᵢⱼᵢⱼ
                # σzx - σzx
                k[6I  , 6J  ] += N[i]*N[j]*w * C⁻¹ᵢⱼᵢⱼ

                # 正应力与剪应力之间的项在各向同性材料的 C⁻¹ 中为 0，所以这里不需要加
            end
        end
    end
end


function ∫∫δSᵢⱼXᵢⱼdxdy_HR_SaintVenantKirchhoff(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    
    for ξ in 𝓖
        E = ξ.E
        ν = ξ.ν
        𝑤 = ξ.𝑤
        C⁻¹ᵢᵢᵢᵢ = 1/E
        C⁻¹ᵢᵢⱼⱼ = -ν/E
        C⁻¹ᵢⱼᵢⱼ = 2*(1+ν)/E
       
        N = ξ[:𝝭]
 
       
        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒)
           S₁₁ += N[i]*xᵢ.dₛ₁₁
           S₂₂ += N[i]*xᵢ.dₛ₂₂
           S₁₂ += N[i]*xᵢ.dₛ₁₂
        end
      

        
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[3*I-2] += (N[i]*C⁻¹ᵢᵢᵢᵢ*S₁₁+N[i]*C⁻¹ᵢᵢⱼⱼ*S₂₂)*𝑤
            f[3*I-1] += (N[i]*C⁻¹ᵢᵢⱼⱼ*S₁₁+N[i]*C⁻¹ᵢᵢᵢᵢ*S₂₂)*𝑤
            f[3*I]   += (N[i]*C⁻¹ᵢⱼᵢⱼ*S₁₂ )*𝑤
        end
    end
end

function ∫∫δSᵢⱼEᵢⱼdxdy_HR(aₛ::T,aᵤ::S, f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
   
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        𝑤 = ξₛ.𝑤

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

        # 2. 计算 Green-Lagrange 应变 E = 0.5 * (F'F - I)
        E₁₁ = 0.5 * (F₁₁*F₁₁ + F₂₁*F₂₁ - 1.0)
        E₂₂ = 0.5 * (F₁₂*F₁₂ + F₂₂*F₂₂ - 1.0)
        E₁₂ = 0.5 * (F₁₁*F₁₂ + F₂₁*F₂₂) # 张量剪切应变

        # 3. 组装残差 (注意剪切项的系数 2.0，对应双点积 S:E)
        for (i, xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            # 对应 δS₁₁ 的方程
            f[3*I-2] += N[i] * E₁₁ * 𝑤
            # 对应 δS₂₂ 的方程
            f[3*I-1] += N[i] * E₂₂ * 𝑤
            # 对应 δS₁₂ 的方程 (能量共轭量是 2*E₁₂)
            f[3*I]   += N[i] * 2.0 * E₁₂ * 𝑤
        end
    end
end


function ∫∫δSᵢⱼΔEᵢⱼdxdy_HR(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
   
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        𝑤 = ξₛ.𝑤

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

        # 2. 计算 Green-Lagrange 应变 E = 0.5 * (F'F - I)
        E₁₁ = 0.5 * (F₁₁*F₁₁ + F₂₁*F₂₁ - 1.0)
        E₂₂ = 0.5 * (F₁₂*F₁₂ + F₂₂*F₂₂ - 1.0)
        E₁₂ = 0.5 * (F₁₁*F₁₂ + F₂₁*F₂₂) # 张量剪切应变

       
        for (i, xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[3*I-2, 2*J-1] += N[i] * F₁₁ * B₁[j] * 𝑤  # Col: u
                k[3*I-2, 2*J]   += N[i] * F₂₁ * B₁[j] * 𝑤  # Col: v

                k[3*I-1, 2*J-1] += N[i] * F₁₂ * B₂[j] * 𝑤
                k[3*I-1, 2*J]   += N[i] * F₂₂ * B₂[j] * 𝑤

                k[3*I, 2*J-1] += N[i] * (F₁₂* B₁[j]+F₁₁ * B₂[j])* 𝑤
                k[3*I, 2*J]   += N[i] * (F₂₂* B₁[j]+F₂₁ * B₂[j]) * 𝑤
               



            end
        end
    end
end



function ∫∫SᵢⱼδEᵢⱼdxdy_HR(aₛ::T,aᵤ::S, f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
   
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        𝑤 = ξₛ.𝑤

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

        # 2. 计算 Green-Lagrange 应变 E = 0.5 * (F'F - I)
        E₁₁ = 0.5 * (F₁₁*F₁₁ + F₂₁*F₂₁ - 1.0)
        E₂₂ = 0.5 * (F₁₂*F₁₂ + F₂₂*F₂₂ - 1.0)
        E₁₂ = 0.5 * (F₁₁*F₁₂ + F₂₁*F₂₂) # 张量剪切应变


        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
        end

     # First Piola-Kirchhoff
        P₁₁ = F₁₁ * S₁₁ + F₁₂ * S₁₂
        P₁₂ = F₁₁ * S₁₂ + F₁₂ * S₂₂
        P₂₁ = F₂₁ * S₁₁ + F₂₂ * S₁₂
        P₂₂ = F₂₁ * S₁₂ + F₂₂ * S₂₂

        # 3. 组装残差 (注意剪切项的系数 2.0，对应双点积 S:E)
        for (i, xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
           f[2*I-1] += (P₁₁ * B₁[i] + P₁₂ * B₂[i]) * 𝑤
           f[2*I]   += (P₂₁ * B₁[i] + P₂₂ * B₂[i]) * 𝑤

        end
    end
end


function ∫∫SᵢⱼδΔEᵢⱼdxdy_HR(aₛ::T,aᵤ::S, f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
   
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        
        𝑤 = ξₛ.𝑤

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

        # 2. 计算 Green-Lagrange 应变 E = 0.5 * (F'F - I)
        E₁₁ = 0.5 * (F₁₁*F₁₁ + F₂₁*F₂₁ - 1.0)
        E₂₂ = 0.5 * (F₁₂*F₁₂ + F₂₂*F₂₂ - 1.0)
        E₁₂ = 0.5 * (F₁₁*F₁₂ + F₂₁*F₂₂) # 张量剪切应变


        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
        end

     # First Piola-Kirchhoff
        P₁₁ = F₁₁ * S₁₁ + F₁₂ * S₁₂
        P₁₂ = F₁₁ * S₁₂ + F₁₂ * S₂₂
        P₂₁ = F₂₁ * S₁₁ + F₂₂ * S₁₂
        P₂₂ = F₂₁ * S₁₂ + F₂₂ * S₂₂

        # 3. 组装残差 (注意剪切项的系数 2.0，对应双点积 S:E)
        for (i, xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
            b1_i = B₁[i]
            b2_i = B₂[i]
            for (j, xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                b1_j = B₁[j]
                b2_j = B₂[j]
                term1 = S₁₁ * b1_i * b1_j
                term2 = S₁₂ * (b1_i * b2_j + b2_i * b1_j) # 对称剪切项
                term3 = S₂₂ * b2_i * b2_j
                g = (term1 + term2 + term3) * 𝑤
                k[2*I-1, 2*J-1] += g
                k[2*I,   2*J]   += g


            end

        end
    end
end




function Δ∫∫δSFnudxdy_HR_SaintVenantKirchhoff(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
   
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
       
       
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

        
      for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                # ∫∫δSΔFnudxdy
                k[3*I-2,2*J-1] -= Nₛ[i]*n₁*u₁*B₁[j]*𝑤
                k[3*I-2,2*J]   -= Nₛ[i]*n₁*u₂*B₁[j]*𝑤
                k[3*I-1,2*J-1] -= Nₛ[i]*n₂*u₁*B₂[j]*𝑤
                k[3*I-1,2*J]   -= Nₛ[i]*n₂*u₂*B₂[j]*𝑤
                k[3*I,2*J-1]   -= (Nₛ[i]*n₂*u₁*B₁[j]+Nₛ[i]*n₁*u₁*B₂[j])*𝑤
                k[3*I,2*J]     -= (Nₛ[i]*n₁*u₂*B₂[j]+Nₛ[i]*n₂*u₂*B₁[j])*𝑤
                # ∫∫δSΔFnudxdy
                k[3*I-2,2*J-1] -= Nₛ[i]*n₁*F₁₁*N[j]*𝑤
                k[3*I-2,2*J]   -= Nₛ[i]*n₁*F₂₁*N[j]*𝑤
                k[3*I-1,2*J-1] -= Nₛ[i]*n₂*F₁₂*N[j]*𝑤
                k[3*I-1,2*J]   -= Nₛ[i]*n₂*F₂₂*N[j]*𝑤
                k[3*I,2*J-1]   -= (Nₛ[i]*n₂*F₁₁*N[j]+Nₛ[i]*n₁*F₁₂*N[j])*𝑤
                k[3*I,2*J]     -= (Nₛ[i]*n₁*F₂₂*N[j]+Nₛ[i]*n₂*F₂₁*N[j])*𝑤


            end
        end
    end
end

function ∫∫δSFnudxdy_HR_SaintVenantKirchhoff(aₛ::T,aᵤ::S,f::AbstractVector{Float64})  where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
   
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
     
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        𝑤 = ξₛ.𝑤
      
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end
   
        for (i,xᵢ) in enumerate(𝓒ₛ)
              I = xᵢ.𝐼
               
               f[3*I-2] -= (Nₛ[i]*n₁*F₁₁*u₁ + Nₛ[i]*n₁*F₂₁*u₂)*𝑤
               f[3*I-1] -= (Nₛ[i]*n₂*F₁₂*u₁ + Nₛ[i]*n₂*F₂₂*u₂)*𝑤
               f[3*I]   -= ((Nₛ[i]*n₂*F₁₁+Nₛ[i]*n₁*F₁₂)*u₁ + (Nₛ[i]*n₁*F₂₂+Nₛ[i]*n₂*F₂₁)*u₂)*𝑤

        end

        
    end
end


function Δ∫∫δ∇SFudxdy_HR_SaintVenantKirchhoff(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
   
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        B₁₁ = ξᵤ[:∂²𝝭∂x²]
        B₁₂ = ξᵤ[:∂²𝝭∂x∂y]
        B₂₂ = ξᵤ[:∂²𝝭∂y²]
        
        Bₛ₁ = ξₛ[:∂𝝭∂x]
        Bₛ₂ = ξₛ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
       
       
       
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0

        F₁₁_₁ = 0.0
        F₁₁_₂ = 0.0
        F₁₂_₁ = 0.0
        F₁₂_₂ = 0.0
        F₂₁_₁ = 0.0
        F₂₁_₂ = 0.0
        F₂₂_₁ = 0.0
        F₂₂_₂ = 0.0

        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂

            F₁₁_₁  += B₁₁[i]*xᵢ.d₁
            F₁₁_₂  += B₁₂[i]*xᵢ.d₁
            F₁₂_₁  += B₁₂[i]*xᵢ.d₁
            F₁₂_₂  += B₂₂[i]*xᵢ.d₁
            F₂₁_₁  += B₁₁[i]*xᵢ.d₂
            F₂₁_₂  += B₁₂[i]*xᵢ.d₂
            F₂₂_₁  += B₁₂[i]*xᵢ.d₂
            F₂₂_₂  += B₂₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

      
        
        
      for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                # ∫∫δ∇SΔFudxdy
                k[3*I-2,2*J-1] += Bₛ₁[i]*u₁*B₁[j]*𝑤
                k[3*I-2,2*J]   += Bₛ₁[i]*u₂*B₁[j]*𝑤
                k[3*I-1,2*J-1] += Bₛ₂[i]*u₁*B₂[j]*𝑤
                k[3*I-1,2*J]   += Bₛ₂[i]*u₂*B₂[j]*𝑤
                k[3*I,2*J-1]   += (Bₛ₂[i]*u₁*B₁[j]+Bₛ₁[i]*u₁*B₂[j])*𝑤
                k[3*I,2*J]     += (Bₛ₁[i]*u₂*B₂[j]+Bₛ₂[i]*u₂*B₁[j])*𝑤
                # ∫∫δ∇SFΔudxdy
                k[3*I-2,2*J-1] += Bₛ₁[i]*F₁₁*N[j]*𝑤
                k[3*I-2,2*J]   += Bₛ₁[i]*F₂₁*N[j]*𝑤
                k[3*I-1,2*J-1] += Bₛ₂[i]*F₁₂*N[j]*𝑤
                k[3*I-1,2*J]   += Bₛ₂[i]*F₂₂*N[j]*𝑤
                k[3*I,2*J-1]   += (Bₛ₂[i]*F₁₁*N[j]+Bₛ₁[i]*F₁₂*N[j])*𝑤
                k[3*I,2*J]     += (Bₛ₁[i]*F₂₂*N[j]+Bₛ₂[i]*F₂₁*N[j])*𝑤

                # ∫∫δS∇ΔFudxdy
                k[3*I-2,2*J-1] += Nₛ[i]*u₁*B₁₁[j]*𝑤
                k[3*I-2,2*J]   += Nₛ[i]*u₂*B₁₁[j]*𝑤
                k[3*I-1,2*J-1] += Nₛ[i]*u₁*B₂₂[j]*𝑤
                k[3*I-1,2*J]   += Nₛ[i]*u₂*B₂₂[j]*𝑤
                k[3*I,2*J-1]   += (Nₛ[i]*u₁*B₁₂[j]+Nₛ[i]*u₁*B₁₂[j])*𝑤
                k[3*I,2*J]     += (Nₛ[i]*u₂*B₁₂[j]+Nₛ[i]*u₂*B₁₂[j])*𝑤

                #  # ∫∫δS∇FΔudxdy
                k[3*I-2,2*J-1] += Nₛ[i]*F₁₁_₁*N[j]*𝑤
                k[3*I-2,2*J]   += Nₛ[i]*F₂₁_₁*N[j]*𝑤
                k[3*I-1,2*J-1] += Nₛ[i]*F₁₂_₂*N[j]*𝑤
                k[3*I-1,2*J]   += Nₛ[i]*F₂₂_₂*N[j]*𝑤
                k[3*I,2*J-1]   += (Nₛ[i]*F₁₁_₂*N[j]+Nₛ[i]*F₁₂_₁*N[j])*𝑤
                k[3*I,2*J]     += (Nₛ[i]*F₂₂_₁*N[j]+Nₛ[i]*F₂₁_₂*N[j])*𝑤




            end
        end
    end
end


 
function ∫∫δ∇SFudxdy_HR_SaintVenantKirchhoff(aₛ::T,aᵤ::S,f::AbstractVector{Float64})  where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
   
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
       
        Bₛ₁ = ξₛ[:∂𝝭∂x]
        Bₛ₂ = ξₛ[:∂𝝭∂y]
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        B₁₁ = ξᵤ[:∂²𝝭∂x²]
        B₁₂ = ξᵤ[:∂²𝝭∂x∂y]
        B₂₂ = ξᵤ[:∂²𝝭∂y²]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
       
        𝑤 = ξₛ.𝑤
       
        
       
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0

        F₁₁_₁ = 0.0
        F₁₁_₂ = 0.0
        F₁₂_₁ = 0.0
        F₁₂_₂ = 0.0
        F₂₁_₁ = 0.0
        F₂₁_₂ = 0.0
        F₂₂_₁ = 0.0
        F₂₂_₂ = 0.0

        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂

            
            F₁₁_₁  += B₁₁[i]*xᵢ.d₁
            F₁₁_₂  += B₁₂[i]*xᵢ.d₁
            F₁₂_₁  += B₁₂[i]*xᵢ.d₁
            F₁₂_₂  += B₂₂[i]*xᵢ.d₁
            F₂₁_₁  += B₁₁[i]*xᵢ.d₂
            F₂₁_₂  += B₁₂[i]*xᵢ.d₂
            F₂₂_₁  += B₁₂[i]*xᵢ.d₂
            F₂₂_₂  += B₂₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

     
       

        for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
               f[3*I-2] += (Bₛ₁[i]*F₁₁*u₁ + Bₛ₁[i]*F₂₁*u₂)*𝑤
               f[3*I-1] += (Bₛ₂[i]*F₁₂*u₁ + Bₛ₂[i]*F₂₂*u₂)*𝑤
               f[3*I]   += ((Bₛ₂[i]*F₁₁+Bₛ₁[i]*F₁₂)*u₁ + (Bₛ₁[i]*F₂₂+Bₛ₂[i]*F₂₁)*u₂)*𝑤

               f[3*I-2] += (Nₛ[i]*F₁₁_₁*u₁ + Nₛ[i]*F₂₁_₁*u₂)*𝑤
               f[3*I-1] += (Nₛ[i]*F₁₂_₂*u₁ + Nₛ[i]*F₂₂_₂*u₂)*𝑤
               f[3*I]   += ((Nₛ[i]*F₁₁_₂+Nₛ[i]*F₁₂_₁)*u₁ + (Nₛ[i]*F₂₂_₁+Nₛ[i]*F₂₁_₂)*u₂)*𝑤


        end

      
    end
end

function ∫∫SδFnΔudxdy_HR_SaintVenantKirchhoff(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
   
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
       
       
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end
        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
        end
          
      for (i,xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
              # ∫∫SΔFnδudxdy
                k[2*I-1,2*J-1] -= (N[i]*B₁[j]*(S₁₁*n₁+S₁₂*n₂)+N[i]*B₂[j]*(S₁₂*n₁+S₂₂*n₂))*𝑤
                k[2*I,2*J]     -= (N[i]*B₁[j]*(S₁₁*n₁+S₁₂*n₂)+N[i]*B₂[j]*(S₁₂*n₁+S₂₂*n₂))*𝑤

                # ∫∫SδFnΔudxdy
                k[2*I-1,2*J-1] -= (N[j]*B₁[i]*(S₁₁*n₁+S₁₂*n₂)+N[j]*B₂[i]*(S₁₂*n₁+S₂₂*n₂))*𝑤
                k[2*I,2*J]     -= (N[j]*B₁[i]*(S₁₁*n₁+S₁₂*n₂)+N[j]*B₂[i]*(S₁₂*n₁+S₂₂*n₂))*𝑤


            end
        end
    end
end
function ∫∫SδFnδudxdy_HR_SaintVenantKirchhoff(aₛ::T,aᵤ::S,f::AbstractVector{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
   
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
       
       
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end
        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
        end
           
        
        for (i,xᵢ) in enumerate(𝓒ᵤ)
                I = xᵢ.𝐼
            
              # ∫∫SδFnudxdy
               f[2*I-1] -= (u₁*B₁[i]*S₁₁*n₁+u₁*B₁[i]*S₁₂*n₂+u₁*B₂[i]*S₁₂*n₁+u₁*B₂[i]*S₂₂*n₂ )*𝑤
               f[2*I]   -= (u₂*B₁[i]*S₁₁*n₁+u₂*B₁[i]*S₁₂*n₂+u₂*B₂[i]*S₁₂*n₁+u₂*B₂[i]*S₂₂*n₂ )*𝑤
                # ∫∫SFnδudxdy
               f[2*I-1] -= (N[i]*F₁₁*S₁₁*n₁+N[i]*F₁₁*S₁₂*n₂+N[i]*F₁₂*S₁₂*n₁+N[i]*F₁₂*S₂₂*n₂ )*𝑤
               f[2*I]   -= (N[i]*F₂₁*S₁₁*n₁+N[i]*F₂₁*S₁₂*n₂+N[i]*F₂₂*S₁₂*n₁+N[i]*F₂₂*S₂₂*n₂ )*𝑤    
       
        end
    end
end


function ∫∫∇SδFΔudxdy_HR_SaintVenantKirchhoff(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
   
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        B₁₁ = ξᵤ[:∂²𝝭∂x²]
        B₁₂ = ξᵤ[:∂²𝝭∂x∂y]
        B₂₂ = ξᵤ[:∂²𝝭∂y²]
        Bₛ₁ = ξₛ[:∂𝝭∂x]
        Bₛ₂ = ξₛ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
       

        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        S₁₁_₁ = 0.0
        S₂₂_₂ = 0.0
        S₁₂_₁ = 0.0
        S₁₂_₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂

           S₁₁_₁ += Bₛ₁[i]*xᵢ.dₛ₁₁
           S₂₂_₂ += Bₛ₂[i]*xᵢ.dₛ₂₂
           S₁₂_₁ += Bₛ₁[i]*xᵢ.dₛ₁₂
           S₁₂_₂ += Bₛ₂[i]*xᵢ.dₛ₁₂

        end
          
      for (i,xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
              # ∫∫∇SΔFδudxdy
                k[2*I-1,2*J-1] += (N[i]*B₁[j]*(S₁₁_₁+S₁₂_₂)+N[i]*B₂[j]*(S₁₂_₁+S₂₂_₂))*𝑤
                k[2*I,2*J]     += (N[i]*B₁[j]*(S₁₁_₁+S₁₂_₂)+N[i]*B₂[j]*(S₁₂_₁+S₂₂_₂))*𝑤

                # ∫∫∇SδFΔudxdy
                k[2*I-1,2*J-1] += (N[j]*B₁[i]*(S₁₁_₁+S₁₂_₂)+N[j]*B₂[i]*(S₁₂_₁+S₂₂_₂))*𝑤
                k[2*I,2*J]     += (N[j]*B₁[i]*(S₁₁_₁+S₁₂_₂)+N[j]*B₂[i]*(S₁₂_₁+S₂₂_₂))*𝑤



                # ∫∫S∇ΔFnδudxdy补充

                k[2*I-1,2*J-1] += (N[i]*(B₁₁[j]*S₁₁+B₁₂[j]*S₁₂)+N[i]*(B₁₂[j]*S₁₂+B₂₂[j]*S₂₂))*𝑤
                k[2*I,2*J]     += (N[i]*(B₁₁[j]*S₁₁+B₁₂[j]*S₁₂)+N[i]*(B₁₂[j]*S₁₂+B₂₂[j]*S₂₂))*𝑤

                k[2*I-1,2*J-1] += (N[j]*(B₁₁[i]*S₁₁+B₁₂[i]*S₁₂)+N[j]*(B₁₂[i]*S₁₂+B₂₂[i]*S₂₂))*𝑤
                k[2*I,2*J]     += (N[j]*(B₁₁[i]*S₁₁+B₁₂[i]*S₁₂)+N[j]*(B₁₂[i]*S₁₂+B₂₂[i]*S₂₂))*𝑤






            end
        end
    end
end

function ∫∫∇SδFδudxdy_HR_SaintVenantKirchhoff(aₛ::T,aᵤ::S,f::AbstractVector{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
   
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
         B₁₁ = ξᵤ[:∂²𝝭∂x²]
        B₁₂ = ξᵤ[:∂²𝝭∂x∂y]
        B₂₂ = ξᵤ[:∂²𝝭∂y²]
        Bₛ₁ = ξₛ[:∂𝝭∂x]
        Bₛ₂ = ξₛ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
     
       
       
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0

        F₁₁_₁ = 0.0
        F₁₁_₂ = 0.0
        F₁₂_₁ = 0.0
        F₁₂_₂ = 0.0
        F₂₁_₁ = 0.0
        F₂₁_₂ = 0.0
        F₂₂_₁ = 0.0
        F₂₂_₂ = 0.0

        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂

           
          
            F₁₁_₁  += B₁₁[i]*xᵢ.d₁
            F₁₁_₂  += B₁₂[i]*xᵢ.d₁
            F₁₂_₁  += B₁₂[i]*xᵢ.d₁
            F₁₂_₂  += B₂₂[i]*xᵢ.d₁
            F₂₁_₁  += B₁₁[i]*xᵢ.d₂
            F₂₁_₂  += B₁₂[i]*xᵢ.d₂
            F₂₂_₁  += B₁₂[i]*xᵢ.d₂
            F₂₂_₂  += B₂₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        S₁₁_₁ = 0.0
        S₂₂_₂ = 0.0
        S₁₂_₁ = 0.0
        S₁₂_₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂

           S₁₁_₁ += Bₛ₁[i]*xᵢ.dₛ₁₁
           S₂₂_₂ += Bₛ₂[i]*xᵢ.dₛ₂₂
           S₁₂_₁ += Bₛ₁[i]*xᵢ.dₛ₁₂
           S₁₂_₂ += Bₛ₂[i]*xᵢ.dₛ₁₂

        end
            
        
      for (i,xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
            
              # ∫∫SδFnudxdy
               f[2*I-1] += (u₁*B₁[i]*S₁₁_₁+u₁*B₁[i]*S₁₂_₂+u₁*B₂[i]*S₁₂_₁+u₁*B₂[i]*S₂₂_₂ )*𝑤
               f[2*I]   += (u₂*B₁[i]*S₁₁_₁+u₂*B₁[i]*S₁₂_₂+u₂*B₂[i]*S₁₂_₁+u₂*B₂[i]*S₂₂_₂ )*𝑤
                # ∫∫SFnδudxdy
               f[2*I-1] += (N[i]*F₁₁*S₁₁_₁+N[i]*F₁₁*S₁₂_₂+N[i]*F₁₂*S₁₂_₁+N[i]*F₁₂*S₂₂_₂ )*𝑤
               f[2*I]   += (N[i]*F₂₁*S₁₁_₁+N[i]*F₂₁*S₁₂_₂+N[i]*F₂₂*S₁₂_₁+N[i]*F₂₂*S₂₂_₂ )*𝑤
               
               f[2*I-1] += (u₁*B₁₁[i]*S₁₁+u₁*B₁₂[i]*S₁₂+u₁*B₁₂[i]*S₁₂+u₁*B₂₂[i]*S₂₂ )*𝑤
               f[2*I]   += (u₂*B₁₁[i]*S₁₁+u₂*B₁₂[i]*S₁₂+u₂*B₁₂[i]*S₁₂+u₂*B₂₂[i]*S₂₂ )*𝑤

               f[2*I-1] += (N[i]*F₁₁_₁*S₁₁+N[i]*F₁₁_₂*S₁₂+N[i]*F₁₂_₁*S₁₂+N[i]*F₁₂_₂*S₂₂ )*𝑤
               f[2*I]   += (N[i]*F₂₁_₁*S₁₁+N[i]*F₂₁_₂*S₁₂+N[i]*F₂₂_₁*S₁₂+N[i]*F₂₂_₂*S₂₂ )*𝑤

        end
    end
end

function ∫∫δSΔFngdxdy_HR_SaintVenantKirchhoff(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
   
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        n₁₁ = ξᵤ.n₁₁
        n₁₂ = ξᵤ.n₁₂
        n₂₂ = ξᵤ.n₂₂
        g₁ = ξᵤ.g₁
        g₂ = ξᵤ.g₂

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end
         
        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
        end
         
        Δu₁ = g₁-u₁
        Δu₂ = g₂-u₂
        
      for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                # ∫∫δSΔFnudxdy
                k[3*I-2,2*J-1] += Nₛ[i]*(n₁*n₁₁*u₁+n₁*n₁₂*u₂)*B₁[j]*𝑤
                k[3*I-2,2*J]   += Nₛ[i]*(n₁*n₁₂*u₁+n₁*n₂₂*u₂)*B₁[j]*𝑤
                k[3*I-1,2*J-1] += Nₛ[i]*(n₂*n₁₁*u₁+n₂*n₁₂*u₂)*B₂[j]*𝑤
                k[3*I-1,2*J]   += Nₛ[i]*(n₂*n₁₂*u₁+n₂*n₂₂*u₂)*B₂[j]*𝑤
                k[3*I,2*J-1]   += (Nₛ[i]*(n₂*n₁₁*u₁+n₂*n₁₂*u₂)*B₁[j]+Nₛ[i]*(n₁*n₁₁*u₁+n₁*n₁₂*u₂)*B₂[j])*𝑤
                k[3*I,2*J]     += (Nₛ[i]*(n₂*n₁₂*u₁+n₂*n₂₂*u₂)*B₁[j]+Nₛ[i]*(n₁*n₁₂*u₁+n₁*n₂₂*u₂)*B₂[j])*𝑤
                

                # ∫∫δSFnΔudxdy
                k[3*I-2,2*J-1] += Nₛ[i]*(n₁*n₁₁*F₁₁+n₁*n₁₂*F₂₁)*N[j]*𝑤
                k[3*I-2,2*J]   += Nₛ[i]*(n₁*n₁₂*F₁₁+n₁*n₂₂*F₂₁)*N[j]*𝑤
                k[3*I-1,2*J-1] += Nₛ[i]*(n₂*n₁₁*F₁₂+n₂*n₁₂*F₂₂)*N[j]*𝑤
                k[3*I-1,2*J]   += Nₛ[i]*(n₂*n₁₂*F₁₂+n₂*n₂₂*F₂₂)*N[j]*𝑤
                k[3*I,2*J-1]   += (Nₛ[i]*(n₂*n₁₁*F₁₁+n₂*n₁₂*F₂₁)*N[j]+Nₛ[i]*(n₁*n₁₁*F₁₂+n₁*n₁₂*F₂₂)*N[j])*𝑤
                k[3*I,2*J]     += (Nₛ[i]*(n₂*n₁₂*F₁₁+n₂*n₂₂*F₂₁)*N[j]+Nₛ[i]*(n₁*n₁₂*F₁₂+n₁*n₂₂*F₂₂)*N[j])*𝑤
                # ∫∫δSΔFngdxdy
                k[3*I-2,2*J-1] -= Nₛ[i]*(n₁*n₁₁*g₁+n₁*n₁₂*g₂)*B₁[j]*𝑤
                k[3*I-2,2*J]   -= Nₛ[i]*(n₁*n₁₂*g₁+n₁*n₂₂*g₂)*B₁[j]*𝑤
                k[3*I-1,2*J-1] -= Nₛ[i]*(n₂*n₁₁*g₁+n₂*n₁₂*g₂)*B₂[j]*𝑤
                k[3*I-1,2*J]   -= Nₛ[i]*(n₂*n₁₂*g₁+n₂*n₂₂*g₂)*B₂[j]*𝑤
                k[3*I,2*J-1]   -= (Nₛ[i]*(n₂*n₁₁*g₁+n₂*n₁₂*g₂)*B₁[j]+Nₛ[i]*(n₁*n₁₁*g₁+n₁*n₁₂*g₂)*B₂[j])*𝑤
                k[3*I,2*J]     -= (Nₛ[i]*(n₂*n₁₂*g₁+n₂*n₂₂*g₂)*B₁[j]+Nₛ[i]*(n₁*n₁₂*g₁+n₁*n₂₂*g₂)*B₂[j])*𝑤
            


            end
            #  ∫∫δSFngdxdy
            f[3*I-2] += Nₛ[i]*((n₁*n₁₁*F₁₁+n₁*n₁₂*F₂₁)*Δu₁ + (n₁*n₁₂*F₁₁+n₁*n₂₂*F₂₁)*Δu₂)*𝑤
            f[3*I-1] += Nₛ[i]*((n₂*n₁₁*F₁₂+n₂*n₁₂*F₂₂)*Δu₁ + (n₂*n₁₂*F₁₂+n₂*n₂₂*F₂₂)*Δu₂)*𝑤
            f[3*I]   += Nₛ[i]*(((n₂*n₁₁*F₁₁+n₂*n₁₂*F₂₁)+(n₁*n₁₁*F₁₂+n₁*n₁₂*F₂₂))*Δu₁ + ((n₂*n₁₂*F₁₁+n₂*n₂₂*F₂₁)+(n₁*n₁₂*F₁₂+n₁*n₂₂*F₂₂))*Δu₂)*𝑤
               
        end
                
    end
end

function ∫∫ΔSδFngdxdy_HR_SaintVenantKirchhoff(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
   
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        n₁₁ = ξᵤ.n₁₁
        n₁₂ = ξᵤ.n₁₂
        n₂₂ = ξᵤ.n₂₂
        g₁ = ξᵤ.g₁
        g₂ = ξᵤ.g₂

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end
         
        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
        end
         
        Δu₁ = g₁-u₁
        Δu₂ = g₂-u₂
        
      for (i,xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                # ∫∫SΔFnδudxdy
                k[2*I-1,2*J-1] += (N[i]*B₁[j]*(S₁₁*n₁+S₁₂*n₂)+N[i]*B₂[j]*(S₁₂*n₁+S₂₂*n₂))*n₁₁*𝑤
                k[2*I-1,2*J]   += (N[i]*B₁[j]*(S₁₁*n₁+S₁₂*n₂)+N[i]*B₂[j]*(S₁₂*n₁+S₂₂*n₂))*n₁₂*𝑤
                k[2*I,2*J-1]   += (N[i]*B₁[j]*(S₁₁*n₁+S₁₂*n₂)+N[i]*B₂[j]*(S₁₂*n₁+S₂₂*n₂))*n₁₂*𝑤
                k[2*I,2*J]     += (N[i]*B₁[j]*(S₁₁*n₁+S₁₂*n₂)+N[i]*B₂[j]*(S₁₂*n₁+S₂₂*n₂))*n₂₂*𝑤

                # ∫∫SδFnΔudxdy
              

                k[2*I-1,2*J-1] += (N[j]*B₁[i]*(S₁₁*n₁+S₁₂*n₂)+N[j]*B₂[i]*(S₁₂*n₁+S₂₂*n₂))*n₁₁*𝑤
                k[2*I-1,2*J]   += (N[j]*B₁[i]*(S₁₁*n₁+S₁₂*n₂)+N[j]*B₂[i]*(S₁₂*n₁+S₂₂*n₂))*n₁₂*𝑤
                k[2*I,2*J-1]   += (N[j]*B₁[i]*(S₁₁*n₁+S₁₂*n₂)+N[j]*B₂[i]*(S₁₂*n₁+S₂₂*n₂))*n₁₂*𝑤
                k[2*I,2*J]     += (N[j]*B₁[i]*(S₁₁*n₁+S₁₂*n₂)+N[j]*B₂[i]*(S₁₂*n₁+S₂₂*n₂))*n₂₂*𝑤
            end
           
        end
                
                for (i,xᵢ) in enumerate(𝓒ᵤ)
                     I = xᵢ.𝐼
                     #  ∫∫SδFngdxdy
                        f[2*I-1] += ((Δu₁*B₁[i]*(S₁₁*n₁+S₁₂*n₂)+Δu₁*B₂[i]*(S₁₂*n₁+S₂₂*n₂))*n₁₁ + (Δu₂*B₁[i]*(S₁₁*n₁+S₁₂*n₂)+Δu₂*B₂[i]*(S₁₂*n₁+S₂₂*n₂))*n₁₂)*𝑤
                        f[2*I]   += ((Δu₁*B₁[i]*(S₁₁*n₁+S₁₂*n₂)+Δu₁*B₂[i]*(S₁₂*n₁+S₂₂*n₂))*n₁₂ + (Δu₂*B₁[i]*(S₁₁*n₁+S₁₂*n₂)+Δu₂*B₂[i]*(S₁₂*n₁+S₂₂*n₂))*n₂₂)*𝑤
                                      
                        # # ∫∫SFnδudxdy
                        # f[2*I-1] -= ((N[i]*F₁₁*(S₁₁*n₁+S₁₂*n₂)+N[i]*F₁₂*(S₁₂*n₁+S₂₂*n₂))*n₁₁ + (N[i]*F₂₁*(S₁₁*n₁+S₁₂*n₂)+N[i]*F₂₂*(S₁₂*n₁+S₂₂*n₂))*n₁₂)*𝑤
                                    
                        # f[2*I]   -= ((N[i]*F₁₁*(S₁₁*n₁+S₁₂*n₂)+N[i]*F₁₂*(S₁₂*n₁+S₂₂*n₂))*n₁₂ + (N[i]*F₂₁*(S₁₁*n₁+S₁₂*n₂)+N[i]*F₂₂*(S₁₂*n₁+S₂₂*n₂))*n₂₂)*𝑤
                         
                 end
    end
end
function ∫∫SFnδudxdy_HR_SaintVenantKirchhoff(aₛ::T,aᵤ::S,f::AbstractVector{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
   
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        n₁₁ = ξᵤ.n₁₁
        n₁₂ = ξᵤ.n₁₂
        n₂₂ = ξᵤ.n₂₂
        g₁ = ξᵤ.g₁
        g₂ = ξᵤ.g₂

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end
         
        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
        end
         
                for (i,xᵢ) in enumerate(𝓒ᵤ)
                     I = xᵢ.𝐼
                    # ∫∫SFnδudxdy
                        f[2*I-1] += ((N[i]*F₁₁*(S₁₁*n₁+S₁₂*n₂)+N[i]*F₁₂*(S₁₂*n₁+S₂₂*n₂))*n₁₁ + (N[i]*F₂₁*(S₁₁*n₁+S₁₂*n₂)+N[i]*F₂₂*(S₁₂*n₁+S₂₂*n₂))*n₁₂)*𝑤
                                    
                        f[2*I]   += ((N[i]*F₁₁*(S₁₁*n₁+S₁₂*n₂)+N[i]*F₁₂*(S₁₂*n₁+S₂₂*n₂))*n₁₂ + (N[i]*F₂₁*(S₁₁*n₁+S₁₂*n₂)+N[i]*F₂₂*(S₁₂*n₁+S₂₂*n₂))*n₂₂)*𝑤
                              
                 end
    end
end

# function Δ∫∫EᵢⱼSᵢⱼdxdy_SaintVenantKirchhoff(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
#     𝓒 = ap.𝓒; 𝓖 = ap.𝓖
   
#     for ξ in 𝓖
#         E = ξ.E
#         ν = ξ.ν
#         B₁ = ξ[:∂𝝭∂x]
#         B₂ = ξ[:∂𝝭∂y]
#         𝑤 = ξ.𝑤
#         Cᵢᵢᵢᵢ = E*(1-ν)/(1+ν)/(1-2*ν)
#         Cᵢᵢⱼⱼ = E*ν/(1+ν)/(1-2*ν)
#         Cᵢⱼᵢⱼ = E/(1+ν)/2
       
#         F₁₁ = 1.0
#         F₁₂ = 0.0
#         F₂₁ = 0.0
#         F₂₂ = 1.0
#         for (i,xᵢ) in  enumerate(𝓒)
#             F₁₁ += B₁[i]*xᵢ.d₁
#             F₁₂ += B₂[i]*xᵢ.d₁
#             F₂₁ += B₁[i]*xᵢ.d₂
#             F₂₂ += B₂[i]*xᵢ.d₂
#         end
        
#         E₁₁ = 0.5*(F₁₁*F₁₁+F₂₁*F₂₁-1.0)
#         E₁₂ = 0.5*(F₁₁*F₁₂+F₂₁*F₂₂)
#         E₂₂ = 0.5*(F₁₂*F₁₂+F₂₂*F₂₂-1.0)
#         S₁₁ = Cᵢᵢᵢᵢ*E₁₁+Cᵢᵢⱼⱼ*E₂₂
#         S₂₂ = Cᵢᵢⱼⱼ*E₁₁+Cᵢᵢᵢᵢ*E₂₂
#         S₁₂ = 2.0*Cᵢⱼᵢⱼ*E₁₂
        
#         for (i,xᵢ) in enumerate(𝓒)
#             I = xᵢ.𝐼
#             for (j,xⱼ) in enumerate(𝓒)
#                 J = xⱼ.𝐼
#                 k[2*I-1,2*J-1] +=((B₁[i]*F₁₁*B₁[j]*F₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂)*Cᵢᵢᵢᵢ
#                                +  (B₁[i]*F₁₁*B₂[j]*F₁₂ + B₂[i]*F₁₂*B₁[j]*F₁₁)*Cᵢᵢⱼⱼ
#                                +  (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*Cᵢⱼᵢⱼ
#                                +   B₁[i]*B₁[j]*S₁₁+B₂[i]*B₂[j]*S₂₂+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂)*𝑤
                              
#                 k[2*I-1,2*J]   += ((B₁[i]*F₁₁*B₁[j]*F₂₁ + B₂[i]*F₁₂*B₂[j]*F₂₂)*Cᵢᵢᵢᵢ
#                                +   (B₁[i]*F₁₁*B₂[j]*F₂₂ + B₂[i]*F₁₂*B₁[j]*F₂₁)*Cᵢᵢⱼⱼ
#                                +   (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*Cᵢⱼᵢⱼ)*𝑤
                               
#                 k[2*I,2*J-1]   += ((B₁[i]*F₂₁*B₁[j]*F₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂)*Cᵢᵢᵢᵢ
#                                +   (B₁[i]*F₂₁*B₂[j]*F₁₂ + B₂[i]*F₂₂*B₁[j]*F₁₁)*Cᵢᵢⱼⱼ
#                                +   (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*Cᵢⱼᵢⱼ)*𝑤

#                 k[2*I,2*J]     += ((B₁[i]*F₂₁*B₁[j]*F₂₁ + B₂[i]*F₂₂*B₂[j]*F₂₂)*Cᵢᵢᵢᵢ 
#                                +   (B₁[i]*F₂₁*B₂[j]*F₂₂ + B₂[i]*F₂₂*B₁[j]*F₂₁)*Cᵢᵢⱼⱼ
#                                +   (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*Cᵢⱼᵢⱼ
#                                +    B₁[i]*B₁[j]*S₁₁+B₂[i]*B₂[j]*S₂₂+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂)*𝑤

#             end
#         end
#     end
# end

function Δ∫∫EᵢⱼSᵢⱼdxdy_SaintVenantKirchhoff(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    
    # 预分配 2x2 矩阵
    F = zeros(2,2)
    S = zeros(2,2)
    
    for ξ in 𝓖
        E_mod = ξ.Ē
        ν = ξ.ν̄ 
        
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        
        # 将原代码中的 C常数 转换为 Lamé 常数
        # 原代码使用的是平面应变(Plane Strain)或三维形式的系数
        # λ = E*ν / ((1+ν)(1-2ν))
        # μ = E / (2(1+ν))
        λ = E_mod * ν / ((1 + ν) * (1 - 2 * ν))
        μ = E_mod / (2 * (1 + ν))
        
        # 1. 计算变形梯度 F (2x2)
        fill!(F, 0.0)
        F[1,1] = 1.0; F[2,2] = 1.0
        
        for (i, xᵢ) in enumerate(𝓒)
            d = [xᵢ.d₁, xᵢ.d₂]
            # F_ij = δ_ij + Σ (u_i * ∂N/∂X_j)
            F[1,1] += d[1] * B₁[i]; F[1,2] += d[1] * B₂[i]
            F[2,1] += d[2] * B₁[i]; F[2,2] += d[2] * B₂[i]
        end
        
        # 2. 计算 Green-Lagrange 应变 E = 0.5 * (F'F - I)
        C = F' * F
        E_strain = 0.5 * (C - I)
        
        # 3. 计算 Second Piola-Kirchhoff 应力 S
        # S = λ*tr(E)*I + 2*μ*E
        trE = tr(E_strain)
        S .= λ * trE * I + 2 * μ * E_strain
        
        # 4. 组装刚度矩阵 (循环自由度 a, b)
        for (i, xᵢ) in enumerate(𝓒)
            I_idx = xᵢ.𝐼
            ∇N_i = [B₁[i], B₂[i]] # 节点 i 的梯度向量
            
            for (j, xⱼ) in enumerate(𝓒)
                J_idx = xⱼ.𝐼
                ∇N_j = [B₁[j], B₂[j]] # 节点 j 的梯度向量
                
                # --- 几何刚度 (Geometric Stiffness) ---
                # scalar = ∇N_i · S · ∇N_j
                s_geo = dot(∇N_i, S, ∇N_j) * 𝑤
                
                # 预计算梯度点积
                grad_dot = dot(∇N_i, ∇N_j)
                
                # 循环自由度: a, b 取 1 (x方向) 或 2 (y方向)
                for a in 1:2 
                    F_a = view(F, a, :) # F 的第 a 行
                    Fa_dot_Ni = dot(F_a, ∇N_i)
                    Fa_dot_Nj = dot(F_a, ∇N_j)
                    
                    for b in 1:2
                        F_b = view(F, b, :) # F 的第 b 行
                        Fb_dot_Nj = dot(F_b, ∇N_j)
                        Fb_dot_Ni = dot(F_b, ∇N_i)
                        Fa_dot_Fb = dot(F_a, F_b)
                        
                        # --- 材料刚度 (Material Stiffness) ---
                        # 使用各向同性张量公式
                        k_val = (λ * Fa_dot_Ni * Fb_dot_Nj + 
                                 μ * Fa_dot_Fb * grad_dot + 
                                 μ * Fb_dot_Ni * Fa_dot_Nj) * 𝑤
                        
                        # 加上几何刚度 (仅对角线项)
                        if a == b
                            k_val += s_geo
                        end
                        
                        # 填入全局刚度矩阵
                        # 索引映射: 1->2I-1, 2->2I
                        row = 2 * I_idx - (2 - a)
                        col = 2 * J_idx - (2 - b)
                        k[row, col] += k_val
                    end
                end
            end
        end
    end
end

function Δ∫∫∫EᵢⱼSᵢⱼdxdydz_SaintVenantKirchhoff(ap::T, k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
   
    for ξ in 𝓖
        E = ξ.E
        ν = ξ.ν
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z]  # 3D: Z方向形函数梯度
        𝑤 = ξ.𝑤
        
        # --- 1. Lamé常数 (用于SVK应力计算和K_M) ---
        λ = E*ν/((1+ν)*(1-2*ν))  # Cᵢᵢⱼⱼ
        μ = E/(2*(1+ν))          # Cᵢⱼᵢⱼ
        
        
       Cᵢᵢᵢᵢ = λ + 2*μ
       Cᵢᵢⱼⱼ = λ
       Cᵢⱼᵢⱼ = μ
        # --- 2. 计算 3x3 变形梯度 F ---
        F₁₁ = 1.0; F₁₂ = 0.0; F₁₃ = 0.0
        F₂₁ = 0.0; F₂₂ = 1.0; F₂₃ = 0.0
        F₃₁ = 0.0; F₃₂ = 0.0; F₃₃ = 1.0
        
        for (i,xᵢ) in  enumerate(𝓒)
            # F_ij = δ_ij + (∂u_i / ∂X_j) = δ_ij + (∂N_k / ∂X_j) * d_i^k
            F₁₁ += B₁[i]*xᵢ.d₁; F₁₂ += B₂[i]*xᵢ.d₁; F₁₃ += B₃[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂; F₂₂ += B₂[i]*xᵢ.d₂; F₂₃ += B₃[i]*xᵢ.d₂
            F₃₁ += B₁[i]*xᵢ.d₃; F₃₂ += B₂[i]*xᵢ.d₃; F₃₃ += B₃[i]*xᵢ.d₃
        end
        
        # --- 3. 计算 Green-Lagrange 应变 E = 0.5*(FᵀF - I) ---
        E₁₁ = 0.5*(F₁₁*F₁₁+F₂₁*F₂₁+F₃₁*F₃₁-1.0)
        E₂₂ = 0.5*(F₁₂*F₁₂+F₂₂*F₂₂+F₃₂*F₃₂-1.0)
        E₃₃ = 0.5*(F₁₃*F₁₃+F₂₃*F₂₃+F₃₃*F₃₃-1.0) # 3D 新增
        E₁₂ = 0.5*(F₁₁*F₁₂+F₂₁*F₂₂+F₃₁*F₃₂)
        E₁₃ = 0.5*(F₁₁*F₁₃+F₂₁*F₂₃+F₃₁*F₃₃)    # 3D 新增
        E₂₃ = 0.5*(F₁₂*F₁₃+F₂₂*F₂₃+F₃₂*F₃₃)    # 3D 新增

        # --- 4. 计算第二 Piola-Kirchhoff 应力 S = C:E (SVK模型) ---
       S₁₁ = λ*(E₁₁ + E₂₂ + E₃₃) + 2*μ*E₁₁
       S₂₂ = λ*(E₁₁ + E₂₂ + E₃₃) + 2*μ*E₂₂
       S₃₃ = λ*(E₁₁ + E₂₂ + E₃₃) + 2*μ*E₃₃
       S₁₂ = 2*μ*E₁₂
       S₁₃ = 2*μ*E₁₃
       S₂₃ = 2*μ*E₂₃ 
        
       
        
       for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] +=((B₁[i]*F₁₁*B₁[j]*F₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂ + B₃[i]*F₁₃*B₃[j]*F₁₃)*Cᵢᵢᵢᵢ
                               +  (B₂[i]*F₁₂*B₁[j]*F₁₁ + B₃[i]*F₁₃*B₁[j]*F₁₁ + B₁[i]*F₁₁*B₂[j]*F₁₂ + B₃[i]*F₁₃*B₂[j]*F₁₂ + B₁[i]*F₁₁*B₃[j]*F₁₃ + B₂[i]*F₁₂*B₃[j]*F₁₃)*Cᵢᵢⱼⱼ
                               +  ((B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)+(B₃[i]*F₁₁+B₁[i]*F₁₃)*(B₃[j]*F₁₁+B₁[j]*F₁₃) + (B₃[i]*F₁₂+B₂[i]*F₁₃)*(B₃[j]*F₁₂+B₂[j]*F₁₃) )*Cᵢⱼᵢⱼ
                               +   B₁[i]*B₁[j]*S₁₁+B₂[i]*B₂[j]*S₂₂+B₃[i]*B₃[j]*S₃₃ + (B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂ + (B₁[i]*B₃[j]+B₃[i]*B₁[j])*S₁₃ + (B₃[i]*B₂[j]+B₂[i]*B₃[j])*S₂₃)*𝑤
                 #  1,2             
                k[3*I-2,3*J-1] += ((B₁[i]*F₁₁*B₁[j]*F₂₁ + B₂[i]*F₁₂*B₂[j]*F₂₂ + B₃[i]*F₁₃*B₃[j]*F₂₃)*Cᵢᵢᵢᵢ
                               +   (B₂[i]*F₁₂*B₁[j]*F₂₁ + B₃[i]*F₁₃*B₁[j]*F₂₁ + B₁[i]*F₁₁*B₂[j]*F₂₂ + B₃[i]*F₁₃*B₂[j]*F₂₂ + B₁[i]*F₁₁*B₃[j]*F₂₃ + B₂[i]*F₁₂*B₃[j]*F₂₃)*Cᵢᵢⱼⱼ
                               +   ((B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁) + (B₃[i]*F₁₁+B₁[i]*F₁₃)*(B₃[j]*F₂₁+B₁[j]*F₂₃) + (B₃[i]*F₁₂+B₂[i]*F₁₃)*(B₃[j]*F₂₂+B₂[j]*F₂₃))*Cᵢⱼᵢⱼ)*𝑤
                 #  1,3                 
                k[3*I-2,3*J]   += ((B₁[i]*F₁₁*B₁[j]*F₃₁ + B₂[i]*F₁₂*B₂[j]*F₃₂ + B₃[i]*F₁₃*B₃[j]*F₃₃)*Cᵢᵢᵢᵢ
                               +   (B₂[i]*F₁₂*B₁[j]*F₃₁ + B₃[i]*F₁₃*B₁[j]*F₃₁ + B₁[i]*F₁₁*B₂[j]*F₃₂ + B₃[i]*F₁₃*B₂[j]*F₃₂ + B₁[i]*F₁₁*B₃[j]*F₃₃ + B₂[i]*F₁₂*B₃[j]*F₃₃)*Cᵢᵢⱼⱼ
                               +   ((B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₃₂+B₂[j]*F₃₁) + (B₃[i]*F₁₁+B₁[i]*F₁₃)*(B₃[j]*F₃₁+B₁[j]*F₃₃) + (B₃[i]*F₁₂+B₂[i]*F₁₃)*(B₃[j]*F₃₂+B₂[j]*F₃₃))*Cᵢⱼᵢⱼ)*𝑤
                #  2,1
                k[3*I-1,3*J-2] +=((B₁[i]*F₂₁*B₁[j]*F₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂ + B₃[i]*F₂₃*B₃[j]*F₁₃)*Cᵢᵢᵢᵢ
                               +  (B₂[i]*F₂₂*B₁[j]*F₁₁ + B₃[i]*F₂₃*B₁[j]*F₁₁ + B₁[i]*F₂₁*B₂[j]*F₁₂ + B₃[i]*F₂₃*B₂[j]*F₁₂ + B₁[i]*F₂₁*B₃[j]*F₁₃ + B₂[i]*F₂₂*B₃[j]*F₁₃)*Cᵢᵢⱼⱼ
                               +  ((B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)+(B₃[j]*F₁₁+B₁[j]*F₁₃)*(B₃[i]*F₂₁+B₁[i]*F₂₃) + (B₃[j]*F₁₂+B₂[j]*F₁₃)*(B₃[i]*F₂₂+B₂[i]*F₂₃) )*Cᵢⱼᵢⱼ)*𝑤
                                 
                
            #    2,2
                k[3*I-1,3*J-1] += ((B₁[i]*F₂₁*B₁[j]*F₂₁ + B₂[i]*F₂₂*B₂[j]*F₂₂ + B₃[i]*F₂₃*B₃[j]*F₂₃)*Cᵢᵢᵢᵢ
                               +  (B₂[i]*F₂₂*B₁[j]*F₂₁ + B₃[i]*F₂₃*B₁[j]*F₂₁ + B₁[i]*F₂₁*B₂[j]*F₂₂ + B₃[i]*F₂₃*B₂[j]*F₂₂ + B₁[i]*F₂₁*B₃[j]*F₂₃ + B₂[i]*F₂₂*B₃[j]*F₂₃)*Cᵢᵢⱼⱼ
                               +  ((B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁) + (B₃[i]*F₂₁+B₁[i]*F₂₃)*(B₃[j]*F₂₁+B₁[j]*F₂₃) + (B₃[i]*F₂₂+B₂[i]*F₂₃)*(B₃[j]*F₂₂+B₂[j]*F₂₃))*Cᵢⱼᵢⱼ
                               +  B₁[i]*B₁[j]*S₁₁+B₂[i]*B₂[j]*S₂₂+B₃[i]*B₃[j]*S₃₃ + (B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂ + (B₁[i]*B₃[j]+B₃[i]*B₁[j])*S₁₃ + (B₃[i]*B₂[j]+B₂[i]*B₃[j])*S₂₃)*𝑤
            #   2,3
                k[3*I-1,3*J]   += ((B₁[i]*F₂₁*B₁[j]*F₃₁ + B₂[i]*F₂₂*B₂[j]*F₃₂ + B₃[i]*F₂₃*B₃[j]*F₃₃)*Cᵢᵢᵢᵢ
                               +   (B₂[i]*F₂₂*B₁[j]*F₃₁ + B₃[i]*F₂₃*B₁[j]*F₃₁ + B₁[i]*F₂₁*B₂[j]*F₃₂ + B₃[i]*F₂₃*B₂[j]*F₃₂ + B₁[i]*F₂₁*B₃[j]*F₃₃ + B₂[i]*F₂₂*B₃[j]*F₃₃)*Cᵢᵢⱼⱼ
                               +   ((B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₃₂+B₂[j]*F₃₁) + (B₃[i]*F₂₁+B₁[i]*F₂₃)*(B₃[j]*F₃₁+B₁[j]*F₃₃) + (B₃[i]*F₂₂+B₂[i]*F₂₃)*(B₃[j]*F₃₂+B₂[j]*F₃₃))*Cᵢⱼᵢⱼ)*𝑤
                             

            # 3,1
                k[3*I,3*J-2]   += ((B₁[i]*F₃₁*B₁[j]*F₁₁ + B₂[i]*F₃₂*B₂[j]*F₁₂ + B₃[i]*F₃₃*B₃[j]*F₁₃)*Cᵢᵢᵢᵢ
                               +   (B₂[i]*F₃₂*B₁[j]*F₁₁ + B₃[i]*F₃₃*B₁[j]*F₁₁ + B₁[i]*F₃₁*B₂[j]*F₁₂ + B₃[i]*F₃₃*B₂[j]*F₁₂ + B₁[i]*F₃₁*B₃[j]*F₁₃ + B₂[i]*F₃₂*B₃[j]*F₁₃)*Cᵢᵢⱼⱼ
                               +   ((B₁[i]*F₃₂+B₂[i]*F₃₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁) + (B₃[i]*F₃₁+B₁[i]*F₃₃)*(B₃[j]*F₁₁+B₁[j]*F₁₃) + (B₃[i]*F₃₂+B₂[i]*F₃₃)*(B₃[j]*F₁₂+B₂[j]*F₁₃))*Cᵢⱼᵢⱼ)*𝑤
                             

                k[3*I,3*J-1]   += ((B₁[i]*F₃₁*B₁[j]*F₂₁ + B₂[i]*F₃₂*B₂[j]*F₂₂ + B₃[i]*F₃₃*B₃[j]*F₂₃)*Cᵢᵢᵢᵢ
                               +   (B₂[i]*F₃₂*B₁[j]*F₂₁ + B₃[i]*F₃₃*B₁[j]*F₂₁ + B₁[i]*F₃₁*B₂[j]*F₂₂ + B₃[i]*F₃₃*B₂[j]*F₂₂ + B₁[i]*F₃₁*B₃[j]*F₂₃ + B₂[i]*F₃₂*B₃[j]*F₂₃)*Cᵢᵢⱼⱼ
                               +   ((B₁[i]*F₃₂+B₂[i]*F₃₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁) + (B₃[i]*F₃₁+B₁[i]*F₃₃)*(B₃[j]*F₂₁+B₁[j]*F₂₃) + (B₃[i]*F₃₂+B₂[i]*F₃₃)*(B₃[j]*F₂₂+B₂[j]*F₂₃))*Cᵢⱼᵢⱼ)*𝑤
                            

                k[3*I,3*J]   += ((B₁[i]*F₃₁*B₁[j]*F₃₁ + B₂[i]*F₃₂*B₂[j]*F₃₂ + B₃[i]*F₃₃*B₃[j]*F₃₃)*Cᵢᵢᵢᵢ
                             +   (B₂[i]*F₃₂*B₁[j]*F₃₁ + B₃[i]*F₃₃*B₁[j]*F₃₁ + B₁[i]*F₃₁*B₂[j]*F₃₂ + B₃[i]*F₃₃*B₂[j]*F₃₂ + B₁[i]*F₃₁*B₃[j]*F₃₃ + B₂[i]*F₃₂*B₃[j]*F₃₃)*Cᵢᵢⱼⱼ
                             +   ((B₁[i]*F₃₂+B₂[i]*F₃₁)*(B₁[j]*F₃₂+B₂[j]*F₃₁) + (B₃[i]*F₃₁+B₁[i]*F₃₃)*(B₃[j]*F₃₁+B₁[j]*F₃₃) + (B₃[i]*F₃₂+B₂[i]*F₃₃)*(B₃[j]*F₃₂+B₂[j]*F₃₃))*Cᵢⱼᵢⱼ
                             +   B₁[i]*B₁[j]*S₁₁+B₂[i]*B₂[j]*S₂₂+B₃[i]*B₃[j]*S₃₃ + (B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂ + (B₁[i]*B₃[j]+B₃[i]*B₁[j])*S₁₃ + (B₃[i]*B₂[j]+B₂[i]*B₃[j])*S₂₃)*𝑤




            end
        end
    end
end





function ∫∫EᵢⱼSᵢⱼdxdy_SaintVenantKirchhoff(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    
    for ξ in 𝓖
        E = ξ.Ē
        ν = ξ.ν̄ 
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        Cᵢᵢᵢᵢ = E*(1-ν)/(1+ν)/(1-2*ν)
        Cᵢᵢⱼⱼ = E*ν/(1+ν)/(1-2*ν)
        Cᵢⱼᵢⱼ = E/(1+ν)/2
       
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        
        E₁₁ = 0.5*(F₁₁*F₁₁+F₂₁*F₂₁-1.0)
        E₁₂ = 0.5*(F₁₁*F₁₂+F₂₁*F₂₂)
        E₂₂ = 0.5*(F₁₂*F₁₂+F₂₂*F₂₂-1.0)
        S₁₁ = Cᵢᵢᵢᵢ*E₁₁+Cᵢᵢⱼⱼ*E₂₂
        S₂₂ = Cᵢᵢⱼⱼ*E₁₁+Cᵢᵢᵢᵢ*E₂₂
        S₁₂ = 2.0*Cᵢⱼᵢⱼ*E₁₂
        
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += (B₁[i]*F₁₁*S₁₁+B₂[i]*F₁₂*S₂₂+(B₁[i]*F₁₂+B₂[i]*F₁₁)*S₁₂)*𝑤
            f[2*I]   += (B₁[i]*F₂₁*S₁₁+B₂[i]*F₂₂*S₂₂+(B₁[i]*F₂₂+B₂[i]*F₂₁)*S₁₂)*𝑤
        end
    end
end

function ∫∫∫EᵢⱼSᵢⱼdxdydz_SaintVenantKirchhoff(ap::T, f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    
    for ξ in 𝓖
        E = ξ.E
        ν = ξ.ν
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z]  # 3D: Z方向形函数梯度
        𝑤 = ξ.𝑤
        
        # --- 1. Lamé常数 ---
        λ = E*ν/((1+ν)*(1-2*ν)) 
        μ = E/(2*(1+ν))           
        Cᵢᵢᵢᵢ = λ + 2*μ           
        
        # --- 2. 计算 3x3 变形梯度 F ---
        F₁₁ = 1.0; F₁₂ = 0.0; F₁₃ = 0.0
        F₂₁ = 0.0; F₂₂ = 1.0; F₂₃ = 0.0
        F₃₁ = 0.0; F₃₂ = 0.0; F₃₃ = 1.0
        
        for (i,xᵢ) in enumerate(𝓒)
            F₁₁ += B₁[i]*xᵢ.d₁; F₁₂ += B₂[i]*xᵢ.d₁; F₁₃ += B₃[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂; F₂₂ += B₂[i]*xᵢ.d₂; F₂₃ += B₃[i]*xᵢ.d₂
            F₃₁ += B₁[i]*xᵢ.d₃; F₃₂ += B₂[i]*xᵢ.d₃; F₃₃ += B₃[i]*xᵢ.d₃
        end
        
        # --- 3. 计算 Green-Lagrange 应变 E = 0.5*(FᵀF - I) ---
        E₁₁ = 0.5*(F₁₁*F₁₁+F₂₁*F₂₁+F₃₁*F₃₁-1.0)
        E₂₂ = 0.5*(F₁₂*F₁₂+F₂₂*F₂₂+F₃₂*F₃₂-1.0)
        E₃₃ = 0.5*(F₁₃*F₁₃+F₂₃*F₂₃+F₃₃*F₃₃-1.0)
        E₁₂ = 0.5*(F₁₁*F₁₂+F₂₁*F₂₂+F₃₁*F₃₂)
        E₁₃ = 0.5*(F₁₁*F₁₃+F₂₁*F₂₃+F₃₁*F₃₃)
        E₂₃ = 0.5*(F₁₂*F₁₃+F₂₂*F₂₃+F₃₂*F₃₃)

        # --- 4. 计算第二 Piola-Kirchhoff 应力 S ---
        S₁₁ = Cᵢᵢᵢᵢ*E₁₁ + λ*E₂₂ + λ*E₃₃
        S₂₂ = λ*E₁₁ + Cᵢᵢᵢᵢ*E₂₂ + λ*E₃₃
        S₃₃ = λ*E₁₁ + λ*E₂₂ + Cᵢᵢᵢᵢ*E₃₃
        S₁₂ = 2.0*μ*E₁₂ 
        S₁₃ = 2.0*μ*E₁₃
        S₂₃ = 2.0*μ*E₂₃
        # S₂₁=S₁₂, S₃₁=S₁₃, S₃₂=S₂₃

        # --- 5. 内力向量 f 组装 ---
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            
            f[3*I-2] += (
                B₁[i] * (F₁₁*S₁₁ + F₁₂*S₁₂ + F₁₃*S₁₃) +
                B₂[i] * (F₁₁*S₁₂ + F₁₂*S₂₂ + F₁₃*S₂₃) +
               B₃[i] * (F₁₁*S₁₃ + F₁₂*S₂₃ + F₁₃*S₃₃)
            ) * 𝑤
         
            f[3*I-1] += (
                B₁[i] * (F₂₁*S₁₁ + F₂₂*S₁₂ + F₂₃*S₁₃) +
                B₂[i] * (F₂₁*S₁₂ + F₂₂*S₂₂ + F₂₃*S₂₃) +
                B₃[i] * (F₂₁*S₁₃ + F₂₂*S₂₃ + F₂₃*S₃₃)
            ) * 𝑤

            f[3*I] += (
               B₁[i] * (F₃₁*S₁₁ + F₃₂*S₁₂ + F₃₃*S₁₃) +
               B₂[i] * (F₃₁*S₁₂ + F₃₂*S₂₂ + F₃₃*S₂₃) +
               B₃[i] * (F₃₁*S₁₃ + F₃₂*S₂₃ + F₃₃*S₃₃)
            ) * 𝑤
        end
    end
end

function Δ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
   
    for ξ in 𝓖
         E=ξ.Ē
         ν=ξ.ν̄ 
         λ=E*ν/(1+ν)/(1-2*ν)
         μ=E/(1+ν)/2
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        C₁₁₁₁=λ*J*(2*J-1.0)*C⁻¹₁₁*C⁻¹₁₁+(μ- λ*J*(J-1.0))*2*C⁻¹₁₁*C⁻¹₁₁
        C₂₂₂₂=λ*J*(2*J-1.0)*C⁻¹₂₂*C⁻¹₂₂+(μ- λ*J*(J-1.0))*2*C⁻¹₂₂*C⁻¹₂₂
        C₁₁₂₂=λ*J*(2*J-1.0)*C⁻¹₁₁*C⁻¹₂₂+(μ- λ*J*(J-1.0))*2*C⁻¹₁₂*C⁻¹₁₂
        C₁₁₁₂=λ*J*(2*J-1.0)*C⁻¹₁₁*C⁻¹₁₂+(μ- λ*J*(J-1.0))*2*C⁻¹₁₁*C⁻¹₁₂ 
        C₂₂₁₂=λ*J*(2*J-1.0)*C⁻¹₂₂*C⁻¹₁₂+(μ- λ*J*(J-1.0))*2*C⁻¹₁₂*C⁻¹₂₂     
        C₁₂₁₂=λ*J*(2*J-1.0)*C⁻¹₁₂*C⁻¹₁₂+(μ- λ*J*(J-1.0))*(C⁻¹₁₁*C⁻¹₂₂+C⁻¹₁₂*C⁻¹₁₂)   

        S₁₁ = λ*J*(J-1)*C⁻¹₁₁ + μ*(1-C⁻¹₁₁)
        S₂₂ = λ*J*(J-1)*C⁻¹₂₂ + μ*(1-C⁻¹₂₂)
        S₁₂ = λ*J*(J-1)*C⁻¹₁₂ - μ*C⁻¹₁₂
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] +=(B₁[i]*F₁₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               + (B₁[i]*F₁₁*B₂[j]*F₁₂ + B₂[i]*F₁₂*B₁[j]*F₁₁)*C₁₁₂₂
                               + (B₁[i]*F₁₁*B₂[j]*F₁₁ + B₁[i]*F₁₁*B₁[j]*F₁₂)*C₁₁₁₂
                               + (B₁[i]*F₁₂*B₁[j]*F₁₁ + B₂[i]*F₁₁*B₁[j]*F₁₁)*C₁₁₁₂
                               + (B₂[i]*F₁₂*B₂[j]*F₁₁ + B₂[i]*F₁₂*B₁[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂*B₂[j]*F₁₂ + B₂[i]*F₁₁*B₂[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂
                               +  B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤
                              
                k[2*I-1,2*J]   += (B₁[i]*F₁₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₂₂*C₂₂₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₂ + B₂[i]*F₁₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₁ + B₁[i]*F₁₁*B₁[j]*F₂₂)*C₁₁₁₂ 
                               +  (B₁[i]*F₁₂*B₁[j]*F₂₁ + B₂[i]*F₁₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₁₂*B₂[j]*F₂₁ + B₂[i]*F₁₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂*B₂[j]*F₂₂ + B₂[i]*F₁₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂)*𝑤
                               
                k[2*I,2*J-1]   += (B₁[i]*F₂₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₂ + B₂[i]*F₂₂*B₁[j]*F₁₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₁ + B₁[i]*F₂₁*B₁[j]*F₁₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₁₁ + B₂[i]*F₂₁*B₁[j]*F₁₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₁₁ + B₂[i]*F₂₂*B₁[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₁₂ + B₂[i]*F₂₁*B₂[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂)*𝑤

                k[2*I,2*J]     += (B₁[i]*F₂₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₂₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₂ + B₂[i]*F₂₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₁ + B₁[i]*F₂₁*B₁[j]*F₂₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₂₁ + B₂[i]*F₂₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₂₁ + B₂[i]*F₂₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₂₂ + B₂[i]*F₂₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂
                               +   B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤

            end
        end
    end
end

# function Δ∫∫∫EᵢⱼSᵢⱼdxdydz_NeoHookean(ap::T, k::AbstractMatrix{Float64}) where T<:AbstractElement
#     𝓒 = ap.𝓒; 𝓖 = ap.𝓖
   
#     for ξ in 𝓖
#         E = ξ.E
#         ν = ξ.ν
#         B₁ = ξ[:∂𝝭∂x]
#         B₂ = ξ[:∂𝝭∂y]
#         B₃ = ξ[:∂𝝭∂z]  # 3D: Z方向形函数梯度
#         𝑤 = ξ.𝑤
        
#         # --- 1. Lamé常数 (用于SVK应力计算和K_M) ---
#         λ = E*ν/((1+ν)*(1-2*ν))  # Cᵢᵢⱼⱼ
#         μ = E/(2*(1+ν))          # Cᵢⱼᵢⱼ
        
        
#        Cᵢᵢᵢᵢ = λ + 2*μ
#        Cᵢᵢⱼⱼ = λ
#        Cᵢⱼᵢⱼ = μ
#         # --- 2. 计算 3x3 变形梯度 F ---
#         F₁₁ = 1.0; F₁₂ = 0.0; F₁₃ = 0.0
#         F₂₁ = 0.0; F₂₂ = 1.0; F₂₃ = 0.0
#         F₃₁ = 0.0; F₃₂ = 0.0; F₃₃ = 1.0
        
#         for (i,xᵢ) in  enumerate(𝓒)
#             # F_ij = δ_ij + (∂u_i / ∂X_j) = δ_ij + (∂N_k / ∂X_j) * d_i^k
#             F₁₁ += B₁[i]*xᵢ.d₁; F₁₂ += B₂[i]*xᵢ.d₁; F₁₃ += B₃[i]*xᵢ.d₁
#             F₂₁ += B₁[i]*xᵢ.d₂; F₂₂ += B₂[i]*xᵢ.d₂; F₂₃ += B₃[i]*xᵢ.d₂
#             F₃₁ += B₁[i]*xᵢ.d₃; F₃₂ += B₂[i]*xᵢ.d₃; F₃₃ += B₃[i]*xᵢ.d₃
#         end
        
#         # --- 3. 计算 Green-Lagrange 应变 E = 0.5*(FᵀF - I) ---
#         E₁₁ = 0.5*(F₁₁*F₁₁+F₂₁*F₂₁+F₃₁*F₃₁-1.0)
#         E₂₂ = 0.5*(F₁₂*F₁₂+F₂₂*F₂₂+F₃₂*F₃₂-1.0)
#         E₃₃ = 0.5*(F₁₃*F₁₃+F₂₃*F₂₃+F₃₃*F₃₃-1.0) # 3D 新增
#         E₁₂ = 0.5*(F₁₁*F₁₂+F₂₁*F₂₂+F₃₁*F₃₂)
#         E₁₃ = 0.5*(F₁₁*F₁₃+F₂₁*F₂₃+F₃₁*F₃₃)    # 3D 新增
#         E₂₃ = 0.5*(F₁₂*F₁₃+F₂₂*F₂₃+F₃₂*F₃₃)    # 3D 新增

#         # --- 4. 计算第二 Piola-Kirchhoff 应力 S = C:E (SVK模型) ---
#        S₁₁ = λ*(E₁₁ + E₂₂ + E₃₃) + 2*μ*E₁₁
#        S₂₂ = λ*(E₁₁ + E₂₂ + E₃₃) + 2*μ*E₂₂
#        S₃₃ = λ*(E₁₁ + E₂₂ + E₃₃) + 2*μ*E₃₃
#        S₁₂ = 2*μ*E₁₂
#        S₁₃ = 2*μ*E₁₃
#        S₂₃ = 2*μ*E₂₃ 
        
       
        
#        for (i,xᵢ) in enumerate(𝓒)
#             I = xᵢ.𝐼
#             for (j,xⱼ) in enumerate(𝓒)
#                 J = xⱼ.𝐼
#                 k[3*I-2,3*J-2] +=((B₁[i]*F₁₁*B₁[j]*F₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂ + B₃[i]*F₁₃*B₃[j]*F₁₃)*Cᵢᵢᵢᵢ
#                                +  (B₂[i]*F₁₂*B₁[j]*F₁₁ + B₃[i]*F₁₃*B₁[j]*F₁₁ + B₁[i]*F₁₁*B₂[j]*F₁₂ + B₃[i]*F₁₃*B₂[j]*F₁₂ + B₁[i]*F₁₁*B₃[j]*F₁₃ + B₂[i]*F₁₂*B₃[j]*F₁₃)*Cᵢᵢⱼⱼ
#                                +  ((B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)+(B₃[i]*F₁₁+B₁[i]*F₁₃)*(B₃[j]*F₁₁+B₁[j]*F₁₃) + (B₃[i]*F₁₂+B₂[i]*F₁₃)*(B₃[j]*F₁₂+B₂[j]*F₁₃) )*Cᵢⱼᵢⱼ
#                                +   B₁[i]*B₁[j]*S₁₁+B₂[i]*B₂[j]*S₂₂+B₃[i]*B₃[j]*S₃₃ + (B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂ + (B₁[i]*B₃[j]+B₃[i]*B₁[j])*S₁₃ + (B₃[i]*B₂[j]+B₂[i]*B₃[j])*S₂₃)*𝑤
#                  #  1,2             
#                 k[3*I-2,3*J-1] += ((B₁[i]*F₁₁*B₁[j]*F₂₁ + B₂[i]*F₁₂*B₂[j]*F₂₂ + B₃[i]*F₁₃*B₃[j]*F₂₃)*Cᵢᵢᵢᵢ
#                                +   (B₂[i]*F₁₂*B₁[j]*F₂₁ + B₃[i]*F₁₃*B₁[j]*F₂₁ + B₁[i]*F₁₁*B₂[j]*F₂₂ + B₃[i]*F₁₃*B₂[j]*F₂₂ + B₁[i]*F₁₁*B₃[j]*F₂₃ + B₂[i]*F₁₂*B₃[j]*F₂₃)*Cᵢᵢⱼⱼ
#                                +   ((B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁) + (B₃[i]*F₁₁+B₁[i]*F₁₃)*(B₃[j]*F₂₁+B₁[j]*F₂₃) + (B₃[i]*F₁₂+B₂[i]*F₁₃)*(B₃[j]*F₂₂+B₂[j]*F₂₃))*Cᵢⱼᵢⱼ)*𝑤
#                  #  1,3                 
#                 k[3*I-2,3*J]   += ((B₁[i]*F₁₁*B₁[j]*F₃₁ + B₂[i]*F₁₂*B₂[j]*F₃₂ + B₃[i]*F₁₃*B₃[j]*F₃₃)*Cᵢᵢᵢᵢ
#                                +   (B₂[i]*F₁₂*B₁[j]*F₃₁ + B₃[i]*F₁₃*B₁[j]*F₃₁ + B₁[i]*F₁₁*B₂[j]*F₃₂ + B₃[i]*F₁₃*B₂[j]*F₃₂ + B₁[i]*F₁₁*B₃[j]*F₃₃ + B₂[i]*F₁₂*B₃[j]*F₃₃)*Cᵢᵢⱼⱼ
#                                +   ((B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₃₂+B₂[j]*F₃₁) + (B₃[i]*F₁₁+B₁[i]*F₁₃)*(B₃[j]*F₃₁+B₁[j]*F₃₃) + (B₃[i]*F₁₂+B₂[i]*F₁₃)*(B₃[j]*F₃₂+B₂[j]*F₃₃))*Cᵢⱼᵢⱼ)*𝑤
#                 #  2,1
#                 k[3*I-1,3*J-2] +=((B₁[i]*F₂₁*B₁[j]*F₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂ + B₃[i]*F₂₃*B₃[j]*F₁₃)*Cᵢᵢᵢᵢ
#                                +  (B₂[i]*F₂₂*B₁[j]*F₁₁ + B₃[i]*F₂₃*B₁[j]*F₁₁ + B₁[i]*F₂₁*B₂[j]*F₁₂ + B₃[i]*F₂₃*B₂[j]*F₁₂ + B₁[i]*F₂₁*B₃[j]*F₁₃ + B₂[i]*F₂₂*B₃[j]*F₁₃)*Cᵢᵢⱼⱼ
#                                +  ((B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)+(B₃[j]*F₁₁+B₁[j]*F₁₃)*(B₃[i]*F₂₁+B₁[i]*F₂₃) + (B₃[j]*F₁₂+B₂[j]*F₁₃)*(B₃[i]*F₂₂+B₂[i]*F₂₃) )*Cᵢⱼᵢⱼ)*𝑤
                                 
                
#             #    2,2
#                 k[3*I-1,3*J-1] += ((B₁[i]*F₂₁*B₁[j]*F₂₁ + B₂[i]*F₂₂*B₂[j]*F₂₂ + B₃[i]*F₂₃*B₃[j]*F₂₃)*Cᵢᵢᵢᵢ
#                                +  (B₂[i]*F₂₂*B₁[j]*F₂₁ + B₃[i]*F₂₃*B₁[j]*F₂₁ + B₁[i]*F₂₁*B₂[j]*F₂₂ + B₃[i]*F₂₃*B₂[j]*F₂₂ + B₁[i]*F₂₁*B₃[j]*F₂₃ + B₂[i]*F₂₂*B₃[j]*F₂₃)*Cᵢᵢⱼⱼ
#                                +  ((B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁) + (B₃[i]*F₂₁+B₁[i]*F₂₃)*(B₃[j]*F₂₁+B₁[j]*F₂₃) + (B₃[i]*F₂₂+B₂[i]*F₂₃)*(B₃[j]*F₂₂+B₂[j]*F₂₃))*Cᵢⱼᵢⱼ
#                                +  B₁[i]*B₁[j]*S₁₁+B₂[i]*B₂[j]*S₂₂+B₃[i]*B₃[j]*S₃₃ + (B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂ + (B₁[i]*B₃[j]+B₃[i]*B₁[j])*S₁₃ + (B₃[i]*B₂[j]+B₂[i]*B₃[j])*S₂₃)*𝑤
#             #   2,3
#                 k[3*I-1,3*J]   += ((B₁[i]*F₂₁*B₁[j]*F₃₁ + B₂[i]*F₂₂*B₂[j]*F₃₂ + B₃[i]*F₂₃*B₃[j]*F₃₃)*Cᵢᵢᵢᵢ
#                                +   (B₂[i]*F₂₂*B₁[j]*F₃₁ + B₃[i]*F₂₃*B₁[j]*F₃₁ + B₁[i]*F₂₁*B₂[j]*F₃₂ + B₃[i]*F₂₃*B₂[j]*F₃₂ + B₁[i]*F₂₁*B₃[j]*F₃₃ + B₂[i]*F₂₂*B₃[j]*F₃₃)*Cᵢᵢⱼⱼ
#                                +   ((B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₃₂+B₂[j]*F₃₁) + (B₃[i]*F₂₁+B₁[i]*F₂₃)*(B₃[j]*F₃₁+B₁[j]*F₃₃) + (B₃[i]*F₂₂+B₂[i]*F₂₃)*(B₃[j]*F₃₂+B₂[j]*F₃₃))*Cᵢⱼᵢⱼ)*𝑤
                             

#             # 3,1
#                 k[3*I,3*J-2]   += ((B₁[i]*F₃₁*B₁[j]*F₁₁ + B₂[i]*F₃₂*B₂[j]*F₁₂ + B₃[i]*F₃₃*B₃[j]*F₁₃)*Cᵢᵢᵢᵢ
#                                +   (B₂[i]*F₃₂*B₁[j]*F₁₁ + B₃[i]*F₃₃*B₁[j]*F₁₁ + B₁[i]*F₃₁*B₂[j]*F₁₂ + B₃[i]*F₃₃*B₂[j]*F₁₂ + B₁[i]*F₃₁*B₃[j]*F₁₃ + B₂[i]*F₃₂*B₃[j]*F₁₃)*Cᵢᵢⱼⱼ
#                                +   ((B₁[i]*F₃₂+B₂[i]*F₃₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁) + (B₃[i]*F₃₁+B₁[i]*F₃₃)*(B₃[j]*F₁₁+B₁[j]*F₁₃) + (B₃[i]*F₃₂+B₂[i]*F₃₃)*(B₃[j]*F₁₂+B₂[j]*F₁₃))*Cᵢⱼᵢⱼ)*𝑤
                             

#                 k[3*I,3*J-1]   += ((B₁[i]*F₃₁*B₁[j]*F₂₁ + B₂[i]*F₃₂*B₂[j]*F₂₂ + B₃[i]*F₃₃*B₃[j]*F₂₃)*Cᵢᵢᵢᵢ
#                                +   (B₂[i]*F₃₂*B₁[j]*F₂₁ + B₃[i]*F₃₃*B₁[j]*F₂₁ + B₁[i]*F₃₁*B₂[j]*F₂₂ + B₃[i]*F₃₃*B₂[j]*F₂₂ + B₁[i]*F₃₁*B₃[j]*F₂₃ + B₂[i]*F₃₂*B₃[j]*F₂₃)*Cᵢᵢⱼⱼ
#                                +   ((B₁[i]*F₃₂+B₂[i]*F₃₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁) + (B₃[i]*F₃₁+B₁[i]*F₃₃)*(B₃[j]*F₂₁+B₁[j]*F₂₃) + (B₃[i]*F₃₂+B₂[i]*F₃₃)*(B₃[j]*F₂₂+B₂[j]*F₂₃))*Cᵢⱼᵢⱼ)*𝑤
                            

#                 k[3*I,3*J]   += ((B₁[i]*F₃₁*B₁[j]*F₃₁ + B₂[i]*F₃₂*B₂[j]*F₃₂ + B₃[i]*F₃₃*B₃[j]*F₃₃)*Cᵢᵢᵢᵢ
#                              +   (B₂[i]*F₃₂*B₁[j]*F₃₁ + B₃[i]*F₃₃*B₁[j]*F₃₁ + B₁[i]*F₃₁*B₂[j]*F₃₂ + B₃[i]*F₃₃*B₂[j]*F₃₂ + B₁[i]*F₃₁*B₃[j]*F₃₃ + B₂[i]*F₃₂*B₃[j]*F₃₃)*Cᵢᵢⱼⱼ
#                              +   ((B₁[i]*F₃₂+B₂[i]*F₃₁)*(B₁[j]*F₃₂+B₂[j]*F₃₁) + (B₃[i]*F₃₁+B₁[i]*F₃₃)*(B₃[j]*F₃₁+B₁[j]*F₃₃) + (B₃[i]*F₃₂+B₂[i]*F₃₃)*(B₃[j]*F₃₂+B₂[j]*F₃₃))*Cᵢⱼᵢⱼ
#                              +   B₁[i]*B₁[j]*S₁₁+B₂[i]*B₂[j]*S₂₂+B₃[i]*B₃[j]*S₃₃ + (B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂ + (B₁[i]*B₃[j]+B₃[i]*B₁[j])*S₁₃ + (B₃[i]*B₂[j]+B₂[i]*B₃[j])*S₂₃)*𝑤




#             end
#         end
#     end
# end

function Δ∫∫∫EᵢⱼSᵢⱼdxdydz_NeoHookean(ap::T, k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    
    # 预先分配临时的 3x3 矩阵以避免循环中的 GC 分配
    F = zeros(3,3)
    C = zeros(3,3)
    Cinv = zeros(3,3)
    S = zeros(3,3)
    
    for ξ in 𝓖
        E = ξ.E
        ν = ξ.ν
        # Neo-Hookean 参数转换
        λ = E * ν / ((1 + ν) * (1 - 2 * ν))
        μ = E / ((1 + ν) * 2)

        # 获取形函数对 x, y, z 的导数
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z] # 需要单元提供 z 方向导数
        𝑤 = ξ.𝑤

        # 1. 计算变形梯度 F (3x3)
        # F_ij = δ_ij + Σ (∂N_k/∂X_j * u_k_i)
        fill!(F, 0.0)
        F[1,1] = 1.0; F[2,2] = 1.0; F[3,3] = 1.0
        
        for (i, xᵢ) in enumerate(𝓒)
            # 节点 i 的位移 d₁, d₂, d₃
            d₁ = xᵢ.d₁
            d₂ = xᵢ.d₂
            d₃ = xᵢ.d₃ # 假设节点数据结构中有 d₃
            
            F[1,1] += B₁[i] * d₁;  F[1,2] += B₂[i] * d₁;  F[1,3] += B₃[i] * d₁
            F[2,1] += B₁[i] * d₂;  F[2,2] += B₂[i] * d₂;  F[2,3] += B₃[i] * d₂
            F[3,1] += B₁[i] * d₃;  F[3,2] += B₂[i] * d₃;  F[3,3] += B₃[i] * d₃
        end

        # 2. 计算右 Cauchy-Green 张量 C = F^T * F
        mul!(C, F', F) # 使用 Julia 内置乘法，或手动展开

        # 3. 计算不变量 J 和 C 的逆
        # 为了性能和数值稳定性，建议使用 StaticArrays 或手动计算 3x3 逆
        detC = det(C)
        J = sqrt(detC) # J = det(F)
        
        # 手动计算 C 的逆 (Cinv) 或者使用 inv(C)
        # 这里使用 inv(C) 简化代码，高性能场景建议手动展开 Cramer 法则
        Cinv .= inv(C)

        # 4. 计算第二 Piola-Kirchhoff 应力 S
        # 对应二维公式: S = λJ(J-1)C⁻¹ + μ(I - C⁻¹)
        # 系数定义
        coeff_pres = λ * J * (J - 1.0)
        
        for i=1:3, j=1:3
            δ = (i == j) ? 1.0 : 0.0
            S[i,j] = coeff_pres * Cinv[i,j] + μ * (δ - Cinv[i,j])
        end

        # 5. 定义材料切线模量系数
        # 对应二维公式中的系数:
        # term1 coeff = λ*J*(2*J-1.0)
        # term2 coeff = 2*(μ - λ*J*(J-1.0))  <-- 注意二维代码中是2倍，对应 3D 张量导数
        
        c_vol = λ * J * (2.0 * J - 1.0)
        c_iso = 2.0 * (μ - λ * J * (J - 1.0)) 

        # 6. 组装刚度矩阵 k
        # 为了避免写出数千行的展开式，这里使用循环进行张量缩并
        # K_AB = K_mat + K_geo
        
        for (I, xI) in enumerate(𝓒)
            row_idx = xI.𝐼 # 节点全局编号
            # 形状函数导数向量 (B_I)
            bI = [B₁[I], B₂[I], B₃[I]] 
            
            for (J, xJ) in enumerate(𝓒)
                col_idx = xJ.𝐼
                bJ = [B₁[J], B₂[J], B₃[J]]

                # 计算几何刚度项 scalar (Geometric Stiffness)
                # K_geo = (B_I · S · B_J) * Identity
                s_geo = 0.0
                for p=1:3, q=1:3
                    s_geo += bI[p] * S[p,q] * bJ[q]
                end
                
                # 遍历自由度 (m: 行自由度 1..3, n: 列自由度 1..3)
                for m = 1:3
                    row = 3 * row_idx - 3 + m
                    for n = 1:3
                        col = 3 * col_idx - 3 + n
                        
                        # --- 几何刚度贡献 ---
                        val = (m == n) ? s_geo : 0.0
                        
                        # --- 材料刚度贡献 ---
                        # K_mat_mn = Σ (∂N_I/∂X_k * F_mk * C_klpq * F_nq * ∂N_J/∂X_p)
                        # C_klpq = c_vol * C⁻¹_kl * C⁻¹_pq + c_iso/2 * (C⁻¹_kq C⁻¹_lp + C⁻¹_kp C⁻¹_lq)
                        
                        term_mat = 0.0
                        for k=1:3, l=1:3, p=1:3, q=1:3
                            # 计算 C_klpq
                            # 注意：二维代码中的系数处理暗示了某种对称性假设
                            # 标准 3D 压缩 Neo-Hookean 导数如下：
                            
                            Isym = 0.5 * (Cinv[k,p]*Cinv[l,q] + Cinv[k,q]*Cinv[l,p])
                            C_ijkl = c_vol * Cinv[k,l] * Cinv[p,q] + c_iso * Isym # c_iso 包含了前面的系数
                            
                            # 缩并: bI[k] * F[m,l] * C_klpq * F[n,p] * bJ[q]
                            # 注意指标对应关系:
                            # 变形梯度的项是 F_mi (spatial m, material i)
                            # 链式法则: ∂x_m/∂X_k = F_mk
                            # 公式: B_I_k * F_mk * C_klpq * F_np * B_J_q
                            
                            term_mat += bI[k] * F[m,l] * C_ijkl * F[n,p] * bJ[q]
                        end
                        
                        val += term_mat
                        
                        # 累加到全局刚度矩阵
                        k[row, col] += val * 𝑤
                    end
                end
            end
        end
    end
end

function Δ∫∫∫_NeoHookean_Dev(ap::T, k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    F = zeros(3,3); C = zeros(3,3); Cinv = zeros(3,3); S_dev = zeros(3,3)

    for ξ in 𝓖
        E, ν, 𝑤 = ξ.E, ξ.ν, ξ.𝑤
        μ = E / (2 * (1 + ν))
        B = [ξ[:∂𝝭∂x] ξ[:∂𝝭∂y] ξ[:∂𝝭∂z]] # 假设为 N_node x 3 矩阵

        # 1. 计算 F, C, J
        fill!(F, 0.0); for i=1:3; F[i,i]=1.0; end
        for (i, xᵢ) in enumerate(𝓒)
            d = [xᵢ.d₁, xᵢ.d₂, xᵢ.d₃]
            F += d * B[i, :]' 
        end
        mul!(C, F', F)
        J = sqrt(det(C))
        Cinv .= inv(C)

        # 2. 偏量第二 PK 应力: S_dev = μ * J^(-2/3) * (I - 1/3 * tr(C)*Cinv)
        # 注意：这里采用了标准 Neo-Hookean 偏量简化形式
        trC = tr(C)
        J23inv = J^(-2/3)
        for i=1:3, j=1:3
            δ = (i == j) ? 1.0 : 0.0
            S_dev[i,j] = μ * J23inv * (δ - (1/3) * trC * Cinv[i,j])
        end

        # 3. 偏量切线模量系数 c_iso
        c_iso = 2.0 * μ * J23inv 

        # 4. 组装 (简化循环逻辑)
        assemble_stiffness!(k, 𝓒, B, F, Cinv, S_dev, c_iso, 0.0, 𝑤) 
    end
end

function Δ∫∫∫_NeoHookean_Vol(ap::T, k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖_low = ap.𝓖_reduced # 使用降阶积分点
    F = zeros(3,3); C = zeros(3,3); Cinv = zeros(3,3); S_vol = zeros(3,3)

    for ξ in 𝓖_low
        E, ν, 𝑤 = ξ.E, ξ.ν, ξ.𝑤
        λ = E * ν / ((1 + ν) * (1 - 2 * ν))
        B = [ξ[:∂𝝭∂x] ξ[:∂𝝭∂y] ξ[:∂𝝭∂z]]

        # 1. 计算 F, J, Cinv
        fill!(F, 0.0); for i=1:3; F[i,i]=1.0; end
        for (i, xᵢ) in enumerate(𝓒)
            d = [xᵢ.d₁, xᵢ.d₂, xᵢ.d₃]
            F += d * B[i, :]'
        end
        detF = det(F)
        Cinv .= inv(F' * F)

        # 2. 体积 PK 应力: S_vol = p * J * Cinv (其中 p = dU/dJ)
        # Neo-Hookean 常用形式: U(J) = 1/2 * λ * (J-1)²
        p = λ * (detF - 1.0) 
        S_vol .= (p * detF) .* Cinv

        # 3. 体积切线系数 c_vol
        # c_vol = J * (p + J * p') = λ * J * (2J - 1)
        c_vol = λ * detF * (2.0 * detF - 1.0)

        # 4. 组装
        assemble_stiffness!(k, 𝓒, B, F, Cinv, S_vol, 0.0, c_vol, 𝑤)
    end
end

function Δ∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E=op.E
    ν=op.ν
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        C₁₁₁₁=λ*J*(2*J-1.0)*C⁻¹₁₁*C⁻¹₁₁- λ*J*(J-1.0)*2*C⁻¹₁₁*C⁻¹₁₁
        C₂₂₂₂=λ*J*(2*J-1.0)*C⁻¹₂₂*C⁻¹₂₂- λ*J*(J-1.0)*2*C⁻¹₂₂*C⁻¹₂₂
        C₁₁₂₂=λ*J*(2*J-1.0)*C⁻¹₁₁*C⁻¹₂₂- λ*J*(J-1.0)*2*C⁻¹₁₂*C⁻¹₁₂
        C₁₁₁₂=λ*J*(2*J-1.0)*C⁻¹₁₁*C⁻¹₁₂- λ*J*(J-1.0)*2*C⁻¹₁₁*C⁻¹₁₂ 
        C₂₂₁₂=λ*J*(2*J-1.0)*C⁻¹₂₂*C⁻¹₁₂- λ*J*(J-1.0)*2*C⁻¹₁₂*C⁻¹₂₂     
        C₁₂₁₂=λ*J*(2*J-1.0)*C⁻¹₁₂*C⁻¹₁₂- λ*J*(J-1.0)*(C⁻¹₁₁*C⁻¹₂₂+C⁻¹₁₂*C⁻¹₁₂)   

        S₁₁ = λ*J*(J-1)*C⁻¹₁₁
        S₂₂ = λ*J*(J-1)*C⁻¹₂₂
        S₁₂ = λ*J*(J-1)*C⁻¹₁₂
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] +=(B₁[i]*F₁₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               + (B₁[i]*F₁₁*B₂[j]*F₁₂ + B₂[i]*F₁₂*B₁[j]*F₁₁)*C₁₁₂₂
                               + (B₁[i]*F₁₁*B₂[j]*F₁₁ + B₁[i]*F₁₁*B₁[j]*F₁₂)*C₁₁₁₂
                               + (B₁[i]*F₁₂*B₁[j]*F₁₁ + B₂[i]*F₁₁*B₁[j]*F₁₁)*C₁₁₁₂
                               + (B₂[i]*F₁₂*B₂[j]*F₁₁ + B₂[i]*F₁₂*B₁[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂*B₂[j]*F₁₂ + B₂[i]*F₁₁*B₂[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂
                               +  B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤
                              
                k[2*I-1,2*J]   += (B₁[i]*F₁₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₂₂*C₂₂₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₂ + B₂[i]*F₁₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₁ + B₁[i]*F₁₁*B₁[j]*F₂₂)*C₁₁₁₂ 
                               +  (B₁[i]*F₁₂*B₁[j]*F₂₁ + B₂[i]*F₁₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₁₂*B₂[j]*F₂₁ + B₂[i]*F₁₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂*B₂[j]*F₂₂ + B₂[i]*F₁₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂)*𝑤
                               
                k[2*I,2*J-1]   += (B₁[i]*F₂₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₂ + B₂[i]*F₂₂*B₁[j]*F₁₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₁ + B₁[i]*F₂₁*B₁[j]*F₁₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₁₁ + B₂[i]*F₂₁*B₁[j]*F₁₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₁₁ + B₂[i]*F₂₂*B₁[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₁₂ + B₂[i]*F₂₁*B₂[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂)*𝑤

                k[2*I,2*J]     += (B₁[i]*F₂₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₂₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₂ + B₂[i]*F₂₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₁ + B₁[i]*F₂₁*B₁[j]*F₂₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₂₁ + B₂[i]*F₂₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₂₁ + B₂[i]*F₂₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₂₂ + B₂[i]*F₂₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂
                               +   B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤

            end
        end
    end
end

function Δ∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E=op.E
    ν=op.ν
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        C₁₁₁₁=μ*2*C⁻¹₁₁*C⁻¹₁₁
        C₂₂₂₂=μ*2*C⁻¹₂₂*C⁻¹₂₂
        C₁₁₂₂=μ*2*C⁻¹₁₂*C⁻¹₁₂
        C₁₁₁₂=μ*2*C⁻¹₁₁*C⁻¹₁₂ 
        C₂₂₁₂=μ*2*C⁻¹₁₂*C⁻¹₂₂     
        C₁₂₁₂=μ*(C⁻¹₁₁*C⁻¹₂₂+C⁻¹₁₂*C⁻¹₁₂)   

        S₁₁ = μ*(1-C⁻¹₁₁)
        S₂₂ = μ*(1-C⁻¹₂₂)
        S₁₂ = - μ*C⁻¹₁₂
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] +=(B₁[i]*F₁₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               + (B₁[i]*F₁₁*B₂[j]*F₁₂ + B₂[i]*F₁₂*B₁[j]*F₁₁)*C₁₁₂₂
                               + (B₁[i]*F₁₁*B₂[j]*F₁₁ + B₁[i]*F₁₁*B₁[j]*F₁₂)*C₁₁₁₂
                               + (B₁[i]*F₁₂*B₁[j]*F₁₁ + B₂[i]*F₁₁*B₁[j]*F₁₁)*C₁₁₁₂
                               + (B₂[i]*F₁₂*B₂[j]*F₁₁ + B₂[i]*F₁₂*B₁[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂*B₂[j]*F₁₂ + B₂[i]*F₁₁*B₂[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂
                               +  B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤
                              
                k[2*I-1,2*J]   += (B₁[i]*F₁₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₂₂*C₂₂₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₂ + B₂[i]*F₁₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₁ + B₁[i]*F₁₁*B₁[j]*F₂₂)*C₁₁₁₂ 
                               +  (B₁[i]*F₁₂*B₁[j]*F₂₁ + B₂[i]*F₁₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₁₂*B₂[j]*F₂₁ + B₂[i]*F₁₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂*B₂[j]*F₂₂ + B₂[i]*F₁₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂)*𝑤
                               
                k[2*I,2*J-1]   += (B₁[i]*F₂₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₂ + B₂[i]*F₂₂*B₁[j]*F₁₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₁ + B₁[i]*F₂₁*B₁[j]*F₁₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₁₁ + B₂[i]*F₂₁*B₁[j]*F₁₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₁₁ + B₂[i]*F₂₂*B₁[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₁₂ + B₂[i]*F₂₁*B₂[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂)*𝑤

                k[2*I,2*J]     += (B₁[i]*F₂₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₂₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₂ + B₂[i]*F₂₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₁ + B₁[i]*F₂₁*B₁[j]*F₂₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₂₁ + B₂[i]*F₂₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₂₁ + B₂[i]*F₂₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₂₂ + B₂[i]*F₂₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂
                               +   B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤

            end
        end
    end
end

function ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
  
  
    for ξ in 𝓖
        E=ξ.Ē
        ν=ξ.ν̄ 
        λ=E*ν/(1+ν)/(1-2*ν)
        μ=E/(1+ν)/2
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        S₁₁ = λ*J*(J-1)*C⁻¹₁₁ + μ*(1-C⁻¹₁₁)
        S₂₂ = λ*J*(J-1)*C⁻¹₂₂ + μ*(1-C⁻¹₂₂)
        S₁₂ = λ*J*(J-1)*C⁻¹₁₂ - μ*C⁻¹₁₂
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += (B₁[i]*F₁₁*S₁₁+B₂[i]*F₁₂*S₂₂+(B₁[i]*F₁₂+B₂[i]*F₁₁)*S₁₂)*𝑤
            f[2*I]   += (B₁[i]*F₂₁*S₁₁+B₂[i]*F₂₂*S₂₂+(B₁[i]*F₂₂+B₂[i]*F₂₁)*S₁₂)*𝑤
        end
    end
end



function ∫∫∫EᵢⱼSᵢⱼdxdydz_NeoHookean(ap::T, f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    
    # 预分配 3x3 矩阵以提升性能 (避免循环内频繁内存分配)
    F = zeros(3,3)
    C = zeros(3,3)
    S = zeros(3,3)
    P = zeros(3,3) # 第一 Piola-Kirchhoff 应力 P = FS

    for ξ in 𝓖
        E = ξ.E
        ν = ξ.ν
        # Neo-Hookean 材料参数
        λ = E * ν / ((1 + ν) * (1 - 2 * ν))
        μ = E / ((1 + ν) * 2)

        # 获取形函数导数 (增加 z 方向 B₃)
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z] 
        𝑤 = ξ.𝑤

        # 1. 计算变形梯度 F (3x3)
        # F = I + ∇u
        fill!(F, 0.0)
        F[1,1] = 1.0; F[2,2] = 1.0; F[3,3] = 1.0
        
        for (i, xᵢ) in enumerate(𝓒)
            d₁ = xᵢ.d₁
            d₂ = xᵢ.d₂
            d₃ = xᵢ.d₃ # 假设节点包含 d₃
            
            # F_ij += d_i * ∂N/∂X_j
            F[1,1] += d₁ * B₁[i]; F[1,2] += d₁ * B₂[i]; F[1,3] += d₁ * B₃[i]
            F[2,1] += d₂ * B₁[i]; F[2,2] += d₂ * B₂[i]; F[2,3] += d₂ * B₃[i]
            F[3,1] += d₃ * B₁[i]; F[3,2] += d₃ * B₂[i]; F[3,3] += d₃ * B₃[i]
        end

        # 2. 计算右 Cauchy-Green 张量 C = F' * F
        mul!(C, F', F) 

        # 3. 计算不变量和逆矩阵
        J = det(F)       # 3D Jacobian
        detC = det(C)    # det(C) = J^2
        Cinv = inv(C)    # 3x3 逆矩阵

        # 4. 计算第二 Piola-Kirchhoff 应力 S
        # 公式: S = λJ(J-1)C⁻¹ + μ(I - C⁻¹)
        coeff_1 = λ * J * (J - 1.0)
        
        for i in 1:3, j in 1:3
            δ = (i == j) ? 1.0 : 0.0
            S[i,j] = coeff_1 * Cinv[i,j] + μ * (δ - Cinv[i,j])
        end

        # 5. 计算第一 Piola-Kirchhoff 应力 P
        # P = F * S (用于将应力推回当前构型，简化力向量组装)
        mul!(P, F, S)

        # 6. 组装内力向量 f
        # f_node = ∫ P ⋅ ∇N dV
        # f_{m,I} += P_{mj} * B_{j,I} * w
        
        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼 # 节点全局编号
            
            # 提取该节点的形函数导数向量 ∇N_i
            dN_dX = B₁[i]
            dN_dY = B₂[i]
            dN_dZ = B₃[i]

            # 自由度 1 (x方向力): P₁₁*Nx + P₁₂*Ny + P₁₃*Nz
            val_x = P[1,1]*dN_dX + P[1,2]*dN_dY + P[1,3]*dN_dZ
            
            # 自由度 2 (y方向力): P₂₁*Nx + P₂₂*Ny + P₂₃*Nz
            val_y = P[2,1]*dN_dX + P[2,2]*dN_dY + P[2,3]*dN_dZ
            
            # 自由度 3 (z方向力): P₃₁*Nx + P₃₂*Ny + P₃₃*Nz
            val_z = P[3,1]*dN_dX + P[3,2]*dN_dY + P[3,3]*dN_dZ

            # 累加到全局力向量 (注意索引: 3*I-2, 3*I-1, 3*I)
            f[3*I-2] += val_x * 𝑤
            f[3*I-1] += val_y * 𝑤
            f[3*I]   += val_z * 𝑤
        end
    end
end

function ∫∫∫_NeoHookean_Force_Dev(ap::T, f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    F = zeros(3,3); C = zeros(3,3); S_dev = zeros(3,3); P_dev = zeros(3,3)

    for ξ in 𝓖
        E, ν, 𝑤 = ξ.E, ξ.ν, ξ.𝑤
        μ = E / (2 * (1 + ν))
        B₁ = ξ[:∂𝝭∂x]; B₂ = ξ[:∂𝝭∂y]; B₃ = ξ[:∂𝝭∂z]

        # 1. 计算变形梯度 F
        fill!(F, 0.0); for i=1:3; F[i,i]=1.0; end
        for (i, xᵢ) in enumerate(𝓒)
            F[1,1]+=xᵢ.d₁*B₁[i]; F[1,2]+=xᵢ.d₁*B₂[i]; F[1,3]+=xᵢ.d₁*B₃[i]
            F[2,1]+=xᵢ.d₂*B₁[i]; F[2,2]+=xᵢ.d₂*B₂[i]; F[2,3]+=xᵢ.d₂*B₃[i]
            F[3,1]+=xᵢ.d₃*B₁[i]; F[3,2]+=xᵢ.d₃*B₂[i]; F[3,3]+=xᵢ.d₃*B₃[i]
        end

        # 2. 计算偏量 PK2 应力 S_dev
        mul!(C, F', F)
        J = det(F)
        J23inv = J^(-2/3)
        Cinv = inv(C)
        trC = tr(C)
        
        # S_dev = μ * J^(-2/3) * (I - 1/3 * tr(C)*Cinv)
        for i=1:3, j=1:3
            δ = (i == j) ? 1.0 : 0.0
            S_dev[i,j] = μ * J23inv * (δ - (1/3) * trC * Cinv[i,j])
        end

        # 3. 转换为 PK1 并组装
        mul!(P_dev, F, S_dev)
        assemble_force!(f, 𝓒, P_dev, B₁, B₂, B₃, 𝑤)
    end
end

function ∫∫∫_NeoHookean_Force_Vol(ap::T, f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖_low = ap.𝓖_reduced # 使用降阶积分点
    F = zeros(3,3); P_vol = zeros(3,3)

    for ξ in 𝓖_low
        E, ν, 𝑤 = ξ.E, ξ.ν, ξ.𝑤
        λ = E * ν / ((1 + ν) * (1 - 2 * ν))
        B₁ = ξ[:∂𝝭∂x]; B₂ = ξ[:∂𝝭∂y]; B₃ = ξ[:∂𝝭∂z]

        # 1. 计算 F 和 J
        fill!(F, 0.0); for i=1:3; F[i,i]=1.0; end
        for (i, xᵢ) in enumerate(𝓒)
            F[1,1]+=xᵢ.d₁*B₁[i]; F[1,2]+=xᵢ.d₁*B₂[i]; F[1,3]+=xᵢ.d₁*B₃[i]
            F[2,1]+=xᵢ.d₂*B₁[i]; F[2,2]+=xᵢ.d₂*B₂[i]; F[2,3]+=xᵢ.d₂*B₃[i]
            F[3,1]+=xᵢ.d₃*B₁[i]; F[3,2]+=xᵢ.d₃*B₂[i]; F[3,3]+=xᵢ.d₃*B₃[i]
        end
        J = det(F)

        # 2. 计算体积项静压力 p 和 PK1 应力 P_vol
        # p = dU/dJ = λ(J-1)
        # P_vol = p * J * F⁻ᵀ (体积项在当前构型的贡献)
        p = λ * (J - 1.0)
        FinvT = inv(F)' # 3x3 逆转置
        P_vol .= (p * J) .* FinvT

        # 3. 组装
        assemble_force!(f, 𝓒, P_vol, B₁, B₂, B₃, 𝑤)
    end
end

function assemble_force!(f, 𝓒, P, B₁, B₂, B₃, 𝑤)
    for (i, xᵢ) in enumerate(𝓒)
        I = xᵢ.𝐼
        # 计算 P ⋅ ∇N (即 P_mj * N_i,j)
        val_x = P[1,1]*B₁[i] + P[1,2]*B₂[i] + P[1,3]*B₃[i]
        val_y = P[2,1]*B₁[i] + P[2,2]*B₂[i] + P[2,3]*B₃[i]
        val_z = P[3,1]*B₁[i] + P[3,2]*B₂[i] + P[3,3]*B₃[i]

        f[3*I-2] += val_x * 𝑤
        f[3*I-1] += val_y * 𝑤
        f[3*I]   += val_z * 𝑤
    end
end

function ∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E=op.E
    ν=op.ν
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        S₁₁ = λ*J*(J-1)*C⁻¹₁₁
        S₂₂ = λ*J*(J-1)*C⁻¹₂₂
        S₁₂ = λ*J*(J-1)*C⁻¹₁₂
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += (B₁[i]*F₁₁*S₁₁+B₂[i]*F₁₂*S₂₂+(B₁[i]*F₁₂+B₂[i]*F₁₁)*S₁₂)*𝑤
            f[2*I]   += (B₁[i]*F₂₁*S₁₁+B₂[i]*F₂₂*S₂₂+(B₁[i]*F₂₂+B₂[i]*F₂₁)*S₁₂)*𝑤
        end
    end
end

function ∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E=op.E
    ν=op.ν
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        S₁₁ = μ*(1-C⁻¹₁₁)
        S₂₂ = μ*(1-C⁻¹₂₂)
        S₁₂ = - μ*C⁻¹₁₂
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += (B₁[i]*F₁₁*S₁₁+B₂[i]*F₁₂*S₂₂+(B₁[i]*F₁₂+B₂[i]*F₁₁)*S₁₂)*𝑤
            f[2*I]   += (B₁[i]*F₂₁*S₁₁+B₂[i]*F₂₂*S₂₂+(B₁[i]*F₂₂+B₂[i]*F₂₁)*S₁₂)*𝑤
        end
    end
end

function ∫vᵢuᵢds(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
  
    for ξ in 𝓖
        α = ξ.α
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        n₁₁ = ξ.n₁₁
        n₂₂ = ξ.n₂₂
        n₁₂ = ξ.n₁₂
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
        end
        Δu₁ = g₁-u₁
        Δu₂ = g₂-u₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += α*N[i]*n₁₁*N[j]*𝑤
                k[2*I,2*J-1]   += α*N[i]*n₁₂*N[j]*𝑤
                k[2*I-1,2*J]   += α*N[i]*n₁₂*N[j]*𝑤
                k[2*I,2*J]     += α*N[i]*n₂₂*N[j]*𝑤
            end
            f[2*I-1] += α*N[i]*(n₁₁*Δu₁+n₁₂*Δu₂)*𝑤
            f[2*I]   += α*N[i]*(n₁₂*Δu₁+n₂₂*Δu₂)*𝑤
        end
    end
end

function ∫vᵢuᵢdΓ(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
  
    for ξ in 𝓖
        α   = ξ.α
        𝑤   = ξ.𝑤
        N   = ξ[:𝝭]

        # 3D 法向张量分量 n_a n_b
        n₁₁ = ξ.n₁₁
        n₂₂ = ξ.n₂₂
        n₃₃ = ξ.n₃₃
        n₁₂ = ξ.n₁₂
        n₁₃ = ξ.n₁₃
        n₂₃ = ξ.n₂₃

        g₁ = ξ.g₁
        g₂ = ξ.g₂
        g₃ = ξ.g₃

        u₁ = 0.0
        u₂ = 0.0
        u₃ = 0.0
        for (i, xᵢ) in enumerate(𝓒)
            u₁ += N[i] * xᵢ.d₁
            u₂ += N[i] * xᵢ.d₂
            u₃ += N[i] * xᵢ.d₃
        end

        Δu₁ = g₁ - u₁
        Δu₂ = g₂ - u₂
        Δu₃ = g₃ - u₃

        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
          
            for (j, xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼

                k[3*I-2,3*J-2] += α*N[i]*n₁₁*N[j]*𝑤
                k[3*I-2,3*J-1] += α*N[i]*n₁₂*N[j]*𝑤
                k[3*I-2,3*J]   += α*N[i]*n₁₃*N[j]*𝑤
                k[3*I-1,3*J-2] += α*N[i]*n₁₂*N[j]*𝑤
                k[3*I-1,3*J-1] += α*N[i]*n₂₂*N[j]*𝑤
                k[3*I-1,3*J]   += α*N[i]*n₂₃*N[j]*𝑤
                k[3*I,3*J-2]   += α*N[i]*n₁₃*N[j]*𝑤
                k[3*I,3*J-1]   += α*N[i]*n₂₃*N[j]*𝑤
                k[3*I,3*J]     += α*N[i]*n₃₃*N[j]*𝑤
            end

            f[3*I-2] += α * N[i] * (n₁₁*Δu₁ + n₁₂*Δu₂ + n₁₃*Δu₃) * 𝑤
            f[3*I-1] += α * N[i] * (n₁₂*Δu₁ + n₂₂*Δu₂ + n₂₃*Δu₃) * 𝑤
            f[3*I]   += α * N[i] * (n₁₃*Δu₁ + n₂₃*Δu₂ + n₃₃*Δu₃) * 𝑤
        end
    end
end


function ∫vᵢuᵢdΓ_center(ap::T, kα::AbstractMatrix{Float64}, fα::AbstractVector{Float64};
                          xc::Float64=0.5, yc::Float64=0.5, σ::Float64=0.02) where T<:AbstractElement
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖

    for ξ in 𝓖
        N = ξ[:𝝭]
        wΓ = ξ.𝑤
        α   = ξ.α
        # midpoint-localized weight (use reference coords from ξ.x, ξ.y)
        dx = ξ.x - xc
        dy = ξ.y - yc
        wloc = exp(-(dx*dx + dy*dy) / (σ*σ))

        # current displacement at this GP (ux, uy) from nodal d1,d2
        ux = 0.0
        uy = 0.0
        for (i, xᵢ) in enumerate(𝓒)
            ux += N[i] * xᵢ.d₁
            uy += N[i] * xᵢ.d₂
        end

        # penalty residual corresponds to: r += -α * ∫ w N^T u dΓ
        # (since you solve (K)Δd = (fext - fint + fα), and want fα = -Kα*u)
        for (i, xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            Ni = N[i]
            coef = α * wloc * Ni * wΓ

            # residual part (only ux, uy)
            fα[3*I-2] += -coef * ux
            fα[3*I-1] += -coef * uy
            # fα[3*I] += -coef * uy
            # stiffness part (only ux, uy blocks)
            for (j, xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                Nj = N[j]
                kij = α * wloc * Ni * Nj * wΓ

                # x-x
                kα[3*I-2, 3*J-2] += kij
                # y-y
                kα[3*I-1, 3*J-1] += kij
                # kα[3*I, 3*J] += kij
                # (no coupling terms, no z terms)
            end
        end
    end
end




function Δ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean2(ap::T, k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
   
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        E=ξ.Ē
        ν=ξ.ν̄ 
        K=E/(1-2*ν)/3
        λ=E*ν/(1+ν)/(1-2*ν)
        μ=E/(1+ν)/2
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        J⁻²=(1.0/J)^2
        # println(J⁻²)
        J⁻²³=cbrt(J⁻²)
        # println(J⁻²³)
        C₁₁₁₁=(-2/3*C⁻¹₁₁-2/3*C⁻¹₁₁)*J⁻²³*μ+(J^2*K+2/9*J⁻²³*I₁*μ)*C⁻¹₁₁*C⁻¹₁₁+(1/3*μ*J⁻²³*I₁-0.5*K*(J^2-1.0))*2*C⁻¹₁₁*C⁻¹₁₁
            #   println(C₁₁₁₁)
        C₂₂₂₂=(-2/3*C⁻¹₂₂-2/3*C⁻¹₂₂)*J⁻²³*μ+(J^2*K+2/9*J⁻²³*I₁*μ)*C⁻¹₂₂*C⁻¹₂₂+(1/3*μ*J⁻²³*I₁-0.5*K*(J^2-1.0))*2*C⁻¹₂₂*C⁻¹₂₂
        C₁₁₂₂=(-2/3*C⁻¹₂₂-2/3*C⁻¹₁₁)*J⁻²³*μ+(J^2*K+2/9*J⁻²³*I₁*μ)*C⁻¹₁₁*C⁻¹₂₂+(1/3*μ*J⁻²³*I₁-0.5*K*(J^2-1.0))*2*C⁻¹₁₂*C⁻¹₁₂
        C₁₁₁₂=(-2/3*C⁻¹₁₂)*J⁻²³*μ+(J^2*K+2/9*J⁻²³*I₁*μ)*C⁻¹₁₁*C⁻¹₁₂+(1/3*μ*J⁻²³*I₁-0.5*K*(J^2-1.0))*2*C⁻¹₁₁*C⁻¹₁₂ 
        C₂₂₁₂=(-2/3*C⁻¹₁₂)*J⁻²³*μ+(J^2*K+2/9*J⁻²³*I₁*μ)*C⁻¹₂₂*C⁻¹₁₂+(1/3*μ*J⁻²³*I₁-0.5*K*(J^2-1.0))*2*C⁻¹₂₂*C⁻¹₁₂     
        C₁₂₁₂=(J^2*K+2/9*J⁻²³*I₁*μ)*C⁻¹₁₂*C⁻¹₁₂+(1/3*μ*J⁻²³*I₁-0.5*K*(J^2-1.0))*(C⁻¹₁₁*C⁻¹₂₂+C⁻¹₁₂*C⁻¹₁₂)   
              
        S₁₁ = μ*J⁻²³*(1.0-1/3*I₁*C⁻¹₁₁)+0.5*K*(J^2-1.0)*C⁻¹₁₁
        S₂₂ = μ*J⁻²³*(1.0-1/3*I₁*C⁻¹₂₂)+0.5*K*(J^2-1.0)*C⁻¹₂₂
        S₁₂ = μ*J⁻²³*(-1/3*I₁*C⁻¹₁₂)+0.5*K*(J^2-1.0)*C⁻¹₁₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] +=(B₁[i]*F₁₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               + (B₁[i]*F₁₁*B₂[j]*F₁₂ + B₂[i]*F₁₂*B₁[j]*F₁₁)*C₁₁₂₂
                               + (B₁[i]*F₁₁*B₂[j]*F₁₁ + B₁[i]*F₁₁*B₁[j]*F₁₂)*C₁₁₁₂
                               + (B₁[i]*F₁₂*B₁[j]*F₁₁ + B₂[i]*F₁₁*B₁[j]*F₁₁)*C₁₁₁₂
                               + (B₂[i]*F₁₂*B₂[j]*F₁₁ + B₂[i]*F₁₂*B₁[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂*B₂[j]*F₁₂ + B₂[i]*F₁₁*B₂[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂
                               +  B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤
                              
                k[2*I-1,2*J]   += (B₁[i]*F₁₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₂₂*C₂₂₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₂ + B₂[i]*F₁₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₁ + B₁[i]*F₁₁*B₁[j]*F₂₂)*C₁₁₁₂ 
                               +  (B₁[i]*F₁₂*B₁[j]*F₂₁ + B₂[i]*F₁₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₁₂*B₂[j]*F₂₁ + B₂[i]*F₁₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂*B₂[j]*F₂₂ + B₂[i]*F₁₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂)*𝑤
                               
                k[2*I,2*J-1]   += (B₁[i]*F₂₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₂ + B₂[i]*F₂₂*B₁[j]*F₁₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₁ + B₁[i]*F₂₁*B₁[j]*F₁₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₁₁ + B₂[i]*F₂₁*B₁[j]*F₁₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₁₁ + B₂[i]*F₂₂*B₁[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₁₂ + B₂[i]*F₂₁*B₂[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂)*𝑤

                k[2*I,2*J]     += (B₁[i]*F₂₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₂₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₂ + B₂[i]*F₂₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₁ + B₁[i]*F₂₁*B₁[j]*F₂₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₂₁ + B₂[i]*F₂₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₂₁ + B₂[i]*F₂₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₂₂ + B₂[i]*F₂₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂
                               +   B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤

            end
        end
    end
end

function Δ∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean2(ap::T, k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        E=ξ.Ē
        ν=ξ.ν̄ 
        K=E/(1-2*ν)/3
        λ=E*ν/(1+ν)/(1-2*ν)
        μ=E/(1+ν)/2
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        J⁻²=(1.0/J)^2
        J⁻²³=cbrt(J⁻²)

        C₁₁₁₁=J^2*K*C⁻¹₁₁*C⁻¹₁₁+0.5*K*(J^2-1.0)*2*C⁻¹₁₁*C⁻¹₁₁
        C₂₂₂₂=J^2*K*C⁻¹₂₂*C⁻¹₂₂+0.5*K*(J^2-1.0)*2*C⁻¹₂₂*C⁻¹₂₂
        C₁₁₂₂=J^2*K*C⁻¹₁₁*C⁻¹₂₂+0.5*K*(J^2-1.0)*2*C⁻¹₁₂*C⁻¹₁₂
        C₁₁₁₂=J^2*K*C⁻¹₁₁*C⁻¹₁₂+0.5*K*(J^2-1.0)*2*C⁻¹₁₁*C⁻¹₁₂ 
        C₂₂₁₂=J^2*K*C⁻¹₂₂*C⁻¹₁₂+0.5*K*(J^2-1.0)*2*C⁻¹₂₂*C⁻¹₁₂     
        C₁₂₁₂=J^2*K*C⁻¹₁₂*C⁻¹₁₂+0.5*K*(J^2-1.0)*(C⁻¹₁₁*C⁻¹₂₂+C⁻¹₁₂*C⁻¹₁₂)   
              
        S₁₁ = 0.5*K*(J^2-1.0)*C⁻¹₁₁
        S₂₂ = 0.5*K*(J^2-1.0)*C⁻¹₂₂
        S₁₂ = 0.5*K*(J^2-1.0)*C⁻¹₁₂
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] +=(B₁[i]*F₁₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               + (B₁[i]*F₁₁*B₂[j]*F₁₂ + B₂[i]*F₁₂*B₁[j]*F₁₁)*C₁₁₂₂
                               + (B₁[i]*F₁₁*B₂[j]*F₁₁ + B₁[i]*F₁₁*B₁[j]*F₁₂)*C₁₁₁₂
                               + (B₁[i]*F₁₂*B₁[j]*F₁₁ + B₂[i]*F₁₁*B₁[j]*F₁₁)*C₁₁₁₂
                               + (B₂[i]*F₁₂*B₂[j]*F₁₁ + B₂[i]*F₁₂*B₁[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂*B₂[j]*F₁₂ + B₂[i]*F₁₁*B₂[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂
                               +  B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤
                              
                k[2*I-1,2*J]   += (B₁[i]*F₁₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₂₂*C₂₂₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₂ + B₂[i]*F₁₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₁ + B₁[i]*F₁₁*B₁[j]*F₂₂)*C₁₁₁₂ 
                               +  (B₁[i]*F₁₂*B₁[j]*F₂₁ + B₂[i]*F₁₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₁₂*B₂[j]*F₂₁ + B₂[i]*F₁₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂*B₂[j]*F₂₂ + B₂[i]*F₁₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂)*𝑤
                               
                k[2*I,2*J-1]   += (B₁[i]*F₂₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₂ + B₂[i]*F₂₂*B₁[j]*F₁₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₁ + B₁[i]*F₂₁*B₁[j]*F₁₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₁₁ + B₂[i]*F₂₁*B₁[j]*F₁₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₁₁ + B₂[i]*F₂₂*B₁[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₁₂ + B₂[i]*F₂₁*B₂[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂)*𝑤

                k[2*I,2*J]     += (B₁[i]*F₂₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₂₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₂ + B₂[i]*F₂₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₁ + B₁[i]*F₂₁*B₁[j]*F₂₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₂₁ + B₂[i]*F₂₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₂₁ + B₂[i]*F₂₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₂₂ + B₂[i]*F₂₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂
                               +   B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤

            end
        end
    end
end
function Δ∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean2(ap::T, k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
   
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        E=ξ.Ē
        ν=ξ.ν̄ 
        K=E/(1-2*ν)/3
        λ=E*ν/(1+ν)/(1-2*ν)
        μ=E/(1+ν)/2
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        J⁻²=(1.0/J)^2
        J⁻²³=cbrt(J⁻²)
        C₁₁₁₁=(-2/3*C⁻¹₁₁-2/3*C⁻¹₁₁)*J⁻²³*μ + (2/9*J⁻²³*I₁*μ)*C⁻¹₁₁*C⁻¹₁₁ + (1/3*μ*J⁻²³*I₁)*2*C⁻¹₁₁*C⁻¹₁₁
        C₂₂₂₂=(-2/3*C⁻¹₂₂-2/3*C⁻¹₂₂)*J⁻²³*μ + (2/9*J⁻²³*I₁*μ)*C⁻¹₂₂*C⁻¹₂₂ + (1/3*μ*J⁻²³*I₁)*2*C⁻¹₂₂*C⁻¹₂₂
        C₁₁₂₂=(-2/3*C⁻¹₂₂-2/3*C⁻¹₁₁)*J⁻²³*μ + (2/9*J⁻²³*I₁*μ)*C⁻¹₁₁*C⁻¹₂₂ + (1/3*μ*J⁻²³*I₁)*2*C⁻¹₁₂*C⁻¹₁₂
        C₁₁₁₂=(-2/3*C⁻¹₁₂)*J⁻²³*μ + (2/9*J⁻²³*I₁*μ)*C⁻¹₁₁*C⁻¹₁₂ + (1/3*μ*J⁻²³*I₁)*2*C⁻¹₁₁*C⁻¹₁₂ 
        C₂₂₁₂=(-2/3*C⁻¹₁₂)*J⁻²³*μ + (2/9*J⁻²³*I₁*μ)*C⁻¹₂₂*C⁻¹₁₂ + (1/3*μ*J⁻²³*I₁)*2*C⁻¹₂₂*C⁻¹₁₂     
        C₁₂₁₂=(2/9*J⁻²³*I₁*μ)*C⁻¹₁₂*C⁻¹₁₂ + (1/3*μ*J⁻²³*I₁)*(C⁻¹₁₁*C⁻¹₂₂+C⁻¹₁₂*C⁻¹₁₂)   
              
        S₁₁ = μ*J⁻²³*(1.0-1/3*I₁*C⁻¹₁₁)
        S₂₂ = μ*J⁻²³*(1.0-1/3*I₁*C⁻¹₂₂)
        S₁₂ = μ*J⁻²³*(-1/3*I₁*C⁻¹₁₂)
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] +=(B₁[i]*F₁₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               + (B₁[i]*F₁₁*B₂[j]*F₁₂ + B₂[i]*F₁₂*B₁[j]*F₁₁)*C₁₁₂₂
                               + (B₁[i]*F₁₁*B₂[j]*F₁₁ + B₁[i]*F₁₁*B₁[j]*F₁₂)*C₁₁₁₂
                               + (B₁[i]*F₁₂*B₁[j]*F₁₁ + B₂[i]*F₁₁*B₁[j]*F₁₁)*C₁₁₁₂
                               + (B₂[i]*F₁₂*B₂[j]*F₁₁ + B₂[i]*F₁₂*B₁[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂*B₂[j]*F₁₂ + B₂[i]*F₁₁*B₂[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂
                               +  B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤
                              
                k[2*I-1,2*J]   += (B₁[i]*F₁₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₂₂*C₂₂₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₂ + B₂[i]*F₁₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₁ + B₁[i]*F₁₁*B₁[j]*F₂₂)*C₁₁₁₂ 
                               +  (B₁[i]*F₁₂*B₁[j]*F₂₁ + B₂[i]*F₁₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₁₂*B₂[j]*F₂₁ + B₂[i]*F₁₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂*B₂[j]*F₂₂ + B₂[i]*F₁₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂)*𝑤
                               
                k[2*I,2*J-1]   += (B₁[i]*F₂₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₂ + B₂[i]*F₂₂*B₁[j]*F₁₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₁ + B₁[i]*F₂₁*B₁[j]*F₁₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₁₁ + B₂[i]*F₂₁*B₁[j]*F₁₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₁₁ + B₂[i]*F₂₂*B₁[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₁₂ + B₂[i]*F₂₁*B₂[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂)*𝑤

                k[2*I,2*J]     += (B₁[i]*F₂₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₂₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₂ + B₂[i]*F₂₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₁ + B₁[i]*F₂₁*B₁[j]*F₂₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₂₁ + B₂[i]*F₂₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₂₁ + B₂[i]*F₂₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₂₂ + B₂[i]*F₂₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂
                               +   B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤

            end
        end
    end
end

function ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean2(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
  
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        E=ξ.Ē
        ν=ξ.ν̄ 
        K=E/(1-2*ν)/3
        λ=E*ν/(1+ν)/(1-2*ν)
        μ=E/(1+ν)/2
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)
        
        J⁻²=(1.0/J)^2
        J⁻²³=cbrt(J⁻²)
        S₁₁ = μ*J⁻²³*(1.0-1/3*I₁*C⁻¹₁₁) + 0.5*K*(J^2-1.0)*C⁻¹₁₁
        S₂₂ = μ*J⁻²³*(1.0-1/3*I₁*C⁻¹₂₂) + 0.5*K*(J^2-1.0)*C⁻¹₂₂
        S₁₂ = μ*J⁻²³*(-1/3*I₁*C⁻¹₁₂) + 0.5*K*(J^2-1.0)*C⁻¹₁₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += (B₁[i]*F₁₁*S₁₁+B₂[i]*F₁₂*S₂₂+(B₁[i]*F₁₂+B₂[i]*F₁₁)*S₁₂)*𝑤
            f[2*I]   += (B₁[i]*F₂₁*S₁₁+B₂[i]*F₂₂*S₂₂+(B₁[i]*F₂₂+B₂[i]*F₂₁)*S₁₂)*𝑤
        end
    end
end
function ∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean2(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
  
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        E=ξ.Ē
        ν=ξ.ν̄ 
        K=E/(1-2*ν)/3
        λ=E*ν/(1+ν)/(1-2*ν)
        μ=E/(1+ν)/2
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        
        S₁₁ = 0.5*K*(J^2-1.0)*C⁻¹₁₁
        S₂₂ = 0.5*K*(J^2-1.0)*C⁻¹₂₂
        S₁₂ = 0.5*K*(J^2-1.0)*C⁻¹₁₂
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += (B₁[i]*F₁₁*S₁₁+B₂[i]*F₁₂*S₂₂+(B₁[i]*F₁₂+B₂[i]*F₁₁)*S₁₂)*𝑤
            f[2*I]   += (B₁[i]*F₂₁*S₁₁+B₂[i]*F₂₂*S₂₂+(B₁[i]*F₂₂+B₂[i]*F₂₁)*S₁₂)*𝑤
        end
    end
end
function ∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean2(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
   
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        E=ξ.Ē
        ν=ξ.ν̄ 
        K=E/(1-2*ν)/3
        λ=E*ν/(1+ν)/(1-2*ν)
        μ=E/(1+ν)/2
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        J⁻²=(1.0/J)^2
        J⁻²³=cbrt(J⁻²)
        S₁₁ = μ*J⁻²³*(1.0-1/3*I₁*C⁻¹₁₁)
        S₂₂ = μ*J⁻²³*(1.0-1/3*I₁*C⁻¹₂₂)
        S₁₂ = μ*J⁻²³*(-1/3*I₁*C⁻¹₁₂)
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += (B₁[i]*F₁₁*S₁₁+B₂[i]*F₁₂*S₂₂+(B₁[i]*F₁₂+B₂[i]*F₁₁)*S₁₂)*𝑤
            f[2*I]   += (B₁[i]*F₂₁*S₁₁+B₂[i]*F₂₂*S₂₂+(B₁[i]*F₂₂+B₂[i]*F₂₁)*S₁₂)*𝑤
        end
    end
end

function get_neoHookean_S_and_H_components(E_vec::Vector{Float64},
                                           E::Float64, ν::Float64)
    E11, E22, E12 = E_vec

    # Lamé 参数
    λ = E*ν/((1+ν)*(1-2ν))
    μ = E/(2*(1+ν))

    # 由 E 得到 C
    C11 = 2*E11 + 1.0
    C22 = 2*E22 + 1.0
    C12 = 2*E12
    C33 = 1.0

    detC = C11*C22 - C12*C12
    if detC <= 0
        detC = 1e-12       # 防守一下
    end
    J = sqrt(detC)

    I1 = C11 + C22 + C33
    I2 = 0.5*(I1^2 - C11^2 - C22^2 - C33^2 - 2*C12^2)

    # C^{-1}
    Cinv11 = 1.0/detC*(C11*C11 + C12*C12 - I1*C11 + I2)
    Cinv12 = 1.0/detC*(C11*C12 + C12*C22 - I1*C12)
    Cinv22 = 1.0/detC*(C12*C12 + C22*C22 - I1*C22 + I2)

    # 应力分量 S_ij(E)
    S11 = λ*J*(J-1.0)*Cinv11 + μ*(1.0 - Cinv11)
    S22 = λ*J*(J-1.0)*Cinv22 + μ*(1.0 - Cinv22)
    S12 = λ*J*(J-1.0)*Cinv12 - μ*Cinv12

    # 四阶切线的分量（对应 H 的 6 个独立值）
    C1111 = λ*J*(2*J-1.0)*Cinv11*Cinv11 +
            (μ-λ*J*(J-1.0))*2*Cinv11*Cinv11

    C2222 = λ*J*(2*J-1.0)*Cinv22*Cinv22 +
            (μ-λ*J*(J-1.0))*2*Cinv22*Cinv22

    C1122 = λ*J*(2*J-1.0)*Cinv11*Cinv22 +
            (μ-λ*J*(J-1.0))*2*Cinv12*Cinv12

    C1112 = λ*J*(2*J-1.0)*Cinv11*Cinv12 +
            (μ-λ*J*(J-1.0))*2*Cinv11*Cinv12

    C2212 = λ*J*(2*J-1.0)*Cinv22*Cinv12 +
            (μ-λ*J*(J-1.0))*2*Cinv12*Cinv22

    C1212 = λ*J*(2*J-1.0)*Cinv12*Cinv12 +
            (μ-λ*J*(J-1.0))*(Cinv11*Cinv22 + Cinv12*Cinv12)

    # 返回应力张量分量 + 切线的 6 个独立分量
    S_vec = [S11, S22, S12]
    H_comp = (C1111, C2222, C1122, C1112, C2212, C1212)

    return S_vec, H_comp
end





# function get_neoHookean_S_and_H(E_vec::Vector{Float64}, E::Float64, ν::Float64)

#     E11, E22, E12 = E_vec

#     # Lamé 参数
#     λ = E*ν/((1+ν)*(1-2ν))
#     μ = E/(2*(1+ν))

#     # 由 E 得到 C
#     C11 = 2*E11 + 1.0
#     C22 = 2*E22 + 1.0
#     C12 = 2*E12
#     C33 = 1.0

#     detC = C11*C22 - C12*C12
#     J    = sqrt(detC)

#     I1 = C11 + C22 + C33
#     I2 = 0.5*(I1^2 - C11^2 - C22^2 - C33^2 - 2*C12^2)

#     # C^{-1}
#     Cinv11 = 1.0/detC*(C11*C11 + C12*C12 - I1*C11 + I2)
#     Cinv12 = 1.0/detC*(C11*C12 + C12*C22 - I1*C12)
#     Cinv22 = 1.0/detC*(C12*C12 + C22*C22 - I1*C22 + I2)

#     # 应力 S
#     S11 = λ*J*(J-1.0)*Cinv11 + μ*(1.0 - Cinv11)
#     S22 = λ*J*(J-1.0)*Cinv22 + μ*(1.0 - Cinv22)
#     S12 = λ*J*(J-1.0)*Cinv12 - μ*Cinv12

#     S_vec = [S11, S22, S12]

#     # 四阶切线分量（完全沿用你现在的表达式）
#     C1111 = λ*J*(2*J-1.0)*Cinv11*Cinv11 + (μ-λ*J*(J-1.0))*2*Cinv11*Cinv11
#     C2222 = λ*J*(2*J-1.0)*Cinv22*Cinv22 + (μ-λ*J*(J-1.0))*2*Cinv22*Cinv22
#     C1122 = λ*J*(2*J-1.0)*Cinv11*Cinv22 + (μ-λ*J*(J-1.0))*2*Cinv12*Cinv12
#     C1112 = λ*J*(2*J-1.0)*Cinv11*Cinv12 + (μ-λ*J*(J-1.0))*2*Cinv11*Cinv12
#     C2212 = λ*J*(2*J-1.0)*Cinv22*Cinv12 + (μ-λ*J*(J-1.0))*2*Cinv12*Cinv22
#     C1212 = λ*J*(2*J-1.0)*Cinv12*Cinv12 + (μ-λ*J*(J-1.0))*(Cinv11*Cinv22 + Cinv12*Cinv12)

#     H = zeros(3,3)
#     H[1,1] = C1111
#     H[1,2] = C1122
#     H[1,3] = C1112

#     H[2,1] = C1122
#     H[2,2] = C2222
#     H[2,3] = C2212

#     H[3,1] = C1112
#     H[3,2] = C2212
#     H[3,3] = C1212

#     return S_vec, H
# end




function get_neoHookean_S_and_H(E_vec::Vector{Float64}, E::Float64, ν::Float64)

    E11, E22, E12 = E_vec

    # Lamé 参数
    λ = E * ν / ((1 + ν) * (1 - 2*ν))
    μ = E / (2 * (1 + ν))

    # 由 E 得到 C (plane strain, C33 = 1)
    C11 = 2.0*E11 + 1.0
    C22 = 2.0*E22 + 1.0
    C12 = 2.0*E12
    C33 = 1.0

    C = @inbounds [
        C11  C12  0.0
        C12  C22  0.0
        0.0  0.0  C33
    ]

    detC = det(C)
    if detC <= 1e-14
        throw(DomainError(detC, "Non-physical state: det(C) <= 0"))
    end
    J = sqrt(detC)

    Cinv = inv(C)
    Cinv11 = Cinv[1,1]
    Cinv22 = Cinv[2,2]
    Cinv12 = Cinv[1,2]  # = Cinv[2,1]

    # ===========================
    # 第二 Piola 应力 S
    # (保持你原来的形式: λ*J*(J-1)*Cinv + μ*(I - Cinv))
    # ===========================
    vol = λ * J * (J - 1.0)

    S11 = vol*Cinv11 + μ*(1.0 - Cinv11)
    S22 = vol*Cinv22 + μ*(1.0 - Cinv22)
    S12 = vol*Cinv12 - μ*Cinv12

    S_vec = [S11, S22, S12]

    # ===========================
    # 一致切线 H (Voigt: [11,22,12])
    # 你原先的写法可以整理为:
    #   A = λ*J*(2J-1)
    #   B = μ - λ*J*(J-1)
    # 然后套用 neo-Hookean 的标准形式:
    #   C_{ijkl} = A*Cinv_ij*Cinv_kl + B*(Cinv_ik*Cinv_jl + Cinv_il*Cinv_jk)
    #
    # 这与您原先逐项写出来的表达式一致，但更清晰且避免抄写错误。
    # ===========================
    A = λ * J * (2.0*J - 1.0)
    B = μ - λ * J * (J - 1.0)

    C1111 = A*Cinv11*Cinv11 + 2.0*B*Cinv11*Cinv11
    C2222 = A*Cinv22*Cinv22 + 2.0*B*Cinv22*Cinv22
    C1122 = A*Cinv11*Cinv22 + 2.0*B*Cinv12*Cinv12
    C1112 = A*Cinv11*Cinv12 + 2.0*B*Cinv11*Cinv12
    C2212 = A*Cinv22*Cinv12 + 2.0*B*Cinv12*Cinv22
    C1212 = A*Cinv12*Cinv12 + B*(Cinv11*Cinv22 + Cinv12*Cinv12)

    H = zeros(3,3)
    H[1,1] = C1111
    H[1,2] = C1122
    H[1,3] = C1112

    H[2,1] = C1122
    H[2,2] = C2222
    H[2,3] = C2212

    H[3,1] = C1112
    H[3,2] = C2212
    H[3,3] = C1212

    return S_vec, H
end




# function get_neoHookean2_S_and_H(E_vec::Vector{Float64}, E::Float64, ν::Float64)

#     E11, E22, E12 = E_vec

#     # Lamé 参数
#     λ = E*ν/((1+ν)*(1-2ν))
#     μ = E/(2*(1+ν))

#     # 由 E 得到 C
#     C11 = 2*E11 + 1.0
#     C22 = 2*E22 + 1.0
#     C12 = 2*E12
#     C33 = 1.0

#     detC = C11*C22 - C12*C12
#     J    = sqrt(detC)

#     I1 = C11 + C22 + C33
#     I2 = 0.5*(I1^2 - C11^2 - C22^2 - C33^2 - 2*C12^2)

#     # C^{-1}
#     Cinv11 = 1.0/detC*(C11*C11 + C12*C12 - I1*C11 + I2)
#     Cinv12 = 1.0/detC*(C11*C12 + C12*C22 - I1*C12)
#     Cinv22 = 1.0/detC*(C12*C12 + C22*C22 - I1*C22 + I2)

#     # 应力 S
#     S11 = λ*J*(J-1.0)*Cinv11 + μ*(1.0 - Cinv11)
#     S22 = λ*J*(J-1.0)*Cinv22 + μ*(1.0 - Cinv22)
#     S12 = λ*J*(J-1.0)*Cinv12 - μ*Cinv12

#     S_vec = [S11, S22, S12]

#     # 四阶切线分量（完全沿用你现在的表达式）
#     C1111 = λ*J*(2*J-1.0)*Cinv11*Cinv11 + (μ-λ*J*(J-1.0))*2*Cinv11*Cinv11
#     C2222 = λ*J*(2*J-1.0)*Cinv22*Cinv22 + (μ-λ*J*(J-1.0))*2*Cinv22*Cinv22
#     C1122 = λ*J*(2*J-1.0)*Cinv11*Cinv22 + (μ-λ*J*(J-1.0))*2*Cinv12*Cinv12
#     C1112 = λ*J*(2*J-1.0)*Cinv11*Cinv12 + (μ-λ*J*(J-1.0))*2*Cinv11*Cinv12
#     C2212 = λ*J*(2*J-1.0)*Cinv22*Cinv12 + (μ-λ*J*(J-1.0))*2*Cinv12*Cinv22
#     C1212 = λ*J*(2*J-1.0)*Cinv12*Cinv12 + (μ-λ*J*(J-1.0))*(Cinv11*Cinv22 + Cinv12*Cinv12)

#     H = zeros(3,3)
#     H[1,1] = C1111
#     H[1,2] = C1122
#     H[1,3] = C1112

#     H[2,1] = C1122
#     H[2,2] = C2222
#     H[2,3] = C2212

#     H[3,1] = C1112
#     H[3,2] = C2212
#     H[3,3] = C1212

#     return S_vec, H
# end

# # 原始可用版本
# function get_neoHookean_lnJ_S_and_H(E_vec::Vector{Float64}, E::Float64, ν::Float64;
#                                     detC_min::Float64=1e-14,
#                                     fd_relstep::Float64=1e-7,
#                                     fd_absstep::Float64=1e-10)

#     E11, E22, E12 = E_vec

#     μ = E / (2.0 * (1.0 + ν))
#     κ = E / (3.0 * (1.0 - 2.0*ν))   # bulk modulus

#     # build C (plane strain: C33=1)
#     C11 = 2.0*E11 + 1.0
#     C22 = 2.0*E22 + 1.0
#     C12 = 2.0*E12
#     C33 = 1.0

#     C = @inbounds [
#         C11  C12  0.0
#         C12  C22  0.0
#         0.0  0.0  C33
#     ]

#     detC = det(C)
#     if !(isfinite(detC) && detC > detC_min)
#         throw(DomainError(detC, "Non-physical state: det(C) <= detC_min"))
#     end

#     J = sqrt(detC)
#     lnJ = log(J)

#     Cinv = inv(C)
#     I1 = tr(C)

#     # --------
#     # stress S = S_iso + S_vol
#     # --------
#     Jm23 = J^(-2.0/3.0)

#     # S_iso = μ J^{-2/3} ( I - 1/3 I1 C^{-1} )
#     S_iso = μ * Jm23 .* (I(3) .- (I1/3.0).*Cinv)

#     # S_vol = κ lnJ C^{-1}
#     S_vol = κ * lnJ .* Cinv

#     S = S_iso .+ S_vol

#     # Voigt output [11,22,12]
#     S_vec = [S[1,1], S[2,2], S[1,2]]

#     # --------
#     # consistent tangent via finite difference on E_vec
#     # H_{ij} = ∂S_i / ∂E_j in Voigt([11,22,12]) space
#     # --------
#     function S_voigt(Ev::Vector{Float64})
#         E11t, E22t, E12t = Ev
#         C11t = 2.0*E11t + 1.0
#         C22t = 2.0*E22t + 1.0
#         C12t = 2.0*E12t
#         Ct = @inbounds [
#             C11t  C12t  0.0
#             C12t  C22t  0.0
#             0.0   0.0   1.0
#         ]
#         detCt = det(Ct)
#         if !(isfinite(detCt) && detCt > detC_min)
#             throw(DomainError(detCt, "Non-physical in fd"))
#         end
#         Jt = sqrt(detCt)
#         lnJt = log(Jt)
#         Cinv_t = inv(Ct)
#         I1t = tr(Ct)
#         Jm23t = Jt^(-2.0/3.0)
#         S_iso_t = μ * Jm23t .* (I(3) .- (I1t/3.0).*Cinv_t)
#         S_vol_t = κ * lnJt .* Cinv_t
#         St = S_iso_t .+ S_vol_t
#         return [St[1,1], St[2,2], St[1,2]]
#     end

#     H = zeros(3,3)
#     for j in 1:3
#         # step size: relative + absolute guard
#         hj = max(fd_absstep, fd_relstep*max(1.0, abs(E_vec[j])))

#         Ep = copy(E_vec); Ep[j] += hj
#         Em = copy(E_vec); Em[j] -= hj

#         Sp = S_voigt(Ep)
#         Sm = S_voigt(Em)

#         @inbounds H[:,j] .= (Sp .- Sm) ./ (2.0*hj)
#     end

#     return S_vec, H
# end



# function get_neoHookean3_S_and_H(E_vec::Vector{Float64}, E::Float64, ν::Float64)

#     E11, E22, E12 = E_vec

#     λ = E*ν/((1+ν)*(1-2ν))
#     μ = E/(2*(1+ν))

#     # C tensor (plane strain: C₃₃ = 1)
#     C₁₁ = 2*E11 + 1
#     C₂₂ = 2*E22 + 1
#     C₁₂ = 2*E12
#     C₃₃ = 1.0

#     detC = C₁₁*C₂₂ - C₁₂*C₁₂
#     J = sqrt(detC)

#     # C inverse (2×2 子块，C33=1 ⇒ Cinv33=1)
#     Cinv11 =  C₂₂/detC
#     Cinv22 =  C₁₁/detC
#     Cinv12 = -C₁₂/detC   # = Cinv21
#     I₁ = C₁₁+C₂₂+C₃₃
#     I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
#     Cinv11=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
#     Cinv12=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
#     Cinv22=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)
#     # ===========================
#     # 第二 Piola 应力 S
#     # S = μ*(I - C^{-1}) + (λ/2)*(J^2 - 1)*C^{-1}
#     # ===========================
#     factor = 0.5 * λ * (J^2 - 1.0)

#     S11 = μ*(1.0 - Cinv11) + factor*Cinv11
#     S22 = μ*(1.0 - Cinv22) + factor*Cinv22
#     S12 =      - μ*Cinv12 + factor*Cinv12

#     S_vec = [S11, S22, S12]

#     # ===========================
#     # 一致切线 H (Voigt: [11,22,12])
#     #
#     # C_{ijkl} = λ J^2 Cinv_ij Cinv_kl
#     #          + B (Cinv_ik Cinv_jl + Cinv_il Cinv_jk)
#     # 其中 B = μ - 0.5*λ*(J^2 - 1)
#     # ===========================
#     A = λ * J^2
#     B = μ - 0.5*λ*(J^2 - 1.0)

#     C1111 = A*Cinv11*Cinv11 + 2.0*B*Cinv11*Cinv11
#     C2222 = A*Cinv22*Cinv22 + 2.0*B*Cinv22*Cinv22
#     C1122 = A*Cinv11*Cinv22 + 2.0*B*Cinv12*Cinv12
#     C1112 = A*Cinv11*Cinv12 + 2.0*B*Cinv11*Cinv12
#     C2212 = A*Cinv22*Cinv12 + 2.0*B*Cinv12*Cinv22
#     C1212 = A*Cinv12*Cinv12 + B*(Cinv11*Cinv22 + Cinv12*Cinv12)

#     H = zeros(3,3)
#     H[1,1] = C1111
#     H[1,2] = C1122
#     H[1,3] = C1112

#     H[2,1] = C1122
#     H[2,2] = C2222
#     H[2,3] = C2212

#     H[3,1] = C1112
#     H[3,2] = C2212
#     H[3,3] = C1212

#     return S_vec, H
# end




function get_neoHookean2_S_and_H(E_vec::Vector{Float64}, E::Float64, ν::Float64)

    E11, E22, E12 = E_vec

    # Lamé 参数
    λ = E * ν / ((1 + ν) * (1 - 2*ν))
    μ = E / (2 * (1 + ν))

    # 由 E 得到 C (plane strain: C33 = 1)
    C11 = 2.0*E11 + 1.0
    C22 = 2.0*E22 + 1.0
    C12 = 2.0*E12/2.0
    C33 = 1.0

    C = @inbounds [
        C11  C12  0.0
        C12  C22  0.0
        0.0  0.0  C33
    ]

    detC = det(C)
    if detC <= 1e-14
        throw(DomainError(detC, "Non-physical state: det(C) <= 0"))
    end
    J = sqrt(detC)

    Cinv = inv(C)
    Cinv11 = Cinv[1,1]
    Cinv22 = Cinv[2,2]
    Cinv12 = Cinv[1,2]  # = Cinv[2,1]

    # ===========================
    # 第二 Piola 应力 S (保持你原形式)
    # ===========================
    vol = λ * J * (J - 1.0)

    S11 = vol*Cinv11 + μ*(1.0 - Cinv11)
    S22 = vol*Cinv22 + μ*(1.0 - Cinv22)
    S12 = vol*Cinv12 - μ*Cinv12

    S_vec = [S11, S22, S12]

    # ===========================
    # 一致切线 H (Voigt: [11,22,12]) - 保持你原表达式结构
    #
    # 你原代码中的两组系数可写成：
    #   A = λ*J*(2J-1)
    #   B = μ - λ*J*(J-1)
    #
    # 然后
    #   C_{ijkl} = A*Cinv_ij*Cinv_kl + B*(Cinv_ik*Cinv_jl + Cinv_il*Cinv_jk)
    # ===========================
    A = λ * J * (2.0*J - 1.0)
    B = μ - λ * J * (J - 1.0)

    C1111 = A*Cinv11*Cinv11 + 2.0*B*Cinv11*Cinv11
    C2222 = A*Cinv22*Cinv22 + 2.0*B*Cinv22*Cinv22
    C1122 = A*Cinv11*Cinv22 + 2.0*B*Cinv12*Cinv12
    C1112 = A*Cinv11*Cinv12 + 2.0*B*Cinv11*Cinv12
    C2212 = A*Cinv22*Cinv12 + 2.0*B*Cinv12*Cinv22
    C1212 = A*Cinv12*Cinv12 + B*(Cinv11*Cinv22 + Cinv12*Cinv12)

    H = zeros(3,3)
    H[1,1] = C1111
    H[1,2] = C1122
    H[1,3] = C1112

    H[2,1] = C1122
    H[2,2] = C2222
    H[2,3] = C2212

    H[3,1] = C1112
    H[3,2] = C2212
    H[3,3] = C1212

    return S_vec, H
end




function get_neoHookean3_S_and_H(E_vec::Vector{Float64}, E::Float64, ν::Float64)
    E11, E22, E12 = E_vec

    # Lamé parameters
    λ = E * ν / ((1 + ν) * (1 - 2*ν))
    μ = E / (2 * (1 + ν))

    # Right Cauchy-Green tensor C (plane strain with C33=1)
    C11 = 2.0*E11 + 1.0
    C22 = 2.0*E22 + 1.0
    C12 = 2.0*E12
    C33 = 1.0

    C = @inbounds [
        C11  C12  0.0
        C12  C22  0.0
        0.0  0.0  C33
    ]

    detC = det(C)
    if detC <= 1e-14
        throw(DomainError(detC, "Non-physical state: det(C) <= 0"))
    end
    J = sqrt(detC)

    Cinv = inv(C)
    Cinv11 = Cinv[1,1]
    Cinv22 = Cinv[2,2]
    Cinv12 = Cinv[1,2]  # = Cinv[2,1]

    # ===========================
    # Second Piola stress S:
    # S = μ*(I - C^{-1}) + (λ/2)*(J^2 - 1)*C^{-1}
    # ===========================
    factor = 0.5 * λ * (J^2 - 1.0)

    S11 = μ*(1.0 - Cinv11) + factor*Cinv11
    S22 = μ*(1.0 - Cinv22) + factor*Cinv22
    S12 =        - μ*Cinv12 + factor*Cinv12

    S_vec = [S11, S22, S12]

    # ===========================
    # Consistent tangent H in Voigt [11,22,12] with tensor shear component (E12)
    #
    # C_{ijkl} = λ J^2 Cinv_ij Cinv_kl
    #          + B (Cinv_ik Cinv_jl + Cinv_il Cinv_jk)
    # where B = μ - 0.5*λ*(J^2 - 1)
    # ===========================
    A = λ * J^2
    B = μ - 0.5*λ*(J^2 - 1.0)

    # Components needed for [11,22,12]
    C1111 = A*Cinv11*Cinv11 + 2.0*B*Cinv11*Cinv11
    C2222 = A*Cinv22*Cinv22 + 2.0*B*Cinv22*Cinv22
    C1122 = A*Cinv11*Cinv22 + 2.0*B*Cinv12*Cinv12
    C1112 = A*Cinv11*Cinv12 + 2.0*B*Cinv11*Cinv12
    C2212 = A*Cinv22*Cinv12 + 2.0*B*Cinv12*Cinv22
    C1212 = A*Cinv12*Cinv12 + B*(Cinv11*Cinv22 + Cinv12*Cinv12)

    H = zeros(3,3)
    H[1,1] = C1111
    H[1,2] = C1122
    H[1,3] = C1112*2.0

    H[2,1] = C1122
    H[2,2] = C2222
    H[2,3] = C2212*2.0

    H[3,1] = C1112
    H[3,2] = C2212
    H[3,3] = C1212*2.0

    return S_vec, H
end





# function update_Econs!(ξ, S_vec; tol=1e-10, maxiter=20)
#     E_mod = ξ.Ē
#     ν     = ξ.ν̄ 

#     # 线弹性刚度矩阵 C_elastic
#     λ = E_mod*ν / ((1+ν)*(1-2ν))
#     μ = E_mod / (2*(1+ν)) 

#     C_elastic = zeros(3,3)
#     C_elastic[1,1] = λ + 2μ
#     C_elastic[1,2] = λ
#     C_elastic[2,1] = λ
#     C_elastic[2,2] = λ + 2μ
#     C_elastic[3,3] = μ 

#     # 用线弹性反算初值 E_vec
#     E_vec = C_elastic \ S_vec

#     # 初始 D：先用线弹性，收敛后再改成 H^{-1}
#     D = copy(C_elastic)
#     converged = false
    
#     for it = 1:maxiter
#         # 一阶导 ∂ψ/∂E (= 应力) 和二阶导 ∂²ψ/∂E²
#         dpsi_dE, H = get_neoHookean_S_and_H(E_vec, E_mod, ν)

#         # 残差 r = S - ∂ψ/∂E
#         r = S_vec .- dpsi_dE
#         # r = dpsi_dE. - S_vec 
#         resnorm = norm(r)

#         if resnorm < tol
#             converged = true
#             D .= inv(H)  # D = H^{-1}
#             @printf("it = %d, res = %.6e (converged)\n", it, resnorm)
#             break
#         end

       
#         # ΔE = inv(H)*r
#         ΔE = H \ r
#         stepnorm = norm(ΔE)
#         @printf("it = %d, res = %.6e, step = %.6e\n", it, resnorm, stepnorm)

#         E_vec .+= ΔE
#     end

#     if !converged
#         @warn "Material iteration not converged "
       
#     end

#     return E_vec, D
# end




#  12.26.可以正确求解
function update_Econs!(ξ, S_vec , E_init::AbstractVector{<:Real}; tol_rel=1e-4, maxiter=30, alpha_max=1.0)
    E_mod = ξ.Ē
    ν     = ξ.ν̄ 
    
    # 1. 计算 Lamé 参数和线弹性刚度矩阵 C_elastic
    λ = E_mod * ν / ((1 + ν) * (1 - 2 * ν))
    μ = E_mod / (2 * (1 + ν))

    C_elastic = zeros(3, 3)
    C_elastic[1, 1] = λ + 2 * μ
    C_elastic[1, 2] = λ
    C_elastic[2, 1] = λ
    C_elastic[2, 2] = λ + 2 * μ
    C_elastic[3, 3] = μ 

    # 2. **关键修改：使用零应变初始化**
    # E_vec = C_elastic \ S_vec  # 原来的线弹性初始化（已注释或删除）
    E_vec = zeros(3)              # 新的零应变初始化
    E_vec_old = copy(E_vec)

    D = copy(C_elastic)
    converged = false

    # --- 关键修改 1: 计算归一化因子 ---
    S_norm = norm(S_vec)
    tol_abs = 1e-10 

    # 预声明变量
    local dpsi_dE_old = zeros(3)
    local H_old = zeros(3, 3)
    local r_old = zeros(3)
    local resnorm_old = 0.0

    # 第一次计算残差 r_old 和切线模量 H_old
    try
      
        dpsi_dE_old, H_old = get_neoHookean2_S_and_H(E_vec_old, E_mod, ν)
    catch e
        @warn "Initial guess failed in get_neoHookean_S_and_H: $(e)"
        return E_vec, D 
    end
    
    r_old .= S_vec .- dpsi_dE_old
    resnorm_old = norm(r_old)

    # 计算目标收敛阈值
    if S_norm < 1e-12 
        tol_check = tol_abs
    else
        tol_check = max(tol_rel * S_norm, tol_abs)
    end


    for it = 1:maxiter
      
        # 3. 计算牛顿步长 ΔE (使用稳定的 H \ r)
        ΔE = H_old \ r_old
        stepnorm = norm(ΔE)
        
        # 4. 回溯线搜索 (Backtracking Line Search)
        alpha = alpha_max
        max_ls_iter = 10 
        
        local r_trial = zeros(3)
        local H_trial = zeros(3, 3)
        
        for ls_it = 1:max_ls_iter
            E_vec_trial = E_vec_old + alpha * ΔE
            
            # --- 4.1 检查非物理状态 (det(C) > 0) ---
            C11_trial = 2.0 * E_vec_trial[1] + 1.0
            C22_trial = 2.0 * E_vec_trial[2] + 1.0
            C12_trial = 2.0 * E_vec_trial[3]
            detC_trial = C11_trial * C22_trial - C12_trial * C12_trial
            
            is_physical = (detC_trial > 1e-12)

            if !is_physical
                alpha /= 2.0
                if alpha < 1e-5
                    @warn "Line search failed: Reached non-physical state (det(C)<=0) and minimum step size."
                    break
                end
                continue
            end
            
            # --- 4.2 计算试验步的残差 ---
            dpsi_dE_trial, H_trial = get_neoHookean2_S_and_H(E_vec_trial, E_mod, ν) 
            
            r_trial .= S_vec .- dpsi_dE_trial
            resnorm_trial = norm(r_trial)
            
            # --- 4.3 检查残差是否显著减小 ---
            if resnorm_trial < resnorm_old
                E_vec .= E_vec_trial
                r_old .= r_trial
                H_old .= H_trial
                resnorm_old = resnorm_trial
                break
            else
                alpha /= 2.0
            end

            if ls_it == max_ls_iter
                # @warn "Line search failed to find a reducing step. (Final alpha=$(alpha))"
            end
        end # end line search
        
        E_vec_old .= E_vec
        
        # @printf("it = %d, res = %.6e, step = %.6e, alpha = %.3f, Rel. Err: %.2e\n", 
        #          it, resnorm_old, norm(alpha * ΔE), alpha, resnorm_old / S_norm)

        # 检查收敛
        if resnorm_old < tol_check
            converged = true

            D .= inv(H_old)
            # @printf("it = %d, res = %.6e (converged, Target: %.2e)\n", it, resnorm_old, tol_check)
            break
        end

    end # end main iteration loop

    if !converged
        # @warn "Material iteration not converged (Final res = $(resnorm_old), Target: $(tol_check))"
    end

    return E_vec, D
end






function Δ∫∫EᵢⱼSᵢⱼdxdy_HR_NeoHookean(e::Int,ap::T,k::AbstractMatrix{Float64},D_hist) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖

    for (g, ξ) in enumerate(𝓖)
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]

      
        S₁₁ = 0.0; S₂₂ = 0.0; S₁₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            S₁₁ += N[i]*xᵢ.dₛ₁₁
            S₂₂ += N[i]*xᵢ.dₛ₂₂
            S₁₂ += N[i]*xᵢ.dₛ₁₂
        end
        S_vec = [S₁₁, S₂₂, S₁₂]

       
        D = D_hist[e,g]
        D11 = D[1,1]; D12 = D[1,2]; D13 = D[1,3]
        D21 = D[2,1]; D22 = D[2,2]; D23 = D[2,3]
        D31 = D[3,1]; D32 = D[3,2]; D33 = D[3,3]
        


      
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[3*I-2,3*J-2] += N[i]*D11*N[j]*𝑤
                k[3*I-2,3*J-1] += N[i]*D12*N[j]*𝑤
                 k[3*I-2, 3*J] += N[i]*D13*N[j]*𝑤
                k[3*I-1,3*J-2] += N[i]*D21*N[j]*𝑤
                k[3*I-1,3*J-1] += N[i]*D22*N[j]*𝑤
                k[3*I-1,3*J] += N[i]*D23*N[j]*𝑤

                k[3*I  ,3*J-2 ] += N[i]*D31*N[j]*𝑤
                k[3*I  ,3*J-1 ] += N[i]*D32*N[j]*𝑤
                k[3*I  ,3*J  ] += N[i]*D33*N[j]*𝑤
                
               
            end
        end
    end
end






function ∫∫EᵢⱼSᵢⱼdxdy_HR_NeoHookean(e::Int, ap::T,f::AbstractVector{Float64},E_cons_hist, D_hist) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    
    for (g, ξ) in enumerate(𝓖)
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]

       
        S₁₁ = 0.0; S₂₂ = 0.0; S₁₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            S₁₁ += N[i]*xᵢ.dₛ₁₁
            S₂₂ += N[i]*xᵢ.dₛ₂₂
            S₁₂ += N[i]*xᵢ.dₛ₁₂
        end
        S_vec = [S₁₁, S₂₂, S₁₂]

      
        # E_vec, D = update_Econs!(ξ, S_vec)
        E_vec, D = update_Econs!(ξ, S_vec, E_cons_hist[e,g])
        E₁₁, E₂₂, E₁₂ = E_vec
        E_cons_hist[e,g] .= E_vec
        D_hist[e,g]      .= D
       
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[3*I-2] += N[i]*E₁₁*𝑤
            f[3*I-1] += N[i]*E₂₂*𝑤
            # f[3*I]   += 2.0*N[i]*E₁₂*𝑤    
            f[3*I]   += N[i]*E₁₂*𝑤 
            

    end
end 
end













function ∫∫δSᵢⱼEᵢⱼdxdy_HR(aₛ::T,aᵤ::S, f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
   
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
       
        𝑤 = ξₛ.𝑤

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

        # 2. 计算 Green-Lagrange 应变 E = 0.5 * (F'F - I)
        E₁₁ = 0.5 * (F₁₁*F₁₁ + F₂₁*F₂₁ - 1.0)
        E₂₂ = 0.5 * (F₁₂*F₁₂ + F₂₂*F₂₂ - 1.0)
        E₁₂ = 0.5 * (F₁₁*F₁₂ + F₂₁*F₂₂) # 张量剪切应变

        # 3. 组装残差 (注意剪切项的系数 2.0，对应双点积 S:E)
        for (i, xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            # 对应 δS₁₁ 的方程
            f[3*I-2] -= Nₛ[i] * E₁₁ * 𝑤
            # 对应 δS₂₂ 的方程
            f[3*I-1] -= Nₛ[i] * E₂₂ * 𝑤
            # 对应 δS₁₂ 的方程 (能量共轭量是 2*E₁₂)
            f[3*I]   -= 2.0*Nₛ[i]  * E₁₂ * 𝑤
        end
    end
end


function ∫∫δSᵢⱼΔEᵢⱼdxdy_HR(aₛ::T,aᵤ::S,k::AbstractMatrix{Float64}) where  {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
   𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
       
        𝑤 = ξₛ.𝑤

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

        # 2. 计算 Green-Lagrange 应变 E = 0.5 * (F'F - I)
        E₁₁ = 0.5 * (F₁₁*F₁₁ + F₂₁*F₂₁ - 1.0)
        E₂₂ = 0.5 * (F₁₂*F₁₂ + F₂₂*F₂₂ - 1.0)
        E₁₂ = 0.5 * (F₁₁*F₁₂ + F₂₁*F₂₂) # 张量剪切应变

       
        for (i, xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[3*I-2, 2*J-1] -= Nₛ[i] * F₁₁ * B₁[j] * 𝑤  # Col: u
                k[3*I-2, 2*J]   -= Nₛ[i] * F₂₁ * B₁[j] * 𝑤  # Col: v

                k[3*I-1, 2*J-1] -= Nₛ[i] * F₁₂ * B₂[j] * 𝑤
                k[3*I-1, 2*J]   -= Nₛ[i] * F₂₂ * B₂[j] * 𝑤

                k[3*I, 2*J-1] -= Nₛ[i] * (F₁₂* B₁[j]+F₁₁ * B₂[j])* 𝑤
                k[3*I, 2*J]   -= Nₛ[i] * (F₂₂* B₁[j]+F₂₁ * B₂[j]) * 𝑤



            end
        end
    end
end



function ∫∫SᵢⱼδEᵢⱼdxdy_HR(aₛ::T,aᵤ::S, f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
   
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
      
        𝑤 = ξₛ.𝑤

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

        # 2. 计算 Green-Lagrange 应变 E = 0.5 * (F'F - I)
        E₁₁ = 0.5 * (F₁₁*F₁₁ + F₂₁*F₂₁ - 1.0)
        E₂₂ = 0.5 * (F₁₂*F₁₂ + F₂₂*F₂₂ - 1.0)
        E₁₂ = 0.5 * (F₁₁*F₁₂ + F₂₁*F₂₂) # 张量剪切应变


        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
        end

     # First Piola-Kirchhoff
        P₁₁ = F₁₁ * S₁₁ + F₁₂ * S₁₂
        P₁₂ = F₁₁ * S₁₂ + F₁₂ * S₂₂
        P₂₁ = F₂₁ * S₁₁ + F₂₂ * S₁₂
        P₂₂ = F₂₁ * S₁₂ + F₂₂ * S₂₂

        # 3. 组装残差 (注意剪切项的系数 2.0，对应双点积 S:E)
        for (i, xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
           f[2*I-1] -= (P₁₁ * B₁[i] + P₁₂ * B₂[i]) * 𝑤
           f[2*I]   -= (P₂₁ * B₁[i] + P₂₂ * B₂[i]) * 𝑤

        end
    end
end


function ∫∫SᵢⱼδΔEᵢⱼdxdy_HR(aₛ::T,aᵤ::S, k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₛ = aₛ.𝓒;𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒;𝓖ᵤ = aᵤ.𝓖
   
    for (ξₛ,ξᵤ) in zip(𝓖ₛ,𝓖ᵤ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        Nₛ = ξₛ[:𝝭]
        N = ξᵤ[:𝝭]
        
        𝑤 = ξₛ.𝑤

        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
            u₁  += N[i]*xᵢ.d₁
            u₂  += N[i]*xᵢ.d₂
        end

        # 2. 计算 Green-Lagrange 应变 E = 0.5 * (F'F - I)
        E₁₁ = 0.5 * (F₁₁*F₁₁ + F₂₁*F₂₁ - 1.0)
        E₂₂ = 0.5 * (F₁₂*F₁₂ + F₂₂*F₂₂ - 1.0)
        E₁₂ = 0.5 * (F₁₁*F₁₂ + F₂₁*F₂₂) # 张量剪切应变


        S₁₁ = 0.0
        S₂₂ = 0.0
        S₁₂ = 0.0
        for (i,xᵢ) in  enumerate(𝓒ₛ)
           S₁₁ += Nₛ[i]*xᵢ.dₛ₁₁
           S₂₂ += Nₛ[i]*xᵢ.dₛ₂₂
           S₁₂ += Nₛ[i]*xᵢ.dₛ₁₂
        end

     # First Piola-Kirchhoff
        P₁₁ = F₁₁ * S₁₁ + F₁₂ * S₁₂
        P₁₂ = F₁₁ * S₁₂ + F₁₂ * S₂₂
        P₂₁ = F₂₁ * S₁₁ + F₂₂ * S₁₂
        P₂₂ = F₂₁ * S₁₂ + F₂₂ * S₂₂

        # 3. 组装残差 (注意剪切项的系数 2.0，对应双点积 S:E)
        for (i, xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
            b1_i = B₁[i]
            b2_i = B₂[i]
            for (j, xⱼ) in enumerate(𝓒ᵤ)
                 J = xⱼ.𝐼
                b1_j = B₁[j]
                b2_j = B₂[j]
                term1 = S₁₁ * b1_i * b1_j
                term2 = S₁₂ * (b1_i * b2_j + b2_i * b1_j) # 对称剪切项
                term3 = S₂₂ * b2_i * b2_j
                g = (term1 + term2 + term3) * 𝑤
                k[2*I-1, 2*J-1] -= g
                k[2*I,   2*J]   -= g


            end

        end
    end
end

function ∫∫Stabilization_Operator_HR(
    aₛ::T, aᵤ::S,
    f_u::AbstractVector{Float64}, f_S::AbstractVector{Float64},
    k_uu::AbstractMatrix{Float64}, # 加入了缺失的 k_uS
    k_Su::AbstractMatrix{Float64}, k_SS::AbstractMatrix{Float64}
) where {T<:AbstractElement, S<:AbstractElement}
    
    𝓒ₛ = aₛ.𝓒; 𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒; 𝓖ᵤ = aᵤ.𝓖
    b₀ = [0.0, 0.0]
    
    for (ξₛ, ξᵤ) in zip(𝓖ₛ, 𝓖ᵤ)
        τ = ξₛ.τ
        B₁ = ξᵤ[:∂𝝭∂x]; B₂ = ξᵤ[:∂𝝭∂y]
        B₁₁ = ξᵤ[:∂²𝝭∂x²]; B₂₂ = ξᵤ[:∂²𝝭∂y²]; B₁₂ = ξᵤ[:∂²𝝭∂x∂y]
        
        Nₛ = ξₛ[:𝝭]; Bₛ₁ = ξₛ[:∂𝝭∂x]; Bₛ₂ = ξₛ[:∂𝝭∂y]
        𝑤 = ξₛ.𝑤
        
        # -------------------------------------------------------------
        # 3. 计算当前积分点的物理场 (F, S, divS, 以及二阶补偿项)
        # -------------------------------------------------------------
        F₁₁=1.0; F₁₂=0.0; F₂₁=0.0; F₂₂=1.0
        u1_11=0.0; u1_12=0.0; u1_22=0.0
        u2_11=0.0; u2_12=0.0; u2_22=0.0
        
        for (m, xₘ) in enumerate(𝓒ᵤ)
            u1 = xₘ.d₁; u2 = xₘ.d₂
            F₁₁ += B₁[m]*u1; F₁₂ += B₂[m]*u1
            F₂₁ += B₁[m]*u2; F₂₂ += B₂[m]*u2
            
            u1_11 += B₁₁[m]*u1; u1_12 += B₁₂[m]*u1; u1_22 += B₂₂[m]*u1
            u2_11 += B₁₁[m]*u2; u2_12 += B₁₂[m]*u2; u2_22 += B₂₂[m]*u2
        end

        S₁₁=0.0; S₂₂=0.0; S₁₂=0.0
        divS_1=0.0; divS_2=0.0
        for (m, xₘ) in enumerate(𝓒ₛ)
            s11=xₘ.dₛ₁₁; s22=xₘ.dₛ₂₂; s12=xₘ.dₛ₁₂
            S₁₁ += Nₛ[m]*s11; S₂₂ += Nₛ[m]*s22; S₁₂ += Nₛ[m]*s12
            divS_1 += Bₛ₁[m]*s11 + Bₛ₂[m]*s12
            divS_2 += Bₛ₁[m]*s12 + Bₛ₂[m]*s22
        end

        # -------------------------------------------------------------
        # 4. 计算严格的强形式残差 R1, R2
        # -------------------------------------------------------------
        comp_1 = S₁₁*u1_11 + 2.0*S₁₂*u1_12 + S₂₂*u1_22
        comp_2 = S₁₁*u2_11 + 2.0*S₁₂*u2_12 + S₂₂*u2_22
        
        R1 = F₁₁*divS_1 + F₁₂*divS_2 + comp_1 + b₀[1]
        R2 = F₂₁*divS_1 + F₂₂*divS_2 + comp_2 + b₀[2]

        # -------------------------------------------------------------
        # 5. 组装残差 f 和刚度矩阵 K
        # -------------------------------------------------------------
        for (i, xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            # 应力测试函数的变分向量: δR^S
            vS11_1 = F₁₁*Bₛ₁[i] + u1_11*Nₛ[i];  vS11_2 = F₂₁*Bₛ₁[i] + u2_11*Nₛ[i]
            vS22_1 = F₁₂*Bₛ₂[i] + u1_22*Nₛ[i];  vS22_2 = F₂₂*Bₛ₂[i] + u2_22*Nₛ[i]
            vS12_1 = F₁₁*Bₛ₂[i] + F₁₂*Bₛ₁[i] + 2.0*u1_12*Nₛ[i]
            vS12_2 = F₂₁*Bₛ₂[i] + F₂₂*Bₛ₁[i] + 2.0*u2_12*Nₛ[i]

            # (A) 组装 f_S
            f_S[3*I-2] -= τ * (vS11_1*R1 + vS11_2*R2) * 𝑤
            f_S[3*I-1] -= τ * (vS22_1*R1 + vS22_2*R2) * 𝑤
            f_S[3*I]   -= τ * (vS12_1*R1 + vS12_2*R2) * 𝑤

            # (B) 组装 K_SS (保持不变)
            for (j, xⱼ) in enumerate(𝓒ₛ)
                J = xⱼ.𝐼
                rS11_1 = F₁₁*Bₛ₁[j] + u1_11*Nₛ[j];  rS11_2 = F₂₁*Bₛ₁[j] + u2_11*Nₛ[j]
                rS22_1 = F₁₂*Bₛ₂[j] + u1_22*Nₛ[j];  rS22_2 = F₂₂*Bₛ₂[j] + u2_22*Nₛ[j]
                rS12_1 = F₁₁*Bₛ₂[j] + F₁₂*Bₛ₁[j] + 2.0*u1_12*Nₛ[j]
                rS12_2 = F₂₁*Bₛ₂[j] + F₂₂*Bₛ₁[j] + 2.0*u2_12*Nₛ[j]

                k_SS[3*I-2, 3*J-2] += τ * (vS11_1*rS11_1 + vS11_2*rS11_2) * 𝑤
                k_SS[3*I-2, 3*J-1] += τ * (vS11_1*rS22_1 + vS11_2*rS22_2) * 𝑤
                k_SS[3*I-2, 3*J]   += τ * (vS11_1*rS12_1 + vS11_2*rS12_2) * 𝑤
                k_SS[3*I-1, 3*J-2] += τ * (vS22_1*rS11_1 + vS22_2*rS11_2) * 𝑤
                k_SS[3*I-1, 3*J-1] += τ * (vS22_1*rS22_1 + vS22_2*rS22_2) * 𝑤
                k_SS[3*I-1, 3*J]   += τ * (vS22_1*rS12_1 + vS22_2*rS12_2) * 𝑤
                k_SS[3*I,   3*J-2] += τ * (vS12_1*rS11_1 + vS12_2*rS11_2) * 𝑤
                k_SS[3*I,   3*J-1] += τ * (vS12_1*rS22_1 + vS12_2*rS22_2) * 𝑤
                k_SS[3*I,   3*J]   += τ * (vS12_1*rS12_1 + vS12_2*rS12_2) * 𝑤
            end
        end

        for (i, xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
            vu_i = B₁[i]*divS_1 + B₂[i]*divS_2 + S₁₁*B₁₁[i] + 2.0*S₁₂*B₁₂[i] + S₂₂*B₂₂[i]

            # (C) 组装 f_u
            f_u[2*I-1] -= τ * (vu_i * R1) * 𝑤
            f_u[2*I]   -= τ * (vu_i * R2) * 𝑤

            # (D) 组装 K_uu (保持不变)
            for (j, xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                ru_j = B₁[j]*divS_1 + B₂[j]*divS_2 + S₁₁*B₁₁[j] + 2.0*S₁₂*B₁₂[j] + S₂₂*B₂₂[j]
                k_uu[2*I-1, 2*J-1] += τ * (vu_i * ru_j) * 𝑤
                k_uu[2*I,   2*J]   += τ * (vu_i * ru_j) * 𝑤
            end

            # (E) 组装 K_uS 和 K_Su 【核心修正区域】
            for (j, xⱼ) in enumerate(𝓒ₛ)
                J = xⱼ.𝐼
                
                # [部分 1] 材料刚度贡献 (来自 δR * ΔR)
                rS11_1 = F₁₁*Bₛ₁[j] + u1_11*Nₛ[j];  rS11_2 = F₂₁*Bₛ₁[j] + u2_11*Nₛ[j]
                rS22_1 = F₁₂*Bₛ₂[j] + u1_22*Nₛ[j];  rS22_2 = F₂₂*Bₛ₂[j] + u2_22*Nₛ[j]
                rS12_1 = F₁₁*Bₛ₂[j] + F₁₂*Bₛ₁[j] + 2.0*u1_12*Nₛ[j]
                rS12_2 = F₂₁*Bₛ₂[j] + F₂₂*Bₛ₁[j] + 2.0*u2_12*Nₛ[j]

                mat_x11 = τ * vu_i * rS11_1 * 𝑤
                mat_x22 = τ * vu_i * rS22_1 * 𝑤
                mat_x12 = τ * vu_i * rS12_1 * 𝑤
                
                mat_y11 = τ * vu_i * rS11_2 * 𝑤
                mat_y22 = τ * vu_i * rS22_2 * 𝑤
                mat_y12 = τ * vu_i * rS12_2 * 𝑤

                # [部分 2] 几何刚度贡献 (来自 Δ(δR) * R) => 补全二次收敛的钥匙
                G11 = B₁[i]*Bₛ₁[j] + B₁₁[i]*Nₛ[j]
                G22 = B₂[i]*Bₛ₂[j] + B₂₂[i]*Nₛ[j]
                G12 = B₁[i]*Bₛ₂[j] + B₂[i]*Bₛ₁[j] + 2.0*B₁₂[i]*Nₛ[j]

                geo_x11 = τ * G11 * R1 * 𝑤
                geo_x22 = τ * G22 * R1 * 𝑤
                geo_x12 = τ * G12 * R1 * 𝑤
                
                geo_y11 = τ * G11 * R2 * 𝑤
                geo_y22 = τ * G22 * R2 * 𝑤
                geo_y12 = τ * G12 * R2 * 𝑤

                # 合并两个贡献
                term_x11 = mat_x11 + geo_x11
                term_x22 = mat_x22 + geo_x22
                term_x12 = mat_x12 + geo_x12
                
                term_y11 = mat_y11 + geo_y11
                term_y22 = mat_y22 + geo_y22
                term_y12 = mat_y12 + geo_y12


                # 填充 K_Su (严格对称位置)
                k_Su[3*J-2, 2*I-1] += term_x11
                k_Su[3*J-1, 2*I-1] += term_x22
                k_Su[3*J,   2*I-1] += term_x12
                
                k_Su[3*J-2, 2*I]   += term_y11
                k_Su[3*J-1, 2*I]   += term_y22
                k_Su[3*J,   2*I]   += term_y12
            end
        end
    end
end

function ∫∫Stabilization_Operator_HR_New(
    aₛ::T, aᵤ::S,
    f_u::AbstractVector{Float64}, f_S::AbstractVector{Float64},
    k_uu::AbstractMatrix{Float64}, 
    k_Su::AbstractMatrix{Float64}, # 只需传入 k_Su，外部组装时转置 k_uS = k_Su'
    k_SS::AbstractMatrix{Float64}
) where {T<:AbstractElement, S<:AbstractElement}
    
    𝓒ₛ = aₛ.𝓒; 𝓖ₛ = aₛ.𝓖
    𝓒ᵤ = aᵤ.𝓒; 𝓖ᵤ = aᵤ.𝓖
    
    for (ξₛ, ξᵤ) in zip(𝓖ₛ, 𝓖ᵤ)
        τ = ξₛ.τ
        ℎ = ξₛ.ℎ
        a=0.9/ℎ^2
       
        B₁ = ξᵤ[:∂𝝭∂x]; B₂ = ξᵤ[:∂𝝭∂y]
        Nₛ = ξₛ[:𝝭]
        𝑤 = ξₛ.𝑤
        
        # -------------------------------------------------------------
        # 1. 计算当前积分点的物理场 F 和 S
        # -------------------------------------------------------------
        F₁₁=1.0; F₁₂=0.0; F₂₁=0.0; F₂₂=1.0
        for (m, xₘ) in enumerate(𝓒ᵤ)
            u1 = xₘ.d₁; u2 = xₘ.d₂
            F₁₁ += B₁[m]*u1; F₁₂ += B₂[m]*u1
            F₂₁ += B₁[m]*u2; F₂₂ += B₂[m]*u2
        end

        S₁₁=0.0; S₂₂=0.0; S₁₂=0.0
        for (m, xₘ) in enumerate(𝓒ₛ)
            s11=xₘ.dₛ₁₁; s22=xₘ.dₛ₂₂; s12=xₘ.dₛ₁₂
            S₁₁ += Nₛ[m]*s11; S₂₂ += Nₛ[m]*s22; S₁₂ += Nₛ[m]*s12
        end

        # -------------------------------------------------------------
        # 2. 计算残差张量 R = FS
        # -------------------------------------------------------------
        R₁₁ = F₁₁*S₁₁ + F₁₂*S₁₂
        R₁₂ = F₁₁*S₁₂ + F₁₂*S₂₂
        R₂₁ = F₂₁*S₁₁ + F₂₂*S₁₂
        R₂₂ = F₂₁*S₁₂ + F₂₂*S₂₂

        # -------------------------------------------------------------
        # 3. 组装残差 f 和刚度矩阵 K
        # -------------------------------------------------------------
        
        # [A] 组装应力部分 f_S 和 k_SS
        for (i, xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            N_i = Nₛ[i]
            
            # 残差 f_S
            f_S[3*I-2] -= a*τ * N_i * (F₁₁*R₁₁ + F₂₁*R₂₁) * 𝑤
            f_S[3*I-1] -= a*τ * N_i * (F₁₂*R₁₂ + F₂₂*R₂₂) * 𝑤
            f_S[3*I]   -= a*τ * N_i * (F₁₂*R₁₁ + F₁₁*R₁₂ + F₂₂*R₂₁ + F₂₁*R₂₂) * 𝑤

            # 刚度 k_SS
            for (j, xⱼ) in enumerate(𝓒ₛ)
                J = xⱼ.𝐼
                NN = N_i * Nₛ[j]
                
                term_11_11 = F₁₁^2 + F₂₁^2
                term_22_22 = F₁₂^2 + F₂₂^2
                term_cross = F₁₁*F₁₂ + F₂₁*F₂₂
                term_12_12 = F₁₂^2 + F₁₁^2 + F₂₂^2 + F₂₁^2

                k_SS[3*I-2, 3*J-2] += a*τ * NN * term_11_11 * 𝑤
                # k_SS[3*I-2, 3*J-1] 是 0
                k_SS[3*I-2, 3*J]   += a*τ * NN * term_cross * 𝑤
                
                # k_SS[3*I-1, 3*J-2] 是 0
                k_SS[3*I-1, 3*J-1] += a*τ * NN * term_22_22 * 𝑤
                k_SS[3*I-1, 3*J]   += a*τ * NN * term_cross * 𝑤
                
                k_SS[3*I,   3*J-2] += a*τ * NN * term_cross * 𝑤
                k_SS[3*I,   3*J-1] += a*τ * NN * term_cross * 𝑤
                k_SS[3*I,   3*J]   += a*τ * NN * term_12_12 * 𝑤
            end
        end

        # [B] 组装位移部分 f_u 和 k_uu
        for (i, xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
            V1_i = B₁[i]*S₁₁ + B₂[i]*S₁₂
            V2_i = B₁[i]*S₁₂ + B₂[i]*S₂₂

            # 残差 f_u
            f_u[2*I-1] -= a*τ * (V1_i*R₁₁ + V2_i*R₁₂) * 𝑤
            f_u[2*I]   -= a*τ * (V1_i*R₂₁ + V2_i*R₂₂) * 𝑤

            # 刚度 k_uu
            for (j, xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                V1_j = B₁[j]*S₁₁ + B₂[j]*S₁₂
                V2_j = B₁[j]*S₁₂ + B₂[j]*S₂₂
                
                uu_val = a*τ * (V1_i*V1_j + V2_i*V2_j) * 𝑤
                
                k_uu[2*I-1, 2*J-1] += uu_val
                k_uu[2*I,   2*J]   += uu_val
                # 交叉项 ux-uy 为 0
            end

            # [C] 组装交叉刚度 k_Su (行:应力, 列:位移)
            for (j, xⱼ) in enumerate(𝓒ₛ)
                J = xⱼ.𝐼
                N_j = Nₛ[j]
                
                # 1. 材料刚度贡献 (来自 δF S : F ΔS)
                mat_x11 = V1_i * F₁₁ * N_j
                mat_x22 = V2_i * F₁₂ * N_j
                mat_x12 = V1_i * F₁₂ * N_j + V2_i * F₁₁ * N_j
                
                mat_y11 = V1_i * F₂₁ * N_j
                mat_y22 = V2_i * F₂₂ * N_j
                mat_y12 = V1_i * F₂₂ * N_j + V2_i * F₂₁ * N_j

                # 2. 几何刚度贡献 (来自 δF ΔS : R)
                geo_x11 = B₁[i] * N_j * R₁₁
                geo_x22 = B₂[i] * N_j * R₁₂
                geo_x12 = B₂[i] * N_j * R₁₁ + B₁[i] * N_j * R₁₂
                
                geo_y11 = B₁[i] * N_j * R₂₁
                geo_y22 = B₂[i] * N_j * R₂₂
                geo_y12 = B₂[i] * N_j * R₂₁ + B₁[i] * N_j * R₂₂

                # 填充 k_Su (利用 k_uS = k_Su'，所以在这里直接给 k_Su 赋值即可)
                k_Su[3*J-2, 2*I-1] += a*τ * (mat_x11 + geo_x11) * 𝑤
                k_Su[3*J-1, 2*I-1] += a*τ * (mat_x22 + geo_x22) * 𝑤
                k_Su[3*J,   2*I-1] += a*τ * (mat_x12 + geo_x12) * 𝑤
                
                k_Su[3*J-2, 2*I]   += a*τ * (mat_y11 + geo_y11) * 𝑤
                k_Su[3*J-1, 2*I]   += a*τ * (mat_y22 + geo_y22) * 𝑤
                k_Su[3*J,   2*I]   += a*τ * (mat_y12 + geo_y12) * 𝑤
            end
        end
    end
end








end