## cal𝗚
function cal𝗚!(ap::ReproducingKernel)
    𝓖 = ap.𝓖
    𝗚 = ap.𝗠[:∇̃]
    n = get𝑛𝒒(ap)
    fill!(𝗚,0.0)
    for ξ in 𝓖
        w = get𝑤(ap,ξ)
        𝒒 = get𝒒(ap,ξ)
        for I in 1:n
            for J in I:n
                𝗚[I,J] += w*𝒒[I]*𝒒[J]
            end
        end
    end
    cholesky!(𝗚)
    U⁻¹ = inverse!(𝗚)
    𝗚⁻¹ = UUᵀ!(U⁻¹)
    return 𝗚⁻¹
end

function cal𝗚!(ap::ReproducingKernel{𝝃,:Quadratic1D,𝑠,𝜙,:Seg2}) where {𝝃<:AbstractNode,𝑠,𝜙}
    𝗚⁻¹ = ap.𝗠[:∇̃]
    fill!(𝗚⁻¹,0.0)
    𝐿 = get𝐿(ap)
    𝗚⁻¹[1] =  4.0/𝐿
    𝗚⁻¹[2] = -6.0/𝐿
    𝗚⁻¹[3] = 12.0/𝐿
    return 𝗚⁻¹
end

function cal𝗚!(ap::ReproducingKernel{𝝃,:Cubic1D,𝑠,𝜙,:Seg2}) where {𝝃<:AbstractNode,𝑠,𝜙}
    𝗚⁻¹ = ap.𝗠[:∇̃]
    fill!(𝗚⁻¹,0.0)
    𝐿 = get𝐿(ap)
    𝗚⁻¹[1] =    9.0/𝐿
    𝗚⁻¹[2] =  -36.0/𝐿
    𝗚⁻¹[3] =  192.0/𝐿
    𝗚⁻¹[4] =   30.0/𝐿
    𝗚⁻¹[5] = -180.0/𝐿
    𝗚⁻¹[6] =  180.0/𝐿
    return 𝗚⁻¹
end

function cal𝗚!(ap::ReproducingKernel{𝝃,:Linear2D,𝑠,𝜙,:Tri3}) where {𝝃<:AbstractNode,𝑠,𝜙}
    𝗚⁻¹ = ap.𝗠[:∇̃]
    fill!(𝗚⁻¹,0.0)
    𝐴 = get𝐴(ap)
    𝗚⁻¹[1] = 1.0/𝐴
    return 𝗚⁻¹
end

function cal𝗚!(ap::ReproducingKernel{𝝃,:Quadratic2D,𝑠,𝜙,:Tri3}) where {𝝃<:AbstractNode,𝑠,𝜙}
    𝗚⁻¹ = ap.𝗠[:∇̃]
    fill!(𝗚⁻¹,0.0)
    𝐴 = get𝐴(ap)
    𝗚⁻¹[1] =   9.0/𝐴
    𝗚⁻¹[2] = -12.0/𝐴
    𝗚⁻¹[3] =  24.0/𝐴
    𝗚⁻¹[4] = -12.0/𝐴
    𝗚⁻¹[5] =  12.0/𝐴
    𝗚⁻¹[6] =  24.0/𝐴
    return 𝗚⁻¹
end

function cal𝗚!(ap::ReproducingKernel{𝝃,:Cubic2D,𝑠,𝜙,:Tri3}) where {𝝃<:AbstractNode,𝑠,𝜙}
    𝗚⁻¹ = ap.𝗠[:∇̃]
    fill!(𝗚⁻¹,0.0)
    𝐴 = get𝐴(ap)
    𝗚⁻¹[1] =   36.0/𝐴
    𝗚⁻¹[2] = -120.0/𝐴
    𝗚⁻¹[3] =  600.0/𝐴
    𝗚⁻¹[4] = -120.0/𝐴
    𝗚⁻¹[5] =  300.0/𝐴
    𝗚⁻¹[6] =  600.0/𝐴
    𝗚⁻¹[7] =   90.0/𝐴
    𝗚⁻¹[8] = -540.0/𝐴
    𝗚⁻¹[9] = -180.0/𝐴
    𝗚⁻¹[10] =  540.0/𝐴
    𝗚⁻¹[11] =  180.0/𝐴
    𝗚⁻¹[12] = -720.0/𝐴
    𝗚⁻¹[13] = -720.0/𝐴
    𝗚⁻¹[14] =  540.0/𝐴
    𝗚⁻¹[15] = 1440.0/𝐴
    𝗚⁻¹[16] =   90.0/𝐴
    𝗚⁻¹[17] = -180.0/𝐴
    𝗚⁻¹[18] = -540.0/𝐴
    𝗚⁻¹[19] =   90.0/𝐴
    𝗚⁻¹[20] =  540.0/𝐴
    𝗚⁻¹[21] =  540.0/𝐴
    return 𝗚⁻¹
end

function cal𝗚₂!(ap::ReproducingKernel{𝝃,:Quadratic2D,𝑠,𝜙,:Tri3}) where {𝝃<:AbstractNode,𝑠,𝜙}
    𝗚⁻¹ = ap.𝗠[:∇̃]
    fill!(𝗚⁻¹,0.0)
    𝐴 = get𝐴(ap)
    𝗚⁻¹[1] = 1.0/𝐴
    return 𝗚⁻¹
end
function cal𝗚₂!(ap::ReproducingKernel{𝝃,:Cubic2D,𝑠,𝜙,:Tri3}) where {𝝃<:AbstractNode,𝑠,𝜙}
    𝗚⁻¹ = ap.𝗠[:∇̃]
    fill!(𝗚⁻¹,0.0)
    𝐴 = get𝐴(ap)
    𝗚⁻¹[1] =   9.0/𝐴
    𝗚⁻¹[2] = -12.0/𝐴
    𝗚⁻¹[3] =  24.0/𝐴
    𝗚⁻¹[4] = -12.0/𝐴
    𝗚⁻¹[5] =  12.0/𝐴
    𝗚⁻¹[6] =  24.0/𝐴
    return 𝗚⁻¹
end

function cal𝗚₂!(ap::ReproducingKernel{𝝃,:Quartic2D,𝑠,𝜙,:Tri3}) where {𝝃<:AbstractNode,𝑠,𝜙}
    𝗚⁻¹ = ap.𝗠[:∇̃]
    fill!(𝗚⁻¹,0.0)
    𝐴 = get𝐴(ap)
    𝗚⁻¹[1] =   36.0/𝐴
    𝗚⁻¹[2] = -120.0/𝐴
    𝗚⁻¹[3] =  600.0/𝐴
    𝗚⁻¹[4] = -120.0/𝐴
    𝗚⁻¹[5] =  300.0/𝐴
    𝗚⁻¹[6] =  600.0/𝐴
    𝗚⁻¹[7] =   90.0/𝐴
    𝗚⁻¹[8] = -540.0/𝐴
    𝗚⁻¹[9] = -180.0/𝐴
    𝗚⁻¹[10] =  540.0/𝐴
    𝗚⁻¹[11] =  180.0/𝐴
    𝗚⁻¹[12] = -720.0/𝐴
    𝗚⁻¹[13] = -720.0/𝐴
    𝗚⁻¹[14] =  540.0/𝐴
    𝗚⁻¹[15] = 1440.0/𝐴
    𝗚⁻¹[16] =   90.0/𝐴
    𝗚⁻¹[17] = -180.0/𝐴
    𝗚⁻¹[18] = -540.0/𝐴
    𝗚⁻¹[19] =   90.0/𝐴
    𝗚⁻¹[20] =  540.0/𝐴
    𝗚⁻¹[21] =  540.0/𝐴
    return 𝗚⁻¹
end

## get shape for SNode
function get∇∇̃²𝝭(ap::ReproducingKernel,ξ::SNode)
    𝝭 = ap.𝝭[:∂1]
    ∂𝝭∂x = ap.𝝭[:∂x]
    ∂𝝭∂y = ap.𝝭[:∂y]
    ∂²𝝭∂x² = ap.𝝭[:∂x²]
    ∂²𝝭∂x∂y = ap.𝝭[:∂x∂y]
    ∂²𝝭∂y² = ap.𝝭[:∂y²]
    ∂³𝝭∂x³ = ap.𝝭[:∂x³]
    ∂³𝝭∂x²∂y = ap.𝝭[:∂x²∂y]
    ∂³𝝭∂x∂y² = ap.𝝭[:∂x∂y²]
    ∂³𝝭∂y³ = ap.𝝭[:∂y³]
    i = ξ.id
    index = ξ.index
    for j in 1:length(ap.𝓒)
        𝝭[j] = ξ.𝝭[:∂1][index[i]+j]
        ∂𝝭∂x[j] = ξ.𝝭[:∂x][index[i]+j]
        ∂𝝭∂y[j] = ξ.𝝭[:∂y][index[i]+j]
        ∂²𝝭∂x²[j] = ξ.𝝭[:∂x²][index[i]+j]
        ∂²𝝭∂x∂y[j] = ξ.𝝭[:∂x∂y][index[i]+j]
        ∂²𝝭∂y²[j] = ξ.𝝭[:∂y²][index[i]+j]
        ∂³𝝭∂x³[j] = ξ.𝝭[:∂x³][index[i]+j]
        ∂³𝝭∂x²∂y[j] = ξ.𝝭[:∂x²∂y][index[i]+j]
        ∂³𝝭∂x∂y²[j] = ξ.𝝭[:∂x∂y²][index[i]+j]
        ∂³𝝭∂y³[j] = ξ.𝝭[:∂y³][index[i]+j]
    end
    return 𝝭, ∂𝝭∂x, ∂𝝭∂y, ∂²𝝭∂x², ∂²𝝭∂x∂y, ∂²𝝭∂y², ∂³𝝭∂x³, ∂³𝝭∂x²∂y, ∂³𝝭∂x∂y², ∂³𝝭∂y³
end

function get∇̄𝝭(ap::ReproducingKernel,ξ::SNode)
    𝝭 = ap.𝝭[:∂1]
    ∂𝝭∂x = ap.𝝭[:∂x]
    ∂𝝭∂y = ap.𝝭[:∂y]
    ∂𝝭∂z = ap.𝝭[:∂z]
    i = ξ.id
    index = ξ.index
    for j in 1:length(ap.𝓒)
        𝝭[j] = ξ.𝝭[:∂1][index[i]+j]
        ∂𝝭∂x[j] = ξ.𝝭[:∂̄x][index[i]+j]
        ∂𝝭∂y[j] = haskey(ξ.𝝭,:∂̄y) ? ξ.𝝭[:∂̄y][index[i]+j] : 0.0
        ∂𝝭∂z[j] = haskey(ξ.𝝭,:∂̄z) ? ξ.𝝭[:∂̄z][index[i]+j] : 0.0
    end
    return 𝝭, ∂𝝭∂x, ∂𝝭∂y, ∂𝝭∂z
end
## RK gradient smoothing
function set∇̃𝝭!(aps::Vector{T}) where T<:ReproducingKernel
    for ap in aps
        set∇̃𝝭!(ap)
    end
end
set∇̃𝝭!(ap::T) where T<:ReproducingKernel{SNode} = set∇̃𝝭!(ap,ap)

function set∇̃²𝝭!(aps::Vector{T}) where T<:ReproducingKernel
    for ap in aps
        set∇̃²𝝭!(ap)
    end
end
set∇̃²𝝭!(ap::T) where T<:ReproducingKernel{SNode} = set∇̃²𝝭!(ap,ap)

function set∇̃𝝭!(gps::Vector{T},aps::Vector{S}) where {T<:ReproducingKernel,S<:ReproducingKernel}
    if length(gps) ≠ length(aps)
        error("Miss match element numbers")
    else
        for i in 1:length(gps)
            set∇̃𝝭!(gps[i],aps[i])
        end
    end
end

function set∇̃²𝝭!(gps::Vector{T},aps::Vector{S}) where {T<:ReproducingKernel,S<:ReproducingKernel}
    if length(gps) ≠ length(aps)
        error("Miss match element numbers")
    else
        for i in 1:length(gps)
            set∇̃²𝝭!(gps[i],aps[i])
        end
    end
end

function set∇∇̃²𝝭!(gps::Vector{T},aps::Vector{S}) where {T<:ReproducingKernel,S<:ReproducingKernel}
    if length(gps) ≠ length(aps)
        error("Miss match element numbers")
    else
        for i in 1:length(gps)
            set∇∇̃²𝝭!(gps[i],aps[i])
        end
    end
end

function set∇̃𝝭!(as::Vector{T},bs::Vector{S},cs::Vector{R}) where {T<:ReproducingKernel,S<:ReproducingKernel,R<:ReproducingKernel}
    if length(as) ≠ length(bs) || length(bs) ≠ length(cs)
        error("Miss match element numbers")
    else
        for i in 1:length(as)
            set∇̃𝝭!(as[i],bs[i],cs[i])
        end
    end
end
function set∇̃𝝭!(gp::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Seg2},ap::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Seg2}) where {𝒑,𝑠,𝜙}
    n₁ =  1.0
    n₂ = -1.0
    𝓒 = gp.𝓒
    𝓖 = gp.𝓖
    for ξ̂ in 𝓖
        𝒒̂ = get𝒒(gp,ξ̂)
        𝗚⁻¹ = cal𝗚!(gp)
        𝒒̂ᵀ𝗚⁻¹ = 𝒒̂*𝗚⁻¹
        ∂𝝭∂x = gp.𝝭[:∂x]
        fill!(∂𝝭∂x,0.0)
        for ξ in ap.𝓖
            w = ξ.w/2
            wᵇ = ξ.wᵇ
            nᵇ₁ = 0.0
            nᵇ₁ += ξ.ξ ==  1.0 ? n₁ : 0.0
            nᵇ₁ += ξ.ξ == -1.0 ? n₂ : 0.0
            𝝭 = get𝝭(ap,ξ)
            𝒒, ∂𝒒∂ξ = get∇𝒒(gp,ξ)
            W₁ = 𝒒̂ᵀ𝗚⁻¹*𝒒*nᵇ₁*wᵇ + 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂ξ*n₁*w
            for i in 1:length(𝓒)
                ∂𝝭∂x[i] += 𝝭[i]*W₁
            end
        end
        for i in 1:length(𝓒)
            ξ̂.𝝭[:∂x][ξ̂.index[ξ̂.id]+i] = ∂𝝭∂x[i]
        end
    end
end

function set∇̃𝝭!(gp::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Tri3},ap::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Tri3}) where {𝒑,𝑠,𝜙}
    x₁ = gp.𝓒[1].x;y₁ = gp.𝓒[1].y
    x₂ = gp.𝓒[2].x;y₂ = gp.𝓒[2].y
    x₃ = gp.𝓒[3].x;y₃ = gp.𝓒[3].y
    n₁₁ = y₃-y₂;n₂₁ = y₁-y₃;n₃₁ = y₂-y₁
    n₁₂ = x₂-x₃;n₂₂ = x₃-x₁;n₃₂ = x₁-x₂
    𝓒 = gp.𝓒
    𝓖 = gp.𝓖
    for ξ̂ in 𝓖
        𝒒̂ = get𝒑₁(gp,ξ̂)
        𝗚⁻¹ = cal𝗚!(gp)
        𝒒̂ᵀ𝗚⁻¹ = 𝒒̂*𝗚⁻¹
        ∂𝝭∂x = gp.𝝭[:∂x]
        ∂𝝭∂y = gp.𝝭[:∂y]
        fill!(∂𝝭∂x,0.0)
        fill!(∂𝝭∂y,0.0)
        for ξ in ap.𝓖
            w = ξ.w
            wᵇ = ξ.wᵇ
            𝝭 = get𝝭(ap,ξ)
            𝒒, ∂𝒒∂ξ, ∂𝒒∂η = get∇𝒑₁(ap,ξ)
            𝒒̂ᵀ𝗚⁻¹𝒒 =  𝒒̂ᵀ𝗚⁻¹*𝒒
            𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ = 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂ξ
            𝒒̂ᵀ𝗚⁻¹∂𝒒∂η = 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂η
            nᵇ₁ = 0.0;nᵇ₂ = 0.0
            ξ.ξ == 0.0 ? (nᵇ₁ += n₁₁;nᵇ₂ += n₁₂) : nothing
            ξ.η == 0.0 ? (nᵇ₁ += n₂₁;nᵇ₂ += n₂₂) : nothing
            ξ.ξ+ξ.η ≈ 1.0 ? (nᵇ₁ += n₃₁;nᵇ₂ += n₃₂) : nothing
            b₁ = 𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ*n₁₁ + 𝒒̂ᵀ𝗚⁻¹∂𝒒∂η*n₂₁
            b₂ = 𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ*n₁₂ + 𝒒̂ᵀ𝗚⁻¹∂𝒒∂η*n₂₂
            W₁ = 𝒒̂ᵀ𝗚⁻¹𝒒*nᵇ₁*wᵇ + b₁*w/2
            W₂ = 𝒒̂ᵀ𝗚⁻¹𝒒*nᵇ₂*wᵇ + b₂*w/2
            for i in 1:length(𝓒)
                ∂𝝭∂x[i] += 𝝭[i]*W₁
                ∂𝝭∂y[i] += 𝝭[i]*W₂
            end
        end
        for i in 1:length(𝓒)
            ξ̂.𝝭[:∂x][ξ̂.index[ξ̂.id]+i] = ∂𝝭∂x[i]
            ξ̂.𝝭[:∂y][ξ̂.index[ξ̂.id]+i] = ∂𝝭∂y[i]
        end
    end
end

function set∇̃²𝝭!(gp::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Tri3},ap::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Tri3}) where {𝒑,𝑠,𝜙}
    x₁ = gp.𝓒[1].x;y₁ = gp.𝓒[1].y
    x₂ = gp.𝓒[2].x;y₂ = gp.𝓒[2].y
    x₃ = gp.𝓒[3].x;y₃ = gp.𝓒[3].y
    𝐴 = get𝐴(gp)
    n₁₁ = y₃-y₂;n₂₁ = y₁-y₃;n₃₁ = y₂-y₁
    n₁₂ = x₂-x₃;n₂₂ = x₃-x₁;n₃₂ = x₁-x₂
    s₁₁ = -n₁₂;s₂₁ = -n₂₂;s₃₁ = -n₃₂
    s₁₂ =  n₁₁;s₂₂ =  n₂₁;s₃₂ =  n₃₁
    𝓒 = gp.𝓒
    𝓖 = gp.𝓖
    for ξ̂ in 𝓖
        𝒒̂ = get𝒑₂(gp,ξ̂)
        𝗚⁻¹ = cal𝗚₂!(gp)
        𝒒̂ᵀ𝗚⁻¹ = 𝒒̂*𝗚⁻¹
        ∂²𝝭∂x² = gp.𝝭[:∂x²]
        ∂²𝝭∂x∂y = gp.𝝭[:∂x∂y]
        ∂²𝝭∂y² = gp.𝝭[:∂y²]
        fill!(∂²𝝭∂x²,0.0)
        fill!(∂²𝝭∂x∂y,0.0)
        fill!(∂²𝝭∂y²,0.0)
        for ξ in ap.𝓖
            w = ξ.w
            wᵇ = ξ.wᵇ
            𝝭,∂𝝭∂x,∂𝝭∂y = get∇𝝭(ap,ξ)
            𝒒, ∂𝒒∂ξ, ∂𝒒∂η, ∂²𝒒∂ξ², ∂²𝒒∂ξ∂η, ∂²𝒒∂η² = get∇²𝒑₂(ap,ξ)
            𝒒̂ᵀ𝗚⁻¹𝒒 =  𝒒̂ᵀ𝗚⁻¹*𝒒
            𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ = 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂ξ
            𝒒̂ᵀ𝗚⁻¹∂𝒒∂η = 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂η
            𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ² = 𝒒̂ᵀ𝗚⁻¹*∂²𝒒∂ξ²
            𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ∂η = 𝒒̂ᵀ𝗚⁻¹*∂²𝒒∂ξ∂η
            𝒒̂ᵀ𝗚⁻¹∂²𝒒∂η² = 𝒒̂ᵀ𝗚⁻¹*∂²𝒒∂η²

            q₁₁ = 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ²*n₁₁*n₁₁ + 2*𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ∂η*n₁₁*n₂₁ + 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂η²*n₂₁*n₂₁
            q₁₂ = 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ²*n₁₁*n₁₂ + 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ∂η*(n₁₁*n₂₂+n₁₂*n₂₁) + 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂η²*n₂₁*n₂₂
            q₂₂ = 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ²*n₁₂*n₁₂ + 2*𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ∂η*n₁₂*n₂₂ + 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂η²*n₂₂*n₂₂

            q₁ = 𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ*n₁₁ + 𝒒̂ᵀ𝗚⁻¹∂𝒒∂η*n₂₁
            q₂ = 𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ*n₁₂ + 𝒒̂ᵀ𝗚⁻¹∂𝒒∂η*n₂₂

            q₁n₁ = 0.0;q₂n₂ = 0.0;q₁n₂ = 0.0;q₂n₁ = 0.0
            mn₁₁n₁ = 0.0;mn₁₁n₂ = 0.0;mn₁₂n₁ = 0.0;mn₁₂n₂ = 0.0;mn₂₂n₁ = 0.0;mn₂₂n₂ = 0.0
            ms₁₁ = 0.0;ms₁₂ = 0.0;ms₂₂ = 0.0
            Δms₁₁ = 0.0;Δms₁₂ = 0.0;Δms₂₂ = 0.0
            𝐿₁² = n₁₁^2+n₁₂^2
            𝐿₂² = n₂₁^2+n₂₂^2
            𝐿₃² = n₃₁^2+n₃₂^2
            if ξ.ξ == 0.0
                q₁n₁ += q₁*n₁₁
                q₁n₂ += q₁*n₁₂
                q₂n₁ += q₂*n₁₁
                q₂n₂ += q₂*n₁₂
                mn₁₁n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₁₁*n₁₁*n₁₁/𝐿₁²
                mn₁₁n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₁₁*n₁₁*n₁₂/𝐿₁²
                mn₁₂n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₁₁*n₁₂*n₁₁/𝐿₁²
                mn₁₂n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₁₁*n₁₂*n₁₂/𝐿₁²
                mn₂₂n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₁₂*n₁₂*n₁₁/𝐿₁²
                mn₂₂n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₁₂*n₁₂*n₁₂/𝐿₁²
                ms₁₁ += (q₁*s₁₁+q₂*s₁₂)*n₁₁*s₁₁/𝐿₁²
                ms₁₂ += (q₁*s₁₁+q₂*s₁₂)*0.5*(n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²
                ms₂₂ += (q₁*s₁₁+q₂*s₁₂)*n₁₂*s₁₂/𝐿₁²
                if ξ.η == 0.0
                    Δms₁₁ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₁₁*s₁₁/𝐿₁²-n₂₁*s₂₁/𝐿₂²)
                    Δms₂₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₁₂*s₁₂/𝐿₁²-n₂₂*s₂₂/𝐿₂²)
                    Δms₁₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*0.5*((n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²-(n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²)
                end
            end
            if  ξ.η == 0.0
                q₁n₁ += q₁*n₂₁
                q₁n₂ += q₁*n₂₂
                q₂n₁ += q₂*n₂₁
                q₂n₂ += q₂*n₂₂
                mn₁₁n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₂₁*n₂₁*n₂₁/𝐿₂²
                mn₁₁n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₂₁*n₂₁*n₂₂/𝐿₂²
                mn₁₂n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₂₁*n₂₂*n₂₁/𝐿₂²
                mn₁₂n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₂₁*n₂₂*n₂₂/𝐿₂²
                mn₂₂n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₂₂*n₂₂*n₂₁/𝐿₂²
                mn₂₂n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₂₂*n₂₂*n₂₂/𝐿₂²
                ms₁₁ += (q₁*s₂₁+q₂*s₂₂)*n₂₁*s₂₁/𝐿₂²
                ms₁₂ += (q₁*s₂₁+q₂*s₂₂)*0.5*(n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²
                ms₂₂ += (q₁*s₂₁+q₂*s₂₂)*n₂₂*s₂₂/𝐿₂²
                if ξ.ξ+ξ.η ≈ 1.0
                    Δms₁₁ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₂₁*s₂₁/𝐿₂²-n₃₁*s₃₁/𝐿₃²)
                    Δms₂₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₂₂*s₂₂/𝐿₂²-n₃₂*s₃₂/𝐿₃²)
                    Δms₁₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*0.5*((n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²-(n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²)
                end
            end
            if ξ.ξ+ξ.η ≈ 1.0
                q₁n₁ += q₁*n₃₁
                q₁n₂ += q₁*n₃₂
                q₂n₁ += q₂*n₃₁
                q₂n₂ += q₂*n₃₂
                mn₁₁n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₃₁*n₃₁*n₃₁/𝐿₃²
                mn₁₁n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₃₁*n₃₁*n₃₂/𝐿₃²
                mn₁₂n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₃₁*n₃₂*n₃₁/𝐿₃²
                mn₁₂n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₃₁*n₃₂*n₃₂/𝐿₃²
                mn₂₂n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₃₂*n₃₂*n₃₁/𝐿₃²
                mn₂₂n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₃₂*n₃₂*n₃₂/𝐿₃²
                ms₁₁ += (q₁*s₃₁+q₂*s₃₂)*n₃₁*s₃₁/𝐿₃²
                ms₁₂ += (q₁*s₃₁+q₂*s₃₂)*0.5*(n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²
                ms₂₂ += (q₁*s₃₁+q₂*s₃₂)*n₃₂*s₃₂/𝐿₃²
                if ξ.ξ == 0.0
                    Δms₁₁ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₃₁*s₃₁/𝐿₃²-n₁₁*s₁₁/𝐿₁²)
                    Δms₂₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₃₂*s₃₂/𝐿₃²-n₁₂*s₁₂/𝐿₁²)
                    Δms₁₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*0.5*((n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²-(n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²)
                end
            end

            W₁₁₁ = mn₁₁n₁*wᵇ
            W₁₁₂ = mn₁₁n₂*wᵇ
            W₁₂₁ = mn₁₂n₁*wᵇ
            W₁₂₂ = mn₁₂n₂*wᵇ
            W₂₂₁ = mn₂₂n₁*wᵇ
            W₂₂₂ = mn₂₂n₂*wᵇ
            W₁₁ = (q₁₁*w + 2*(q₁n₁+ms₁₁)*wᵇ)/4/𝐴 + Δms₁₁
            W₁₂ = (q₁₂*w + (q₁n₂+q₂n₁+2*ms₁₂)*wᵇ)/4/𝐴 + Δms₁₂
            W₂₂ = (q₂₂*w + 2*(q₂n₂+ms₂₂)*wᵇ)/4/𝐴 + Δms₂₂
            for i in 1:length(𝓒)
                ∂²𝝭∂x²[i] += 𝝭[i]*W₁₁ + ∂𝝭∂x[i]*W₁₁₁ + ∂𝝭∂y[i]*W₁₁₂
                ∂²𝝭∂x∂y[i] += 𝝭[i]*W₁₂ + ∂𝝭∂x[i]*W₁₂₁ + ∂𝝭∂y[i]*W₁₂₂
                ∂²𝝭∂y²[i] += 𝝭[i]*W₂₂ + ∂𝝭∂x[i]*W₂₂₁ + ∂𝝭∂y[i]*W₂₂₂
            end
        end
        for i in 1:length(𝓒)
            ξ̂.𝝭[:∂x²][ξ̂.index[ξ̂.id]+i] = ∂²𝝭∂x²[i]
            ξ̂.𝝭[:∂x∂y][ξ̂.index[ξ̂.id]+i] = ∂²𝝭∂x∂y[i]
            ξ̂.𝝭[:∂y²][ξ̂.index[ξ̂.id]+i] = ∂²𝝭∂y²[i]
        end
    end
end

function set∇∇̃²𝝭!(gp::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Tri3},ap::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Tri3}) where {𝒑,𝑠,𝜙}
    x₁ = ap.𝓒[1].x;y₁ = ap.𝓒[1].y
    x₂ = ap.𝓒[2].x;y₂ = ap.𝓒[2].y
    x₃ = ap.𝓒[3].x;y₃ = ap.𝓒[3].y
    𝐴 = get𝐴(ap)
    n₁₁ = y₃-y₂;n₂₁ = y₁-y₃;n₃₁ = y₂-y₁
    n₁₂ = x₂-x₃;n₂₂ = x₃-x₁;n₃₂ = x₁-x₂
    s₁₁ = -n₁₂;s₂₁ = -n₂₂;s₃₁ = -n₃₂
    s₁₂ =  n₁₁;s₂₂ =  n₂₁;s₃₂ =  n₃₁
    𝓒 = gp.𝓒
    𝓖 = gp.𝓖
    for ξ̂ in 𝓖
        𝒒̂,∂𝒒̂∂ξ,∂𝒒̂∂η = get∇𝒑₂(gp,ξ̂)
        𝗚⁻¹ = cal𝗚₂!(gp)
        𝒒̂ᵀ𝗚⁻¹ = 𝒒̂*𝗚⁻¹
        ∂𝒒̂ᵀ∂ξ𝗚⁻¹ = ∂𝒒̂∂ξ*𝗚⁻¹
        ∂𝒒̂ᵀ∂η𝗚⁻¹ = ∂𝒒̂∂η*𝗚⁻¹
        ∂²𝝭∂x² = gp.𝝭[:∂x²]
        ∂²𝝭∂x∂y = gp.𝝭[:∂x∂y]
        ∂²𝝭∂y² = gp.𝝭[:∂y²]
        ∂∂²𝝭∂x²∂x = gp.𝝭[:∂x²∂x]
        ∂∂²𝝭∂x²∂y = gp.𝝭[:∂x²∂y]
        ∂∂²𝝭∂x∂y∂x = gp.𝝭[:∂x∂y∂x]
        ∂∂²𝝭∂x∂y∂y = gp.𝝭[:∂x∂y∂y]
        ∂∂²𝝭∂y²∂x = gp.𝝭[:∂y²∂x]
        ∂∂²𝝭∂y²∂y = gp.𝝭[:∂y²∂y]
        fill!(∂²𝝭∂x²,0.0)
        fill!(∂²𝝭∂x∂y,0.0)
        fill!(∂²𝝭∂y²,0.0)
        fill!(∂∂²𝝭∂x²∂x,0.0)
        fill!(∂∂²𝝭∂x²∂y,0.0)
        fill!(∂∂²𝝭∂x∂y∂x,0.0)
        fill!(∂∂²𝝭∂x∂y∂y,0.0)
        fill!(∂∂²𝝭∂y²∂x,0.0)
        fill!(∂∂²𝝭∂y²∂y,0.0)
        for ξ in ap.𝓖
            w = ξ.w
            wᵇ = ξ.wᵇ
            𝝭,∂𝝭∂x,∂𝝭∂y = get∇𝝭(ap,ξ)
            𝒒, ∂𝒒∂ξ, ∂𝒒∂η, ∂²𝒒∂ξ², ∂²𝒒∂ξ∂η, ∂²𝒒∂η² = get∇²𝒑₂(ap,ξ)

            𝒒̂ᵀ𝗚⁻¹𝒒 =  𝒒̂ᵀ𝗚⁻¹*𝒒
            ∂𝒒̂ᵀ∂ξ𝗚⁻¹𝒒 =  ∂𝒒̂ᵀ∂ξ𝗚⁻¹*𝒒
            ∂𝒒̂ᵀ∂η𝗚⁻¹𝒒 =  ∂𝒒̂ᵀ∂η𝗚⁻¹*𝒒
            ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒 =  ∂𝒒̂ᵀ∂ξ𝗚⁻¹𝒒*n₁₁+∂𝒒̂ᵀ∂η𝗚⁻¹𝒒*n₂₁
            ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒 =  ∂𝒒̂ᵀ∂ξ𝗚⁻¹𝒒*n₁₂+∂𝒒̂ᵀ∂η𝗚⁻¹𝒒*n₂₂

            𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ = 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂ξ
            ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂ξ = ∂𝒒̂ᵀ∂ξ𝗚⁻¹*∂𝒒∂ξ
            ∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂ξ = ∂𝒒̂ᵀ∂η𝗚⁻¹*∂𝒒∂ξ
            ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂ξ =  ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂ξ*n₁₁+∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂ξ*n₂₁
            ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂ξ =  ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂ξ*n₁₂+∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂ξ*n₂₂

            𝒒̂ᵀ𝗚⁻¹∂𝒒∂η = 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂η
            ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂η = ∂𝒒̂ᵀ∂ξ𝗚⁻¹*∂𝒒∂η
            ∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂η = ∂𝒒̂ᵀ∂η𝗚⁻¹*∂𝒒∂η
            ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂η =  ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂η*n₁₁+∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂η*n₂₁
            ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂η =  ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂η*n₁₂+∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂η*n₂₂

            𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ² = 𝒒̂ᵀ𝗚⁻¹*∂²𝒒∂ξ²
            ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂²𝒒∂ξ² = ∂𝒒̂ᵀ∂ξ𝗚⁻¹*∂²𝒒∂ξ²
            ∂𝒒̂ᵀ∂η𝗚⁻¹∂²𝒒∂ξ² = ∂𝒒̂ᵀ∂η𝗚⁻¹*∂²𝒒∂ξ²
            ∂𝒒̂ᵀ∂x𝗚⁻¹∂²𝒒∂ξ² =  ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂²𝒒∂ξ²*n₁₁+∂𝒒̂ᵀ∂η𝗚⁻¹∂²𝒒∂ξ²*n₂₁
            ∂𝒒̂ᵀ∂y𝗚⁻¹∂²𝒒∂ξ² =  ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂²𝒒∂ξ²*n₁₂+∂𝒒̂ᵀ∂η𝗚⁻¹∂²𝒒∂ξ²*n₂₂

            𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ∂η = 𝒒̂ᵀ𝗚⁻¹*∂²𝒒∂ξ∂η
            ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂²𝒒∂ξ∂η = ∂𝒒̂ᵀ∂ξ𝗚⁻¹*∂²𝒒∂ξ∂η
            ∂𝒒̂ᵀ∂η𝗚⁻¹∂²𝒒∂ξ∂η = ∂𝒒̂ᵀ∂η𝗚⁻¹*∂²𝒒∂ξ∂η
            ∂𝒒̂ᵀ∂x𝗚⁻¹∂²𝒒∂ξ∂η =  ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂²𝒒∂ξ∂η*n₁₁+∂𝒒̂ᵀ∂η𝗚⁻¹∂²𝒒∂ξ∂η*n₂₁
            ∂𝒒̂ᵀ∂y𝗚⁻¹∂²𝒒∂ξ∂η =  ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂²𝒒∂ξ∂η*n₁₂+∂𝒒̂ᵀ∂η𝗚⁻¹∂²𝒒∂ξ∂η*n₂₂

            𝒒̂ᵀ𝗚⁻¹∂²𝒒∂η² = 𝒒̂ᵀ𝗚⁻¹*∂²𝒒∂η²
            ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂²𝒒∂η² = ∂𝒒̂ᵀ∂ξ𝗚⁻¹*∂²𝒒∂η²
            ∂𝒒̂ᵀ∂η𝗚⁻¹∂²𝒒∂η² = ∂𝒒̂ᵀ∂η𝗚⁻¹*∂²𝒒∂η²
            ∂𝒒̂ᵀ∂x𝗚⁻¹∂²𝒒∂η² =  ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂²𝒒∂η²*n₁₁+∂𝒒̂ᵀ∂η𝗚⁻¹∂²𝒒∂η²*n₂₁
            ∂𝒒̂ᵀ∂y𝗚⁻¹∂²𝒒∂η² =  ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂²𝒒∂η²*n₁₂+∂𝒒̂ᵀ∂η𝗚⁻¹∂²𝒒∂η²*n₂₂

            q₁₁ = 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ²*n₁₁*n₁₁ + 2*𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ∂η*n₁₁*n₂₁ + 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂η²*n₂₁*n₂₁
            ∂q₁₁∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹∂²𝒒∂ξ²*n₁₁*n₁₁ + 2*∂𝒒̂ᵀ∂x𝗚⁻¹∂²𝒒∂ξ∂η*n₁₁*n₂₁ + ∂𝒒̂ᵀ∂x𝗚⁻¹∂²𝒒∂η²*n₂₁*n₂₁
            ∂q₁₁∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹∂²𝒒∂ξ²*n₁₁*n₁₁ + 2*∂𝒒̂ᵀ∂y𝗚⁻¹∂²𝒒∂ξ∂η*n₁₁*n₂₁ + ∂𝒒̂ᵀ∂y𝗚⁻¹∂²𝒒∂η²*n₂₁*n₂₁

            q₁₂ = 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ²*n₁₁*n₁₂ + 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ∂η*(n₁₁*n₂₂+n₁₂*n₂₁) + 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂η²*n₂₁*n₂₂
            ∂q₁₂∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹∂²𝒒∂ξ²*n₁₁*n₁₂ + ∂𝒒̂ᵀ∂x𝗚⁻¹∂²𝒒∂ξ∂η*(n₁₁*n₂₂+n₁₂*n₂₁) + ∂𝒒̂ᵀ∂x𝗚⁻¹∂²𝒒∂η²*n₂₁*n₂₂
            ∂q₁₂∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹∂²𝒒∂ξ²*n₁₁*n₁₂ + ∂𝒒̂ᵀ∂y𝗚⁻¹∂²𝒒∂ξ∂η*(n₁₁*n₂₂+n₁₂*n₂₁) + ∂𝒒̂ᵀ∂y𝗚⁻¹∂²𝒒∂η²*n₂₁*n₂₂

            q₂₂ = 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ²*n₁₂*n₁₂ + 2*𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ∂η*n₁₂*n₂₂ + 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂η²*n₂₂*n₂₂
            ∂q₂₂∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹∂²𝒒∂ξ²*n₁₂*n₁₂ + 2*∂𝒒̂ᵀ∂x𝗚⁻¹∂²𝒒∂ξ∂η*n₁₂*n₂₂ + ∂𝒒̂ᵀ∂x𝗚⁻¹∂²𝒒∂η²*n₂₂*n₂₂
            ∂q₂₂∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹∂²𝒒∂ξ²*n₁₂*n₁₂ + 2*∂𝒒̂ᵀ∂y𝗚⁻¹∂²𝒒∂ξ∂η*n₁₂*n₂₂ + ∂𝒒̂ᵀ∂y𝗚⁻¹∂²𝒒∂η²*n₂₂*n₂₂

            q₁ = 𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ*n₁₁ + 𝒒̂ᵀ𝗚⁻¹∂𝒒∂η*n₂₁
            ∂q₁∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂ξ*n₁₁ + ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂η*n₂₁
            ∂q₁∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂ξ*n₁₁ + ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂η*n₂₁
            q₂ = 𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ*n₁₂ + 𝒒̂ᵀ𝗚⁻¹∂𝒒∂η*n₂₂
            ∂q₂∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂ξ*n₁₂ + ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂η*n₂₂
            ∂q₂∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂ξ*n₁₂ + ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂η*n₂₂

            q₁n₁ = 0.0;q₂n₂ = 0.0;q₁n₂ = 0.0;q₂n₁ = 0.0
            ∂q₁∂xn₁ = 0.0;∂q₂∂xn₂ = 0.0;∂q₁∂xn₂ = 0.0;∂q₂∂xn₁ = 0.0
            ∂q₁∂yn₁ = 0.0;∂q₂∂yn₂ = 0.0;∂q₁∂yn₂ = 0.0;∂q₂∂yn₁ = 0.0
            mn₁₁n₁ = 0.0;mn₁₁n₂ = 0.0;mn₁₂n₁ = 0.0;mn₁₂n₂ = 0.0;mn₂₂n₁ = 0.0;mn₂₂n₂ = 0.0
            ∂mn₁₁∂xn₁ = 0.0;∂mn₁₁∂xn₂ = 0.0;∂mn₁₂∂xn₁ = 0.0;∂mn₁₂∂xn₂ = 0.0;∂mn₂₂∂xn₁ = 0.0;∂mn₂₂∂xn₂ = 0.0
            ∂mn₁₁∂yn₁ = 0.0;∂mn₁₁∂yn₂ = 0.0;∂mn₁₂∂yn₁ = 0.0;∂mn₁₂∂yn₂ = 0.0;∂mn₂₂∂yn₁ = 0.0;∂mn₂₂∂yn₂ = 0.0
            ms₁₁ = 0.0;ms₁₂ = 0.0;ms₂₂ = 0.0
            ∂ms₁₁∂x = 0.0;∂ms₁₂∂x = 0.0;∂ms₂₂∂x = 0.0
            ∂ms₁₁∂y = 0.0;∂ms₁₂∂y = 0.0;∂ms₂₂∂y = 0.0
            Δms₁₁ = 0.0;Δms₁₂ = 0.0;Δms₂₂ = 0.0
            Δ∂ms₁₁∂x = 0.0;Δ∂ms₁₂∂x = 0.0;Δ∂ms₂₂∂x = 0.0
            Δ∂ms₁₁∂y = 0.0;Δ∂ms₁₂∂y = 0.0;Δ∂ms₂₂∂y = 0.0
            𝐿₁² = n₁₁^2+n₁₂^2
            𝐿₂² = n₂₁^2+n₂₂^2
            𝐿₃² = n₃₁^2+n₃₂^2
            if ξ.ξ == 0.0
                q₁n₁ += q₁*n₁₁
                ∂q₁∂xn₁ += ∂q₁∂x*n₁₁
                ∂q₁∂yn₁ += ∂q₁∂y*n₁₁
                q₁n₂ += q₁*n₁₂
                ∂q₁∂xn₂ += ∂q₁∂x*n₁₂
                ∂q₁∂yn₂ += ∂q₁∂y*n₁₂
                q₂n₁ += q₂*n₁₁
                ∂q₂∂xn₁ += ∂q₂∂x*n₁₁
                ∂q₂∂yn₁ += ∂q₂∂y*n₁₁
                q₂n₂ += q₂*n₁₂
                ∂q₂∂xn₂ += ∂q₂∂x*n₁₂
                ∂q₂∂yn₂ += ∂q₂∂y*n₁₂
                mn₁₁n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₁₁*n₁₁*n₁₁/𝐿₁²
                ∂mn₁₁∂xn₁ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₁₁*n₁₁*n₁₁/𝐿₁²
                ∂mn₁₁∂yn₁ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₁₁*n₁₁*n₁₁/𝐿₁²
                mn₁₁n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₁₁*n₁₁*n₁₂/𝐿₁²
                ∂mn₁₁∂xn₂ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₁₁*n₁₁*n₁₂/𝐿₁²
                ∂mn₁₁∂yn₂ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₁₁*n₁₁*n₁₂/𝐿₁²
                mn₁₂n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₁₁*n₁₂*n₁₁/𝐿₁²
                ∂mn₁₂∂xn₁ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₁₁*n₁₂*n₁₁/𝐿₁²
                ∂mn₁₂∂yn₁ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₁₁*n₁₂*n₁₁/𝐿₁²
                mn₁₂n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₁₁*n₁₂*n₁₂/𝐿₁²
                ∂mn₁₂∂xn₂ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₁₁*n₁₂*n₁₂/𝐿₁²
                ∂mn₁₂∂yn₂ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₁₁*n₁₂*n₁₂/𝐿₁²
                mn₂₂n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₁₂*n₁₂*n₁₁/𝐿₁²
                ∂mn₂₂∂xn₁ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₁₂*n₁₂*n₁₁/𝐿₁²
                ∂mn₂₂∂yn₁ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₁₂*n₁₂*n₁₁/𝐿₁²
                mn₂₂n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₁₂*n₁₂*n₁₂/𝐿₁²
                ∂mn₂₂∂xn₂ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₁₂*n₁₂*n₁₂/𝐿₁²
                ∂mn₂₂∂yn₂ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₁₂*n₁₂*n₁₂/𝐿₁²
                ms₁₁ += (q₁*s₁₁+q₂*s₁₂)*n₁₁*s₁₁/𝐿₁²
                ∂ms₁₁∂x += (∂q₁∂x*s₁₁+∂q₂∂x*s₁₂)*n₁₁*s₁₁/𝐿₁²
                ∂ms₁₁∂y += (∂q₁∂y*s₁₁+∂q₂∂y*s₁₂)*n₁₁*s₁₁/𝐿₁²
                ms₁₂ += (q₁*s₁₁+q₂*s₁₂)*0.5*(n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²
                ∂ms₁₂∂x += (∂q₁∂x*s₁₁+∂q₂∂x*s₁₂)*0.5*(n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²
                ∂ms₁₂∂y += (∂q₁∂y*s₁₁+∂q₂∂y*s₁₂)*0.5*(n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²
                ms₂₂ += (q₁*s₁₁+q₂*s₁₂)*n₁₂*s₁₂/𝐿₁²
                ∂ms₂₂∂x += (∂q₁∂x*s₁₁+∂q₂∂x*s₁₂)*n₁₂*s₁₂/𝐿₁²
                ∂ms₂₂∂y += (∂q₁∂y*s₁₁+∂q₂∂y*s₁₂)*n₁₂*s₁₂/𝐿₁²
                if ξ.η == 0.0
                    Δms₁₁ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₁₁*s₁₁/𝐿₁²-n₂₁*s₂₁/𝐿₂²)
                    Δ∂ms₁₁∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*(n₁₁*s₁₁/𝐿₁²-n₂₁*s₂₁/𝐿₂²)
                    Δ∂ms₁₁∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*(n₁₁*s₁₁/𝐿₁²-n₂₁*s₂₁/𝐿₂²)
                    Δms₂₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₁₂*s₁₂/𝐿₁²-n₂₂*s₂₂/𝐿₂²)
                    Δ∂ms₂₂∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*(n₁₂*s₁₂/𝐿₁²-n₂₂*s₂₂/𝐿₂²)
                    Δ∂ms₂₂∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*(n₁₂*s₁₂/𝐿₁²-n₂₂*s₂₂/𝐿₂²)
                    Δms₁₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*0.5*((n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²-(n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²)
                    Δ∂ms₁₂∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*0.5*((n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²-(n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²)
                    Δ∂ms₁₂∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*0.5*((n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²-(n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²)
                end
            end
            if  ξ.η == 0.0
                q₁n₁ += q₁*n₂₁
                ∂q₁∂xn₁ += ∂q₁∂x*n₂₁
                ∂q₁∂yn₁ += ∂q₁∂y*n₂₁
                q₁n₂ += q₁*n₂₂
                ∂q₁∂xn₂ += ∂q₁∂x*n₂₂
                ∂q₁∂yn₂ += ∂q₁∂y*n₂₂
                q₂n₁ += q₂*n₂₁
                ∂q₂∂xn₁ += ∂q₂∂x*n₂₁
                ∂q₂∂yn₁ += ∂q₂∂y*n₂₁
                q₂n₂ += q₂*n₂₂
                ∂q₂∂xn₂ += ∂q₂∂x*n₂₂
                ∂q₂∂yn₂ += ∂q₂∂y*n₂₂
                mn₁₁n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₂₁*n₂₁*n₂₁/𝐿₂²
                ∂mn₁₁∂xn₁ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₂₁*n₂₁*n₂₁/𝐿₂²
                ∂mn₁₁∂yn₁ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₂₁*n₂₁*n₂₁/𝐿₂²
                mn₁₁n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₂₁*n₂₁*n₂₂/𝐿₂²
                ∂mn₁₁∂xn₂ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₂₁*n₂₁*n₂₂/𝐿₂²
                ∂mn₁₁∂yn₂ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₂₁*n₂₁*n₂₂/𝐿₂²
                mn₁₂n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₂₁*n₂₂*n₂₁/𝐿₂²
                ∂mn₁₂∂xn₁ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₂₁*n₂₂*n₂₁/𝐿₂²
                ∂mn₁₂∂yn₁ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₂₁*n₂₂*n₂₁/𝐿₂²
                mn₁₂n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₂₁*n₂₂*n₂₂/𝐿₂²
                ∂mn₁₂∂xn₂ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₂₁*n₂₂*n₂₂/𝐿₂²
                ∂mn₁₂∂yn₂ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₂₁*n₂₂*n₂₂/𝐿₂²
                mn₂₂n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₂₂*n₂₂*n₂₁/𝐿₂²
                ∂mn₂₂∂xn₁ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₂₂*n₂₂*n₂₁/𝐿₂²
                ∂mn₂₂∂yn₁ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₂₂*n₂₂*n₂₁/𝐿₂²
                mn₂₂n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₂₂*n₂₂*n₂₂/𝐿₂²
                ∂mn₂₂∂xn₂ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₂₂*n₂₂*n₂₂/𝐿₂²
                ∂mn₂₂∂yn₂ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₂₂*n₂₂*n₂₂/𝐿₂²
                ms₁₁ += (q₁*s₂₁+q₂*s₂₂)*n₂₁*s₂₁/𝐿₂²
                ∂ms₁₁∂x += (∂q₁∂x*s₂₁+∂q₂∂x*s₂₂)*n₂₁*s₂₁/𝐿₂²
                ∂ms₁₁∂y += (∂q₁∂y*s₂₁+∂q₂∂y*s₂₂)*n₂₁*s₂₁/𝐿₂²
                ms₁₂ += (q₁*s₂₁+q₂*s₂₂)*0.5*(n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²
                ∂ms₁₂∂x += (∂q₁∂x*s₂₁+∂q₂∂x*s₂₂)*0.5*(n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²
                ∂ms₁₂∂y += (∂q₁∂y*s₂₁+∂q₂∂y*s₂₂)*0.5*(n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²
                ms₂₂ += (q₁*s₂₁+q₂*s₂₂)*n₂₂*s₂₂/𝐿₂²
                ∂ms₂₂∂x += (∂q₁∂x*s₂₁+∂q₂∂x*s₂₂)*n₂₂*s₂₂/𝐿₂²
                ∂ms₂₂∂y += (∂q₁∂y*s₂₁+∂q₂∂y*s₂₂)*n₂₂*s₂₂/𝐿₂²
                if ξ.ξ+ξ.η ≈ 1.0
                    Δms₁₁ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₂₁*s₂₁/𝐿₂²-n₃₁*s₃₁/𝐿₃²)
                    Δ∂ms₁₁∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*(n₂₁*s₂₁/𝐿₂²-n₃₁*s₃₁/𝐿₃²)
                    Δ∂ms₁₁∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*(n₂₁*s₂₁/𝐿₂²-n₃₁*s₃₁/𝐿₃²)
                    Δms₂₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₂₂*s₂₂/𝐿₂²-n₃₂*s₃₂/𝐿₃²)
                    Δ∂ms₂₂∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*(n₂₂*s₂₂/𝐿₂²-n₃₂*s₃₂/𝐿₃²)
                    Δ∂ms₂₂∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*(n₂₂*s₂₂/𝐿₂²-n₃₂*s₃₂/𝐿₃²)
                    Δms₁₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*0.5*((n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²-(n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²)
                    Δ∂ms₁₂∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*0.5*((n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²-(n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²)
                    Δ∂ms₁₂∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*0.5*((n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²-(n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²)
                end
            end
            if ξ.ξ+ξ.η ≈ 1.0
                q₁n₁ += q₁*n₃₁
                ∂q₁∂xn₁ += ∂q₁∂x*n₃₁
                ∂q₁∂yn₁ += ∂q₁∂y*n₃₁
                q₁n₂ += q₁*n₃₂
                ∂q₁∂xn₂ += ∂q₁∂x*n₃₂
                ∂q₁∂yn₂ += ∂q₁∂y*n₃₂
                q₂n₁ += q₂*n₃₁
                ∂q₂∂xn₁ += ∂q₂∂x*n₃₁
                ∂q₂∂yn₁ += ∂q₂∂y*n₃₁
                q₂n₂ += q₂*n₃₂
                ∂q₂∂xn₂ += ∂q₂∂x*n₃₂
                ∂q₂∂yn₂ += ∂q₂∂y*n₃₂
                mn₁₁n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₃₁*n₃₁*n₃₁/𝐿₃²
                ∂mn₁₁∂xn₁ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₃₁*n₃₁*n₃₁/𝐿₃²
                ∂mn₁₁∂yn₁ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₃₁*n₃₁*n₃₁/𝐿₃²
                mn₁₁n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₃₁*n₃₁*n₃₂/𝐿₃²
                ∂mn₁₁∂xn₂ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₃₁*n₃₁*n₃₂/𝐿₃²
                ∂mn₁₁∂yn₂ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₃₁*n₃₁*n₃₂/𝐿₃²
                mn₁₂n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₃₁*n₃₂*n₃₁/𝐿₃²
                ∂mn₁₂∂xn₁ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₃₁*n₃₂*n₃₁/𝐿₃²
                ∂mn₁₂∂yn₁ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₃₁*n₃₂*n₃₁/𝐿₃²
                mn₁₂n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₃₁*n₃₂*n₃₂/𝐿₃²
                ∂mn₁₂∂xn₂ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₃₁*n₃₂*n₃₂/𝐿₃²
                ∂mn₁₂∂yn₂ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₃₁*n₃₂*n₃₂/𝐿₃²
                mn₂₂n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₃₂*n₃₂*n₃₁/𝐿₃²
                ∂mn₂₂∂xn₁ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₃₂*n₃₂*n₃₁/𝐿₃²
                ∂mn₂₂∂yn₁ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₃₂*n₃₂*n₃₁/𝐿₃²
                mn₂₂n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₃₂*n₃₂*n₃₂/𝐿₃²
                ∂mn₂₂∂xn₂ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₃₂*n₃₂*n₃₂/𝐿₃²
                ∂mn₂₂∂yn₂ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₃₂*n₃₂*n₃₂/𝐿₃²
                ms₁₁ += (q₁*s₃₁+q₂*s₃₂)*n₃₁*s₃₁/𝐿₃²
                ∂ms₁₁∂x += (∂q₁∂x*s₃₁+∂q₂∂x*s₃₂)*n₃₁*s₃₁/𝐿₃²
                ∂ms₁₁∂y += (∂q₁∂y*s₃₁+∂q₂∂y*s₃₂)*n₃₁*s₃₁/𝐿₃²
                ms₁₂ += (q₁*s₃₁+q₂*s₃₂)*0.5*(n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²
                ∂ms₁₂∂x += (∂q₁∂x*s₃₁+∂q₂∂x*s₃₂)*0.5*(n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²
                ∂ms₁₂∂y += (∂q₁∂y*s₃₁+∂q₂∂y*s₃₂)*0.5*(n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²
                ms₂₂ += (q₁*s₃₁+q₂*s₃₂)*n₃₂*s₃₂/𝐿₃²
                ∂ms₂₂∂x += (∂q₁∂x*s₃₁+∂q₂∂x*s₃₂)*n₃₂*s₃₂/𝐿₃²
                ∂ms₂₂∂y += (∂q₁∂y*s₃₁+∂q₂∂y*s₃₂)*n₃₂*s₃₂/𝐿₃²
                if ξ.ξ == 0.0
                    Δms₁₁ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₃₁*s₃₁/𝐿₃²-n₁₁*s₁₁/𝐿₃²)
                    Δ∂ms₁₁∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*(n₃₁*s₃₁/𝐿₃²-n₁₁*s₁₁/𝐿₃²)
                    Δ∂ms₁₁∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*(n₃₁*s₃₁/𝐿₃²-n₁₁*s₁₁/𝐿₃²)
                    Δms₂₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₃₂*s₃₂/𝐿₃²-n₁₂*s₁₂/𝐿₁²)
                    Δ∂ms₂₂∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*(n₃₂*s₃₂/𝐿₃²-n₁₂*s₁₂/𝐿₁²)
                    Δ∂ms₂₂∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*(n₃₂*s₃₂/𝐿₃²-n₁₂*s₁₂/𝐿₁²)
                    Δms₁₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*0.5*((n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²-(n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²)
                    Δ∂ms₁₂∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*0.5*((n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²-(n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²)
                    Δ∂ms₁₂∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*0.5*((n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²-(n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²)
                end
            end

            W₁₁₁ = mn₁₁n₁*wᵇ
            ∂W₁₁₁∂x = ∂mn₁₁∂xn₁*wᵇ
            ∂W₁₁₁∂y = ∂mn₁₁∂yn₁*wᵇ
            W₁₁₂ = mn₁₁n₂*wᵇ
            ∂W₁₁₂∂x = ∂mn₁₁∂xn₂*wᵇ
            ∂W₁₁₂∂y = ∂mn₁₁∂yn₂*wᵇ
            W₁₂₁ = mn₁₂n₁*wᵇ
            ∂W₁₂₁∂x = ∂mn₁₂∂xn₁*wᵇ
            ∂W₁₂₁∂y = ∂mn₁₂∂yn₁*wᵇ
            W₁₂₂ = mn₁₂n₂*wᵇ
            ∂W₁₂₂∂x = ∂mn₁₂∂xn₂*wᵇ
            ∂W₁₂₂∂y = ∂mn₁₂∂yn₂*wᵇ
            W₂₂₁ = mn₂₂n₁*wᵇ
            ∂W₂₂₁∂x = ∂mn₂₂∂xn₁*wᵇ
            ∂W₂₂₁∂y = ∂mn₂₂∂yn₁*wᵇ
            W₂₂₂ = mn₂₂n₂*wᵇ
            ∂W₂₂₂∂x = ∂mn₂₂∂xn₂*wᵇ
            ∂W₂₂₂∂y = ∂mn₂₂∂yn₂*wᵇ
            W₁₁ = (q₁₁*w + 2*(q₁n₁+ms₁₁)*wᵇ)/4/𝐴 + Δms₁₁
            ∂W₁₁∂x = (∂q₁₁∂x*w + 2*(∂q₁∂xn₁+∂ms₁₁∂x)*wᵇ)/4/𝐴 + Δ∂ms₁₁∂x
            ∂W₁₁∂y = (∂q₁₁∂y*w + 2*(∂q₁∂yn₁+∂ms₁₁∂y)*wᵇ)/4/𝐴 + Δ∂ms₁₁∂y
            W₁₂ = (q₁₂*w + (q₁n₂+q₂n₁+2*ms₁₂)*wᵇ)/4/𝐴 + Δms₁₂
            ∂W₁₂∂x = (∂q₁₂∂x*w + (∂q₁∂xn₂+∂q₂∂xn₁+2*∂ms₁₂∂x)*wᵇ)/4/𝐴 + Δ∂ms₁₂∂x
            ∂W₁₂∂y = (∂q₁₂∂y*w + (∂q₁∂yn₂+∂q₂∂yn₁+2*∂ms₁₂∂y)*wᵇ)/4/𝐴 + Δ∂ms₁₂∂y
            W₂₂ = (q₂₂*w + 2*(q₂n₂+ms₂₂)*wᵇ)/4/𝐴 + Δms₂₂
            ∂W₂₂∂x = (∂q₂₂∂x*w + 2*(∂q₂∂xn₂+∂ms₂₂∂x)*wᵇ)/4/𝐴 + Δ∂ms₂₂∂x
            ∂W₂₂∂y = (∂q₂₂∂y*w + 2*(∂q₂∂yn₂+∂ms₂₂∂y)*wᵇ)/4/𝐴 + Δ∂ms₂₂∂y
            for i in 1:length(𝓒)
                ∂²𝝭∂x²[i] += 𝝭[i]*W₁₁ + ∂𝝭∂x[i]*W₁₁₁ + ∂𝝭∂y[i]*W₁₁₂
                ∂∂²𝝭∂x²∂x[i] += 𝝭[i]*∂W₁₁∂x + ∂𝝭∂x[i]*∂W₁₁₁∂x + ∂𝝭∂y[i]*∂W₁₁₂∂x
                ∂∂²𝝭∂x²∂y[i] += 𝝭[i]*∂W₁₁∂y + ∂𝝭∂x[i]*∂W₁₁₁∂y + ∂𝝭∂y[i]*∂W₁₁₂∂y
                ∂²𝝭∂x∂y[i] += 𝝭[i]*W₁₂ + ∂𝝭∂x[i]*W₁₂₁ + ∂𝝭∂y[i]*W₁₂₂
                ∂∂²𝝭∂x∂y∂x[i] += 𝝭[i]*∂W₁₂∂x + ∂𝝭∂x[i]*∂W₁₂₁∂x + ∂𝝭∂y[i]*∂W₁₂₂∂x
                ∂∂²𝝭∂x∂y∂y[i] += 𝝭[i]*∂W₁₂∂y + ∂𝝭∂x[i]*∂W₁₂₁∂y + ∂𝝭∂y[i]*∂W₁₂₂∂y
                ∂²𝝭∂y²[i] += 𝝭[i]*W₂₂ + ∂𝝭∂x[i]*W₂₂₁ + ∂𝝭∂y[i]*W₂₂₂
                ∂∂²𝝭∂y²∂x[i] += 𝝭[i]*∂W₂₂∂x + ∂𝝭∂x[i]*∂W₂₂₁∂x + ∂𝝭∂y[i]*∂W₂₂₂∂x
                ∂∂²𝝭∂y²∂y[i] += 𝝭[i]*∂W₂₂∂y + ∂𝝭∂x[i]*∂W₂₂₁∂y + ∂𝝭∂y[i]*∂W₂₂₂∂y
            end
        end
        for i in 1:length(𝓒)
            ξ̂.𝝭[:∂x²][ξ̂.index[ξ̂.id]+i] = ∂²𝝭∂x²[i]
            ξ̂.𝝭[:∂x²∂x][ξ̂.index[ξ̂.id]+i] = ∂∂²𝝭∂x²∂x[i]
            ξ̂.𝝭[:∂x²∂y][ξ̂.index[ξ̂.id]+i] = ∂∂²𝝭∂x²∂y[i]
            ξ̂.𝝭[:∂x∂y][ξ̂.index[ξ̂.id]+i] = ∂²𝝭∂x∂y[i]
            ξ̂.𝝭[:∂x∂y∂x][ξ̂.index[ξ̂.id]+i] = ∂∂²𝝭∂x∂y∂x[i]
            ξ̂.𝝭[:∂x∂y∂y][ξ̂.index[ξ̂.id]+i] = ∂∂²𝝭∂x∂y∂y[i]
            ξ̂.𝝭[:∂y²][ξ̂.index[ξ̂.id]+i] = ∂²𝝭∂y²[i]
            ξ̂.𝝭[:∂y²∂x][ξ̂.index[ξ̂.id]+i] = ∂∂²𝝭∂y²∂x[i]
            ξ̂.𝝭[:∂y²∂y][ξ̂.index[ξ̂.id]+i] = ∂∂²𝝭∂y²∂y[i]
        end
    end
end

function set∇̄𝝭!(aps::Vector{T}) where T<:AbstractElement
    set𝝭!(aps)
    for ap in aps
        set∇̄𝝭!(ap)
    end
end

function set∇̄𝝭!(ap::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Seg2}) where {𝒑,𝑠,𝜙}
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ̂ in 𝓖
        𝒒̂ = get𝒑₁(ap,ξ̂)
        𝗚⁻¹ = cal𝗚!(ap)
        𝒒̂ᵀ𝗚⁻¹ = 𝒒̂*𝗚⁻¹
        ∂𝝭∂x = ap.𝝭[:∂x]
        fill!(∂𝝭∂x,0.0)
        for ξ in ap.𝓖
            w = ξ.w
            n = ξ.n₁
            𝝭 = get𝝭(ap,ξ)
            𝒒 = get𝒑₁(ap,ξ)
            W₁ = 𝒒̂ᵀ𝗚⁻¹*𝒒*n*w
            for i in 1:length(𝓒)
                ∂𝝭∂x[i] += 𝝭[i]*W₁
            end
        end
        for i in 1:length(𝓒)
            ξ̂.𝝭[:∂̄x][ξ̂.index[ξ̂.id]+i] = ∂𝝭∂x[i]
        end
    end
end

function set∇̄𝝭!(ap::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Tri3}) where {𝒑,𝑠,𝜙}
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ̂ in 𝓖
        𝒒̂ = get𝒑₁(ap,ξ̂)
        𝗚⁻¹ = cal𝗚!(ap)
        𝒒̂ᵀ𝗚⁻¹ = 𝒒̂*𝗚⁻¹
        ∂𝝭∂x = ap.𝝭[:∂x]
        ∂𝝭∂y = ap.𝝭[:∂y]
        fill!(∂𝝭∂x,0.0)
        fill!(∂𝝭∂y,0.0)
        for ξ in ap.𝓖
            w = ξ.w
            n₁ = ξ.n₁
            n₂ = ξ.n₂
            𝝭 = get𝝭(ap,ξ)
            𝒒 = get𝒑₁(ap,ξ)
            𝒒̂ᵀ𝗚⁻¹𝒒 = 𝒒̂ᵀ𝗚⁻¹*𝒒
            W₁ = 𝒒̂ᵀ𝗚⁻¹𝒒*n₁*w
            W₂ = 𝒒̂ᵀ𝗚⁻¹𝒒*n₂*w
            for i in 1:length(𝓒)
                ∂𝝭∂x[i] += 𝝭[i]*W₁
                ∂𝝭∂y[i] += 𝝭[i]*W₂
            end
        end
        for i in 1:length(𝓒)
            ξ̂.𝝭[:∂̄x][ξ̂.index[ξ̂.id]+i] = ∂𝝭∂x[i]
            ξ̂.𝝭[:∂̄y][ξ̂.index[ξ̂.id]+i] = ∂𝝭∂y[i]
        end
    end
end

function set∇̄²𝝭!(ap::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Tri3}) where {𝒑,𝑠,𝜙}
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ̂ in 𝓖
        𝒒̂ = get𝒑₂(ap,ξ̂)
        𝗚⁻¹ = cal𝗚₂!(ap)
        𝒒̂ᵀ𝗚⁻¹ = 𝒒̂*𝗚⁻¹
        M₁₁ = ap.𝝭[:M₁₁]
        M₁₂ = ap.𝝭[:M₁₂]
        M₂₂ = ap.𝝭[:M₂₂]
        fill!(M₁₁,0.0)
        fill!(M₁₂,0.0)
        fill!(M₂₂,0.0)
        for ξ in 𝓖
            w = ξ.w
            n₁ = ξ.n₁
            n₂ = ξ.n₂
            𝐿² = ξ.𝐿^2
            _,∂𝝭∂x,∂𝝭∂y = get∇𝝭(ap,ξ)
            𝒒 = get𝒑₂(ap,ξ)
            𝒒̂ᵀ𝗚⁻¹𝒒 = 𝒒̂ᵀ𝗚⁻¹*𝒒
            W₁₁₁ = 𝒒̂ᵀ𝗚⁻¹𝒒*n₁*n₁*n₁*w/𝐿²
            W₁₁₂ = 𝒒̂ᵀ𝗚⁻¹𝒒*n₁*n₁*n₂*w/𝐿²
            W₂₂₁ = 𝒒̂ᵀ𝗚⁻¹𝒒*n₂*n₂*n₁*w/𝐿²
            W₂₂₂ = 𝒒̂ᵀ𝗚⁻¹𝒒*n₂*n₂*n₂*w/𝐿²
            W₁₂₁ = 𝒒̂ᵀ𝗚⁻¹𝒒*n₁*n₂*n₁*w/𝐿²
            W₁₂₂ = 𝒒̂ᵀ𝗚⁻¹𝒒*n₁*n₂*n₂*w/𝐿²
            for i in 1:length(𝓒)
                M₁₁[i] += ∂𝝭∂x[i]*W₁₁₁ + ∂𝝭∂y[i]*W₁₁₂
                M₁₂[i] += ∂𝝭∂x[i]*W₁₂₁ + ∂𝝭∂y[i]*W₁₂₂
                M₂₂[i] += ∂𝝭∂x[i]*W₂₂₁ + ∂𝝭∂y[i]*W₂₂₂
            end
        end
        n₁ = ξ̂.n₁
        n₂ = ξ̂.n₂
        for i in 1:length(𝓒)
            ξ̂.𝝭[:Mₙₙ][ξ̂.index[ξ̂.id]+i] = M₁₁[i]*n₁*n₁ + 2*M₁₂[i]*n₁*n₂ + M₂₂[i]*n₂*n₂
        end
    end
end

function setV̄ₙ!(ap::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Tri3}) where {𝒑,𝑠,𝜙}
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    x₁ = 𝓒[1].x;y₁ = 𝓒[1].y
    x₂ = 𝓒[2].x;y₂ = 𝓒[2].y
    x₃ = 𝓒[3].x;y₃ = 𝓒[3].y
    𝐴 = get𝐴(gp)
    n₁₁ = y₃-y₂;n₂₁ = y₁-y₃;n₃₁ = y₂-y₁
    n₁₂ = x₂-x₃;n₂₂ = x₃-x₁;n₃₂ = x₁-x₂

    for ξ̂ in 𝓖
        _,∂𝒒̂∂ξ,∂𝒒̂∂η = get∇𝒑₂(ap,ξ̂)
        𝗚⁻¹ = cal𝗚₂!(ap)
        ∂𝒒̂ᵀ∂ξ𝗚⁻¹ = ∂𝒒̂∂ξ*𝗚⁻¹
        ∂𝒒̂ᵀ∂η𝗚⁻¹ = ∂𝒒̂∂η*𝗚⁻¹
        V₁₁₁ = ap.𝝭[:V₁₁₁]
        V₁₁₂ = ap.𝝭[:V₁₁₂]
        V₁₂₁ = ap.𝝭[:V₁₂₁]
        V₁₂₂ = ap.𝝭[:V₁₂₂]
        V₂₁₁ = ap.𝝭[:V₂₁₁]
        V₂₁₂ = ap.𝝭[:V₂₁₂]
        V₂₂₁ = ap.𝝭[:V₂₂₁]
        V₂₂₂ = ap.𝝭[:V₂₂₂]
        fill!(V₁₁₁,0.0)
        fill!(V₁₁₂,0.0)
        fill!(V₁₂₁,0.0)
        fill!(V₁₂₂,0.0)
        fill!(V₂₁₁,0.0)
        fill!(V₂₁₂,0.0)
        fill!(V₂₂₁,0.0)
        fill!(V₂₂₂,0.0)
        for ξ in ap.𝓖
            w = ξ.w
            n₁ = ξ.n₁
            n₂ = ξ.n₂
            s₁ = -ξ.n₂
            s₂ = ξ.n₁
            𝐿² = ξ.𝐿^2
            𝝭 = get𝝭(ap,ξ)
            _,∂𝒒∂ξ,∂𝒒∂η = get∇𝒑₂(ap,ξ)
            ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂ξ = ∂𝒒̂ᵀ∂ξ𝗚⁻¹*∂𝒒∂ξ
            ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂η = ∂𝒒̂ᵀ∂ξ𝗚⁻¹*∂𝒒∂η
            ∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂ξ = ∂𝒒̂ᵀ∂η𝗚⁻¹*∂𝒒∂ξ
            ∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂η = ∂𝒒̂ᵀ∂η𝗚⁻¹*∂𝒒∂η
            ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂x = ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂ξ*n₁₁*n₁₁ + (∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂η + ∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂ξ)*n₁₁*n₂₁ + ∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂η*n₂₁*n₂₁
            ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂y = ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂ξ*n₁₁*n₁₂ + ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂η*n₁₁*n₂₂ + ∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂ξ*n₁₂*n₂₁ + ∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂η*n₂₁*n₂₂
            ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂x = ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂ξ*n₁₁*n₁₂ + ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂η*n₁₂*n₂₁ + ∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂ξ*n₁₁*n₂₂ + ∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂η*n₂₁*n₂₂
            ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂y = ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂ξ*n₁₂*n₁₂ + (∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂η + ∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂ξ)*n₁₂*n₂₂ + ∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂η*n₂₂*n₂₂
            W₁₁₁ = ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂x*n₁ + (∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂x*s₁ + ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂y*s₂)*s₁*n₁*w/𝐿²
            W₁₁₂ = ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂x*n₁ + (∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂x*s₁ + ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂y*s₂)*s₁*n₁*w/𝐿²
            W₁₂₁ = ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂x*n₂ + (∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂x*s₁ + ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂y*s₂)*s₁*n₂*w/𝐿²
            W₁₂₂ = ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂x*n₂ + (∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂x*s₁ + ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂y*s₂)*s₁*n₂*w/𝐿²
            W₂₁₁ = ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂y*n₁ + (∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂x*s₁ + ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂y*s₂)*s₂*n₁*w/𝐿²
            W₂₁₂ = ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂y*n₁ + (∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂x*s₁ + ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂y*s₂)*s₂*n₁*w/𝐿²
            W₂₂₁ = ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂y*n₂ + (∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂x*s₁ + ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂y*s₂)*s₂*n₂*w/𝐿²
            W₂₂₂ = ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂y*n₂ + (∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂x*s₁ + ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂y*s₂)*s₂*n₂*w/𝐿²
            for i in 1:length(𝓒)
                V₁₁₁[i] += 𝝭[i]*W₁₁₁
                V₁₁₂[i] += 𝝭[i]*W₁₁₂
                V₁₂₁[i] += 𝝭[i]*W₁₂₁
                V₁₂₂[i] += 𝝭[i]*W₁₂₂
                V₂₁₁[i] += 𝝭[i]*W₂₁₁
                V₂₁₂[i] += 𝝭[i]*W₂₁₂
                V₂₂₁[i] += 𝝭[i]*W₂₂₁
                V₂₂₂[i] += 𝝭[i]*W₂₂₂
            end
        end
        n₁ = ξ̂.n₁
        n₂ = ξ̂.n₂
        s₁ = -ξ̂.n₂
        s₂ = ξ̂.n₁
        for i in 1:length(𝓒)
            ξ̂.𝝭[:Vₙ][ξ̂.index[ξ̂.id]+i] = V₁₁₁[i]*n₁ + V₁₂₂[i]*n₁ + V₂₁₁[i]*n₂ + V₂₂₂[i]*n₂ + (V₁₁₁*s₁+V₁₁₂*s₂)*s₁*n₁ + (V₁₂₁*s₁+V₁₂₂*s₂)*s₁*n₂ + (V₂₁₁*s₁+V₂₁₂*s₂)*s₂*n₁ + (V₂₂₁*s₁+V₂₂₂*s₂)*s₂*n₂
        end
    end
end
