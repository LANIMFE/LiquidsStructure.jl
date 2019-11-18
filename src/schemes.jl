### Types
abstract type ApproximationScheme end

struct RosenfeldFMT{T} <: ApproximationScheme
    A::T
    B::T
    G::T
end

struct PercusYevick{T} <: ApproximationScheme
    α::T
    β::T
    δ::T
end

struct VerletWeis{T} <: ApproximationScheme
    coreliquid::HardSpheres{T}
    subscheme::PercusYevick{T}
    α::T
end

struct SharmaSharma{S, T} <: ApproximationScheme
    coreliquid::HardSpheres{T}
    subscheme::S
end

struct MSA{S, T} <: ApproximationScheme
    coreliquid::HardSpheres{T}
    subscheme::S
    α₀::T
    α₁::T
    α₂::T
    α₃::T
    β₀::T
    β₁::T
    β₂::T
    β₃::T
    a₁::T
    a₃::T
    b₁::T
    b₃::T
end


### Constructors
function RosenfeldFMT(liquid::HardDisks)
    η  = volume_fraction(liquid)
    η² = η * η
    η³ = η² * η
    η⁴ = η³ * η
    η₋⁻¹ = 1 / (1 - η)
    η₋⁻² = η₋⁻¹ * η₋⁻¹

    Z  = (1 + 0.128η² + 0.027η³ + 0.06η⁴) * η₋⁻²
    Z′ = ((0.256η + 0.081η² + 0.24η³) * η₋⁻²) + (2Z * η₋⁻¹)
    χ  = Z + η * Z′

    G = sqrt(Z′ / 2)
    A = (1 + (2η - 1) * χ + 2η * G) / η
    B = (-1 + (1 - η) * χ - 3η * G) / η

    return RosenfeldFMT{typeof(A)}(A, B, G)
end

function PercusYevick(liquid::HardSpheres)
    η = volume_fraction(liquid)
    η₀⁴ = @fastmath (1 - η)^4
    η₁² = (1 + 2η)^2
    η₂² = (2 +  η)^2

    α = -η₁² / η₀⁴
    β = 3η * η₂² / 2η₀⁴
    δ = -η * η₁² / 2η₀⁴

    return PercusYevick{typeof(α)}(α, β, δ)
end

function VerletWeis(liquid::HardSpheres)
    η  = volume_fraction(liquid)
    κ  = 1 - η / 16
    η′ = η * κ
    α  = ∛(κ)

    coreliquid = HardSpheres(η′)
    subscheme = PercusYevick(coreliquid)

    return VerletWeis{typeof(α)}(coreliquid, subscheme, α)
end

function SharmaSharma{S}(liquid::AttractiveHardSpheres{TT}) where
    {S<:ApproximationScheme, TT}

    η = volume_fraction(liquid)
    coreliquid = HardSpheres(η)
    subscheme = S(coreliquid)
    SS = typeof(subscheme)

    return SharmaSharma{SS, TT}(coreliquid, subscheme)
end

function MSA{S}(liquid::DipolarHardSpheres{TT}; tol = √(eps(TT))) where
    {S<:ApproximationScheme, TT}

    η = volume_fraction(liquid)
    T = temperature(liquid)
    coreliquid = HardSpheres(η)
    subscheme = S(coreliquid)
    SS = typeof(subscheme)
    κ = dhs_msa_parameter(liquid, tol)

    T⁻¹ = 1 / T
    ξ   = κ * η

    ξ₁² = (1 - 2ξ)^2
    ξ₁⁴ = ξ₁² * ξ₁²
    ξ₂⁴ = (1 +  ξ)^4
    ξ₃² = (1 + 4ξ)^2

    r₁ = ξ₃² / ξ₁⁴
    r₂ = ξ₁² / ξ₂⁴
    r₃ = 8 * (1 + ξ)^2 / ξ₁⁴
    r₄ = (2 - ξ)^2 / ξ₂⁴

    α₀′ = 2κ * (-r₁ + r₂)
    α₁′ = 3ξ * κ * ( r₃ + r₄)
    α₃′ = -ξ * κ * (2r₁ + r₂)
    β₁′ = 3ξ * κ / 8 * (2r₃ - r₄)
    β₃′ = -ξ * κ / 4 * (4r₁ - r₂)

    α₀ = α₀′ / 9  + α₁′ / 12  + α₃′ / 18
    α₁ = α₀′ / 90 + α₁′ / 108 + α₃′ / 144
    α₂ = α₀′ +  α₁′ +  α₃′
    α₃ = α₀′ + 2α₁′ + 4α₃′

    β₀ = -T⁻¹ / 9
    β₁ =  β₁′ / 270 + β₃′ / 360 + β₀ / 10
    β₂ =  β₁′ +  β₃′ - T⁻¹
    β₃ = 4β₁′ + 6β₃′ + β₂

    a₁ = 2 * (α₁′ + 8β₁′) / 3
    a₃ = 8 * (α₃′ + 4β₃′)
    b₁ = 2 * (α₁′ - 4β₁′) / 3
    b₃ = 8 * (α₃′ - 2β₃′)

    return MSA{SS, TT}(coreliquid, subscheme,
                       α₀, α₁, α₂, α₃, β₀, β₁, β₂, β₃, a₁, a₃, b₁, b₃)
end
