abstract type AbstractStructureFactor <: Function end

struct StructureFactor{L, S} <: AbstractStructureFactor
    liquid::L
    scheme::S
end

function StructureFactor(liquid::L, T::Type; kw...) where {L <: Liquid}
    scheme = T(liquid, kw...)
    return StructureFactor{L, typeof(scheme)}(liquid, scheme)
end

(S::StructureFactor)(k) = structure_factor(S.liquid, S.scheme, k)

"""
    structure_factor(liquid, scheme, k)

Returns the static structure factor of a `liquid` using the approximation
defined by `scheme` at the wavevector value `k`.
"""
function structure_factor(liquid, scheme, k)
    ρCk = ρĈ(liquid, scheme, k)
    return 1 / (1 - ρCk)
end

function structure_factor(liquid::DipolarHardSpheres, scheme::MSA, k)
    ρC₀₀k, ρC₁₀k, ρC₁₁k = ρĈ(liquid, scheme, k)

    S₀₀k = 1 / (1 - ρC₀₀k)
    S₁₀k = 1 / (1 - ρC₁₀k)
    S₁₁k = 1 / (1 - ρC₁₁k)

    return (S₀₀k, S₁₀k, S₁₁k)
end

function ρĈ(liquid::HardDisks, scheme::RosenfeldFMT, k)
    η = liquid.η

    A, B, G = scheme.A, scheme.B, scheme.G

    k² = k * k
    J₀ = besselj0(k / 2)
    J₁ = besselj1(k / 2)
    J₁′ = besselj1(k)

    smallk = k < 0.075

    C = smallk ? (A * (1 -  k² / 16) / 4 + B * (1 - 3k² / 32) / 2 +
                  G * (1 -  k² / 8 )) :
                 (A * (2J₁ / k)^2 + B * 2J₀ * J₁ / k + G * 2J₁′ / k)

    return -4η * C
end

"""
    ρĈ(liquid::HardSpheres, scheme::PercusYevick, k)

Returns the product of the bulk density and the fourier transform of the direct
correlation function for a hard-spheres liquid using the Percus-Yevick
approximation.

`k` is the wavevector value
"""
function ρĈ(liquid::HardSpheres, scheme::PercusYevick, k)
    η = liquid.η

    α, β, δ = scheme.α, scheme.β, scheme.δ

    k² = k * k
    k³ = k² * k
    k⁴ = k² * k²
    k⁶ = k³ * k³

    sink, cosk = sin(k), cos(k)

    smallk = k < 0.075

    C₀ = smallk ? (1 - k² / 10) / 3 :
                  (sink - k * cosk) / k³
    C₁ = smallk ? (1 - k² / 9 ) / 4 :
                  (2k * sink - (k² - 2) * cosk - 2) / k⁴
    C₃ = smallk ? (1 - k² / 8 ) / 6 :
                  ((4k³ - 24k) * sink - (k⁴ - 12k² + 24) * cosk + 24) / k⁶

    return 24η * (α * C₀ + β * C₁ + δ * C₃)
end

"""
    ρĈ(::HardSpheres, ::VerletWeis, k)

Returns the product of the bulk density and the fourier transform of the direct
correlation function for a hard-spheres liquid using the Percus-Yevick
approximation with the Verlet-Weis correction.

`k` is the wavevector value
"""
ρĈ(liquid::HardSpheres, scheme::VerletWeis, k) =
    ρĈ(scheme.coreliquid, scheme.subscheme, scheme.α * k)

"""
    ρĈ(liquid::DipolarHardSpheres, scheme::MSA, k)

Returns the product of the bulk density and the fourier transform of
projections of the direct correlation function `(ρC₀₀, ρC₁₀, ρC₁₁)` for a
dipolar hard-spheres liquid using the MSA approximation.

`k` is the wavevector value
"""
function ρĈ(liquid::DipolarHardSpheres, scheme::MSA, k)
    η = liquid.η

    α₀, α₁, α₂, α₃ = scheme.α₀, scheme.α₁, scheme.α₂, scheme.α₃
    β₀, β₁, β₂, β₃ = scheme.β₀, scheme.β₁, scheme.β₂, scheme.β₃
    a₁, a₃, b₁, b₃ = scheme.a₁, scheme.a₃, scheme.b₁, scheme.b₃

    k² = k * k
    k⁶ = k² * k² * k²
    sink, cosk = sin(k), cos(k)

    smallk = k < 0.075

    C₁₀ = smallk ? (α₀ + 2β₀ - (α₁ + 2β₁) * k²) :
                   (  a₃ - a₁ * k² + (-a₃ + (α₃ + 2β₃) * k² / 3) * k * sink +
                    (-a₃ + ((a₁ + a₃ / 2) - (α₂ + 2β₂) * k² / 3) * k²) * cosk
                   ) / k⁶
    C₁₁ = smallk ? (α₀ -  β₀ - (α₁ -  β₁) * k²) :
                   (  b₃ - b₁ * k² + (-b₃ + (α₃ -  β₃) * k² / 3) * k * sink +
                    (-b₃ + ((b₁ + b₃ / 2) - (α₂ -  β₂) * k² / 3) * k²) * cosk
                   ) / k⁶

    ρC₀₀ = ρĈ(scheme.coreliquid, scheme.subscheme, k)
    ρC₁₀ = 24η * C₁₀
    ρC₁₁ = 24η * C₁₁

    return (ρC₀₀, ρC₁₀, ρC₁₁)
end
