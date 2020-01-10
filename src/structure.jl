const SMALL_K = 0.075


abstract type AbstractStructureFactor <: Function end

struct StructureFactor{L, S} <: AbstractStructureFactor
    liquid::L
    scheme::S
end

function StructureFactor(liquid::L, T::Type; kw...) where {L<:AbstractLiquid}
    scheme = T(liquid; kw...)
    return StructureFactor{L, typeof(scheme)}(liquid, scheme)
end

(S::StructureFactor)(k) = structure_factor(S.liquid, S.scheme, k)

"""
    structure_factor(liquid, scheme, k)

Returns the static structure factor of a `liquid` using the approximation
defined by `scheme` at the wavenumber `k`.
"""
function structure_factor(liquid, scheme, k)
    Ck = Ĉ(liquid, scheme, k)
    return 1 / (1 - Ck)
end

function structure_factor(liquid::DipolarHardSpheres, scheme::MSA, k)
    C₀₀k, C₁₀k, C₁₁k = Ĉ(liquid, scheme, k)

    S₀₀k = 1 / (1 - C₀₀k)
    S₁₀k = 1 / (1 - C₁₀k)
    S₁₁k = 1 / (1 - C₁₁k)

    return (S₀₀k, S₁₀k, S₁₁k)
end

function Ĉ(liquid::HardDisks, scheme::RosenfeldFMT, k′)
    η = volume_fraction(liquid)
    σ = length_scale(potential(liquid))

    k  = σ * k′
    k² = k * k
    J₀ = besselj0(k / 2)
    J₁ = besselj1(k / 2)
    J₁′ = besselj1(k)

    A, B, G = scheme.A, scheme.B, scheme.G

    smallk = k < SMALL_K # for small wavenumber use a Taylor series approximation

    C = smallk ? (A * (1 -  k² / 16) / 4 + B * (1 - 3k² / 32) / 2 +
                  G * (1 -  k² / 8 )) :
                 (A * (2J₁ / k)^2 + B * 2J₀ * J₁ / k + G * 2J₁′ / k)

    return -4η * C
end

"""
    Ĉ(liquid::HardSpheres, scheme::PercusYevick, k)

Returns the product of the bulk density and the Fourier transform of the direct
correlation function for a hard-spheres liquid using the Percus-Yevick
approximation.

`k` is the wavenumber.
"""
function Ĉ(liquid::HardSpheres, scheme::PercusYevick, k′)
    η = volume_fraction(liquid)
    σ = length_scale(potential(liquid))

    k  = σ * k′
    k² = k * k
    k³ = k² * k
    k⁴ = k² * k²
    k⁶ = k³ * k³
    sink, cosk = sin(k), cos(k)

    α, β, δ = scheme.α, scheme.β, scheme.δ

    smallk = k < SMALL_K # for small wavenumber use a Taylor series approximation

    C₀ = smallk ? (1 - k² / 10) / 3 :
                  (sink - k * cosk) / k³
    C₁ = smallk ? (1 - k² / 9 ) / 4 :
                  (2k * sink - (k² - 2) * cosk - 2) / k⁴
    C₃ = smallk ? (1 - k² / 8 ) / 6 :
                  ((4k³ - 24k) * sink - (k⁴ - 12k² + 24) * cosk + 24) / k⁶

    return 24η * (α * C₀ + β * C₁ + δ * C₃)
end

"""
    Ĉ(::HardSpheres, ::VerletWeis, k)

Returns the product of the bulk density and the Fourier transform of the direct
correlation function for a hard-spheres liquid using the Percus-Yevick
approximation with the Verlet-Weis correction.

`k` is the wavenumber.
"""
Ĉ(liquid::HardSpheres, scheme::VerletWeis, k) =
    Ĉ(scheme.coreliquid, scheme.subscheme, scheme.α * k)

"""
    Ĉ(::AttractiveHardSpheres{<:YHSPotential}, ::SharmaSharma, k)

Returns the product of the bulk density and the Fourier transform (FT) of the
direct correlation function for an attractive Yukawa hard-spheres liquid using
the Sharma and Sharma approximation, that is, the sum of the FT of the direct
correlation function for the hard spheres core and ``\\beta\\cdot\\hat{u}(k)``,
where ``\\hat{u}(k)`` is the FT of the Yukawa potential:

```math
    u(r) = -\\epsilon\\exp(-zr/\\sigma) / (r/\\sigma)
```

`k` is the wavenumber.
"""
function Ĉ(liquid::AttractiveHardSpheres{P, U}, scheme::SharmaSharma, k′) where
    {P, U<:YHSPotential}

    η = volume_fraction(liquid)
    T = temperature(liquid)
    σ, z = length_scale(potential(liquid))

    k  = σ * k′
    k² = k * k
    k⁴ = k² * k²
    sink, cosk = sin(k), cos(k)

    smallk = k < SMALL_K # for small wavenumber use a Taylor series approximation

    C = smallk ? ((1 + z) - (3 + z) / 6 * k² + (5 + z) / 120 * k⁴) :
                  (k * cosk + z * sink) / k
    C = C / (k² + z * z)

    C₀ = Ĉ(scheme.coreliquid, scheme.subscheme, k′)

    return C₀ + 24η / T * C
end

"""
    Ĉ(::AttractiveHardSpheres{<:SWHSPotential}, ::SharmaSharma, k)

Returns the product of the bulk density and the Fourier transform (FT) of the
direct correlation function for an attractive square-well hard-spheres liquid
using the Sharma and Sharma approximation, that is, the sum of the FT of the
direct correlation function for the hard spheres core and
``\\beta\\cdot\\hat{u}(k)``, where ``\\hat{u}(k)`` is the FT of the square well
potential:

```math
    u(r) = -\\epsilon\\sigma, \\qquad \\sigma < r \\le \\lambda\\sigma
```

`k` is the wavenumber.
"""
function Ĉ(liquid::AttractiveHardSpheres{P, U}, scheme::SharmaSharma, k′) where
    {P, U<:SWHSPotential}

    η = volume_fraction(liquid)
    T = temperature(liquid)
    σ, λ = length_scale(potential(liquid))

    λ³ = λ * λ * λ
    λ⁵ = λ³ * λ * λ

    k  = σ * k′
    k² = k * k
    k³ = k² * k
    sink, cosk = sin(k), cos(k)
    sinλk, cosλk = sin(λ * k), cos(λ * k)

    smallk = k < SMALL_K # for small wavenumber use a Taylor series approximation

    C = smallk ? (λ³ - 1) / 3 - (λ⁵ - 1) / 30 * k² :
                 (cosk - λ * cosλk) / k² + (sinλk - sink) / k³

    C₀ = Ĉ(scheme.coreliquid, scheme.subscheme, k′)

    return C₀ + 24η / T * C
end

"""
    Ĉ(::AttractiveHardSpheres{<:SWHSPotential}, ::NonLinearSharmaSharma, k)

Returns the product of the bulk density and the Fourier transform (FT) of the
direct correlation function for an attractive square-well hard-spheres liquid
using the Non-Linear Sharma and Sharma approximation, that is, the sum of the FT of the
direct correlation function for the hard spheres core and
``exp[\\beta\\cdot\\u(r)]-1``, where ``u(r)`` is the potential:

```math
    u(r) = -\\epsilon\\sigma, \\qquad \\sigma < r \\le \\lambda\\sigma
```

`k` is the wavenumber.
"""
function Ĉ(liquid::AttractiveHardSpheres{P, U}, scheme::NonLinearSharmaSharma, k′) where
    {P, U<:SWHSPotential}

    η = volume_fraction(liquid)
    T = temperature(liquid)
    σ, λ = length_scale(potential(liquid))
    x = exp(1 / T) - 1

    λ³ = λ * λ * λ
    λ⁵ = λ³ * λ * λ

    k = σ * k′
    k² = k * k
    k³ = k² * k
    sink, cosk = sin(k), cos(k)
    sinλk, cosλk = sin(λ * k), cos(λ * k)

    smallk = k < 0.075

    C = smallk ? (λ³ - 1) / 3 - (λ⁵ - 1) / 30 * k² :
                 (cosk - λ * cosλk) / k² + (sinλk - sink) / k³

    C₀ = Ĉ(scheme.coreliquid, scheme.subscheme, k′)

    return C₀ + 24η * x * C
end

"""
    Ĉ(liquid::DipolarHardSpheres, scheme::MSA, k)

Returns the product of the bulk density and the Fourier transform of
projections of the direct correlation function `(C₀₀, C₁₀, C₁₁)` for a
dipolar hard-spheres liquid using the MSA approximation.

`k` is the wavenumber.
"""
function Ĉ(liquid::DipolarHardSpheres, scheme::MSA, k′)
    η = volume_fraction(liquid)
    σ = length_scale(potential(liquid))

    k  = σ * k′
    k² = k * k
    k⁶ = k² * k² * k²
    sink, cosk = sin(k), cos(k)

    α₀, α₁, α₂, α₃ = scheme.α₀, scheme.α₁, scheme.α₂, scheme.α₃
    β₀, β₁, β₂, β₃ = scheme.β₀, scheme.β₁, scheme.β₂, scheme.β₃
    a₁, a₃, b₁, b₃ = scheme.a₁, scheme.a₃, scheme.b₁, scheme.b₃

    smallk = k < SMALL_K # for small wavenumber use a Taylor series approximation

    C₁₀ = smallk ? (α₀ + 2β₀ - (α₁ + 2β₁) * k²) :
                   (  a₃ - a₁ * k² + (-a₃ + (α₃ + 2β₃) * k² / 3) * k * sink +
                    (-a₃ + ((a₁ + a₃ / 2) - (α₂ + 2β₂) * k² / 3) * k²) * cosk
                   ) / k⁶
    C₁₁ = smallk ? (α₀ -  β₀ - (α₁ -  β₁) * k²) :
                   (  b₃ - b₁ * k² + (-b₃ + (α₃ -  β₃) * k² / 3) * k * sink +
                    (-b₃ + ((b₁ + b₃ / 2) - (α₂ -  β₂) * k² / 3) * k²) * cosk
                   ) / k⁶

    C₀₀ = Ĉ(scheme.coreliquid, scheme.subscheme, k)
    C₁₀ = 24η * C₁₀
    C₁₁ = 24η * C₁₁

    return (C₀₀, C₁₀, C₁₁)
end
