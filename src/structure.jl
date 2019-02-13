struct StructureFactor{U, C, T}
    u::U
    c::C
end

function StructureFactor(u::U, c::C) where
         {U <: InteractionPotential, T, C <: ApproximationScheme{T}}
    return StructureFactor{U, C, T}(u, c)
end

function (S::StructureFactor{U, C, T})(k) where {U, C, T}
    return structure_factor(S.u, S.c, k)
end

"""
    structure_factor(::HardSphere, ::PercusYevick, η, k)

Returns the static structure factor for hard-spheres
using the Percus-Yevick approximation.

`η` is the effective volume fraction\n
`k` is the wavevector value
"""
function structure_factor(::HardSphere, c::PercusYevick, k)
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

    Ck = 24c.η * ( c.α * C₀ + c.β * C₁ + c.δ * C₃)

    return Sk = 1 / (1 - Ck)
end

"""
    structure_factor(p::DipolarHardSphere, c::MSA, k)

Returns the the proyections of the static structure factor
`(S₀₀, S₁₀, S₁₁)` for dipolar hard-spheres using the MSA
approximation.

`k` is the wavevector value
"""
function structure_factor(p::DipolarHardSphere, c::MSA, k)
    k  = k
    k² = k * k
    k⁶ = k² * k² * k²
    sink, cosk = sin(k), cos(k)

    smallk = k < 0.075

    α₀, α₁, α₂, α₃ = c.α₀, c.α₁, c.α₂, c.α₃
    β₀, β₁, β₂, β₃ = c.β₀, c.β₁, c.β₂, c.β₃
    a₁, a₃, b₁, b₃ = c.a₁, c.a₃, c.b₁, c.b₃

    C₁₀ = smallk ? (α₀ + 2β₀ - (α₁ + 2β₁) * k²) :
                   (  a₃ - a₁ * k² + (-a₃ + (α₃ + 2β₃) * k² / 3) * k * sink +
                    (-a₃ + ((a₁ + a₃ / 2) - (α₂ + 2β₂) * k² / 3) * k²) * cosk
                   ) / k⁶
    C₁₁ = smallk ? (α₀ -  β₀ - (α₁ -  β₁) * k²) :
                   (  b₃ - b₁ * k² + (-b₃ + (α₃ -  β₃) * k² / 3) * k * sink +
                    (-b₃ + ((b₁ + b₃ / 2) - (α₂ -  β₂) * k² / 3) * k²) * cosk
                   ) / k⁶

    Ck₁₀ = 24c.η * C₁₀
    Ck₁₁ = 24c.η * C₁₁

    Sk₀₀ = structure_factor(HardSphere(), c.c, k)
    Sk₁₀ = 1 / (1 - Ck₁₀)
    Sk₁₁ = 1 / (1 - Ck₁₁)

    return (Sk₀₀, Sk₁₀, Sk₁₁)
end

function structure_factor(::HardDisk, c::RosenfeldFMT, k)
    k² = k * k
    J₀ = besselj0(k / 2)
    J₁ = besselj1(k / 2)
    J₁′ = besselj1(k)

    smallk = k < 0.075

    C = smallk ? (c.A * (1 - k² / 16) / 4 + c.B * (1 - 3k² / 32) / 2 +
                  c.G * (1 - k² / 8 )) :
                 (c.A * (2J₁ / k)^2 + c.B * 2J₀ * J₁ / k + c.G * 2J₁′ / k)

    Ck = -4c.η * C

    return Sk = 1 / (1 - Ck)
end

"""
    structure_factor(::HardSphere, ::VerletWeis, η, k)

Returns the static structure factor for hard-spheres
using the Percus-Yevick approximation with the Verlet-Weis
correction.

`η` is the effective volume fraction\n
`k` is the wavevector value
"""
structure_factor(::HardSphere, c::VerletWeis, k) =
    structure_factor(HardSphere(), c.py, c.α * k)
