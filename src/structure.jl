struct StructureFactor{U, C, T}
    u::U
    c::C
end

StructureFactor(u::U, c::C) where {U<:ApproximationScheme, T, C<:ApproximationScheme{T}} =
    StructureFactor{U, C, T}(u, c)

(S::StructureFactor{U, C, T})(k) where {U, C, T} = structure_factor(S.u, S.c, k)

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
    sink, cosk = sin(k), cos(k)

    T⁻¹ = 1 / c.T′
    α₀, α₁, α₃, β₁, β₃ = c.α₀, c.α₁, c.α₃, c.β₁, c.β₃
    a₁, a₃, b₁, b₃ = c.a₁, c.a₃, c.b₁, c.b₃

    C₁₀ = a₃ - a₁ * k² +
          (-a₃ + (α₀ + 2α₁ + 4α₃ + 2 * (5β₁ + 7β₃ - T⁻¹)) * k² / 3) * k * sink +
          (-a₃ + ((a₁ + a₃ / 2) - (α₀ + α₁ + α₃ + 2 * (β₁ + β₃ - T⁻¹)) * k² / 3) * k²) * cosk

    C₁₁ = b₃ - b₁ * k² +
          (-b₃ + (α₀ + 2α₁ + 4α₃ - (5β₁ + 7β₃ - T⁻¹)) * k² / 3) * k * sink +
          (-b₃ + ((b₁ + b₃ / 2) - (α₀ + α₁ + α₃ - (β₁ + β₃ - T⁻¹)) * k² / 3) * k²) * cosk

    S₀₀ = structure_factor(HardSphere(), c.py, k)
    S₁₀ = 1 / (1 - 24c.η * C₁₀ / k²^3)
    S₁₁ = 1 / (1 - 24c.η * C₁₁ / k²^3)

    return (S₀₀, S₁₀, S₁₁)
end

function structure_factor(::HardDisk, c::RosenfeldFMT, k)
    k² = k * k
    J₀ = besselj0(k / 2)
    J₁ = besselj1(k / 2)
    J₁′ = besselj1(k)

    smallk = k < 0.075

    α = smallk ? c.A * (0.25 - 0.015625k²) + c.B * (0.5 - 0.046875k²) + c.G *  (1 - 0.125k²) :
                 c.A * (2J₁ / k)^2  + (c.B * 2J₀ * J₁ / k) + (c.G * 2J₁′ / k)

    Ck = -4c.η * α

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
