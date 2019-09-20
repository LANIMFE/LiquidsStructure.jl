abstract type PotentialSymmetry end
abstract type Isotropic   <: PotentialSymmetry end
abstract type Anisotropic <: PotentialSymmetry end

abstract type Potential{Symmetry} end

abstract type CorePotential{S}       <: Potential{S} end
abstract type LongRangePotential{S}  <: Potential{S} end
abstract type ShortRangePotential{S} <: Potential{S} end

abstract type HardCorePotential{S} <: CorePotential{S} end
abstract type SoftCorePotential{S} <: CorePotential{S} end

struct Yukawa{T} <: LongRangePotential{Isotropic}
    z::T
    σ::T
    ϵ::T
end

function Yukawa(z; σ = 1.0, ϵ = 1.0)
    @assert z > 0
    @assert σ > 0

    z′, σ′, ϵ′ = promote(z, σ, ϵ)

    return Yukawa{typeof(z′)}(z′, σ′, ϵ′)
end

struct SquareWell{T} <: ShortRangePotential{Isotropic}
    λ::T
    σ::T
    ϵ::T
end

function SquareWell(λ; σ = 1.0, ϵ = 1.0)
    @assert λ > 1
    @assert σ > 0

    λ′, σ′, ϵ′ = promote(λ, σ, ϵ)

    return SquareWell{typeof(λ′)}(λ′, σ′, ϵ′)
end
