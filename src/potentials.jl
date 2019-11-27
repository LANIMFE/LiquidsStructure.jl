abstract type Potential{N} end
abstract type RepulsivePotential{N}  <: Potential{N} end
abstract type AttractivePotential{N} <: Potential{N} end

struct CompositePotential{N, U₁<:Potential{N}, U₂<:Potential{N}} <: Potential{N}
    u₁::U₁
    u₂::U₂

    function CompositePotential{N, U₁, U₂}(u₁::U₁, u₂::U₂) where
        {N, U₁<:Potential{N}, U₂<:Potential{N}}
        return new(u₁, u₂)
    end
end
#
function CompositePotential(u₁::U₁, u₂::U₂) where
    {N, U₁<:Potential{N}, U₂<:Potential{N}}
    return CompositePotential{N, U₁, U₂}(u₁, u₂)
end
#

struct HardCore{T, N} <: RepulsivePotential{N}
    σ::T
end
#
HardCore{N}(σ::T) where {T, N} = (@assert σ > 0; HardCore{T, N}(σ))
HardCore{N}() where {N} = HardCore{N}(1)

struct DipoleDipole{T} <: Potential{3}
    μ::T
end
#
DipoleDipole() = DipoleDipole(1)

struct Yukawa{T, N} <: AttractivePotential{N}
    z::T
    ϵ::T
end
#
function Yukawa{N}(z, ϵ) where {N}
    @assert z > 0
    z′, ϵ′ = promote(z, ϵ)

    return Yukawa{typeof(z′), N}(z′, ϵ′)
end
#
Yukawa{N}(z; ϵ = 1) where {N} = Yukawa{N}(z, ϵ)

struct SquareWell{T, N} <: AttractivePotential{N}
    λ::T
    ϵ::T
end
#
function SquareWell{N}(λ, ϵ) where {N}
    @assert λ > 1
    λ′, ϵ′ = promote(λ, ϵ)

    return SquareWell{typeof(λ′), N}(λ′, ϵ′)
end
#
SquareWell{N}(λ; ϵ = 1) where {N} = SquareWell{N}(λ, ϵ)

# Aliases for composite potentials
const DHSPotential = CompositePotential{3, <:HardCore, <:DipoleDipole}
const YHSPotential = CompositePotential{3, <:HardCore, <:Yukawa}
const SWHSPotential = CompositePotential{3, <:HardCore, <:SquareWell}

energy_scale(potential::Potential) = 1
energy_scale(potential::AttractivePotential) = potential.ϵ

length_scale(potential::Potential) = 1
length_scale(potential::HardCore) = potential.σ
length_scale(potential::SquareWell) = potential.λ
length_scale(potential::Yukawa) = potential.z
length_scale(potential::DHSPotential) = length_scale(potential.u₁)
length_scale(potential::CompositePotential) =
    (length_scale(potential.u₁), length_scale(potential.u₂))
