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
HardCore{N}(σ::T) where {T, N} = HardCore{T, N}(σ)
HardCore{N}() where {N} = HardCore{N}(1)

struct DipoleDipole{T} <: Potential{3}
    μ::T
end
#
DipoleDipole() = DipoleDipole(1)

struct Yukawa{T, N} <: AttractivePotential{N}
    z::T
    σ::T
    ϵ::T
end
#
function Yukawa{N}(z, σ, ϵ) where {N}
    @assert z > 0
    @assert σ > 0

    z′, σ′, ϵ′ = promote(z, σ, ϵ)

    return Yukawa{typeof(z′), N}(z′, σ′, ϵ′)
end
#
Yukawa{N}(z; σ = 1, ϵ = 1) where {N} = Yukawa{N}(z, σ, ϵ)

struct SquareWell{T, N} <: AttractivePotential{N}
    λ::T
    σ::T
    ϵ::T
end
#
function SquareWell{N}(λ, σ, ϵ) where {N}
    @assert λ > 1
    @assert σ > 0

    λ′, σ′, ϵ′ = promote(λ, σ, ϵ)

    return SquareWell{typeof(λ′), N}(λ′, σ′, ϵ′)
end
#
SquareWell{N}(λ; σ = 1, ϵ = 1) where {N} = SquareWell{N}(λ, σ, ϵ)

energy_scale(potential::AttractivePotential) = potential.ϵ
length_scale(potential::CompositePotential) = length_scale(potential.u₁)
length_scale(potential::Union{AttractivePotential, HardCore}) = potential.σ

# Aliases for composite potentials
const DHSPotential = CompositePotential{3, <:HardCore, <:DipoleDipole}
const YHSPotential = CompositePotential{3, <:HardCore, <:Yukawa}
const SWHSPotential = CompositePotential{3, <:HardCore, <:SquareWell}
