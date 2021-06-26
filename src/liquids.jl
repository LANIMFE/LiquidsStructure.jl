struct LiquidProperties{P, U}
    η::P # volume fraction
    T::P # reduced temperature
    potential::U
end

LiquidProperties(η::P, potential::U) where {P, U<:Potential} =
    LiquidProperties{P, U}(η, one(P), potential)

potential(p::LiquidProperties) = p.potential
temperature(p::LiquidProperties) = p.T
volume_fraction(p::LiquidProperties) = p.η


### Liquids types
abstract type AbstractLiquid{U} end

const LIQUIDS_NAMES = (
    :AttractiveHardSpheres,
    :DipolarHardSpheres,
    :HardDisks,
    :HardSpheres
)

for liquid in LIQUIDS_NAMES
    @eval begin
        struct $liquid{P, U} <: AbstractLiquid{U}
            properties::LiquidProperties{P, U}
        end
    end
end


### Liquids constructors
function HardDisks(η::Real)
    properties = LiquidProperties(η, HardCore{2}())
    return HardDisks(properties)
end

function HardSpheres(η::Real)
    properties = LiquidProperties(η, HardCore{3}())
    return HardSpheres(properties)
end

function AttractiveHardSpheres(η::Real, T::Real, potential::CompositePotential)
    η′, T′ = promote(η, T)
    properties = LiquidProperties(η′, T′, potential)
    return AttractiveHardSpheres(properties)
end
#
function AttractiveHardSpheres(η::Real, T::Real, potential::AttractivePotential)
    T′ = T / energy_scale(potential)
    potential′ = CompositePotential(HardCore{3}(), potential)
    return AttractiveHardSpheres(η, T′, potential′)
end

function DipolarHardSpheres(η::Real, T::Real)
    η′, T′ = promote(η, T)
    potential = CompositePotential(HardCore{3}(), DipoleDipole())
    properties = LiquidProperties(η′, T′, potential)
    return DipolarHardSpheres(properties)
end

potential(l::AbstractLiquid) = potential(l.properties)
temperature(l::AbstractLiquid) = temperature(l.properties)
volume_fraction(l::AbstractLiquid) = volume_fraction(l.properties)
