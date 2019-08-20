abstract type Liquid{N} end

struct HardDisks{T} <: Liquid{2}
    η::T
end

struct HardSpheres{T} <: Liquid{3}
    η::T
end

struct AttractiveHardSpheres{U, T} <: Liquid{3}
    η ::T
    T′::T
    potential::U

    function AttractiveHardSpheres(η, T, potential::U) where {U}
        θ = T / potential.ϵ
        η′, T′ = promote(η, θ)

        return new{U, typeof(η′)}(η′, T′, potential)
    end
end

struct DipolarHardSpheres{T} <: Liquid{3}
    η ::T
    T′::T
end
