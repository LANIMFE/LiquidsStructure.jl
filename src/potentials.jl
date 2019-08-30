abstract type Potential end

struct Yukawa{T} <: Potential
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

struct SquareWell{T} <: Potential
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
