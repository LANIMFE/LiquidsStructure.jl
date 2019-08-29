abstract type Potential end

struct Yukawa{T} <: Potential
    Z::T
    σ::T
    ϵ::T
end

function Yukawa(Z; σ = 1.0, ϵ = 1.0)
    @assert Z > 0
    @assert σ > 0

    Z′, σ′, ϵ′ = promote(Z, σ, ϵ)

    return Yukawa{typeof(Z′)}(Z′, σ′, ϵ′)
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
