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
