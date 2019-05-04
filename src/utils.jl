"""
    dimensionality(::Liquid)

Returns the dimensionality of a liquid system.
"""
dimensionality(::L) where {N, L <: Liquid{N}} = N
dimensionality(sk::StructureFactor) = dimensionality(sk.liquid)

"""
    dhs_msa_parameter(T′, η::T, tol)

Returns the solution `κ` to (see [Wertheim1971])

```math
    \\frac{(1 + 4\\kappa\\eta)^2}{(1 - 2\\kappa\\eta)^4} -
    \\frac{(1 - 2\\kappa\\eta)^2}{(1 + \\kappa\\eta)^4} =
    \\frac{8\\eta}{T^{\\prime}}
```

The solution is found by Newton's method
(up to the specified tolerance `tol` for `ξ = κ * η`, when possible).

[Wertheim1971] M. S. Wertheim, J. Chem. Phys. 55, 4291 (1971).
"""
function dhs_msa_parameter(liquid::DipolarHardSpheres, tol)
    y = 8 * liquid.η / liquid.T′

    ξ = 0.499 # 0 < ξ < 1 / 2

    while true
        ξ₋ = ξ

        ξ₁ = 1 - 2ξ
        ξ₂ = 1 +  ξ
        ξ₃ = 1 + 4ξ
        ξ₄ = 2 -  ξ

        ξ₁₂  = ξ₁ * ξ₂
        ξ₁₂⁴ = ξ₁₂^4
        ξ₁⁶  = ξ₁^6
        ξ₂⁴  = ξ₂^4
        ξ₃²  = ξ₃^2
        ξ₂⁶  = ξ₂⁴ * ξ₂^2

        ξ = ξ + (ξ₁₂ / 4) * (ξ₁⁶ - ξ₃² * ξ₂⁴ + y * ξ₁₂⁴) / (4ξ₃ * ξ₂⁶ + ξ₁⁶ * ξ₄)

        r = ξ / ξ₋

        if (1 - r) ≤ tol
            if r > 1
                ξ = ξ₋
            end
            break
        end
    end

    return ξ / liquid.η
end

# Printing methods
Base.show(io::IO, f::StructureFactor) = print(io, f.liquid, ',', f.scheme)
Base.show(io::IO, ::S) where {S <: ApproximationScheme} = print(io, nameof(S))
Base.show(io::IO, ::MSA{SS}) where {SS} = print(io, "MSA{", nameof(SS), '}')
function Base.show(io::IO, liquid::L) where {L <: Liquid}
    names = fieldnames(L)
    params = join((string(s, '=', getfield(liquid, s)) for s in names), ',')
    print(io, nameof(L), '(', params, ')')
end
