abstract type Liquid{N} end

struct HardDisks{T} <: Liquid{2}
    η::T
end

struct HardSpheres{T} <: Liquid{3}
    η::T
end

struct DipolarHardSpheres{T} <: Liquid{3}
    η ::T
    T′::T
end
