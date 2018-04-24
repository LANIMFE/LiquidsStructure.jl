__precompile__()


module LiquidStructure


### Implementation
include("potentials.jl")
include("schemes.jl")
include("structure.jl")


### Exports
export ApproximationScheme, DipolarHardSphere, HardDisk, HardSphere,
       InteractionPotential, MSA, PercusYevick, RosenfeldFMT, StructureFactor,
       VerletWeis,
       structure_factor


end # module
