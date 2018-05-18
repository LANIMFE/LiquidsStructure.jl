__precompile__()


module LiquidsStructure


### Implementation
include("potentials.jl")
include("schemes.jl")
include("structure.jl")
include("utils.jl")


### Exports
export ApproximationScheme, DipolarHardSphere, HardDisk, HardSphere,
       InteractionPotential, MSA, PercusYevick, RosenfeldFMT, StructureFactor,
       VerletWeis,
       dimension, structure_factor


end # module
