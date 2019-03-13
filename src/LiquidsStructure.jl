__precompile__()


module LiquidsStructure


### Imports
using SpecialFunctions


### Implementation
include("potentials.jl")
include("schemes.jl")
include("structure.jl")
include("utils.jl")


### Exports
export DipolarHardSphere, HardDisk, HardSphere,
       MSA, PercusYevick, RosenfeldFMT, StructureFactor, VerletWeis,
       dimension, structure_factor


end # module
