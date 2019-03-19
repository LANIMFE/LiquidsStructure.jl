__precompile__()


module LiquidsStructure


### Imports
using SpecialFunctions


### Implementation
include("liquids.jl")
include("schemes.jl")
include("structure.jl")
include("utils.jl")


### Exports
export StructureFactor, DipolarHardSpheres, HardDisks, HardSpheres,
       MSA, PercusYevick, RosenfeldFMT, VerletWeis,
       dimensionality


end # module
