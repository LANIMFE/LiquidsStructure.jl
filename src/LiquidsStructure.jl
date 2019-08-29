__precompile__()


module LiquidsStructure


### Imports
using SpecialFunctions


### Implementation
include("liquids.jl")
include("potentials.jl")
include("schemes.jl")
include("structure.jl")
include("utils.jl")


### Exports
export StructureFactor, AttractiveHardSpheres, DipolarHardSpheres, HardDisks,
       HardSpheres, MSA, PercusYevick, RosenfeldFMT, SharmaSharma, VerletWeis,
       Yukawa, SquareWell, dimensionality


end # module
