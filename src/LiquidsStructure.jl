module LiquidsStructure


### Imports
using SpecialFunctions


### Exports
export AttractiveHardSpheres, DipolarHardSpheres, HardDisks, HardSpheres, MSA,
       PercusYevick, RosenfeldFMT, SharmaSharma, SquareWell, StructureFactor,
       VerletWeis, Yukawa, dimensionality


### Implementation
include("liquids.jl")
include("potentials.jl")
include("schemes.jl")
include("structure.jl")
include("utils.jl")


end # module
