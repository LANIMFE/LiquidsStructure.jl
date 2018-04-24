abstract type InteractionPotential end

struct HardDisk          <: InteractionPotential end
struct HardSphere        <: InteractionPotential end
struct DipolarHardSphere <: InteractionPotential end
