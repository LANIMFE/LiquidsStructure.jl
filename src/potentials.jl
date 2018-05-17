abstract type InteractionPotential{N} end

struct HardDisk          <: InteractionPotential{2} end
struct HardSphere        <: InteractionPotential{3} end
struct DipolarHardSphere <: InteractionPotential{3} end
