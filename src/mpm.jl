module MPM

include("types.jl")
include("materials.jl")
include("boundary_conditions.jl")
include("helper_functions.jl")
include("shapefunctions.jl")
include("p2g.jl")
include("double_mapping.jl")
include("g2p.jl")
include("timestep.jl")


export timestep!, MaterialPointGroup, Grid, MPMSimulation, LinearElastic, NeoHookean

end