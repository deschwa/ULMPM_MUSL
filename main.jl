include("src/mpm.jl")
using .MPM

include("scripts/read_csv.jl")
include("scripts/write_csv.jl")

input_csv = "input/input.csv"
input_yaml = "input/input.yaml"

particle_path(t) = "output/particles/dump_p.$t.dat"
grid_path(t) = "output/grid/dump_g.$t.dat"

mat_dict, grid_dict, time_dict = get_sim_data(input_yaml)



# Create Simulation from files
mp_groups = create_particle_distribution_from_csv(input_csv, mat_dict)
grid = create_grid_from_yaml(grid_dict)

sim = MPMSimulation(
    mp_groups,
    grid,
    time_dict["dt"],
    0.0
)

set_borders_dirichlet!(sim.grid)

display(sim.mp_groups[1])

for mp_group in sim.mp_groups
    mp_group.ext_force_density[1] .+= @MVector [0.0, -9.81]
end

println("Starting simulation...")
calculated_steps = 0
while sim.t < time_dict["total_time"]
    timestep!(sim, 1.0)
    print("Calculating time: ", round(sim.t, digits=4), " or $(round(sim.t / time_dict["total_time"], digits=4)*100)%\r")
    if calculated_steps % 50 == 0
        write_particle_csv(sim, mat_dict, particle_path(calculated_steps))
        write_grid_csv(sim, grid_path(calculated_steps))
    end
    global calculated_steps += 1
    # println("$(sim.t),$(sim.mp_groups[1].pos[1][2]), $(sim.mp_groups[1].vel[1][2])")
end
println("Simulation complete.                   ")


# write_particle_csv(sim, mat_dict, output_particles)
# write_grid_csv(sim, output_grid)