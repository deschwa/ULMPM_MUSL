using StaticArrays



function reset_grid!(grid::Grid)
    zero_vector = @MVector [0.0, 0.0]

    fill!(grid.v, zero_vector)
    fill!(grid.v_new, zero_vector)
    fill!(grid.momentum, zero_vector)
    fill!(grid.momentum_new, zero_vector)
    fill!(grid.f_ext, zero_vector)
    fill!(grid.f_int, zero_vector)
    fill!(grid.mass, 0.0)
    
end



function create_bin_map_2x2(grid::Grid, mp_group::MaterialPointGroup)
    Nx = size(grid.m, 1)
    Ny = size(grid.m, 2)

    bin_map = [Vector{Int64}() for _ in (1:Nx, 1:Ny)]

    for particle_idx in 1:length(mp_group.mass)
        pos = mp_group.pos[particle_idx]
        
        for (i,j) in get_adjacent_grid_nodes(pos, grid)
            push!(bin_map[i, j], particle_idx)
        end
    end
end


function p2g!(sim::MPMSimulation)
    grid = sim.grid

    for (i,j) in CartesianIndex(grid.mass)
        pos_ij = grid.pos[i,j]

        local_m = 0.0
        local_momentum = @MVector [0.0, 0.0]
        local_f_int = @MVector [0.0, 0.0]
        local_f_ext = @MVector [0.0, 0.0]

        for mp_group in sim.mp_groups
            bin_map = create_bin_map_2x2(grid, mp_group)
            
            particle_indices = bin_map[i,j]
            


            for p_idx in particle_indices
                mass = sim.mp_group.mass[p_idx]

                rel_pos = sim.mp_group.pos[p_idx] - pos_ij
                N_Ip, ∇N_Ip = shape_function(rel_pos, grid.dx, grid.dy)

                local_m += N_Ip * mass
                local_momentum .+= N_Ip * mass * sim.mp_group.vel[p_idx]
                
                local_f_ext .+= N_Ip * sim.mp_group.ext_force_density[p_idx] * mass 
                local_f_int .-= sim.mp_group.volume[p_idx] * (sim.mp_group.σ[p_idx] * ∇N_Ip)
            end
        
        end

        grid.mass[i,j] = local_m
        grid.momentum[i,j] .= local_momentum
        grid.f_ext[i,j] .= local_f_ext
        grid.f_int[i,j] .= local_f_int

        # Momentum update and velocity update
        grid.momentum_new[i,j] .= local_momentum + (local_f_ext + local_f_int) * sim.dt
        grid.v_new[i,j] .= grid.momentum_new[i,j] / local_m
        grid.v[i,j] .= local_momentum / local_m
    end
end