using StaticArrays
using Base.Threads


function double_mapping!(sim::MPMSimulation, alpha::Float64)
    grid = sim.grid
    Nx, Ny = size(grid.mass)

    # Update Particle Positions
    for mp_group in sim.mp_groups
        @threads for p_idx in 1:length(mp_group.mass)
            
            # Update particle position
            pos_p = mp_group.pos[p_idx]
            p_sum = @MVector [0.0, 0.0]

            for (i,j) in get_adjacent_grid_nodes(pos_p, grid)
                pos_ij = grid.pos[i,j]
                rel_pos = pos_p - pos_ij
                N_Ip, ∇N_Ip = shape_function(rel_pos, grid.dx, grid.dy)
                p_sum .+= N_Ip * grid.v_new[i,j]
            end            

            mp_group.pos[p_idx] .+= p_sum * sim.dt


            # Update velocity
            v_sum = @MVector [0.0, 0.0]
            pos_p = mp_group.pos[p_idx]

            for (i,j) in get_adjacent_grid_nodes(pos_p, grid)
                pos_ij = grid.pos[i,j]
                rel_pos = pos_p - pos_ij
                N_Ip, ∇N_Ip = shape_function(rel_pos, grid.dx, grid.dy)
                v_sum .+= N_Ip * (grid.v_new[i,j] - alpha * grid.v[i,j])
            end

            mp_group.vel[p_idx] .*= alpha
            mp_group.vel[p_idx] .+= v_sum
        end
    end

    #Update grid momenta and velocities
    for i in 1:Nx, j in 1:Ny
        pos_ij = grid.pos[i,j]

        local_momentum_new = @MVector [0.0, 0.0]
        local_mass = 0.0
        

        for mp_group in sim.mp_groups
            bin_map = create_bin_map_2x2(grid, mp_group)
            
            particle_indices = bin_map[i,j]

            
            for p_idx in particle_indices
                pos_p = mp_group.pos[p_idx]

                rel_pos = pos_p - pos_ij
                N_Ip, ∇N_Ip = shape_function(rel_pos, grid.dx, grid.dy)

                local_momentum_new .+= N_Ip * mp_group.mass[p_idx] * mp_group.vel[p_idx]
                local_mass += N_Ip * mp_group.mass[p_idx]
            end        
        end

        grid.momentum_new[i,j] .= local_momentum_new
        grid.v_new[i,j] .= local_momentum_new / local_mass
        
    end

    # Fix boundaries
    fix_dirichlet(sim)
    
end