using StaticArrays
using Base.Threads


function double_mapping_depr!(sim::MPMSimulation, alpha::Float64)
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





function double_mapping!(sim::MPMSimulation, alpha::Float64)
    grid = sim.grid
    Nx, Ny = size(grid.mass)

    

    bin_maps = Dict{MaterialPointGroup, Array{Vector{Int64},2}}()
    @threads for mp_group in sim.mp_groups
        bin_maps[mp_group] = create_bin_map_2x2(grid, mp_group)

        # Get shape functions
        cache_shape_functions!(mp_group, grid)
    end


    
    
    
    # update particle positions and velocities
    @threads for mp_group in sim.mp_groups
        for p_idx in 1:length(mp_group.mass)
            pos_p = mp_group.pos[p_idx]

            v_update_temp = @MVector [0.0, 0.0]
            x_update_temp = @MVector [0.0, 0.0]

            for (idx, (i,j)) in enumerate(mp_group.node_cache[p_idx])
                N_Ip = mp_group.N_cache[p_idx][idx]

                x_update_temp .+= N_Ip * grid.v_new[i,j]
                v_update_temp .+= N_Ip * (grid.v_new[i,j] - alpha * grid.v[i,j])
            end

            mp_group.pos[p_idx] .+= sim.dt * x_update_temp

            mp_group.vel[p_idx] .*= alpha
            mp_group.vel[p_idx] .+= v_update_temp

        end
    end 

    # Update grid Momenta and velocities
    @threads for i in 1:Nx
        for j in 1:Ny
            momentum_new_temp = @MVector [0.0, 0.0]
            for mp_group in sim.mp_groups
                bin_map = bin_maps[mp_group]
                
                for p_idx in bin_map[i,j]
                    N_Ip = 0.0
                    for k in 1:4
                        if mp_group.node_cache[p_idx][k] == (i,j)
                            N_Ip = mp_group.N_cache[p_idx][k]
                            break
                        end
                    end
                    @assert N_Ip != 0.0 "Some error in the grid node caching"
                    
                    momentum_new_temp .+= N_Ip * mp_group.mass[p_idx] * mp_group.vel[p_idx]
                
                end

                grid.momentum_new[i,j] .= momentum_new_temp
                grid.v_new[i,j] .= momentum_new_temp / grid.mass[i,j]
            end
        
        end

    end


end