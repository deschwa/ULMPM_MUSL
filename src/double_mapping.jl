using StaticArrays
using Base.Threads



function double_mapping!(sim::MPMSimulation, alpha::Float64)
    grid = sim.grid
    Nx, Ny = size(grid.mass)    
    
    
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
                bin_map = mp_group.bin_map_cache
                
                for p_idx in bin_map[i,j]
                    # Set N_ip via cache
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

            end
            grid.momentum_new[i,j] .= momentum_new_temp
            grid.v_new[i,j] .= momentum_new_temp / grid.mass[i,j]
            
        end

    end


end