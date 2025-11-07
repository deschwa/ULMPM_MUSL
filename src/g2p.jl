using StaticArrays
using Base.Threads

function g2p!(sim::MPMSimulation)
    grid = sim.grid

    
    for mp_group in sim.mp_groups
        @threads for p_idx in 1:length(mp_group.mass)
            pos_p = mp_group.pos[p_idx]

            L_sum = MMatrix{2,2,Float64,4}(zeros(2,2))

            for (i,j) in get_adjacent_grid_nodes(pos_p, grid)
                pos_ij = grid.pos[i,j]
                rel_pos = pos_p - pos_ij
                N_Ip, ∇N_Ip = shape_function(rel_pos, grid.dx, grid.dy)

                L_sum .+=  ∇N_Ip * mp_group.vel[p_idx]'
            end

            mp_group.L[p_idx] .= L_sum

            mp_group.F[p_idx] .= (MMatrix{2,2,Float64,4}(I) + L_sum * sim.dt) * mp_group.F[p_idx]

            mp_group.volume[p_idx] = det(mp_group.F[p_idx]) * mp_group.volume_0[p_idx]


        end

        stress_update!(mp_group)
    end

end