using StaticArrays
using Base.Threads

function g2p!(sim::MPMSimulation)
    grid = sim.grid

    
    for mp_group in sim.mp_groups
        @threads for p_idx in 1:length(mp_group.mass)
            pos_p = mp_group.pos[p_idx]

            L_sum = MMatrix{2,2,Float64,4}(zeros(2,2))

            for (idx, (i,j)) in enumerate(mp_group.node_cache[p_idx])
                ∇N_Ip = mp_group.∇N_cache[p_idx][idx]

                L_sum .+=  grid.v_new[i,j] * ∇N_Ip'
            end

            mp_group.L[p_idx] .= L_sum

            mp_group.F[p_idx] .= (MMatrix{2,2,Float64,4}(I) + L_sum * sim.dt) * mp_group.F[p_idx]

            mp_group.volume[p_idx] = det(mp_group.F[p_idx]) * mp_group.volume_0[p_idx]


        end

        stress_update!(mp_group, sim.dt)
    end

end