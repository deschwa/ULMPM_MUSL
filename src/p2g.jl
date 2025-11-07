using StaticArrays
using Base.Threads


function reset_grid!(grid::Grid)
    for i in 1:size(grid.v, 1), j in 1:size(grid.v, 2)
        grid.v[i,j] = @MVector [0.0, 0.0]
        grid.v_new[i,j] = @MVector [0.0, 0.0]
        grid.momentum[i,j] = @MVector [0.0, 0.0]
        grid.momentum_new[i,j] = @MVector [0.0, 0.0]
        grid.f_ext[i,j] = @MVector [0.0, 0.0]
        grid.f_int[i,j] = @MVector [0.0, 0.0]
    end
    
    fill!(grid.mass, 0.0) 
end



function create_bin_map_2x2(grid::Grid, mp_group::MaterialPointGroup)
    Nx = size(grid.mass, 1)
    Ny = size(grid.mass, 2)

    bin_map = [Vector{Int64}() for i in 1:Nx, j in 1:Ny]

    for particle_idx in 1:length(mp_group.mass)
        pos = mp_group.pos[particle_idx]
        
        for (i,j) in get_adjacent_grid_nodes(pos, grid)
            push!(bin_map[i, j], particle_idx)
        end
    end

    return bin_map
end


function p2g!(sim::MPMSimulation)
    grid = sim.grid

    Nx, Ny = size(grid.mass)

    bin_maps = Dict{MaterialPointGroup, Array{Vector{Int64},2}}()
    @threads for mp_group in sim.mp_groups
        bin_maps[mp_group] = create_bin_map_2x2(grid, mp_group)
    end

    @threads for i in 1:Nx
            for j in 1:Ny
            
            pos_ij = grid.pos[i,j]

            local_m = 0.0
            local_momentum = @MVector [0.0, 0.0]
            local_f_int = @MVector [0.0, 0.0]
            local_f_ext = @MVector [0.0, 0.0]

            for mp_group in sim.mp_groups
                bin_map = bin_maps[mp_group]
                
                particle_indices = bin_map[i,j]
                


                for p_idx in particle_indices
                    mass = mp_group.mass[p_idx]

                    rel_pos = mp_group.pos[p_idx] - pos_ij
                    N_Ip, ∇N_Ip = shape_function(rel_pos, grid.dx, grid.dy)

                    # println("Particle $p_idx at grid node ($i,$j): N_Ip = $N_Ip, ∇N_Ip = $∇N_Ip")

                    local_m += N_Ip * mass
                    local_momentum .+= N_Ip * mass * mp_group.vel[p_idx]
                    
                    local_f_ext .+= N_Ip * mp_group.ext_force_density[p_idx] * mass 
                    local_f_int .-= mp_group.volume[p_idx] * (mp_group.σ[p_idx] * ∇N_Ip)
                end
                # println("f_ext_local = $local_f_ext")
            
            end

            grid.mass[i,j] = local_m
            grid.momentum[i,j] .= local_momentum
            grid.f_ext[i,j] .= local_f_ext
            grid.f_int[i,j] .= local_f_int

            if local_f_ext[2] != 0.0
                # println("Node ($i,$j): local_f_ext = $local_f_ext, new momentum = $(grid.momentum_new[i,j])")
                # println("Node ($i,$j): mass = $(grid.mass[i,j])")
            end

            # Momentum update and velocity update
            grid.momentum_new[i,j] .= local_momentum .+ (local_f_ext .+ local_f_int) .* sim.dt
            if local_m <= 1e-12
                grid.v_new[i,j] .= @MVector [0.0, 0.0]
                grid.v[i,j] .= @MVector [0.0, 0.0]
                # println("Grid node ($i,$j): mass = $(grid.mass[i,j]), local_f_ext = $local_f_ext, local_f_int = $local_f_int")
            else
                grid.v_new[i,j] .= grid.momentum_new[i,j] / local_m
                grid.v[i,j] .= local_momentum / local_m

                # println("Node ($i,$j): momentum new should be $(local_momentum + (local_f_ext + local_f_int) * sim.dt), calculated: $(grid.momentum_new[i,j])")
                # println("f_ext at node ($i,$j): $local_f_ext, new momentum: $(grid.momentum_new[i,j]), new velocity: $(grid.v_new[i,j])")
            end
        end
    end
end