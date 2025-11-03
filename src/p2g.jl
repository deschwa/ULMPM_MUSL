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

        # Cell indices. Ceil is used to handle the 1-based array indices in julia
        cell_idx_x = ceil(Int64, pos[1]/grid.dx)    # <=> floor(...) + 1
        cell_idx_y = ceil(Int64, pos[1]/grid.dy)
        
        # push particle index to 
        push!(bin_map[cell_idx_x, cell_idx_y], particle_idx)
    end
end