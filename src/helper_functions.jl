using StaticArrays

"""
(i  ,j+1)           (i+1,j+1)
        
           Particle

(i  ,j  )           (i+1,j  )
"""
function get_adjacent_grid_nodes(pos::MVector{2, Float64}, grid::Grid)

    i = floor(Int64, (pos[1]-grid.minx) / grid.dx) + 1  # 1-based indexing
    j = floor(Int64, (pos[2]-grid.miny) / grid.dy) + 1  # 1-based indexing

    return [
        (i, j),
        (i+1, j),
        (i, j+1),
        (i+1, j+1)
    ]
end






