using StaticArrays

function fix_dirichlet(sim::MPMSimulation)
    for (i,j) in sim.grid.dirichlet_nodes
        sim.grid.v[i,j] .= @MVector [0.0, 0.0]
        sim.grid.v_new[i,j] .= @MVector [0.0, 0.0]
        sim.grid.momentum[i,j] .= @MVector [0.0, 0.0]
        sim.grid.momentum_new[i,j] .= @MVector [0.0, 0.0]
    end
end


function set_borders_dirichlet!(grid::Grid)
    Nx, Ny = size(grid.mass)

    for i in 1:Nx
        push!(grid.dirichlet_nodes, (i, 1))       # Bottom border
        push!(grid.dirichlet_nodes, (i, Ny))      # Top border
    end

    for j in 1:Ny
        push!(grid.dirichlet_nodes, (1, j))       # Left border
        push!(grid.dirichlet_nodes, (Nx, j))      # Right border
    end

    println("Dirichlet nodes set: ", length(grid.dirichlet_nodes))
end