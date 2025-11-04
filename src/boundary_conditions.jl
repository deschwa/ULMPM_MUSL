using StaticArrays

function fix_dirichlet(sim::MPMSimulation)
    for (i,j) in sim.grid.dirichlet_nodes
        sim.grid.v[i,j] .= @MVector [0.0, 0.0]
        sim.grid.v_new[i,j] .= @MVector [0.0, 0.0]
        sim.grid.momentum[i,j] .= @MVector [0.0, 0.0]
        sim.grid.momentum_new[i,j] .= @MVector [0.0, 0.0]
    end
end