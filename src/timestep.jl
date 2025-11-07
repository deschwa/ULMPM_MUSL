function timestep!(sim::MPMSimulation, alpha::Float64)
    reset_grid!(sim.grid)

    p2g!(sim)

    fix_dirichlet(sim)

    double_mapping!(sim, alpha)

    g2p!(sim)

    sim.t += sim.dt
end