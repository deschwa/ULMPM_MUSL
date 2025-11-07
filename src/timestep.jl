function timestep!(sim::MPMSimulation, alpha::Float64)
    reset_grid!(sim.grid)
    
    # println("Just before P2G: v=")
    # display(sim.grid.v)
    # println("\nJust before P2G: v_new=")
    # display(sim.grid.v_new)
    # println("Just before P2G: mass=")
    # display(sim.grid.mass)

    p2g!(sim)

    # println("Just after P2G: v=")
    # display(sim.grid.v)
    # println("\nJust after P2G: v_new=")
    # display(sim.grid.v_new)
    # println("Just after P2G: mass=")
    # display(sim.grid.mass)
    # println("Just after P2G: momentum=")
    # display(sim.grid.momentum)
    # println("Just after P2G: momentum_new=")
    # display(sim.grid.momentum_new)
    # println("Just after P2G: f_int=")
    # display(sim.grid.f_int)
    # println("Just after P2G: f_ext=")
    # display(sim.grid.f_ext) 

    fix_dirichlet(sim)

    # println("Just after boundaries: v=")
    # display(sim.grid.v)
    # println("\nJust after boundaries: v_new=")
    # display(sim.grid.v_new)
    # println("Just after boundaries: mass=")
    # display(sim.grid.mass)

    double_mapping!(sim, alpha)

    g2p!(sim)

    sim.t += sim.dt
end