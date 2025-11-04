using StaticArrays


"""
Materials
"""
abstract type AbstractMaterial end

struct LinearElastic<:AbstractMaterial
    E::Float64        # Young's Modulus
    ν::Float64        # Poisson's Ratio
    ρ::Float64        # Density
end

struct NeoHookean<:AbstractMaterial
    μ::Float64        # Shear Modulus
    λ::Float64       # Lamé's First Parameter
    ρ::Float64        # Density
end


"""
struct MaterialPointGroup

Contains all information about a group of Material Points, eg. of one body or one material.
"""
struct MaterialPointGroup{MaterialType<:AbstractMaterial}
    pos::Vector{MVector{2, Float64}}        # Positions
    vel::Vector{MVector{2, Float64}}        # Velocities
    ext_force::Vector{MVector{2, Float64}}  # External Forces

    mass::Vector{Float64}                   # Masses
    volume::Vector{Float64}                 # Volumes
    volume_0::Vector{Float64}               # Initial Volumes
    density::Vector{Float64}                # Densities

    F::Vector{MMatrix{2,2,Float64,4}}       # Deformation Gradients
    σ::Vector{MMatrix{2,2,Float64,4}}       # Cauchy Stresses  
    L::Vector{MMatrix{2,2,Float64,4}}       # Velocity Gradients

    material::MaterialType                  # Material Model
    

    # Constructor using only pos, vel, mass, volume, and material
    function MaterialPointGroup(pos::Vector{MVector{2, Float64}}, 
                                vel::Vector{MVector{2, Float64}},
                                mass::Vector{Float64}, 
                                volume::Vector{Float64}, 
                                material::MaterialType) where {MaterialType<:AbstractMaterial}


        N = length(pos)
        density = [m/v for (m,v) in zip(mass, volume)]
        volume_0 = copy(volume)
        F = [MMatrix{2,2,Float64,4}(I(2)) for _ in 1:N]
        σ = [MMatrix{2,2,Float64,4}(zeros(2,2)) for _ in 1:N]
        L = [MMatrix{2,2,Float64,4}(zeros(2,2)) for _ in 1:N]
        ext_force = [MVector(0.0, 0.0) for _ in 1:N]
        new{material}(pos,
            vel,
            ext_force,
            mass,
            volume,
            volume_0,
            density,
            F,
            σ,
            L,
            material)

    end

end


"""
struct Grid

Contains all necessary information about the grid.
"""
struct Grid
    pos::Array{SVector{2, Float64}, 2}              # Positions

    v::Array{MVector{2, Float64}, 2}                # Velocities
    v_new::Array{MVector{2, Float64}, 2}            # New Velocities

    momentum::Array{MVector{2, Float64}, 2}         # Momenta
    momentum_new::Array{MVector{2, Float64}, 2}     # New Momenta

    f_ext::Array{MVector{2, Float64}, 2}            # External Forces
    f_int::Array{MVector{2, Float64}, 2}            # Internal Forces

    mass::Array{Float64, 2}                         # Masses

    dx::Float64                                     # Grid Spacing in x
    dy::Float64                                     # Grid Spacing in y

    dirichlet_nodes::Vector{Tuple{Int64, Int64}}    # List of Dirichlet Node indices

    function Grid(Nx, Ny, minx, maxx, miny, maxy)
        dx = (maxx - minx) / Nx
        dy = (maxy - miny) / Ny

        pos = [@SVector [minx + (i-1)*dx, miny + (j-1)*dy] for (i,j) in (1:Nx,1:Ny)]

        zero_vec_arr = [@MVector [0.0, 0.0] for _ in (1:Nx,1:Ny)]
        
        m = [0.0 for _ in (1:Nx,1:Ny)]
        
        new(Nx, Ny, pos, 
        zero_vec_arr, zero_vec_arr, zero_vec_arr, zero_vec_arr,
        m, dx, dy)
    end
end


"""
MPMSimulation Type. Contains all necessary information about an MPM simulation.
"""
mutable struct MPMSimulation{MPGroupTouple <: Tuple}
    mp_groups::MPGroupTouple        # Tuple of MaterialPointGroups
    grid::Grid
    dt::Float64
    t::Float64
    
end


