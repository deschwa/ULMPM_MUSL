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


# Equality overloading
import Base.==

function ==(m1::LinearElastic, m2::LinearElastic)
    return m1.E == m2.E && m1.ν == m2.ν && m1.ρ == m2.ρ
end

function ==(m1::NeoHookean, m2::NeoHookean)
    return m1.μ == m2.μ && m1.λ == m2.λ && m1.ρ == m2.ρ
end

"""
struct MaterialPointGroup

Contains all information about a group of Material Points, eg. of one body or one material.
"""
struct MaterialPointGroup{MaterialType<:AbstractMaterial}
    pos::Vector{MVector{2, Float64}}        # Positions
    vel::Vector{MVector{2, Float64}}        # Velocities
    ext_force_density::Vector{MVector{2, Float64}}  # External Forces

    mass::Vector{Float64}                   # Masses
    volume::Vector{Float64}                 # Volumes
    volume_0::Vector{Float64}               # Initial Volumes
    density::Vector{Float64}                # Densities

    F::Vector{MMatrix{2,2,Float64,4}}       # Deformation Gradients
    σ::Vector{MMatrix{2,2,Float64,4}}       # Cauchy Stresses  
    L::Vector{MMatrix{2,2,Float64,4}}       # Velocity Gradients

    material::MaterialType                  # Material Model
    type::String                            # Material Type

    node_cache::Vector{MVector{4, Tuple{Int64, Int64}}} # Cache for Indices of adjacent nodes
    N_cache::Vector{MVector{4, Float64}}                # Cache for Shape function values of corresponding nodes
    ∇N_cache::Vector{MVector{4, MVector{2, Float64}}}   # Cache for ∇N    

    # Constructor using only pos, vel, mass, volume, and material. AbstractString is used because sometimes CSVs are read as String7
    function MaterialPointGroup(pos::Vector{MVector{2, Float64}}, 
                                vel::Vector{MVector{2, Float64}},
                                mass::Vector{Float64}, 
                                volume::Vector{Float64}, 
                                material::MaterialType, type::AbstractString) where {MaterialType<:AbstractMaterial}


        N = length(pos)

        density = [m/v for (m,v) in zip(mass, volume)]
        volume_0 = deepcopy(volume)

        F = [MMatrix{2,2,Float64,4}(I(2)) for _ in 1:N]
        σ = [MMatrix{2,2,Float64,4}(zeros(2,2)) for _ in 1:N]
        L = [MMatrix{2,2,Float64,4}(zeros(2,2)) for _ in 1:N]

        ext_force = [MVector(0.0, 0.0) for _ in 1:N]

        node_cache = [@MVector [(0,0) for _ in 1:4] for _ in 1:N]
        N_cache = [@MVector [0.0 for _ in 1:4] for _ in 1:N]
        ∇N_cache = [@MVector [@MVector [0.0, 0.0] for _ in 1:4] for _ in 1:N]


        new{MaterialType}(pos,
            vel,
            ext_force,
            mass,
            volume,
            volume_0,
            density,
            F,
            σ,
            L,
            material,
            type,
            node_cache,
            N_cache,
            ∇N_cache)

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
    minx::Float64                                   # Minimum x-coordinate
    miny::Float64                                   # Minimum y-coordinate

    dirichlet_nodes::Vector{Tuple{Int64, Int64}}    # List of Dirichlet Node indices

    function Grid(Nx::Int64, Ny::Int64, minx::Float64, maxx::Float64, miny::Float64, maxy::Float64)
        dx::Float64 = (maxx - minx) / (Nx-1)
        dy::Float64 = (maxy - miny) / (Ny-1)

        pos = [@SVector [round(minx + (i-1)*dx, digits=15), round(miny + (j-1)*dy, digits=15)] for i in 1:Nx, j in 1:Ny]

        v = [@MVector [0.0, 0.0] for i in 1:Nx, j in 1:Ny]
        v_new = [@MVector [0.0, 0.0] for i in 1:Nx, j in 1:Ny]
        momentum = [@MVector [0.0, 0.0] for i in 1:Nx, j in 1:Ny]
        momentum_new = [@MVector [0.0, 0.0] for i in 1:Nx, j in 1:Ny]
        f_ext = [@MVector [0.0, 0.0] for i in 1:Nx, j in 1:Ny]
        f_int = [@MVector [0.0, 0.0] for i in 1:Nx, j in 1:Ny]
        
        m = [0.0 for i in 1:Nx, j in 1:Ny]

        empty_tuple_vector = Vector{Tuple{Int64, Int64}}()
    

        new(pos, 
        v, v_new, 
        momentum, momentum_new,
        f_ext, f_int,
        m, dx, dy, minx, miny, empty_tuple_vector)
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


