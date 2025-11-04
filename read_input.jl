include("src/mpm.jl")
using CSV
using DataFrames
using StaticArrays
using .MPM
using YAML

input_csv = "input/input.csv"
input_yaml = "input/input.yaml"

function get_sim_data(input_yaml)
    file = YAML.load_file(input_yaml)

    # Material dict
    mat_dict = Dict{String, Any}()
    for material in file["materials"]
        if material["type"] == "LinearElastic"
            mat_dict[material["name"]] = LinearElastic(
                material["E"], material["nu"], material["rho"]
            )
        elseif material["type"] == "NeoHookean"
            mat_dict[material["name"]] = NeoHookean(
                material["E"], material["nu"], material["rho"]
            )
        end
    end

    # Grid Data
    grid_dict =  file["grid"]

    # time Data
    time_dict = file["time"]

    return mat_dict, grid_dict, time_dict
    
    
end

function create_particle_distribution_from_csv(input_file, material_dict)
    file = CSV.read(input_file, DataFrame)

    grouped_df = groupby(file, :type)

    mp_groups = []
    
    for (key, group) in pairs(grouped_df)
        type = key.type
        pos = [@MVector [Float64(group.x[i]), Float64(group.y[i])] for i in 1:length(group.x)]
        vel = [@MVector [Float64(group.vx[i]), Float64(group.vy[i])] for i in 1:length(group.vx)]
        mass = [group.mass[i] for i in 1:length(group.mass)]
        volume = [group.volume[i] for i in 1:length(group.volume)]
        material = material_dict[string(type)]
        mp_group = MaterialPointGroup(pos, vel, mass, volume, material)
        push!(mp_groups, mp_group)
    end

    mp_groups = tuple(mp_groups...)
    return mp_groups
end


function create_grid_from_yaml(grid_dict)
    Nx = grid_dict["Nx"]
    Ny = grid_dict["Ny"]
    minx = grid_dict["minx"]
    maxx = grid_dict["maxx"]
    miny = grid_dict["miny"]
    maxy = grid_dict["maxy"]

    return Grid(Nx, Ny, minx, maxx, miny, maxy)
    
end


function create_sim_from_csv(input_csv, input_yaml)
    material_dict, grid_dict, time_dict = get_sim_data(input_yaml)
    mp_groups = create_particle_distribution_from_csv(input_csv, material_dict)
    grid = create_grid_from_yaml(grid_dict)
    dt = time_dict["dt"]
    t_end = time_dict["t_end"]

    sim = MPMSimulation(mp_groups, grid, dt, t_end)
    return sim
end


create_sim_from_csv(input_csv, input_yaml)