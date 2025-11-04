using CSV
using DataFrames
using .MPM
using YAML

input_csv = "input/input.csv"
input_yaml = "input/input.yaml"

function get_sim_data(input_yaml)
    file = YAML.load_file(input_yaml)

    # Material dict
    mat_dict = {}
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

function read_particle_distribution(input_file, material_dict)
    file = CSV.read(input_file, DataFrames)

    grouped_df = groupby(file, :type)
    
    for
end