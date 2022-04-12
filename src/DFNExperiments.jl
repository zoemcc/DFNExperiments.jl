module DFNExperiments

using NeuralPDE, ModelingToolkit, Symbolics, DomainSets
using JSON
using PyCall

abstract type AbstractPyBaMMModel end

pybamm_func_str(::AbstractPyBaMMModel) = throw("Not implemented")

struct SPMModel <: AbstractPyBaMMModel
end

pybamm_func_str(spm::SPMModel) = "spm"


ty(x) = typeof(x)
fn(x) = fieldnames(x)
fnty = fn ∘ ty

include("generate_py.jl")

function __init__()
    initialize_pybamm_funcs()
end

function generate_sim_model(model::M; current_input=false, output_dir=nothing) where {M <: AbstractPyBaMMModel}
    model_str = pybamm_func_str(model)
    current_input_str = current_input ? "True" : "False" 
    sim, mtk_str, variables = py"solve_plot_generate(*$$(model_str)(), current_input=$$(current_input_str))"

    if typeof(output_dir) <: AbstractString
        if !isdir(output_dir)
            mkpath(output_dir)
        end
        model_filename = joinpath(output_dir, "model.jl")
    else
        model_filename = Base.tempname()
    end

    open(model_filename, "w") do f
        write(f, mtk_str)
    end
    include(model_filename)

    sol_data_json = sim.solution.save_data(variables=variables, to_format="json")
    sol_data = JSON.parse(sol_data_json)

    iv_names = nameof.(pde_system.ivs)
    dv_names = (nameof ∘ operation ∘ Symbolics.value).(pde_system.dvs)

    get_one_of_two(data_dict, key1, key2) = haskey(data_dict, key1) ? data_dict[key1] : data_dict[key2]

    ivs_sim::Vector{Vector{Float64}} = map(iv_names) do iv_name 
        data = get_one_of_two(sol_data, independent_variables_to_pybamm_names[iv_name], string(iv_name))
        first_cols = first.(data)
    end

    dvs_sim::Vector{Matrix{Float64}} = map(dv_names) do dv_name
        data = get_one_of_two(sol_data, dependent_variables_to_pybamm_names[dv_name], string(dv_name))
        data_mat = reduce(hcat, data)
        # make sure time is last axis, this needs the check so we don't flip data that only depends on time
        if size(data_mat, 1) != 1
            data_mat = collect(data_mat')
        end
        data_mat
    end

    sim_data = Dict(
        :ivs => Dict(
            :names => iv_names,
            :data => ivs_sim
        ),
        :dvs => Dict(
            :names => dv_names,
            :data => dvs_sim
        ) 
    )

    if typeof(output_dir) <: AbstractString
        sim_filename = joinpath(output_dir, "sim.json")
        open(sim_filename, "w") do f
            JSON.print(f, sim_data)
        end
    end

    (sim_data=sim_data, pde_system=pde_system)
end

function read_sim_data(output_dir::AbstractString)
    # either the folder that contains the file or the file itself is acceptable input
    if isfile(output_dir)
        sim_filename = output_dir
    elseif isdir(output_dir)
        sim_filename = joinpath(output_dir, "sim.json")
    else
        throw("Invalid path.  Path must be either the simulation json file or the directory that contains the sim.json file.")
    end
    sim_data_strs = JSON.parse(open(f->read(f, String), sim_filename, "r"))
    sim_data = Dict(
        :ivs => Dict(
            :names => Symbol.(sim_data_strs["ivs"]["names"]),
            :data => map(x->Float64.(x), (sim_data_strs["ivs"]["data"]))
        ),
        :dvs => Dict(
            :name => Symbol(sim_data_strs["dvs"]["names"]),
            :data => reduce.((v1, v2) -> hcat(Float64.(v1), Float64.(v2)), sim_data_strs["dvs"]["data"])
        )
    )

end

export ty, fn, fnty
export generate_sim_model, read_sim_data, pybamm_func_str
export AbstractPyBaMMModel, SPMModel


end
