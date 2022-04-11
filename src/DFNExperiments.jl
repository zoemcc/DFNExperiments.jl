module DFNExperiments


using NeuralPDE, ModelingToolkit, Symbolics, DomainSets
using JSON
using PyCall

ty(x) = typeof(x)
fn(x) = fieldnames(x)
fnty = fn ∘ ty

include("generate_py.jl")

function __init__()
    initialize_pybamm_funcs()
end

function generate_sim_model(model_str; current_input=false, model_filename=nothing, sim_filename=nothing)
    current_input_str = current_input ? "True" : "False" 
    sim, mtk_str, variables = py"solve_plot_generate(*$$(model_str)(), current_input=$$(current_input_str))"

    if typeof(model_filename) <: AbstractString
        model_dir, _ = splitdir(model_filename)
        if !isdir(model_dir)
            mkpath(model_dir)
        end
    else
        model_filename = Base.tempname()
    end

    open(model_filename, "w") do f
        write(f, mtk_str)
    end
    include(model_filename)

    sol_data_json = sim.solution.save_data(variables=variables, to_format="json")
    sol_data = JSON.parse(sol_data_json)

    iv_names = [var.val.name for var in pde_system.ivs]
    dv_names = [var.val.f.name for var in pde_system.dvs]

    ivs_sim = Vector{Float64}[]
    for iv_name in iv_names
        full_name = independent_variables_to_pybamm_names[iv_name]
        if haskey(sol_data, full_name)
            data = sol_data[full_name]
        else
            data = sol_data[string(iv_name)]
        end 
        @show length(data)
        first_col = [datum[1] for datum in data]
        push!(ivs_sim, first_col)
    end

    dvs_sim = Matrix{Float64}[]
    for dv_name in dv_names
        full_name = dependent_variables_to_pybamm_names[dv_name]
        if haskey(sol_data, full_name)
            data = sol_data[full_name]
        else
            data = sol_data[string(iv_name)]
        end 
        @show length(data)
        
        data_mat = reduce(hcat, data)
        @show size(data_mat)
        push!(dvs_sim, data_mat)
    end

    sim_data = Dict(
        :ivs => [
            Dict(
                :name => iv_name,
                :data => data
            ) for (iv_name, data) in zip(iv_names, ivs_sim)
        ],
        :dvs => [
            Dict(
                :name => dv_name,
                :data => data
            ) for (dv_name, data) in zip(dv_names, dvs_sim)
        ]
    )

    if typeof(sim_filename) <: AbstractString
        sim_dir, _ = splitdir(sim_filename)
        if !isdir(sim_dir)
            mkpath(sim_dir)
        end
        open(sim_filename, "w") do f
            JSON.print(f, sim_data)
        end
    end

    (sim_data=sim_data, pde_system=pde_system)
end

function read_sim_data(sim_filename)
    sim_data_strs = JSON.parse(open(f->read(f, String), sim_filename, "r"))
    sim_data = Dict(
        :ivs => [
            Dict(
                :name => Symbol(iv["name"]),
                :data => Float64.(iv["data"])
            ) for iv in sim_data_strs["ivs"]
        ],
        :dvs => [
            Dict(
                :name => Symbol(dv["name"]),
                :data => reduce((v1, v2) -> hcat(Float64.(v1), Float64.(v2)), dv["data"])
            ) for dv in sim_data_strs["dvs"]
        ]
    )

end

export ty, fn, fnty
export generate_sim_model, read_sim_data


end
