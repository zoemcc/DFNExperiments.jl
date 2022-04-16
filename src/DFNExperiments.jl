module DFNExperiments

using DiffEqFlux, Flux
using NeuralPDE, ModelingToolkit, Symbolics, DomainSets
import ModelingToolkit: Interval, infimum, supremum
using LabelledArrays
using JSON
using PyCall
using IfElse
using Plots, LinearAlgebra, TensorBoardLogger

include("utils.jl")

abstract type AbstractApplyFuncType end

include("multi_dimensional_function.jl")

abstract type AbstractPyBaMMModel end

pybamm_func_str(::AbstractPyBaMMModel) = throw("Not implemented")

get_model_dir(model::AbstractPyBaMMModel) = abspath(joinpath(@__DIR__, "..", "models", "$(pybamm_func_str(model))"))
function load_model(model::AbstractPyBaMMModel)
    model_dir = get_model_dir(model)
    include(joinpath(model_dir, "model.jl"))
    pde_system
end

struct SPMModel <: AbstractPyBaMMModel
end

pybamm_func_str(spm::SPMModel) = "spm"
struct SPMeModel <: AbstractPyBaMMModel
end

pybamm_func_str(spm::SPMeModel) = "spme"


include("generate_py.jl")

function __init__()
    initialize_pybamm_funcs()
end

include("plot_log.jl")

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
    @show iv_names
    @show independent_variables_to_pybamm_names
    dv_names = (nameof ∘ operation ∘ Symbolics.value).(pde_system.dvs)
    @show dv_names
    @show dependent_variables_to_pybamm_names
    @show keys(sol_data)
    #return sol_data, variables, sim

    get_one_of_two(data_dict, key1, key2) = haskey(data_dict, key1) ? data_dict[key1] : data_dict[key2]

    ivs_sim::Vector{Vector{Float64}} = map(iv_names) do iv_name 
        data = get_one_of_two(sol_data, independent_variables_to_pybamm_names[iv_name], string(iv_name))
        first_cols = first.(data)
    end

    dvs_sim::Vector{Array{Float64}} = map(dv_names) do dv_name
        data = get_one_of_two(sol_data, dependent_variables_to_pybamm_names[dv_name], string(dv_name))
        data_array = recursive_array_to_array(data; do_reverse=true)
    end

    sim_data = Dict(
        :ivs => Dict(
            :names => iv_names,
            :data => ivs_sim,
        ),
        :dvs => Dict(
            :names => dv_names,
            :deps => [collect(dependent_variables_to_dependencies[dv_name]) for dv_name in dv_names], 
            :data => dvs_sim,
        ) 
    )

    if typeof(output_dir) <: AbstractString
        sim_filename = joinpath(output_dir, "sim.json")
    else
        sim_filename = Base.tempname()
    end

    open(sim_filename, "w") do f
        JSON.print(f, sim_data)
    end
    sim_data_nt = read_sim_data(sim_filename)

    (sim_data=sim_data_nt, pde_system=pde_system)
end

function read_sim_data(output_dir_or_file::AbstractString)
    # either the folder that contains the file or the file itself is acceptable input
    if isfile(output_dir_or_file)
        sim_filename = output_dir_or_file
    elseif isdir(output_dir_or_file)
        sim_filename = joinpath(output_dir_or_file, "sim.json")
    else
        throw("Invalid path.  Path must be either the simulation json file or the directory that contains the sim.json file.")
    end
    sim_data_strs = JSON.parse(open(f->read(f, String), sim_filename, "r"))
    iv_syms = tuple(Symbol.(sim_data_strs["ivs"]["names"])...)
    iv_data = tuple(map(x->Float64.(x), (sim_data_strs["ivs"]["data"]))...)
    iv_nt = NamedTuple{iv_syms}(iv_data)
    dv_syms = tuple(Symbol.(sim_data_strs["dvs"]["names"])...)
    dv_deps = tuple(map(x->NamedTuple{tuple(Symbol.(x)...)}(tuple((1:length(x))...)), sim_data_strs["dvs"]["deps"])...)
    dv_data = tuple(map(recursive_array_to_array, sim_data_strs["dvs"]["data"])...) 
    dv_data_nt = NamedTuple{dv_syms}(dv_data)
    dv_deps_nt = NamedTuple{dv_syms}(dv_deps)
    sim_data = (ivs = iv_nt, dvs = dv_data_nt, dv_deps = dv_deps_nt)

end
read_sim_data(model::AbstractPyBaMMModel) = read_sim_data(get_model_dir(model))

# this is required to make the symbolic functions get broadcasted on ifelse functions, which we need for 
# neuralpde loss functions.  SymbolicUtils.Code uses a cond ? true : false operator for some reason which 
# can't get broadcasted the same way 
function SymbolicUtils.Code.function_to_expr(::typeof(IfElse.ifelse), O, st::SymbolicUtils.Code.LazyState)
    #@show O
    args = arguments(O)
    :(IfElse.ifelse($(toexpr(args[1], st)), $(toexpr(args[2], st)), $(toexpr(args[3], st))))
end

export ty, fn, fnty
export generate_sim_model, read_sim_data, pybamm_func_str, load_model, get_model_dir
export AbstractPyBaMMModel, SPMModel, SPMeModel
export AbstractApplyFuncType, ParameterizedMatrixApplyFuncType, VectorOfParameterizedMDFApplyFuncType
export MultiDimensionalFunction, cartesian_product
export get_eval_network_at_sim_data_func, get_cb_func, get_plot_function, do_plot


end
