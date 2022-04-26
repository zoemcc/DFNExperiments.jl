module DFNExperiments

using DiffEqFlux, Flux
using NeuralPDE, ModelingToolkit, Symbolics, DomainSets
import ModelingToolkit: Interval, infimum, supremum
using LabelledArrays
using Interpolations
using JSON
using PyCall
using IfElse
using Plots, LinearAlgebra, TensorBoardLogger

include("utils.jl")

abstract type AbstractApplyFuncType end

include("multi_dimensional_function.jl")

abstract type AbstractPyBaMMModel end

Base.nameof(term::Term{Real, Base.ImmutableDict{DataType, Any}}) = (nameof ∘ operation ∘ Symbolics.value)(term)

pybamm_func_str(::AbstractPyBaMMModel) = throw("Not implemented")

get_model_dir(model::AbstractPyBaMMModel) = abspath(joinpath(@__DIR__, "..", "models", "$(pybamm_func_str(model))"))
function load_model(model::AbstractPyBaMMModel)
    model_dir = get_model_dir(model)
    include(joinpath(model_dir, "model.jl"))
    pde_system
end

struct SPMnoRModel <: AbstractPyBaMMModel
end

pybamm_func_str(spm::SPMnoRModel) = "spm_no_r"

struct SPMModel <: AbstractPyBaMMModel
end

pybamm_func_str(spm::SPMModel) = "spm"


struct ReducedCModel <: AbstractPyBaMMModel
end

pybamm_func_str(spm::ReducedCModel) = "reduced_c"


struct SPMeModel <: AbstractPyBaMMModel
end

pybamm_func_str(spm::SPMeModel) = "spme"

struct ReducedCPhiModel <: AbstractPyBaMMModel
end

pybamm_func_str(spm::ReducedCPhiModel) = "reduced_c_phi"

struct ReducedCPhiJModel <: AbstractPyBaMMModel
end

pybamm_func_str(spm::ReducedCPhiJModel) = "reduced_c_phi_j"

struct DFNnoRModel <: AbstractPyBaMMModel
end

pybamm_func_str(spm::DFNnoRModel) = "dfn_no_r"

struct DFNModel <: AbstractPyBaMMModel
end

pybamm_func_str(spm::DFNModel) = "dfn"


include("generate_py.jl")

function __init__()
    initialize_pybamm_funcs()
end

include("plot_log.jl")

struct FastChainInterpolator{I}
    interpolator::I
end
function (fastchaininterpolator::FastChainInterpolator)(v, θ) 
    println("in fastchaininterpolator")
    @show v
    inputs = [@view v[[i], :] for i in 1:size(v, 1)]
    fastchaininterpolator.interpolator.(inputs...)
    #reshape(fastchaininterpolator.interpolator.(inputs...), 1, size(v, 2))
end
DiffEqFlux.initial_params(::FastChainInterpolator) = Float64[]


function generate_sim_model_and_test(model::M; current_input=false, output_dir=nothing, num_pts=100) where {M <: AbstractPyBaMMModel}
    model_str = pybamm_func_str(model)
    current_input_str = current_input ? "True" : "False" 
    sim, mtk_str, variables = py"solve_plot_generate(*$$(model_str)(), current_input=$$(current_input_str), num_pts=$$(num_pts))"

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
    dv_names = nameof.(pde_system.dvs)
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

    dvs_interpolation_fastchain = map(dep_vars) do dv
        dv_pybamm_name = dependent_variables_to_pybamm_names[nameof(dv)]
        dv_processed = sim.solution.__getitem__(dv_pybamm_name)
        dv_pybamm_interpolation_function = dv_processed._interpolation_function
        deps = dependent_variables_to_dependencies[nameof(dv)]
        num_deps = length(deps)
        py_axis_name = (:x, :y, :z)
        iv_ranges = map(1:num_deps) do i
            iv = deps[i]
            iv_pybamm_name = independent_variables_to_pybamm_names[iv]
            iv_axis_name = py_axis_name[i]
            iv_grid = dv_pybamm_interpolation_function[iv_axis_name]
            # assume time is first
            iv_scale = i == 1 ? dv_processed.timescale : dv_processed.length_scales[iv_pybamm_name]
            range(Interval(iv_grid[1] / iv_scale, iv_grid[end] / iv_scale), length(iv_grid))
        end
        dv_grid = permutedims(reshape(dv_pybamm_interpolation_function[py_axis_name[num_deps + 1]], reverse(length.(iv_ranges))...), reverse(1:num_deps))
        dv_interpolation = FastChainInterpolator(scale(interpolate(dv_grid, BSpline(Cubic(Line(OnGrid())))), iv_ranges...))

        dv_fastchain = FastChain(dv_interpolation)

        (dv_interpolation, dv_fastchain)
    end
    dvs_interpolation = map(first, dvs_interpolation_fastchain)
    dvs_fastchain = map(x->x[2], dvs_interpolation_fastchain)

    @named modded_pde_system = PDESystem(pde_system.eqs[[4]], pde_system.bcs[[1]], pde_system.domain, pde_system.ivs, pde_system.dvs)
    strategy =  NeuralPDE.StochasticTraining(8, 8)
    discretization = NeuralPDE.PhysicsInformedNN(dvs_fastchain, strategy)
    prob = NeuralPDE.discretize(modded_pde_system, discretization)
    symb_modded_pde_system = NeuralPDE.symbolic_discretize(modded_pde_system, discretization)
    total_loss = prob.f(Float64[], Float64[])

    (sim_data=sim_data_nt, pde_system=pde_system, sim=sim, variables=variables, 
    independent_variables_to_pybamm_names=independent_variables_to_pybamm_names, 
    dependent_variables_to_pybamm_names=dependent_variables_to_pybamm_names,
    dependent_variables_to_dependencies=dependent_variables_to_dependencies,
    dvs_interpolation=dvs_interpolation,
    dvs_fastchain=dvs_fastchain,
    prob=prob,
    total_loss=total_loss,
    modded_pde_system=modded_pde_system,
    symb_modded_pde_system=symb_modded_pde_system)
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
export generate_sim_model_and_test, read_sim_data, pybamm_func_str, load_model, get_model_dir
export AbstractPyBaMMModel
export AbstractApplyFuncType, ParameterizedMatrixApplyFuncType, VectorOfParameterizedMDFApplyFuncType
export MultiDimensionalFunction, cartesian_product
export get_eval_network_at_sim_data_func, get_cb_func, get_plot_function, do_plot
export SPMnoRModel, SPMModel, ReducedCModel, SPMeModel, ReducedCPhiModel, ReducedCPhiJModel, DFNnoRModel, DFNModel


end
