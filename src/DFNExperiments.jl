module DFNExperiments

using Requires
using DiffEqFlux, Flux
using NeuralPDE, ModelingToolkit, Symbolics, DomainSets
import ModelingToolkit: Interval, infimum, supremum
using LabelledArrays
using Interpolations
using JSON
using IterTools
using IfElse
using Plots, LinearAlgebra, TensorBoardLogger

include("utils.jl")

abstract type AbstractApplyFuncType end

include("multi_dimensional_function.jl")

abstract type AbstractPyBaMMModel end

Base.nameof(term::Term{Real, Base.ImmutableDict{DataType, Any}}) = (nameof ∘ operation ∘ Symbolics.value)(term)

pybamm_func_str(::AbstractPyBaMMModel) = throw("Not implemented")
include_q_model(::AbstractPyBaMMModel) = throw("Not implemented")

get_model_dir(model::AbstractPyBaMMModel) = abspath(joinpath(@__DIR__, "..", "models", "$(pybamm_func_str(model))"))
function load_model(model::AbstractPyBaMMModel)
    model_dir = get_model_dir(model)
    include(joinpath(model_dir, "model.jl"))
    pde_system
end

struct SPMnoRModel <: AbstractPyBaMMModel
end

pybamm_func_str(::SPMnoRModel) = "spm_no_r"
include_q_model(::SPMnoRModel) = true

struct SPMModel <: AbstractPyBaMMModel
end

pybamm_func_str(::SPMModel) = "spm"
include_q_model(::SPMModel) = true


struct ReducedCModel <: AbstractPyBaMMModel
end

pybamm_func_str(::ReducedCModel) = "reduced_c"
include_q_model(::ReducedCModel) = false


struct SPMeModel <: AbstractPyBaMMModel
end

pybamm_func_str(::SPMeModel) = "spme"
include_q_model(::SPMeModel) = true

struct ReducedCPhiModel <: AbstractPyBaMMModel
end

pybamm_func_str(::ReducedCPhiModel) = "reduced_c_phi"
include_q_model(::ReducedCPhiModel) = false

struct ReducedCPhiJModel <: AbstractPyBaMMModel
end

pybamm_func_str(::ReducedCPhiJModel) = "reduced_c_phi_j"
include_q_model(::ReducedCPhiJModel) = false

struct DFNnoRModel <: AbstractPyBaMMModel
end

pybamm_func_str(::DFNnoRModel) = "dfn_no_r"
include_q_model(::DFNnoRModel) = true

struct DFNModel <: AbstractPyBaMMModel
end

pybamm_func_str(::DFNModel) = "dfn"
include_q_model(::DFNModel) = true


include("generate_py.jl")

function __init__()
    #@require PyCall="438e738f-606a-5dbb-bf0a-cddfbfd45ab0" begin 
    begin
        initialize_pybamm_funcs()
    end
end

include("plot_log.jl")

struct FastChainInterpolator{I}
    interpolator::I
end
function (fastchaininterpolator::FastChainInterpolator)(v, θ) 
    #println("in fastchaininterpolator")
    #@show v
    if size(v, 1) > 2
        # this uses RegularGridInterpolator which is an ND interpolator so it expects a vector in
        inputs = [v[:, i] for i in 1:size(v, 2)]
        output = fastchaininterpolator.interpolator.(inputs)
    else
        # this uses 1-d or 2-d interpolators which expect a tuple of arguments
        inputs = [@view v[[i], :] for i in 1:size(v, 1)]
        output = fastchaininterpolator.interpolator.(inputs...)
    end
    #@show output
    output
end
DiffEqFlux.initial_params(::FastChainInterpolator) = Float64[]

function generate_sim_model_and_test(model::M; current_input=false, include_q=true, output_dir=nothing, num_pts=100, 
        large_interp_grid_length=1000, small_interp_grid_length=100, num_stochastic_samples_from_loss=1024, writemodel=true) where {M <: AbstractPyBaMMModel}
    throw("PyCall and PyBaMM need to be installed for this function to run")
    model_str = pybamm_func_str(model)
    current_input_str = current_input ? "True" : "False" 
    include_q_str = include_q ? "True" : "False"
    #sim, mtk_str, variables = py"solve_plot_generate(*$$(model_str)(), current_input=$$(current_input_str), include_q=$$(include_q_str), num_pts=$$(num_pts))"

    if typeof(output_dir) <: AbstractString
        if !isdir(output_dir)
            mkpath(output_dir)
        end
        model_filename = joinpath(output_dir, "model.jl")
    else
        model_filename = Base.tempname()
    end

    @show writemodel
    if writemodel
        open(model_filename, "w") do f
            write(f, mtk_str)
        end
    end
    include(model_filename)

    #sol_data_json = sim.solution.save_data(variables=variables, to_format="json")
    #sol_data = JSON.parse(sol_data_json)

    iv_names = nameof.(pde_system.ivs)
    @show iv_names
    @show independent_variables_to_pybamm_names
    dv_names = nameof.(pde_system.dvs)
    @show dv_names
    @show dependent_variables_to_pybamm_names
    #@show keys(sol_data)

    dvs_interpolation_fastchain_datablob = map(dep_vars) do dv
        println("interpolating $(nameof(dv))")
        dv_pybamm_name = dependent_variables_to_pybamm_names[nameof(dv)]
        dv_processed = sim.solution.__getitem__(dv_pybamm_name)
        dv_pybamm_interpolation_function = dv_processed._interpolation_function
        deps = dependent_variables_to_dependencies[nameof(dv)]
        num_deps = length(deps)
        py_axis_name = (:x, :y, :z)
        function iv_ranges_from_grid_length(grid_length)
            iv_ranges = map(1:num_deps) do i
                iv = deps[i]
                iv_pybamm_name = independent_variables_to_pybamm_names[iv]
                if num_deps > 2
                    grid_indices = reverse(collect(1:num_deps))
                    iv_grid = dv_pybamm_interpolation_function.grid[grid_indices[i]]
                else
                    iv_axis_name = py_axis_name[i]
                    iv_grid = dv_pybamm_interpolation_function[iv_axis_name]
                end
                # assume time is first
                iv_scale = i == 1 ? dv_processed.timescale : dv_processed.length_scales[iv_pybamm_name]
                unscaled_iv_range = range(Interval(iv_grid[1], iv_grid[end]), grid_length)
                scaled_iv_range = range(Interval(iv_grid[1] / iv_scale, iv_grid[end] / iv_scale), grid_length)
                (unscaled_iv_range, scaled_iv_range)
            end
            unscaled_iv_ranges, scaled_iv_ranges = unzip(iv_ranges)
        end
        unscaled_iv_ranges, scaled_iv_ranges = iv_ranges_from_grid_length(large_interp_grid_length)

        unscaled_ivs_mat = reduce(hcat, map(collect, vec(collect(product(unscaled_iv_ranges...)))))
        dv_pybamm_interpolator = FastChainInterpolator(dv_pybamm_interpolation_function)
        dv_mat = dv_pybamm_interpolator(unscaled_ivs_mat, Float64[])
        @show size(dv_mat), typeof(dv_mat), size(dv_mat[1])
        dv_mat = first.(dv_mat)
        #if num_deps <= 2
            #dv_mat = first.(dv_mat)
        #else
            #dv_mat = first(dv_mat)
        #end
        dv_grid = reshape(dv_mat, fill(large_interp_grid_length, num_deps)...)
        @show size(dv_grid)
        @show size.(scaled_iv_ranges)
        #dv_grid = permutedims(reshape(dv_pybamm_interpolation_function[py_axis_name[num_deps + 1]], reverse(length.(iv_ranges))...), reverse(1:num_deps))
        dv_interpolation = FastChainInterpolator(extrapolate(scale(interpolate(dv_grid, BSpline(Cubic(Line(OnGrid())))), scaled_iv_ranges...), Line()))
        @show typeof(dv_interpolation)

        dv_fastchain = FastChain(dv_interpolation)
        unscaled_iv_ranges_small, scaled_iv_ranges_small = iv_ranges_from_grid_length(small_interp_grid_length)
        scaled_ivs_small_mat = reduce(hcat, map(collect, vec(collect(product(scaled_iv_ranges_small...)))))
        dv_data_blob = reshape(dv_fastchain(scaled_ivs_small_mat, Float64[]), (length.(scaled_iv_ranges_small))...)

        (dv_interpolation, dv_fastchain, dv_data_blob)
    end

    ivs_range = map(pde_system.domain) do var_domain
        domain = var_domain.domain
        iv_range = range(domain, small_interp_grid_length)
    end
    @show ivs_range
    @show typeof(ivs_range)

    dvs_interpolation, dvs_fastchain, dvs_data_blob = unzip(dvs_interpolation_fastchain_datablob)
    @show size(dvs_data_blob)
    @show typeof(dvs_data_blob)

    sim_data = Dict(
        :ivs => Dict(
            :names => iv_names,
            :data => collect.(ivs_range),
        ),
        :dvs => Dict(
            :names => dv_names,
            :deps => [collect(dependent_variables_to_dependencies[dv_name]) for dv_name in dv_names], 
            :data => dvs_data_blob,
        ) 
    )

    if typeof(output_dir) <: AbstractString
        sim_filename = joinpath(output_dir, "sim.json")
    else
        sim_filename = Base.tempname()
    end

    @show writesimdata
    if writesimdata
        open(sim_filename, "w") do f
            JSON.print(f, sim_data)
        end
    end
    sim_data_nt = read_sim_data(sim_filename)

    println("generating PINN discretization from the simulation interpolation to test against")
    strategy =  NeuralPDE.StochasticTraining(num_stochastic_samples_from_loss, num_stochastic_samples_from_loss)
    discretization = NeuralPDE.PhysicsInformedNN(dvs_fastchain, strategy)
    
    if !isdefined(DFNExperiments, :subdomain_relations)
        DFNExperiments.subdomain_relations = nothing
    end
    if !isdefined(DFNExperiments, :eqs_integration_domains)
        DFNExperiments.eqs_integration_domains = nothing
    end
    if !isdefined(DFNExperiments, :ics_bcs_integration_domains)
        DFNExperiments.ics_bcs_integration_domains = nothing
    end

    symb_modded_pde_system = NeuralPDE.symbolic_discretize(pde_system, discretization; 
        subdomain_relations=DFNExperiments.subdomain_relations, 
        eqs_integration_domains=DFNExperiments.eqs_integration_domains,
        ics_bcs_integration_domains=DFNExperiments.ics_bcs_integration_domains,
    )
    prob = NeuralPDE.discretize(pde_system, discretization;
        subdomain_relations=DFNExperiments.subdomain_relations, 
        eqs_integration_domains=DFNExperiments.eqs_integration_domains,
        ics_bcs_integration_domains=DFNExperiments.ics_bcs_integration_domains,
    )
    @show DFNExperiments.subdomain_relations
    @show DFNExperiments.eqs_integration_domains
    @show DFNExperiments.ics_bcs_integration_domains

    total_loss = prob.f(Float64[], Float64[])


    (sim_data=sim_data_nt, pde_system=pde_system, sim=sim, variables=variables, 
    independent_variables_to_pybamm_names=independent_variables_to_pybamm_names, 
    dependent_variables_to_pybamm_names=dependent_variables_to_pybamm_names,
    dependent_variables_to_dependencies=dependent_variables_to_dependencies,
    dvs_interpolation=dvs_interpolation,
    dvs_fastchain=dvs_fastchain,
    prob=prob,
    total_loss=total_loss,
    symb_modded_pde_system=symb_modded_pde_system,
    discretization=discretization)
    #"""
    #(sim_data=sim_data_nt, pde_system=pde_system, sim=sim, variables=variables, 
    #independent_variables_to_pybamm_names=independent_variables_to_pybamm_names, 
    #dependent_variables_to_pybamm_names=dependent_variables_to_pybamm_names,
    #dependent_variables_to_dependencies=dependent_variables_to_dependencies,
    #dvs_interpolation=dvs_interpolation,
    #dvs_fastchain=dvs_fastchain,)
    #"""
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
