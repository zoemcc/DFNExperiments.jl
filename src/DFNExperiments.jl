module DFNExperiments

using Requires
using DiffEqFlux, Flux
using Random
using Lux
using NeuralPDE, ModelingToolkit, Symbolics, DomainSets
import ModelingToolkit: Interval, infimum, supremum
using LabelledArrays
using Interpolations
using JSON
using IterTools
using IfElse
using Plots, LinearAlgebra, TensorBoardLogger
using PyCall

include("utils.jl")

abstract type AbstractApplyFuncType end

include("multi_dimensional_function.jl")

abstract type AbstractPyBaMMModel end

Base.nameof(term::Term{Real, Base.ImmutableDict{DataType, Any}}) = (nameof ∘ operation ∘ Symbolics.value)(term)

pybamm_func_str(::AbstractPyBaMMModel) = throw("Not implemented")
include_q_model(::AbstractPyBaMMModel) = throw("Not implemented")
num_pts_validation(::AbstractPyBaMMModel) = throw("Not implemented")
large_interp_grid_length_validation(::AbstractPyBaMMModel) = throw("Not implemented")
small_interp_grid_length_validation(::AbstractPyBaMMModel) = throw("Not implemented")
num_stochastic_samples_from_loss_validation(::AbstractPyBaMMModel) = throw("Not implemented")

get_model_dir(model::AbstractPyBaMMModel) = abspath(joinpath(@__DIR__, "..", "models", "$(pybamm_func_str(model))"))
function load_model(model::AbstractPyBaMMModel)
    model_dir = get_model_dir(model)
    include(joinpath(model_dir, "model.jl"))
    pde_system
end

struct SPMnoRModel <: AbstractPyBaMMModel
end

pybamm_func_str(::SPMnoRModel) = "spm_no_r"
include_q_model(::SPMnoRModel) = false
num_pts_validation(::SPMnoRModel) = 4000
large_interp_grid_length_validation(::SPMnoRModel) = 1000
small_interp_grid_length_validation(::SPMnoRModel) = 400
num_stochastic_samples_from_loss_validation(::SPMnoRModel) = 4096

struct SPMModel <: AbstractPyBaMMModel
end

pybamm_func_str(::SPMModel) = "spm"
include_q_model(::SPMModel) = false
num_pts_validation(::SPMModel) = 4000
large_interp_grid_length_validation(::SPMModel) = 1000
small_interp_grid_length_validation(::SPMModel) = 400
num_stochastic_samples_from_loss_validation(::SPMModel) = 4096


struct ReducedCModel <: AbstractPyBaMMModel
end

pybamm_func_str(::ReducedCModel) = "reduced_c"
include_q_model(::ReducedCModel) = false
num_pts_validation(::ReducedCModel) = 4000
large_interp_grid_length_validation(::ReducedCModel) = 1000
small_interp_grid_length_validation(::ReducedCModel) = 400
num_stochastic_samples_from_loss_validation(::ReducedCModel) = 4096


struct SPMeModel <: AbstractPyBaMMModel
end

pybamm_func_str(::SPMeModel) = "spme"
include_q_model(::SPMeModel) = false
num_pts_validation(::SPMeModel) = 4000
large_interp_grid_length_validation(::SPMeModel) = 1000
small_interp_grid_length_validation(::SPMeModel) = 400
num_stochastic_samples_from_loss_validation(::SPMeModel) = 4096

struct ReducedCPhiModel <: AbstractPyBaMMModel
end

pybamm_func_str(::ReducedCPhiModel) = "reduced_c_phi"
include_q_model(::ReducedCPhiModel) = false
num_pts_validation(::ReducedCPhiModel) = 4000
large_interp_grid_length_validation(::ReducedCPhiModel) = 1000
small_interp_grid_length_validation(::ReducedCPhiModel) = 400
num_stochastic_samples_from_loss_validation(::ReducedCPhiModel) = 4096

struct ReducedCPhiJModel <: AbstractPyBaMMModel
end

pybamm_func_str(::ReducedCPhiJModel) = "reduced_c_phi_j"
include_q_model(::ReducedCPhiJModel) = false
num_pts_validation(::ReducedCPhiJModel) = 4000
large_interp_grid_length_validation(::ReducedCPhiJModel) = 1000
small_interp_grid_length_validation(::ReducedCPhiJModel) = 400
num_stochastic_samples_from_loss_validation(::ReducedCPhiJModel) = 4096

struct DFNnoRModel <: AbstractPyBaMMModel
end

pybamm_func_str(::DFNnoRModel) = "dfn_no_r"
include_q_model(::DFNnoRModel) = false
num_pts_validation(::DFNnoRModel) = 1000
large_interp_grid_length_validation(::DFNnoRModel) = 600
small_interp_grid_length_validation(::DFNnoRModel) = 400
num_stochastic_samples_from_loss_validation(::DFNnoRModel) = 4096

struct DFNModel <: AbstractPyBaMMModel
end

pybamm_func_str(::DFNModel) = "dfn"
include_q_model(::DFNModel) = false
num_pts_validation(::DFNModel) = 100
large_interp_grid_length_validation(::DFNModel) = 60
small_interp_grid_length_validation(::DFNModel) = 30
num_stochastic_samples_from_loss_validation(::DFNModel) = 4096


include("generate_py.jl")

function __init__()
    #@require PyCall="438e738f-606a-5dbb-bf0a-cddfbfd45ab0" begin 
    begin
        initialize_pybamm_funcs()
    end
end

include("plot_log.jl")

"""
    AbstractExplicitLayer

Abstract Type for all Lux Layers

Users implementing their custom layer, **must** implement

  - `initialparameters(rng::AbstractRNG, layer::CustomAbstractExplicitLayer)` -- This
    returns a `NamedTuple` containing the trainable parameters for the layer.
  - `initialstates(rng::AbstractRNG, layer::CustomAbstractExplicitLayer)` -- This returns a
    NamedTuple containing the current state for the layer. For most layers this is typically
    empty. Layers that would potentially contain this include `BatchNorm`, `LSTM`, `GRU` etc.

Optionally:

  - `parameterlength(layer::CustomAbstractExplicitLayer)` -- These can be automatically
    calculated, but it is recommended that the user defines these.
  - `statelength(layer::CustomAbstractExplicitLayer)` -- These can be automatically
    calculated, but it is recommended that the user defines these.

See also [`AbstractExplicitContainerLayer`](@ref)
"""

struct LuxChainInterpolator{I} <: Lux.AbstractExplicitLayer
    interpolator::I
end

Lux.initialparameters(::Random.AbstractRNG, ::LuxChainInterpolator) = (int=[1.0], )
Lux.initialstates(::Random.AbstractRNG, ::LuxChainInterpolator) = NamedTuple()
Lux.parameterlength(::LuxChainInterpolator) = 1
Lux.statelength(::LuxChainInterpolator) = 0

function (lux_chain_interpolator::LuxChainInterpolator)(x::AbstractArray, ps, st::NamedTuple)
    if !(typeof(lux_chain_interpolator.interpolator) <: AbstractInterpolation) && size(x, 1) > 2 
        # this uses RegularGridInterpolator which is an ND interpolator so it expects a vector in
        #inputs = [@view x[:, i] for i in 1:size(x, 2)]
        inputs = [reverse(x[:, i]) for i in 1:size(x, 2)]
        #@show inputs
        #@show typeof(inputs), size(inputs)
        output = lux_chain_interpolator.interpolator.(inputs)
        #@show typeof(output), size(output)
        #@show output
        #@show lux_chain_interpolator.interpolator.grid
        output
    else
        # this uses 1-d or 2-d interpolators which expect a tuple of arguments
        inputs = [@view x[[i], :] for i in 1:size(x, 1)]
        #@show typeof(inputs), size(inputs)
        output = lux_chain_interpolator.interpolator.(inputs...)
        #@show typeof(output), size(output)
        output
    end
    #@show output
    (output, st)
end

struct FastChainInterpolator{I}
    interpolator::I
end
function (fastchaininterpolator::FastChainInterpolator)(v, θ) 
    #println("in fastchaininterpolator")
    #@show v
    if size(v, 1) > 2
        # this uses RegularGridInterpolator which is an ND interpolator so it expects a vector in
        inputs = [@view v[:, i] for i in 1:size(v, 2)]
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

function generate_sim_model_and_test(model::M; current_input=false, include_q=false, output_dir=nothing, num_pts=100, 
        large_interp_grid_length=1000, small_interp_grid_length=100, num_stochastic_samples_from_loss=1024, writemodel=true, writesimdata=true) where {M <: AbstractPyBaMMModel}
    model_str = pybamm_func_str(model)
    current_input_str = current_input ? "True" : "False" 
    include_q_str = include_q ? "True" : "False"
    sim, mtk_str, variables = py"solve_plot_generate(*$$(model_str)(), current_input=$$(current_input_str), include_q=$$(include_q_str), num_pts=$$(num_pts))"

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

    dvs_interpolation_luxchain_datablob = map(dep_vars) do dv
        println("interpolating $(nameof(dv))")
        extraprints = nameof(dv) == :c_s_n
        @show extraprints
        dv_pybamm_name = dependent_variables_to_pybamm_names[nameof(dv)]
        @show dv_pybamm_name
        dv_processed = sim.solution.__getitem__(dv_pybamm_name)
        dv_pybamm_interpolation_function = dv_processed._interpolation_function
        @show dv_pybamm_interpolation_function
        deps = dependent_variables_to_dependencies[nameof(dv)]
        num_deps = length(deps)
        if extraprints @show num_deps end
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
        @show size(unscaled_iv_ranges)

        unscaled_ivs_mat = reduce(hcat, map(collect, vec(collect(product(unscaled_iv_ranges...)))))
        @show size(unscaled_ivs_mat)
        dv_pybamm_interpolator = LuxChainInterpolator(dv_pybamm_interpolation_function)
        dv_mat = dv_pybamm_interpolator(unscaled_ivs_mat, NamedTuple(), NamedTuple())[1]
        @show size(dv_mat), typeof(dv_mat)
        #@show dv_mat
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
        extrapolate_interpolate_function = extrapolate(scale(interpolate(dv_grid, BSpline(Cubic(Line(OnGrid())))), scaled_iv_ranges...), Line())
        dv_interpolation = LuxChainInterpolator(extrapolate_interpolate_function)
        #@show typeof(dv_interpolation)

        dv_luxchain = Lux.Chain(dv_interpolation)
        unscaled_iv_ranges_small, scaled_iv_ranges_small = iv_ranges_from_grid_length(small_interp_grid_length)
        scaled_ivs_small_mat = reduce(hcat, map(collect, vec(collect(product(scaled_iv_ranges_small...)))))
        dv_data_mat = dv_luxchain(scaled_ivs_small_mat, NamedTuple(), NamedTuple())[1]
        #@show dv_data_mat
        dv_data_blob = reshape(dv_data_mat, (length.(scaled_iv_ranges_small))...)

        (dv_interpolation, dv_luxchain, dv_data_blob)
    end

    ivs_range = map(pde_system.domain) do var_domain
        domain = var_domain.domain
        iv_range = range(domain, small_interp_grid_length)
    end
    @show ivs_range
    @show typeof(ivs_range)

    dvs_interpolation, dvs_luxchain, dvs_data_blob = unzip(dvs_interpolation_luxchain_datablob)
    @show dvs_interpolation
    @show dvs_luxchain
    @show size(dvs_data_blob)
    @show size(dvs_data_blob[1])
    #@show dvs_data_blob[1]
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
    #@show sim_data_nt
    #csn_data = sim_data_nt[:dvs][:c_s_n]
    #@show csn_data
    #@show size(csn_data)

    #cen_data = sim_data_nt[:dvs][:c_e_n]
    #@show cen_data
    #@show size(cen_data)



    println("generating PINN discretization from the simulation interpolation to test against")
    strategy =  NeuralPDE.StochasticTraining(num_stochastic_samples_from_loss, num_stochastic_samples_from_loss)
    discretization = NeuralPDE.PhysicsInformedNN(dvs_luxchain, strategy)

    pinnrep = NeuralPDE.symbolic_discretize(pde_system, discretization)
    flat_init_params = pinnrep.flat_init_params
    @show flat_init_params

    eq_fake = pde_system.dvs[1] ~ pde_system.dvs[1]
    ic_bc_fake = pde_system.dvs[1] ~ pde_system.dvs[1]

    eq_loss_certs = []
    for (i, eq) in enumerate(eqs)
        println("generating loss for equation $i")
        @named single_eq_pde_system = PDESystem([eq], [ic_bc_fake], domains, ind_vars, dep_vars)
        prob = NeuralPDE.discretize(single_eq_pde_system, discretization)
        eq_loss = prob.f(flat_init_params, Float64[])
        push!(eq_loss_certs, (i, eq, eq_loss))
    end

    ic_bc_loss_certs = []
    for (i, ic_bc) in enumerate(ics_bcs)
        println("generating loss for ic_bc $i")
        @named single_ic_bc_pde_system = PDESystem([eq_fake], [ic_bc], domains, ind_vars, dep_vars)
        prob = NeuralPDE.discretize(single_ic_bc_pde_system, discretization)
        ic_bc_loss = prob.f(flat_init_params, Float64[])
        push!(ic_bc_loss_certs, (i, ic_bc, ic_bc_loss))
    end

    prob = NeuralPDE.discretize(pde_system, discretization)
    total_loss = prob.f(flat_init_params, Float64[])
    #@show total_loss

    loss_cert_string = "total_loss: $total_loss\n\n" * "eq losses:\n\n" * 
    join(map(eq_loss_certs) do (i, eq, loss)
        string("eq $i:\n", eq, "\nloss: ", loss)
    end, "\n\n") * "\n\nic bc losses:\n\n" * join(map(ic_bc_loss_certs) do (i, ic_bc, loss)
        string("ic_bc $i:\n", ic_bc, "\nloss: ", loss)
    end, "\n\n") * "\n"

    # write loss_cert_string to file 
    if typeof(output_dir) <: AbstractString
        loss_cert_filename = joinpath(output_dir, "individual_loss_cert.txt")
    else
        loss_cert_filename = Base.tempname()
    end
    open(loss_cert_filename, "w") do f
        print(f, loss_cert_string)
    end


    #@show prob
    #@show typeof(prob)
    #@show prob.f



    (sim_data=sim_data_nt, pde_system=pde_system, sim=sim, variables=variables, 
    independent_variables_to_pybamm_names=independent_variables_to_pybamm_names, 
    dependent_variables_to_pybamm_names=dependent_variables_to_pybamm_names,
    dependent_variables_to_dependencies=dependent_variables_to_dependencies,
    dvs_interpolation=dvs_interpolation,
    dvs_luxchain=dvs_luxchain,
    prob=prob,
    total_loss=total_loss,
    pinnrep=pinnrep,
    discretization=discretization)
    #"""
    #(sim_data=sim_data_nt, pde_system=pde_system, sim=sim, variables=variables, 
    #independent_variables_to_pybamm_names=independent_variables_to_pybamm_names, 
    #dependent_variables_to_pybamm_names=dependent_variables_to_pybamm_names,
    #dependent_variables_to_dependencies=dependent_variables_to_dependencies,
    #dvs_interpolation=dvs_interpolation,
    #dvs_luxchain=dvs_luxchain,)
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
export num_pts_validation, large_interp_grid_length_validation, small_interp_grid_length_validation, num_stochastic_samples_from_loss_validation


end
