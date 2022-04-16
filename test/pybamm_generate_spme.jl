begin
    using Distributed
    NUM_WORKERS = 0
    if NUM_WORKERS > 0
        if nworkers() != NUM_WORKERS
            addprocs(NUM_WORKERS - nworkers())
        end
    end
end
@everywhere begin 
    using Pkg
    Pkg.activate(abspath(joinpath(@__DIR__, "..")))
end
@everywhere begin 
    using DFNExperiments
    using ModelingToolkit, Symbolics, DomainSets
    import ModelingToolkit: Interval, infimum, supremum
    using LabelledArrays
    using NeuralPDE
    using PyCall
    using JSON
    using Random
    using TensorBoardLogger
    using Flux, DiffEqFlux
    using LinearAlgebra
    using Plots
    using GalacticOptim
    using Printf
    using IfElse
end


begin 
    model = SPMeModel()
    seed = 1
    num_experiments = 120
    start_experiment = 1
    plot_frequency = 5_00
    log_frequency = 1_00
    iterations = 1_000

    models_dir = abspath(joinpath(@__DIR__, "..", "models"))
    model_str = pybamm_func_str(model)
    output_dir = joinpath(models_dir, "$(model_str)")
end

@everywhere function eval_network_at_sim_data_func(chains)
    sim_data = read_sim_data(output_dir)
    ivs = keys(sim_data[:ivs])
    dvs = keys(sim_data[:dvs])
    dv_deps = collect(map(keys, sim_data[:dv_deps]))
    num_dvs = length(dvs)

    ground_truths = map(dv->sim_data[:dvs][dv], dvs)
    ground_truths_norm = LinearAlgebra.norm.(ground_truths)

    mdf = MultiDimensionalFunction(chains, VectorOfParameterizedMDFApplyFuncType(), ivs, dvs, dv_deps)

    function eval_network_at_sim_data(p)
        network_evals = mdf(p, values(sim_data[:ivs])...; flat=Val{false})

        errors = map(i -> network_evals[i] .- ground_truths[i], 1:num_dvs)
        errors_norm = LinearAlgebra.norm.(errors)
        errors_rel = map(i -> errors_norm[i] / ground_truths_norm[i], 1:num_dvs)

        return (dvs=dvs, network_evals=network_evals, errors=errors,
                errors_norm=errors_norm, errors_rel=errors_rel, )
    end

end

@everywhere function get_cb(logger, iteration, chains)
    eval_network_at_sim_data = eval_network_at_sim_data_func(chains)
    cb = function (p,l)
        if iteration[1] % log_frequency == 0
            res = eval_network_at_sim_data(p)
            error_rel_max = max(res[:errors_rel]...)
            log_value(logger, "errors/max_error_relative", error_rel_max; step=iteration[1])
            for (i, dv) in enumerate(res[:dvs])
                log_value(logger, "errors/$(string(dv))_error_norm", res[:errors_norm][i]; step=iteration[1])
                log_value(logger, "errors/$(string(dv))_error_relative", res[:errors_rel][i]; step=iteration[1])
            end
        end
        return false
    end
    return cb
end

function do_plot(data_array::AbstractArray{R, N}, iv_data, deps, name) where {R, N}
    used_ivs = map(iv->iv_data[iv], deps)
    if N == 1
        plot(used_ivs..., data_array, title=name)
    elseif N == 2
        plot(reverse(used_ivs)..., data_array, linetype=:contourf, title=name)
    else
        println("no plot defined for arrays with more than 2 axes")
        nothing
    end
end

@everywhere function get_plot_function(logger, iteration, chains)
    sim_data = read_sim_data(output_dir)
    dvs = keys(sim_data[:dvs])
    dv_deps = collect(map(keys, sim_data[:dv_deps]))
    num_dvs = length(dvs)
    iv_data = sim_data[:ivs]
    pybamm_plots = map(i->(name="$(string(dvs[i]))_pybamm", image=do_plot(sim_data[:dvs][i], iv_data, dv_deps[i], "$(string(dvs[i]))_pybamm")), 1:num_dvs)

    have_given_pybamm_plots = [false]

    eval_network_at_sim_data = eval_network_at_sim_data_func(chains)

    function plot_function(θ, adaloss)
        res = eval_network_at_sim_data(θ)
        eval_plots = map(i->(name="$(string(dvs[i]))_nn", image=do_plot(res[:network_evals][i], iv_data, dv_deps[i], "$(string(dvs[i]))_nn")), 1:num_dvs)
        error_plots = map(i->(name="$(string(dvs[i]))_error", image=do_plot(res[:errors][i], iv_data, dv_deps[i], "$(string(dvs[i]))_error")), 1:num_dvs)

        if !(have_given_pybamm_plots[1])
            have_given_pybamm_plots[1] = true
            vcat(pybamm_plots, eval_plots, error_plots)
        else
            vcat(eval_plots, error_plots)
        end
    end
    return plot_function
end

function main(model, seed, num_experiments, start_experiment)
    Random.seed!(seed)

    models_dir = abspath(joinpath(@__DIR__, "..", "models"))
    model_str = pybamm_func_str(model)
    output_dir = joinpath(models_dir, "$(model_str)")
    sim_file = joinpath(output_dir, "sim.json")
    model_file = joinpath(output_dir, "model.jl")
    include(joinpath(output_dir, "model.jl"))
    sim_data = read_sim_data(output_dir)

    sg = StructGenerator(
        :CompositeHyperParameter,
        RandomChoice(1:2^10), # seed
        StructGenerator( # nn
            :SimpleFeedForwardNetwork, # type/constructor name
            RandomChoice(3:5),
            RandomChoice(64, 128),
            RandomChoice(:GELUNonLin, :SigmoidNonLin),
            :GlorotUniformParams
        ),
        StructGenerator( # training
            :StochasticTraining,
            RandomChoice(128, 256)
        ),
        RandomChoice( # adaptive loss
            StructGenerator( 
                :MiniMaxAdaptiveLoss
            ),
            StructGenerator(
                :GradientScaleAdaptiveLoss
            )
        ),
        RandomChoice( # optimizer
            StructGenerator(:ADAMOptimiser, iterations, 3e-4),
            StructGenerator(:ADAMOptimiser, iterations, 3e-3),
            StructGenerator(:ADAMOptimiser, iterations, 1e-3)
        )
    )

    hyperparametersweep = StructGeneratorHyperParameterSweep(seed, num_experiments, sg)
    hyperparameters = generate_hyperparameters(hyperparametersweep)

    log_options = NeuralPDE.LogOptions(;plot_function=get_plot_function, log_dir=abspath(joinpath(@__DIR__, "..", "logs", model_str)), 
        log_frequency=log_frequency, plot_frequency=plot_frequency)
    neuralpde_workers = map(NeuralPDE.NeuralPDEWorker, workers())
    experiment_manager = NeuralPDE.ExperimentManager(pde_system, hyperparameters, get_cb, log_options, neuralpde_workers)
    experiment_index = start_experiment
    NeuralPDE.remote_run_neuralpde_with_logs(pde_system, hyperparameters[experiment_index], get_cb, log_options, experiment_index, false)
    #NeuralPDE.run_experiment_queue(experiment_manager; remote=false, startat=start_experiment)
end

main(model, seed, num_experiments, start_experiment)
nothing
