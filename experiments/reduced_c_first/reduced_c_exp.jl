begin
    using Pkg
    env_path = abspath(joinpath(@__DIR__, "..", ".."))
    Pkg.activate(env_path)
    using Distributed
    NUM_WORKERS = 0
    if NUM_WORKERS > 0
        if nworkers() != NUM_WORKERS
            addprocs(NUM_WORKERS - nworkers())
        end
    end
end
@everywhere begin 
    # https://github.com/JuliaPlots/Plots.jl/issues/1076
    # for headless plots
    ENV["GKSwstype"] = "100"
    using Pkg
    Pkg.activate(env_path)
end
@everywhere begin 
    using DFNExperiments, NeuralPDE
    using ModelingToolkit, Symbolics, DomainSets
    import ModelingToolkit: Interval, infimum, supremum
    using IfElse
end

begin 
    model = ReducedCModel()
    hyperseed = 1
    num_experiments = 64
    if haskey(ENV, "SLURM_ARRAY_TASK_ID")
        task_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
        num_tasks = parse(Int, ENV["SLURM_ARRAY_TASK_COUNT"])
        @show task_id
        @show num_tasks
        num_hyperparameters_per_task = Int(ceil(num_experiments/num_tasks))
        @show num_hyperparameters_per_task
        hyperparameter_indices_to_compute = range((task_id - 1) * num_hyperparameters_per_task + 1, min(task_id * num_hyperparameters_per_task, num_experiments))
        @show hyperparameter_indices_to_compute
        start_experiment = hyperparameter_indices_to_compute[1]
        endat = hyperparameter_indices_to_compute[end]
    else
        start_experiment = 1
        endat = num_experiments
    end
    plot_frequency = 5_000
    log_frequency = 1_000
    checkpoint_frequency = 5_000
    iterations = 200_000
    log_string = "reduced_c"
    if haskey(ENV, "LOG_DIR")
        log_dir = ENV["LOG_DIR"]
    else
        log_dir = abspath(joinpath(@__DIR__, "..", "logs", log_string))
    end
    distributed = false

    hyperparameter_generator = StructGenerator(
        :CompositeHyperParameter,
        RandomChoice(1:2^10), # seed
        StructGenerator( # nn
            :SimpleFeedForwardNetwork, # type/constructor name
            RandomChoice(5:8),
            RandomChoice(64, 128),
            RandomChoice(:GELUNonLin, :SigmoidNonLin),
            :GlorotUniformParams
        ),
        StructGenerator( # training
            :StochasticTraining,
            RandomChoice(256, 128)
        ),
        RandomChoice( # adaptive loss
            StructGenerator( 
                :MiniMaxAdaptiveLoss
            )
            #StructGenerator(
                #:GradientScaleAdaptiveLoss
            #)
        ),
        RandomChoice( # optimizer
            StructGenerator(:ADAMOptimiser, iterations, 3e-3),
            StructGenerator(:ADAMOptimiser, iterations, 1e-3),
            StructGenerator(:ADAMOptimiser, iterations, 3e-4)
        )
    )

end

function main(model, hyperparameter_generator, hyperseed, num_experiments, start_experiment, endat, plot_frequency, log_frequency, checkpoint_frequency, distributed)
    pde_system = load_model(model)
    hyperparametersweep = StructGeneratorHyperParameterSweep(hyperseed, num_experiments, hyperparameter_generator)
    hyperparameters = generate_hyperparameters(hyperparametersweep)
    log_options = NeuralPDE.LogOptions(;plot_function=get_plot_function(model), log_dir=log_dir, 
        log_frequency=log_frequency, plot_frequency=plot_frequency, checkpoint_frequency=checkpoint_frequency)
    neuralpde_workers = map(NeuralPDE.NeuralPDEWorker, workers())
    experiment_manager = NeuralPDE.ExperimentManager(pde_system, hyperparameters, get_cb_func(model, log_frequency), log_options, neuralpde_workers)

    #NeuralPDE.remote_run_neuralpde_with_logs(pde_system, hyperparameters[start_experiment], get_cb_func(model, log_frequency), log_options, start_experiment, false)
    NeuralPDE.run_experiment_queue(experiment_manager; remote=false, startat=start_experiment, endat=endat, distributed=distributed)
end

main(model, hyperparameter_generator, hyperseed, num_experiments, start_experiment, endat, plot_frequency, log_frequency, checkpoint_frequency, distributed)
nothing
