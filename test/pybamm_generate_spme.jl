begin
    using Distributed
    NUM_WORKERS = 13
    if NUM_WORKERS > 0
        if nworkers() != NUM_WORKERS
            addprocs(NUM_WORKERS - nworkers())
        end
    end
end
@everywhere begin 
    ENV["GKSwstype"] = "100"
    using Pkg
    Pkg.activate(abspath(joinpath(@__DIR__, "..")))
end
@everywhere begin 
    using DFNExperiments, NeuralPDE
    using ModelingToolkit, Symbolics, DomainSets
    import ModelingToolkit: Interval, infimum, supremum
end

begin 
    model = SPMeModel()
    hyperseed = 1
    num_experiments = 120
    start_experiment = 27
    plot_frequency = 5_000
    log_frequency = 1_000
    checkpoint_frequency = 5_000
    iterations = 400_000

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
            RandomChoice(128, 256)
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

function main(model, hyperparameter_generator, hyperseed, num_experiments, start_experiment, plot_frequency, log_frequency, checkpoint_frequency)
    model_str = pybamm_func_str(model)
    pde_system = load_model(model)
    hyperparametersweep = StructGeneratorHyperParameterSweep(hyperseed, num_experiments, hyperparameter_generator)
    hyperparameters = generate_hyperparameters(hyperparametersweep)
    log_options = NeuralPDE.LogOptions(;plot_function=get_plot_function(model), log_dir=abspath(joinpath(@__DIR__, "..", "logs", model_str)), 
        log_frequency=log_frequency, plot_frequency=plot_frequency, checkpoint_frequency=checkpoint_frequency)
    neuralpde_workers = map(NeuralPDE.NeuralPDEWorker, workers())
    experiment_manager = NeuralPDE.ExperimentManager(pde_system, hyperparameters, get_cb_func(model, log_frequency), log_options, neuralpde_workers)

    #NeuralPDE.remote_run_neuralpde_with_logs(pde_system, hyperparameters[start_experiment], get_cb_func(model, log_frequency), log_options, start_experiment, false)
    NeuralPDE.run_experiment_queue(experiment_manager; remote=false, startat=start_experiment)
end

main(model, hyperparameter_generator, hyperseed, num_experiments, start_experiment, plot_frequency, log_frequency, checkpoint_frequency)
nothing
