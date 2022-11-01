using DFNExperiments, NeuralPDE, Lux, Interpolations, Random

model = SPMModel()
interpolators = DFNExperiments.get_interpolator_from_model(model)
exp_title = "SPMModel_minimax"
num_cores = 1
start_experiment, end_experiment = NeuralPDE.slurm_find_indices_to_execute(num_cores)

sg = NeuralPDE.StructGenerator(
    :CompositeHyperParameter,
    NeuralPDE.RandomChoice(1:2^10),
    NeuralPDE.
    StructGenerator(:SimpleFeedForwardNetwork,
        NeuralPDE.RandomChoice(4:7),
        NeuralPDE.RandomChoice(10, 20, 30),
        NeuralPDE.RandomChoice(:SigmoidNonLin, :GELUNonLin)),
    NeuralPDE.StructGenerator(:StochasticTraining,
        NeuralPDE.RandomChoice(64, 128, 256),
        NeuralPDE.RandomChoice(32, 64, 128)),
    NeuralPDE.StructGenerator(:MiniMaxAdaptiveLoss, 100),
    NeuralPDE.RandomChoice(
        NeuralPDE.StructGenerator(:ADAMOptimiser, 1000000, 0.001),
        NeuralPDE.StructGenerator(:ADAMOptimiser, 1000000, 0.0003)))
hyperparametersweep = NeuralPDE.StructGeneratorHyperParameterSweep(3, num_cores, sg)
hyperparameters = NeuralPDE.generate_hyperparameters(hyperparametersweep)
hyperparam = hyperparameters[1]

pde_system = load_model(model)
log_dir = joinpath(ENV["HOME"], "logs", "full_length_2/SPMModel", "SPMModel_minimax")
seed = NeuralPDE.getseed(hyperparam)
rng = Random.MersenneTwister(seed)


# Neural network
num_ivs_for_dvs = map(pde_system.dvs) do dv
    # assumes dv is in the form u(t,x) etc 
    return length(dv.val.arguments)
end
#chains, init_params = NeuralPDE.getfunction(rng, hyperparam, num_ivs_for_dvs)
#interpolators
#i = 1
#ci, inti = chains[i], interpolators[i]
#plusi = Lux.Parallel(+, ci, inti)
#thi = Lux.initialstates(rng, plusi)
int_plus_net_spec = DFNExperiments.InterpolationPlusNetwork(model, hyperparam.nn)
int_plus_net_chains, int_plus_net_init_params = NeuralPDE.getfunction(rng, int_plus_net_spec, num_ivs_for_dvs)
i = 1
ipni, ipnith = int_plus_net_chains[i], int_plus_net_init_params[i]
inti = ipni.layers[1](x, NamedTuple(), ps)
ipnith2[:layer_2]
#neti = ipni.layers[2](x, ipnith2[:layer_2], ps)
ipnith2, ps = Lux.setup(rng, ipni)
#ps = NamedTuple()
x0 = [0.0; 0.1]
x = [0.0 0.1 0.15; 0.1 0.2 0.3]
ipni(x, ipnith2, ps)
ipni(x, ipnith, ps)
nothing
#log_frequency, plot_frequency, checkpoint_frequency = 5000, 10000, 10000
#log_options = NeuralPDE.LogOptions(; plot_function=DFNExperiments.get_plot_function(model), log_dir=log_dir,
#log_frequency=log_frequency, plot_frequency=plot_frequency, checkpoint_frequency=checkpoint_frequency)
#neuralpde_workers = map(NeuralPDE.NeuralPDEWorker, workers())
#experiment_manager = NeuralPDE.ExperimentManager(pde_system, hyperparameters, get_cb_func(model, log_frequency), log_options, neuralpde_workers)

#NeuralPDE.run_experiment_queue(experiment_manager; remote=false, startat=start_experiment, endat=end_experiment, distributed=false)
nothing
