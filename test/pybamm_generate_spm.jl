#ENV["PYTHON"] = joinpath(homedir(), "anaconda3", "envs", "pybamm_dev", "bin", "python")
#using Pkg
#Pkg.build("PyCall")
using Distributed
addprocs(20)
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
end
seed = 1
Random.seed!(seed)

@everywhere model = SPMModel()
@everywhere models_dir = abspath(joinpath(@__DIR__, "..", "models"))
@everywhere output_dir = joinpath(models_dir, "$(pybamm_func_str(model))")
regenerate_sim = true
if regenerate_sim
    pybamm = pyimport("pybamm")
    current_input = false
    sim_data, pde_system = generate_sim_model(model; current_input=current_input, output_dir=output_dir)
else 
    include(joinpath(output_dir, "model.jl"))
end

sim_data = read_sim_data(output_dir)
@everywhere sim_data = $sim_data
@everywhere pde_system = $pde_system

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
        StructGenerator(:ADAMOptimiser, 100_000, 3e-4),
        StructGenerator(:ADAMOptimiser, 100_000, 3e-3),
        StructGenerator(:ADAMOptimiser, 100_000, 1e-3)
    )
)


hyperparametersweep = StructGeneratorHyperParameterSweep(1, 120, sg)
hyperparameters = generate_hyperparameters(hyperparametersweep)

#iteration = [0]
@everywhere function get_cb(logger, iteration, chains, param_indices)
    eval_network_at_sim_data = eval_network_at_sim_data_func()
    cb = function (p,l)
        if iteration[1] % 50 == 0
            res = eval_network_at_sim_data(p, chains, param_indices)
            #@info "Iteration: $(iteration[1]), loss: $(l)" 
            #@info "q_error_norm: $(@sprintf("%.4e", res[:error_norms][:q])), csn_error_norm: $(@sprintf("%.4e", res[:error_norms][:csn])), csp_error_norm: $(@sprintf("%.4e", res[:error_norms][:csp])), q_error_relative: $(@sprintf("%.4e", res[:error_rels][:q])), csn_error_relative: $(@sprintf("%.4e", res[:error_rels][:csn])), csp_error_relative: $(@sprintf("%.4e", res[:error_rels][:csp]))"
            log_value(logger, "errors/q_error_norm", res[:error_norms][:q]; step=iteration[1])
            log_value(logger, "errors/csn_error_norm", res[:error_norms][:csn]; step=iteration[1])
            log_value(logger, "errors/csp_error_norm", res[:error_norms][:csp]; step=iteration[1])
            log_value(logger, "errors/q_error_relative", res[:error_rels][:q]; step=iteration[1])
            log_value(logger, "errors/csn_error_relative", res[:error_rels][:csn]; step=iteration[1])
            log_value(logger, "errors/csp_error_relative", res[:error_rels][:csp]; step=iteration[1])
        end
        #iteration[1] += 1

        return false
    end
    return cb
end

@everywhere function get_plot_function()
    sim_data = read_sim_data(output_dir)
    ts = sim_data[:ivs][:t]
    r_ns = sim_data[:ivs][:r_n]
    r_ps = sim_data[:ivs][:r_p]
    qs = sim_data[:dvs][:Q_Ah]
    csns = sim_data[:dvs][:c_s_n_xav]
    csps = sim_data[:dvs][:c_s_p_xav]


    q_pybamm = plot(ts, qs, title = "q_pybamm");
    csn_pybamm = plot(ts, r_ns, csns, linetype=:contourf,title = "csn_pybamm")
    csp_pybamm = plot(ts, r_ps, csps, linetype=:contourf,title = "csp_pybamm")

    pybamm_plots = [(name="q_pybamm", image=q_pybamm), (name="csn_pybamm", image=csn_pybamm), (name="csp_pybamm", image=csp_pybamm)]

    has_given_pybamm_plots = [false]

    eval_network_at_sim_data = eval_network_at_sim_data_func()

    function plot_function(logger, step, phi, θ, indices_in_params, adaloss)
        res = eval_network_at_sim_data(θ, phi, indices_in_params)

        q_nn = plot(ts, res[:evals][:q], title = "q_nn");
        q_error = plot(ts, res[:errors][:q], title = "q_error");

        csn_nn = plot(ts, r_ns, res[:evals][:csn], linetype=:contourf,title = "csn_nn")
        csn_error = plot(ts, r_ns, res[:errors][:csn], linetype=:contourf,title = "csn_error")

        csp_nn = plot(ts, r_ps, res[:evals][:csp], linetype=:contourf,title = "csp_nn")
        csp_error = plot(ts, r_ps, res[:errors][:csp], linetype=:contourf,title = "csp_error")



        eval_compare_plots = [(name="q_nn", image=q_nn), (name="q_error", image=q_error),
                              (name="csn_nn", image=csn_nn), (name="csn_error", image=csn_error),
                              (name="csp_nn", image=csp_nn), (name="csp_error", image=csp_error)]
        if !(has_given_pybamm_plots[1])
            has_given_pybamm_plots[1] = true
            vcat(pybamm_plots, eval_compare_plots)
        else
            eval_compare_plots
        end
    end
    return plot_function
end


#run = 1
#loggerloc = abspath(joinpath(@__DIR__, "..", "logs", "$run"))
#if isdir(loggerloc)
    #rm(loggerloc, recursive=true)
#end
#logger = TBLogger(loggerloc, tb_append) #create tensorboard logger
#hid = 40

#network = SimpleFeedForwardNetwork(3, 64, GELUNonLin(), GlorotUniformParams())
#Random.seed!(seed)
#chains, initialparams = NeuralPDE.getfunction(network, [1, 2, 2])
#flat_initialparams = vcat(initialparams...)
#param_lengths = length.(initialparams)
#param_cumsums = vcat([0], cumsum(param_lengths))
#param_indices = [param_cumsums[i] + 1 : param_cumsums[i+1] for i in 1:length(param_lengths)]


@everywhere function eval_network_at_sim_data_func()
    sim_data = read_sim_data(output_dir)
    ts = sim_data[:ivs][:t]
    r_ns = sim_data[:ivs][:r_n]
    r_ps = sim_data[:ivs][:r_p]
    qs = sim_data[:dvs][:Q_Ah]
    csns = sim_data[:dvs][:c_s_n_xav]
    csps = sim_data[:dvs][:c_s_p_xav]
    qs_norm = LinearAlgebra.norm(qs)
    csns_norm = LinearAlgebra.norm(csns)
    csps_norm = LinearAlgebra.norm(csps)

    function eval_network_at_sim_data(p, chains, param_indices)
        qs_eval = Vector{Float64}(undef, length(ts))
        csns_eval = Matrix{Float64}(undef, length(ts), length(r_ns))
        csps_eval = Matrix{Float64}(undef, length(ts), length(r_ps))

        for (i, ti) in enumerate(ts)
            qs_eval[i] = first(chains[1]([ti], @view p[param_indices[1]]))
        end
        for (i, ti) in enumerate(ts), (j, r_nj) in enumerate(r_ns)
            csns_eval[i, j] = first(chains[2]([ti, r_nj], @view p[param_indices[2]]))
        end
        for (i, ti) in enumerate(ts), (j, r_pj) in enumerate(r_ps)
            csps_eval[i, j] = first(chains[3]([ti, r_pj], @view p[param_indices[3]]))
        end

        qs_error = qs_eval .- qs
        csns_error = csns_eval .- csns
        csps_error = csps_eval .- csps

        qs_error_norm = LinearAlgebra.norm(qs_error)
        csns_error_norm = LinearAlgebra.norm(csns_error)
        csps_error_norm = LinearAlgebra.norm(csps_error)

        qs_error_rel = qs_error_norm / qs_norm
        csns_error_rel = csns_error_norm / csns_norm
        csps_error_rel = csps_error_norm / csps_norm

        return (evals=(q=qs_eval, csn=csns_eval, csp=csps_eval),
                errors=(q=qs_error, csn=csns_error, csp=csps_error),
                error_norms=(q=qs_error_norm, csn=csns_error_norm, csp=csps_error_norm),
                error_rels=(q=qs_error_rel, csn=csns_error_rel, csp=csps_error_rel))
    end

end


#eval_network_at_sim_data = eval_network_at_sim_data_func()
#strategy_ =  NeuralPDE.StochasticTraining(256)
#adaptive_loss = MiniMaxAdaptiveLoss(;reweight_every=50, pde_max_optimiser=Flux.ADAM(3e-3))
log_options = NeuralPDE.LogOptions(;plot_function=get_plot_function(), log_dir=abspath(joinpath(@__DIR__, "..", "logs")))
#discretization = NeuralPDE.PhysicsInformedNN(chains,
                                                #strategy_;
                                                #init_params = initialparams,
                                                #adaptive_loss = adaptive_loss,
                                                #logger = logger,
                                                #log_options = log_options,
                                                #iteration = iteration)



#prob = NeuralPDE.discretize(pde_system,discretization)
#phi = discretization.phi

#maxiters = 100_000
#nn_res = GalacticOptim.solve(prob, ADAM(3e-3); maxiters=maxiters, cb=cb_func)

neuralpde_workers = map(NeuralPDE.NeuralPDEWorker, workers())
experiment_manager = NeuralPDE.ExperimentManager(pde_system, hyperparameters, get_cb, log_options, neuralpde_workers)


NeuralPDE.run_experiment_queue(experiment_manager; remote=false)

nothing
