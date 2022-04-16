rebuild = false
if rebuild
    ENV["PYTHON"] = joinpath(ENV["ANACONDA_LOCATION"], "envs", "pybamm_dev", "bin", "python")
    using Pkg
    Pkg.build("PyCall")
end
begin 
    using DFNExperiments
    using ModelingToolkit, Symbolics, DomainSets
    import ModelingToolkit: Interval, infimum, supremum
    using Infiltrator
end

begin 
    #models = [SPMModel(), SPMeModel()]
    #all_models = [SPMnoRModel(), SPMModel(), ReducedCModel(), SPMeModel(), ReducedCPhiModel(), ReducedCPhiJModel(), DFNnoRModel(), DFNModel()]
    all_models = [ReducedCModel(), SPMeModel(), ReducedCPhiModel(), ReducedCPhiJModel(), DFNnoRModel(), DFNModel()]
    #model = models[2]
    for model in all_models
    #begin
        models_dir = abspath(joinpath(@__DIR__, "..", "models"))
        model_str = pybamm_func_str(model)
        output_dir = joinpath(models_dir, "$(model_str)")
        model_file = joinpath(output_dir, "model.jl")
        current_input = false
        @show model
        sim_data, pde_system = generate_sim_model(model; current_input=current_input, output_dir=output_dir)
        #include(model_file)
    end
end
