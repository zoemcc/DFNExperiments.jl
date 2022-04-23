rebuild = false
if rebuild
    ENV["PYTHON"] = joinpath(ENV["ANACONDA_LOCATION"], "envs", "pybamm_dev", "bin", "python")
    using Pkg
    Pkg.build("PyCall")
end
begin 
    using DFNExperiments
    using FiniteDifferences
    using ModelingToolkit, Symbolics, DomainSets
    import ModelingToolkit: Interval, infimum, supremum
    using IterTools
    using Infiltrator
end

begin 
    #models = [SPMModel(), SPMeModel()]
    all_models = [SPMnoRModel(), SPMModel(), ReducedCModel(), SPMeModel(), ReducedCPhiModel(), ReducedCPhiJModel(), DFNnoRModel(), DFNModel()]
    num_pts = 100
    #all_models = [ReducedCModel(), SPMeModel(), ReducedCPhiModel(), ReducedCPhiJModel(), DFNnoRModel(), DFNModel()]
    model = all_models[2]
    #for model in all_models
    begin
        models_dir = abspath(joinpath(@__DIR__, "..", "models"))
        model_str = pybamm_func_str(model)
        output_dir = joinpath(models_dir, "$(model_str)")
        model_file = joinpath(output_dir, "model.jl")
        current_input = false
        @show model
        #sim_data, pde_system, sim, variables = generate_sim_model(model; current_input=current_input, output_dir=output_dir, num_pts=num_pts)
        sim_data, pde_system, sim, variables = generate_sim_model(model; current_input=current_input, num_pts=num_pts)
        #include(model_file)
        # solvars = [sim.solution.__getitem__(var) for var in variables]
        # solcsn = solvars[end]
        # solcsn(t=0.15930183773127454 * solcsn.timescale, r=1.0*solcsn.length_scales["negative particle size"], x=-100.0)
    end
end
solvars = Dict([(var=>sim.solution.__getitem__(var)) for var in variables]...)
solcsp = solvars["X-averaged positive particle concentration"]
solcsn = solvars["X-averaged negative particle concentration"]
csn(t, r) = solcsn(t=t*0.15930183773127454 * solcsn.timescale, r=r*solcsn.length_scales["negative particle size"], x=0.0)[1]
#csn(t, r) = solcsn(t=t * solcsn.timescale, r=r*solcsn.length_scales["negative particle size"], x=0.0)[1]
csntf(t) = begin csnr(r) = csn(t, r) end
csp(t, r) = solcsp(t=t*0.15930183773127454 * solcsp.timescale, r=r*solcsp.length_scales["positive particle size"], x=1.0)[1]
csptf(t) = begin cspr(r) = csp(t, r) end
@show abs(central_fdm(5, 1, max_range=1e-6)(csntf(1), 1) - (-0.14182855923368468))
@show abs(central_fdm(5, 1, max_range=1e-6)(csntf(1), 0) - 0)
@show abs(forward_fdm(5, 1)(csntf(1), 0) - 0)
@show abs(backward_fdm(5, 1)(csntf(1), 1) - (-0.14182855923368468))
@show abs(central_fdm(5, 1, max_range=1e-6)(csptf(1), 1) - (0.03237700710041634))
@show abs(central_fdm(5, 1, max_range=1e-6)(csptf(1), 0) - 0)
@show abs(forward_fdm(5, 1)(csptf(1), 0) - 0)
@show abs(backward_fdm(5, 1)(csptf(1), 1) - (0.03237700710041634))
@show central_fdm(5, 1, max_range=1e-6)(csptf(1), 1)
ts = range(Interval(0,1), 100)
csn_accurate = map(tr->csn(tr[1], tr[2]), collect(product(ts, ts)))
csp_accurate = map(tr->csp(tr[1], tr[2]), collect(product(ts, ts)))
@show (csn_accurate[end, end] - csn_accurate[end, end-1]) / (ts[end] - ts[end-1])
@show (csp_accurate[end, end] - csp_accurate[end, end-1]) / (ts[end] - ts[end-1])
