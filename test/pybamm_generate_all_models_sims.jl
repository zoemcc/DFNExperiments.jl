rebuild = false
if rebuild
    ENV["PYTHON"] = joinpath(ENV["ANACONDA_LOCATION"], "envs", "pybamm_dev", "bin", "python")
    using Pkg
    Pkg.build("PyCall")
end
begin 
    using DFNExperiments
    using FiniteDifferences
    using DiffEqFlux, NeuralPDE
    using Interpolations
    using ModelingToolkit, Symbolics, DomainSets
    import ModelingToolkit: Interval, infimum, supremum
    using IterTools
    using Infiltrator
end

begin 
    #models = [SPMModel(), SPMeModel()]
    all_models = [SPMnoRModel(), SPMModel(), ReducedCModel(), SPMeModel(), ReducedCPhiModel(), ReducedCPhiJModel(), DFNnoRModel(), DFNModel()]
    num_pts = 4000
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
#@named modded_pde_system = PDESystem(pde_system.eqs, pde_system.bcs, pde_system.domain[[1,3,2]], pde_system.ivs[[1,3,2]], pde_system.dvs)
#@named modded_pde_system = PDESystem([pde_system.dvs[1] ~ pde_system.dvs[1]], pde_system.bcs, pde_system.domain, pde_system.ivs, pde_system.dvs)
#@named modded_pde_system = PDESystem(pde_system.eqs, [pde_system.dvs[1] ~ pde_system.dvs[1]], pde_system.domain, pde_system.ivs, pde_system.dvs)
@named modded_pde_system = PDESystem(pde_system.eqs[[2]], [pde_system.dvs[1] ~ pde_system.dvs[1]], pde_system.domain, pde_system.ivs, pde_system.dvs)
@named modded_pde_system = PDESystem(pde_system.eqs, [pde_system.dvs[1] ~ pde_system.dvs[1]], pde_system.domain, pde_system.ivs, pde_system.dvs)
solvars = Dict([(var=>sim.solution.__getitem__(var)) for var in variables]...)
solq = solvars["Discharge capacity [A.h]"]
solcsp = solvars["X-averaged positive particle concentration"]
solcsn = solvars["X-averaged negative particle concentration"]
interpq_py = solq._interpolation_function
interpcsn_py = solcsn._interpolation_function
interpcsp_py = solcsp._interpolation_function
tgrid = interpq_py.x
trange = range(Interval(tgrid[1] / solq.timescale, tgrid[end] / solq.timescale), length(tgrid))
qgrid = interpq_py.y
tgridn = interpcsn_py.x
rgridn = interpcsn_py.y
rnrange = range(Interval(rgridn[1] / solcsn.length_scales["negative particle size"], rgridn[end] / solcsn.length_scales["negative particle size"]), length(rgridn))
rprange = range(Interval(rgridp[1] / solcsp.length_scales["positive particle size"], rgridp[end] / solcsp.length_scales["positive particle size"]), length(rgridp))
tgridp = interpcsp_py.x
rgridp = interpcsp_py.y
csngrid = collect(reshape(interpcsn_py.z, size(rgridn)..., size(tgridn)...)')
cspgrid = collect(reshape(interpcsp_py.z, size(rgridp)..., size(tgridp)...)')

q_interp = scale(interpolate(qgrid, BSpline(Cubic(Line(OnGrid())))), trange)
csn_interp = scale(interpolate(csngrid, BSpline(Cubic(Line(OnGrid())))), trange, rnrange)
csp_interp = scale(interpolate(cspgrid, BSpline(Cubic(Line(OnGrid())))), trange, rprange)

q(t) = solq(t=t * solcsn.timescale)[1]
csn(t, r) = solcsn(t=t * solcsn.timescale, r=r*solcsn.length_scales["negative particle size"], x=0.0)[1]
#csn(t, r) = solcsn(t=t * solcsn.timescale, r=r*solcsn.length_scales["negative particle size"], x=0.0)[1]
csntf(t) = begin csnr(r) = csn(t, r) end
csp(t, r) = solcsp(t=t * solcsp.timescale, r=r*solcsp.length_scales["positive particle size"], x=1.0)[1]
csptf(t) = begin cspr(r) = csp(t, r) end
#@show abs(central_fdm(5, 1, max_range=1e-6)(csntf(1), 1) - (-0.14182855923368468))
#@show abs(central_fdm(5, 1, max_range=1e-6)(csntf(1), 0) - 0)
#@show abs(forward_fdm(5, 1)(csntf(1), 0) - 0)
#@show abs(backward_fdm(5, 1)(csntf(1), 1) - (-0.14182855923368468))
#@show abs(central_fdm(5, 1, max_range=1e-6)(csptf(1), 1) - (0.03237700710041634))
#@show abs(central_fdm(5, 1, max_range=1e-6)(csptf(1), 0) - 0)
#@show abs(forward_fdm(5, 1)(csptf(1), 0) - 0)
#@show abs(backward_fdm(5, 1)(csptf(1), 1) - (0.03237700710041634))
#@show central_fdm(5, 1, max_range=1e-6)(csptf(1), 1)
ts = range(Interval(0,1), 100)
#csn_accurate = map(tr->csn(tr[1], tr[2]), collect(product(ts, ts)))
#csp_accurate = map(tr->csp(tr[1], tr[2]), collect(product(ts, ts)))
#@show (csn_accurate[end, end] - csn_accurate[end, end-1]) / (ts[end] - ts[end-1])
#@show (csp_accurate[end, end] - csp_accurate[end, end-1]) / (ts[end] - ts[end-1])
#@show abs((csn_accurate[end, end] - csn_accurate[end, end-1]) / (ts[end] - ts[end-1]) - (-0.14182855923368468))
#@show abs((csp_accurate[end, end] - csp_accurate[end, end-1]) / (ts[end] - ts[end-1]) - (0.03237700710041634))
#@show abs((csn_accurate[end, 2] - csn_accurate[end, 1]) / (ts[2] - ts[1]))
#@show abs((csp_accurate[end, 2] - csp_accurate[end, 1]) / (ts[2] - ts[1]))
#disc = sim._disc
#@show disc.y_slices_explicit


function matrixapply_q(v, θ)
    #println("inside q")
    #@show v
    #@show size(v)
    retval = reshape(q_interp.(v), 1, size(v, 2))
    #@show retval
    #@show size(retval)
    retval
end
function matrixapply_csn(v, θ)
    #println("inside csn")
    #@show v
    #@show size(v)
    retval = reshape(csn_interp.(v[1, :], v[2, :]), 1, size(v, 2))
    #@show retval
    #@show size(retval)
    retval
end
function matrixapply_csp(v, θ)
    #println("inside csp")
    #@show v
    #@show size(v)
    retval = reshape(csp_interp.(v[1, :], v[2, :]), 1, size(v, 2))
    #@show retval
    #@show size(retval)
    retval
end
q_fc = FastChain(matrixapply_q)
csn_fc = FastChain(matrixapply_csn)
csp_fc = FastChain(matrixapply_csp)

ts = reduce(hcat, map(collect, vec(collect(range(Interval(0.000,0.159), 100)))))
trs = reduce(hcat, map(collect, vec(collect(product(range(Interval(0.000,0.159), 100), range(Interval(0.0, 1.0), 100))))))
trs_small = trs[:, 1:10]
trs_onedim = trs[1:2, 1]
#matrixapply_csn(trs)
#matrixapply_csp(trs)

spm_fcs = [q_fc, csn_fc, csp_fc]
q_fc(ts, [])
csn_fc(trs, [])

DiffEqFlux.initial_params(::typeof(matrixapply_q)) = Float64[]
DiffEqFlux.initial_params(::typeof(matrixapply_csn)) = Float64[]
DiffEqFlux.initial_params(::typeof(matrixapply_csp)) = Float64[]

strategy =  NeuralPDE.StochasticTraining(1024, 1024)
#derivative =  NeuralPDE.get_numeric_derivative(epsvar=1e-6)
discretization = NeuralPDE.PhysicsInformedNN(spm_fcs, strategy)
prob = NeuralPDE.discretize(pde_system,discretization)
sym_prob = NeuralPDE.symbolic_discretize(modded_pde_system,discretization)
prob.f(Float64[], Float64[])
nothing

