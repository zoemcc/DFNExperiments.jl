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
    using Term
    #using Infiltrator
end

begin 
    ev = Float64[]
    #models = [SPMModel(), SPMeModel()]
    #SPMnoRModel(), 
    all_models = [SPMModel(), ReducedCModel(), SPMeModel(), ReducedCPhiModel(), ReducedCPhiJModel(), DFNnoRModel(), DFNModel()]
    num_pts = 4000
    large_interp_grid_lengths = 1000
    small_interp_grid_length = 100
    num_stochastic_samples_from_loss = 1024
    current_input = false
    #model = all_models[1]
    #for i in 1:5
    j = 5
    i = j
    #for i in j:j
    begin
        model = all_models[i]
        include_q = DFNExperiments.include_q_model(model)
        models_dir = abspath(joinpath(@__DIR__, "..", "models"))
        model_str = pybamm_func_str(model)
        output_dir = joinpath(models_dir, "$(model_str)")
        model_file = joinpath(output_dir, "model.jl")
        loss_file = joinpath(output_dir, "loss_certificate.txt")
        writemodel=false
        writesimdata=true
        writeloss=true

        @show model
        #sim_data, pde_system, sim, variables = generate_sim_model(model; current_input=current_input, output_dir=output_dir, num_pts=num_pts)
        sim_data, pde_system, sim, variables, independent_variables_to_pybamm_names, 
            dependent_variables_to_pybamm_names, dependent_variables_to_dependencies, dvs_interpolation,
            dvs_fastchain, prob, total_loss, symb_modded_pde_system, discretization = 
                generate_sim_model_and_test(model; current_input=current_input, include_q=include_q, output_dir=output_dir, num_pts=num_pts,
                    large_interp_grid_length=large_interp_grid_lengths, small_interp_grid_length=small_interp_grid_length,
                    num_stochastic_samples_from_loss=num_stochastic_samples_from_loss,
                    writemodel=writemodel, writesimdata=writesimdata
                    )
        """
        results_nt = generate_sim_model_and_test(model; current_input=current_input, include_q=include_q, output_dir=output_dir, num_pts=num_pts,
                    large_interp_grid_length=large_interp_grid_lengths, small_interp_grid_length=small_interp_grid_length,
                    num_stochastic_samples_from_loss=num_stochastic_samples_from_loss,
                    writemodel=writemodel,
                    )
        """
        nothing
        if !isnothing(prob)
            example_loss = prob.f(ev,ev)
            @show example_loss
            if writeloss
                open(loss_file, "w") do f
                    print(f, string(example_loss) * "\n")
                end
            end
        end

        #include(model_file)
        #solvars = [sim.solution.__getitem__(var) for var in variables]
        # solcsn = solvars[end]
        # solcsn(t=0.15930183773127454 * solcsn.timescale, r=1.0*solcsn.length_scales["negative particle size"], x=-100.0)
    end
end
#loss_file = joinpath(output_dir, "loss_certificate.txt")
        #solvars = [sim.solution.__getitem__(var) for var in variables]

nothing
"""
#@named modded_pde_system = PDESystem(pde_system.eqs, pde_system.bcs, pde_system.domain[[1,3,2]], pde_system.ivs[[1,3,2]], pde_system.dvs)
#@named modded_pde_system = PDESystem([pde_system.dvs[1] ~ pde_system.dvs[1]], pde_system.bcs, pde_system.domain, pde_system.ivs, pde_system.dvs)
#@named modded_pde_system = PDESystem(pde_system.eqs, [pde_system.dvs[1] ~ pde_system.dvs[1]], pde_system.domain, pde_system.ivs, pde_system.dvs)
#@named modded_pde_system = PDESystem(pde_system.eqs[[2]], [pde_system.dvs[1] ~ pde_system.dvs[1]], pde_system.domain, pde_system.ivs, pde_system.dvs)
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
rgridp = interpcsp_py.y
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


function matrixapply_q(v, ??)
    #println("inside q")
    #@show v
    #@show size(v)
    retval = q_interp.(v[[1], :])
    #@show retval
    #@show size(retval)
    retval
end
function matrixapply_csn(v, ??)
    #println("inside csn")
    #@show v
    #@show size(v)
    retval = reshape(csn_interp.(v[1, :], v[2, :]), 1, size(v, 2))
    #@show retval
    #@show size(retval)
    retval
end
function matrixapply_csp(v, ??)
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

"""

"""
unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))

j = 4
dv = pde_system.dvs[j]
dv_pybamm_name = dependent_variables_to_pybamm_names[nameof(dv)]
dv_processed = sim.solution.__getitem__(dv_pybamm_name)
dv_pybamm_interpolation_function = dv_processed._interpolation_function
deps = dependent_variables_to_dependencies[nameof(dv)]
num_deps = length(deps)
py_axis_name = (:x, :y, :z)
new_iv_grid_length = 100
iv_ranges = map(1:num_deps) do i
    iv = deps[i]
    iv_pybamm_name = independent_variables_to_pybamm_names[iv]
    iv_axis_name = py_axis_name[i]
    iv_grid = dv_pybamm_interpolation_function[iv_axis_name]
    # assume time is first
    iv_scale = i == 1 ? dv_processed.timescale : dv_processed.length_scales[iv_pybamm_name]
    unscaled_iv_range = range(Interval(iv_grid[1], iv_grid[end]), new_iv_grid_length)
    scaled_iv_range = range(Interval(iv_grid[1] / iv_scale, iv_grid[end] / iv_scale), new_iv_grid_length)
    (unscaled_iv_range, scaled_iv_range)
end
unscaled_iv_ranges, scaled_iv_ranges = unzip(iv_ranges)
unscaled_ivs_mat = reduce(hcat, map(collect, vec(collect(product(unscaled_iv_ranges...)))))
dv_pybamm_interpolator = DFNExperiments.FastChainInterpolator(dv_pybamm_interpolation_function)
dv_mat = first.(dv_pybamm_interpolator(unscaled_ivs_mat, Float64[]))
dv_grid = reshape(dv_mat, fill(new_iv_grid_length, num_deps)...)
dv_interpolation = DFNExperiments.FastChainInterpolator(scale(interpolate(dv_grid, BSpline(Cubic(Line(OnGrid())))), scaled_iv_ranges...))

using Plots
t = 0.0
ts = range(pde_system.domain[1].domain, length=100)
xs = range(pde_system.domain[4].domain, length=100)
xs = range(start=0.0,stop=1., length=100)
t = ts[1]
txs = vcat(fill(t, 100)', xs')
c_es = vec(dvs_interpolation[4](txs, Float64[]))
plot(xs, c_es)
anim = @animate for i in 1:length(ts)
    t = ts[i]
    txs = vcat(fill(t, 100)', xs')
    eps_c_es = vec(dvs_interpolation[4](txs, Float64[]))
    plot(xs, eps_c_es)
end
gif(anim, "gifs/spme_interp.gif",fps=30)

anim = @animate for i in 1:10
    t = ts[i]
    txs = vcat(fill(t, 100)', xs')
    eps_c_es = vec(dvs_interpolation[4](txs, Float64[]))
    plot(xs, eps_c_es)
end
gif(anim, "gifs/spme_interp_slow.gif",fps=1)
"""