rebuild = false
if rebuild
    ENV["PYTHON"] = joinpath(ENV["ANACONDA_LOCATION"], "envs", "pybamm_dev", "bin", "python")
    using Pkg
    Pkg.build("PyCall")
end
begin
    using DFNExperiments
    using DFNExperimentsPyBaMM
    using FiniteDifferences
    using DiffEqFlux, NeuralPDE
    using Interpolations
    using ModelingToolkit, Symbolics, DomainSets
    import ModelingToolkit: Interval, infimum, supremum
    using IterTools
    using Term
    using Plots
    using Infiltrator
end

function main()
    ev = Float64[]
    #models = [SPMModel(), SPMeModel()]
    #SPMnoRModel(), 
    all_models = [SPMModel(), ReducedCModel(), SPMeModel(), ReducedCPhiModel(), ReducedCPhiJModel(), DFNnoRModel(), DFNModel()]
    #num_pts = 4000
    #num_pts = 1000
    #num_pts = 100
    #large_interp_grid_length = 1000
    #large_interp_grid_length = 400
    #large_interp_grid_length = 40
    #small_interp_grid_length = 400
    #small_interp_grid_length = 100
    #small_interp_grid_length = 30
    #num_stochastic_samples_from_loss = 1024
    #num_stochastic_samples_from_loss = 4096
    current_input = false
    include_q = false
    #model = all_models[1]
    j = 1
    i = j
    # 5
    bad_models = [3]
    #for i in bad_models
    #for i in 7:7
    begin
        model = all_models[i]
        #include_q = DFNExperiments.include_q_model(model)
        models_dir = abspath(joinpath(@__DIR__, "..", "..", "..", "models"))
        model_str = pybamm_func_str(model)
        num_pts = num_pts_validation(model)
        num_tsteps = num_tsteps_validation(model)
        large_interp_grid_length = large_interp_grid_length_validation(model)
        large_grid_tsteps_length = large_grid_tsteps_validation(model)
        small_interp_grid_length = small_interp_grid_length_validation(model)
        small_grid_tsteps_length = small_grid_tsteps_validation(model)
        num_stochastic_samples_from_loss = num_stochastic_samples_from_loss_validation(model)
        @show num_pts
        @show large_interp_grid_length
        @show large_grid_tsteps_length
        @show small_interp_grid_length
        @show small_grid_tsteps_length
        @show num_stochastic_samples_from_loss
        output_dir = joinpath(models_dir, "$(model_str)")
        model_file = joinpath(output_dir, "model.jl")
        loss_file = joinpath(output_dir, "loss_certificate.txt")
        writemodel = false
        writesimdata = true
        writeloss = true

        @show model
        #sim_data, pde_system, sim, variables = generate_sim_model(model; current_input=current_input, output_dir=output_dir, num_pts=num_pts)
        full_results =
            generate_sim_model_and_test(model; current_input=current_input, include_q=include_q, output_dir=output_dir, num_pts=num_pts,
                num_tsteps=num_tsteps, large_interp_grid_length=large_interp_grid_length, small_interp_grid_length=small_interp_grid_length,
                large_grid_tsteps_length=large_grid_tsteps_length, small_grid_tsteps_length=small_grid_tsteps_length,
                num_stochastic_samples_from_loss=num_stochastic_samples_from_loss,
                writemodel=writemodel, writesimdata=writesimdata
            )
        @unpack sim_data_nt, pde_system, sim, variables, independent_variables_to_pybamm_names,
        dependent_variables_to_pybamm_names, dependent_variables_to_dependencies, dvs_interpolation,
        dvs_luxchain, prob, total_loss, pinnrep, discretization = full_results
        """
        results_nt = generate_sim_model_and_test(model; current_input=current_input, include_q=include_q, output_dir=output_dir, num_pts=num_pts,
                    large_interp_grid_length=large_interp_grid_lengths, small_interp_grid_length=small_interp_grid_length,
                    num_stochastic_samples_from_loss=num_stochastic_samples_from_loss,
                    writemodel=writemodel,
                    )
        """
        #nothing
        if writeloss && !isnothing(total_loss)
            @show total_loss
            if writeloss
                open(loss_file, "w") do f
                    print(f, string(total_loss) * "\n")
                end
            end
        end

        #include(model_file)
        #solvars = [sim.solution.__getitem__(var) for var in variables]
        # solcsn = solvars[end]
        # solcsn(t=0.15930183773127454 * solcsn.timescale, r=1.0*solcsn.length_scales["negative particle size"], x=-100.0)
    end
    return full_results
end
full_results = main()
#loss_file = joinpath(output_dir, "loss_certificate.txt")
#solvars = [sim.solution.__getitem__(var) for var in variables]

c_e_n_lf = pinnrep.loss_functions.datafree_pde_loss_functions[3]
c_e_p_lf = pinnrep.loss_functions.datafree_pde_loss_functions[5]
flat_init_params = pinnrep.flat_init_params
Nt = 400
Nx = 200
ts = range(0.0, 0.159, Nt)
x_ns = range(0.0, 0.4444444444444445, Nx)
x_ns = range(0.0, 0.4444444444444445, Nx)
ts_x_ns = cartesian_product(ts, x_ns; flat=Val{true})
c_e_n_s_raw = reshape(c_e_n_lf(ts_x_ns, flat_init_params), (Nt, Nx))
c_e_n_s_abs = abs.(c_e_n_s_raw)
c_e_n_s_log_abs = log.(c_e_n_s_abs)
p1 = heatmap(ts, x_ns, c_e_n_s_abs', title="abs loss for the c_e_n differential equation")
p2 = heatmap(ts, x_ns, c_e_n_s_log_abs', title="log abs loss for the c_e_n differential equation")

x_ps = range(0.5555555555555556, 1.0, Nx)
ts_x_ps = cartesian_product(ts, x_ps; flat=Val{true})
c_e_p_s_raw = reshape(c_e_p_lf(ts_x_ps, flat_init_params), (Nt, Nx))
c_e_p_s_abs = abs.(c_e_p_s_raw)
c_e_p_s_log_abs = log.(c_e_p_s_abs)
heatmap(ts, x_ps, c_e_p_s_abs', title="abs loss for the c_e_p differential equation")
heatmap(ts, x_ps, c_e_p_s_log_abs', title="log abs loss for the c_e_p differential equation")

@unpack domains, eqs, bcs, dict_indvars, dict_depvars, flat_init_params = pinnrep
strategy = discretization.strategy
eltypeθ = eltype(pinnrep.flat_init_params)
bounds = get_bounds(domains, eqs, bcs, eltypeθ, dict_indvars, dict_depvars, strategy)

dist = NeuralPDE.distance_to_boundary_function(bounds[1][1], 1 / 10)

dist(ts_x_ns)
ts_x = collect(hcat(ts, fill(0.1, 400))')
dists = dist(ts_x)[1, :]
plot(ts, dists, title="distance to boundary")
plot(ts, dists .^ 2, title="distance squared to boundary")
nothing



"""
c_s_n_lf = pinnrep.loss_functions.datafree_pde_loss_functions[1]
c_n_p_lf = pinnrep.loss_functions.datafree_pde_loss_functions[2]
c_s_n_ic_lf = pinnrep.loss_functions.datafree_bc_loss_functions[1]
flat_init_params = pinnrep.flat_init_params
@unpack domains, eqs, bcs, dict_indvars, dict_depvars, flat_init_params = pinnrep
strategy = discretization.strategy
eltypeθ = eltype(pinnrep.flat_init_params)
bounds = get_bounds(domains, eqs, bcs, eltypeθ, dict_indvars, dict_depvars,
                    strategy)
function distance_to_boundary_function(bound, ratio_to_full_weight)
    this_distance = function distance_to_boundary(x)
        distances = map(enumerate(bound)) do (i, variable) 
            if variable isa Array
                bound_length = variable[2] - variable[1]
                scale = 1 / (ratio_to_full_weight * bound_length)
                x_i = @view x[[i], :]
                min.(max.(min.(x_i .- variable[1], variable[2] .- x_i), 0) .* scale, 1)
            else
                1
            end
        end
        min.(distances...)
    end
    this_distance
end

dfs = distance_to_boundary_function(bounds[1][1], 1/10)
variable = bounds[1][1][1]
ratio_to_full_weight = 1/10
bound_length = variable[2] - variable[1]
scale = 1 / (bound_length * ratio_to_full_weight)
x = ts
d1 = x .- variable[1]
d2 = variable[2] .- x
m1 = min.(d1, d2)
m2 = max.(m1, 0) .* scale
m3 = min.(m2, 1)
distances = distance_to_boundary(ts_r_ns, bounds[1][1], 1/10)
#df(x) = max.(min.(x .- variable[1], variable[2] .- x), 0.0)
ts = @view ts_r_ns[[1], :]
dfs[1](ts)
distance_to_boundary_function(bounds[1][1])
Nt = 400
Nr = 200
ts = range(0., 0.159, Nt)
r_ns = range(0.0, 1.0, Nr)
ts_r_ns = cartesian_product(ts, r_ns; flat=Val{true})
c_s_n_s_raw = reshape(c_s_n_lf(ts_r_ns, flat_init_params), (Nt, Nr))
c_s_n_s_abs = abs.(c_s_n_s_raw)
c_s_n_s_log_abs = log.(c_s_n_s_abs)
heatmap(ts, r_ns, c_s_n_s_abs', title="abs loss for the c_s_n differential equation")
heatmap(ts, r_ns, c_s_n_s_log_abs', title="log abs loss for the c_s_n differential equation")

x_ps = range(0.5555555555555556, 1.0, Nx)
ts_x_ps = cartesian_product(ts, x_ps; flat=Val{true})
c_e_p_s_raw = reshape(c_e_p_lf(ts_x_ps, flat_init_params), (Nt, Nx))
c_e_p_s_abs = abs.(c_e_p_s_raw)
c_e_p_s_log_abs = log.(c_e_p_s_abs)
heatmap(ts, x_ps, c_e_p_s_abs', title="abs loss for the c_e_p differential equation")
heatmap(ts, x_ps, c_e_p_s_log_abs', title="log abs loss for the c_e_p differential equation")
#plot(x_ns, ts, c_e_n_s_abs, linetype=:contourf, title="abs loss for the c_e_n differential equation")
"""
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


function matrixapply_q(v, θ)
    #println("inside q")
    #@show v
    #@show size(v)
    retval = q_interp.(v[[1], :])
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