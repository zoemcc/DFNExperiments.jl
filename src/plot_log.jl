norm2(x) = LinearAlgebra.norm(x, 2)
norm1(x) = LinearAlgebra.norm(x, 1)

function get_eval_network_at_sim_data_func(model, chains)
    sim_data = read_sim_data(model)
    ivs = keys(sim_data[:ivs])
    dvs = keys(sim_data[:dvs])
    dv_deps = collect(map(keys, sim_data[:dv_deps]))
    num_dvs = length(dvs)
    ground_truths = map(dv -> sim_data[:dvs][dv], dvs)
    ground_truths_norm = LinearAlgebra.norm.(ground_truths)
    ground_truths_l2 = norm2.(ground_truths)
    ground_truths_l1 = norm1.(ground_truths)
    mdf = MultiDimensionalFunction(chains, VectorOfParameterizedMDFApplyFuncType(), ivs, dvs, dv_deps)

    function eval_network_at_sim_data(p)
        eval_at = values(sim_data[:ivs])
        network_evals = mdf(p, eval_at...; flat=Val{false})
        #println("")
        #@show ivs
        #@show dvs
        #@show dv_deps
        #@show num_dvs
        #println("")
        #@show typeof(eval_at)
        #@show size.(eval_at)
        #@show typeof(network_evals)
        #@show typeof.(network_evals)
        #@show size.(network_evals)
        #@show typeof(ground_truths)
        #@show typeof.(ground_truths)
        #@show size.(ground_truths)

        errors = map(i -> network_evals[i] .- ground_truths[i], 1:num_dvs)
        errors_l2 = norm2.(errors)
        errors_l1 = norm1.(errors)
        errors_rel_l2 = map(i -> errors_l2[i] / ground_truths_l2[i], 1:num_dvs)
        errors_rel_l1 = map(i -> errors_l1[i] / ground_truths_l1[i], 1:num_dvs)
        #@infiltrate

        return (; dvs, network_evals, errors, errors_l2, errors_l1, errors_rel_l2, errors_rel_l1)
    end
end

function get_cb_func(model, log_frequency)
    function get_cb_func_inner(logger, iteration, chains)
        eval_network_at_sim_data = get_eval_network_at_sim_data_func(model, chains)

        cb = function (p, l)
            if iteration[1] % log_frequency == 0
                res = eval_network_at_sim_data(p)
                error_rel_l2_max = max(res[:errors_rel_l2]...)
                error_rel_l1_max = max(res[:errors_rel_l1]...)
                #@show res[:errors_rel]
                #@show typeof(res[:errors_rel])
                #@show error_rel_max
                #@show typeof(error_rel_max)
                NeuralPDE.logscalar(logger, error_rel_l2_max, "errors/max_error_l2_relative", iteration[1])
                NeuralPDE.logscalar(logger, error_rel_l1_max, "errors/max_error_l1_relative", iteration[1])
                for (i, dv) in enumerate(res[:dvs])
                    NeuralPDE.logscalar(logger, res[:errors_l2][i], "errors/$(string(dv))_error_l2", iteration[1])
                    NeuralPDE.logscalar(logger, res[:errors_l1][i], "errors/$(string(dv))_error_l1", iteration[1])
                    NeuralPDE.logscalar(logger, res[:errors_rel_l2][i], "errors/$(string(dv))_error_relative_l2", iteration[1])
                    NeuralPDE.logscalar(logger, res[:errors_rel_l1][i], "errors/$(string(dv))_error_relative_l1", iteration[1])
                end
                println("Iteration $(iteration[1]), loss $(l)")
            end
            return false
        end
        return cb
    end
end

function do_plot(data_array::AbstractArray{R,N}, iv_data, deps, name) where {R,N}
    used_ivs = map(iv -> iv_data[iv], deps)
    if N == 1
        plot(used_ivs..., data_array, title=name)
    elseif N == 2
        plot(reverse(used_ivs)..., data_array, linetype=:contourf, title=name)
    else
        println("no plot defined for arrays with more than 2 axes")
        nothing
    end
end


function get_plot_function(model)
    function get_plot_function_inner(logger, chains)
        eval_network_at_sim_data = get_eval_network_at_sim_data_func(model, chains)
        sim_data = read_sim_data(model)
        dvs = keys(sim_data[:dvs])
        dv_deps = collect(map(keys, sim_data[:dv_deps]))
        num_dvs = length(dvs)
        iv_data = sim_data[:ivs]
        pybamm_plots = map(i -> (name="$(string(dvs[i]))_pybamm", image=do_plot(sim_data[:dvs][i], iv_data, dv_deps[i], "$(string(dvs[i]))_pybamm")), 1:num_dvs)
        have_given_pybamm_plots = [false]

        function plot_function(θ, pinnrep)
            # error plots
            res = eval_network_at_sim_data(θ)
            eval_plots = map(i -> (name="$(string(dvs[i]))_nn", image=do_plot(res[:network_evals][i], iv_data, dv_deps[i], "$(string(dvs[i]))_nn")), 1:num_dvs)
            error_plots = map(i -> (name="$(string(dvs[i]))_error", image=do_plot(res[:errors][i], iv_data, dv_deps[i], "$(string(dvs[i]))_error")), 1:num_dvs)
            abs_error_plots = map(i -> (name="$(string(dvs[i]))_abs_error", image=do_plot(abs.(res[:errors][i]), iv_data, dv_deps[i], "$(string(dvs[i]))_abs_error")), 1:num_dvs)

            # loss function plots
            @unpack domains, eqs, bcs, dict_indvars, dict_depvars, flat_init_params = pinnrep
            eltypeθ = eltype(pinnrep.flat_init_params)
            bounds = get_bounds(domains, eqs, bcs, eltypeθ, dict_indvars, dict_depvars, StochasticTraining(64, 64))
            pde_bounds, bcs_bounds = bounds
            rangetype = typeof(range(eltypeθ(0.0), eltypeθ(1.0), length=100))
            pde_vars = NeuralPDE.get_variables(eqs, dict_indvars, dict_depvars)
            pde_args = NeuralPDE.get_argument(eqs, dict_indvars, dict_depvars)
            bc_vars = NeuralPDE.get_variables(bcs, dict_indvars, dict_depvars)
            bc_args = NeuralPDE.get_argument(bcs, dict_indvars, dict_depvars)

            num_perdim = 101

            # this is the bc_vals array we're looking for. this could be plotted 
            # (using bc_vars as the axes name, the index as the title, 
            # and the bc_vals as the error of the lf to be plotted (should be squared first))
            # should the bc vals also use the interior range? this seems to sample the range at the boundary (and not the interior only)


            pde_loss_plots = map(enumerate(pinnrep.loss_functions.datafree_pde_loss_functions)) do (i, pde_lfi)
                name = "pde_loss_$(string(i))"
                num_dimsi = length(pde_vars[i])
                rangeboundi = map(pde_bounds[i]) do boundj
                    if boundj isa Vector
                        bound_eltype = eltypeθ.(boundj)
                        range(start=bound_eltype[1], stop=bound_eltype[2], length=num_perdim)
                    else
                        eltypeθ(boundj)
                    end
                end
                only_rangesi = filter(x -> x isa AbstractRange, rangeboundi)
                Ni = length(only_rangesi)
                input_cartprodi = cartesian_product(rangeboundi...; flat=Val{true})
                new_sizei = tuple(fill(num_perdim, num_dimsi)...)
                pde_valsi = reshape(pde_lfi(input_cartprodi, θ), new_sizei) .^ 2
                image = if Ni == 1
                    plot(only_rangesi..., pde_valsi, title=name)
                elseif Ni == 2
                    plot(reverse(only_rangesi)..., pde_valsi, linetype=:contourf, title=name)
                else
                    println("no plot defined for arrays with more than 2 axes")
                    nothing
                end

                # this is the bc_vals array we're looking for. this could be plotted 
                # (using bc_vars as the axes name, the index as the title, 
                #image = do_plot(loss_function, title=name)
                (; name, image)
            end

            bc_loss_plots = map(enumerate(pinnrep.loss_functions.datafree_bc_loss_functions)) do (i, bc_lfi)
                name = "bc_loss_$(string(i))"
                num_dimsi = length(bc_vars[i])
                rangeboundi = map(bcs_bounds[i]) do boundj
                    if boundj isa Vector
                        bound_eltype = eltypeθ.(boundj)
                        range(start=bound_eltype[1], stop=bound_eltype[2], length=num_perdim)
                    else
                        eltypeθ(boundj)
                    end
                end
                only_rangesi = filter(x -> x isa AbstractRange, rangeboundi)
                Ni = length(only_rangesi)
                input_cartprodi = cartesian_product(rangeboundi...; flat=Val{true})
                new_sizei = tuple(fill(num_perdim, num_dimsi)...)
                bc_valsi = reshape(bc_lfi(input_cartprodi, θ), new_sizei) .^ 2
                image = if Ni == 1
                    plot(only_rangesi..., bc_valsi, title=name)
                elseif Ni == 2
                    plot(reverse(only_rangesi)..., bc_valsi, linetype=:contourf, title=name)
                else
                    println("no plot defined for arrays with more than 2 axes")
                    nothing
                end

                # this is the bc_vals array we're looking for. this could be plotted 
                # (using bc_vars as the axes name, the index as the title, 
                #image = do_plot(loss_function, title=name)
                (; name, image)
            end

            if !(have_given_pybamm_plots[1])
                have_given_pybamm_plots[1] = true
                vcat(pybamm_plots, eval_plots, error_plots, abs_error_plots, pde_loss_plots, bc_loss_plots)
            else
                vcat(eval_plots, error_plots, abs_error_plots, pde_loss_plots, bc_loss_plots)
            end
        end
        return plot_function
    end
end
