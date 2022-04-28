function get_eval_network_at_sim_data_func(model, chains)
    sim_data = read_sim_data(model)
    ivs = keys(sim_data[:ivs])
    dvs = keys(sim_data[:dvs])
    dv_deps = collect(map(keys, sim_data[:dv_deps]))
    num_dvs = length(dvs)
    ground_truths = map(dv->sim_data[:dvs][dv], dvs)
    ground_truths_norm = LinearAlgebra.norm.(ground_truths)
    mdf = MultiDimensionalFunction(chains, VectorOfParameterizedMDFApplyFuncType(), ivs, dvs, dv_deps)

    function eval_network_at_sim_data(p)
        network_evals = mdf(p, values(sim_data[:ivs])...; flat=Val{false})

        errors = map(i -> network_evals[i] .- ground_truths[i], 1:num_dvs)
        errors_norm = LinearAlgebra.norm.(errors)
        errors_rel = map(i -> errors_norm[i] / ground_truths_norm[i], 1:num_dvs)

        return (dvs=dvs, network_evals=network_evals, errors=errors,
                errors_norm=errors_norm, errors_rel=errors_rel)
    end
end

function get_cb_func(model, log_frequency)
    function get_cb_func_inner(logger, iteration, chains)
        eval_network_at_sim_data = get_eval_network_at_sim_data_func(model, chains)

        cb = function (p,l)
            if iteration[1] % log_frequency == 0
                res = eval_network_at_sim_data(p)
                error_rel_max = max(res[:errors_rel]...)
                log_value(logger, "errors/max_error_relative", error_rel_max; step=iteration[1])
                for (i, dv) in enumerate(res[:dvs])
                    log_value(logger, "errors/$(string(dv))_error_norm", res[:errors_norm][i]; step=iteration[1])
                    log_value(logger, "errors/$(string(dv))_error_relative", res[:errors_rel][i]; step=iteration[1])
                end
            end
            return false
        end
        return cb
    end
end

function do_plot(data_array::AbstractArray{R, N}, iv_data, deps, name) where {R, N}
    used_ivs = map(iv->iv_data[iv], deps)
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
    function get_plot_function_inner(logger, iteration, chains)
        eval_network_at_sim_data = get_eval_network_at_sim_data_func(model, chains)
        sim_data = read_sim_data(model)
        dvs = keys(sim_data[:dvs])
        dv_deps = collect(map(keys, sim_data[:dv_deps]))
        num_dvs = length(dvs)
        iv_data = sim_data[:ivs]
        pybamm_plots = map(i->(name="$(string(dvs[i]))_pybamm", image=do_plot(sim_data[:dvs][i], iv_data, dv_deps[i], "$(string(dvs[i]))_pybamm")), 1:num_dvs)
        have_given_pybamm_plots = [false]

        function plot_function(θ, adaloss)
            res = eval_network_at_sim_data(θ)
            eval_plots = map(i->(name="$(string(dvs[i]))_nn", image=do_plot(res[:network_evals][i], iv_data, dv_deps[i], "$(string(dvs[i]))_nn")), 1:num_dvs)
            error_plots = map(i->(name="$(string(dvs[i]))_error", image=do_plot(res[:errors][i], iv_data, dv_deps[i], "$(string(dvs[i]))_error")), 1:num_dvs)
            abs_error_plots = map(i->(name="$(string(dvs[i]))_abs_error", image=do_plot(abs.(res[:errors][i]), iv_data, dv_deps[i], "$(string(dvs[i]))_abs_error")), 1:num_dvs)

            if !(have_given_pybamm_plots[1])
                have_given_pybamm_plots[1] = true
                vcat(pybamm_plots, eval_plots, error_plots, abs_error_plots)
            else
                vcat(eval_plots, error_plots, abs_error_plots)
            end
        end
        return plot_function
    end
end
