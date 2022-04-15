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
    models = [SPMModel(), SPMeModel()]
    model = models[2]
    for model in models
    #begin
        models_dir = abspath(joinpath(@__DIR__, "..", "models"))
        model_str = pybamm_func_str(model)
        output_dir = joinpath(models_dir, "$(model_str)")
        model_file = joinpath(output_dir, "model.jl")
        current_input = false
        sim_data, pde_system = generate_sim_model(model; current_input=current_input, output_dir=output_dir)
        #sol_data, variables, sim = generate_sim_model(model; current_input=current_input, output_dir=output_dir)
        include(model_file)
    end
end

"""
dvs_sim = map(safehouse.dv_names) do dv_name
    data = safehouse.get_one_of_two(safehouse.sol_data, dependent_variables_to_pybamm_names[dv_name], string(dv_name))
    data_mat = reduce(hcat, data)
    # make sure time is first axis, this needs the check so we don't flip data that only depends on time
    if size(data_mat, 1) == 1
        data_mat = data_mat
    end
    data_mat
end

dv_i = 1
dv_name = safehouse.dv_names[dv_i]
data = safehouse.get_one_of_two(safehouse.sol_data, dependent_variables_to_pybamm_names[dv_name], string(dv_name))
array_size = get_array_size(data)
data_mat = reduce(hcat, data)
# make sure time is first axis, this needs the check so we don't flip data that only depends on time
if size(data_mat, 1) == 1
    data_mat = data_mat
end

data_array = Array{Float64, length(array_size)}(undef, array_size...)
indices = CartesianIndices(axes(data_array))
index = indices[2]

recursive_array_to_array(data)
"""

