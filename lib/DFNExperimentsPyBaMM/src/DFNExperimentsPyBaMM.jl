module DFNExperimentsPyBaMM

using DFNExperiments

include("generate_py.jl")

function __init__()
    #@require PyCall="438e738f-606a-5dbb-bf0a-cddfbfd45ab0" begin 
    begin
        initialize_pybamm_funcs()
    end
end


function DFNExperiments.pybamm_generate(model_str, current_input, include_q, num_pts, num_tsteps)
    current_input_str = current_input ? "True" : "False"
    include_q_str = include_q ? "True" : "False"
    sim, mtk_str, variables = py"solve_plot_generate(*$$(model_str)(), current_input=$$(current_input_str), include_q=$$(include_q_str), num_pts=$$(num_pts), num_tsteps=$$(num_tsteps))"
    return sim, mtk_str, variables
end




end
