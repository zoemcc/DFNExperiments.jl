function initialize_pybamm_funcs()
    @pyinclude(joinpath(@__DIR__, "..", "python", "generate.py"))
end

