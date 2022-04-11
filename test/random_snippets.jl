"""
model = pybamm.lithium_ion.SPM()
sim = pybamm.Simulation(model)
sim.solve([0, 3600])

sim.plot()

model.variable_names()

params = model.default_parameter_values


rhs = (keys=collect(keys(sim.model.rhs)), values=collect(values(sim.model.rhs)))
init_cond = (keys=collect(keys(sim.model.initial_conditions)), values=collect(values(sim.model.initial_conditions)))

first_rhs_val = rhs[:values][1]
first_rhs_key = rhs[:keys][1]
"""


# spm_mtk_str_strip = spm_mtk_str[7:end-14]


