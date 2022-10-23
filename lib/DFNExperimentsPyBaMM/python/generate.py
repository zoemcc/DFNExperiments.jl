import pybamm
import rich
import numpy as np
from icecream import ic
from rich import inspect
from IPython import embed

class BaseModel(pybamm.lithium_ion.BaseModel):
    def __init__(self, name):
        super().__init__({"timescale": 3600}, name)
        #super().__init__({}, name)
        #self.timescale = pybamm.Scalar(3600)
        self.variables = {
            "Time": pybamm.t,
            "x": pybamm.standard_spatial_vars.x,
            "x_n": pybamm.standard_spatial_vars.x_n,
            "x_p": pybamm.standard_spatial_vars.x_p,
            "r_n": pybamm.SpatialVariable("r_n", domain=["negative particle"]),
            "r_p": pybamm.SpatialVariable("r_p", domain=["positive particle"]),
        }

    def new_empty_copy(self):
        return pybamm.BaseModel.new_empty_copy(self)


def spm():
    model = pybamm.lithium_ion.SPM(name="SPM")
    model.variables.update(
        {
            "r_n": pybamm.SpatialVariable("r_n", domain=["negative particle"]),
            "r_p": pybamm.SpatialVariable("r_p", domain=["positive particle"]),
        }
    )
    variables = [
        "X-averaged negative particle concentration",
        "X-averaged positive particle concentration",
    ]
    return model, variables


def spm_no_r():
    model = pybamm.lithium_ion.SPM({"particle": "uniform profile"}, name="SPM_no_r")
    variables = [
        "Average negative particle concentration",
        "Average positive particle concentration",
    ]
    return model, variables


def spme():
    model = pybamm.lithium_ion.SPMe(name="SPMe")
    model.variables.update(
        {
            "r_n": pybamm.SpatialVariable("r_n", domain=["negative particle"]),
            "r_p": pybamm.SpatialVariable("r_p", domain=["positive particle"]),
        }
    )
    variables = [
        "X-averaged negative particle concentration",
        "X-averaged positive particle concentration",
        "Electrolyte concentration",
        #"Porosity times concentration"
    ]
    return model, variables


def dfn_no_r():
    model = pybamm.lithium_ion.DFN({"particle": "uniform profile"}, name="DFN_no_r")
    variables = [
        "R-averaged negative particle concentration",
        "R-averaged positive particle concentration",
        "Electrolyte concentration",
        "Electrolyte potential",
        "Negative electrode potential",
        "Positive electrode potential",
    ]
    return model, variables


def dfn():
    model = pybamm.lithium_ion.DFN(name="DFN")
    variables = [
        "Negative particle concentration",
        "Positive particle concentration",
        "Electrolyte concentration",
        "Electrolyte potential",
        "Negative electrode potential",
        "Positive electrode potential",
    ]
    return model, variables


def reduced_c():
    class Model(BaseModel):
        def __init__(self, name="reduced_c"):
            super().__init__(name)
            param = self.param
            #c_e = pybamm.Variable(
                #"Electrolyte concentration",
                #["negative electrode", "separator", "positive electrode"],
            #)
            #c_e.short_name = "c_e"
            #c_e.print_name = "c_e"
            c_e_n = pybamm.Variable(
                "Electrolyte concentration neg", "negative electrode"
            )
            c_e_s = pybamm.Variable("Electrolyte concentration sep", "separator")
            c_e_p = pybamm.Variable(
                "Electrolyte concentration pos", "positive electrode"
            )
            c_e_n.short_name = "c_e_n"
            c_e_n.print_name = "c_e_n"
            c_e_s.short_name = "c_e_s"
            c_e_s.print_name = "c_e_s"
            c_e_p.short_name = "c_e_p"
            c_e_p.print_name = "c_e_p"
            c_e = pybamm.concatenation(c_e_n, c_e_s, c_e_p)
            c_e.short_name = "c_e"
            c_e.print_name = "c_e"

            j = pybamm.concatenation(
                pybamm.PrimaryBroadcast(1 / param.l_n, "negative electrode"),
                pybamm.PrimaryBroadcast(0, "separator"),
                pybamm.PrimaryBroadcast(-1 / param.l_p, "positive electrode"),
            )
            self.rhs[c_e] = pybamm.div(pybamm.grad(c_e)) + j

            self.initial_conditions[c_e] = pybamm.Scalar(1)
            self.boundary_conditions[c_e] = {
                "left": (pybamm.Scalar(0), "Neumann"),
                "right": (pybamm.Scalar(0), "Neumann"),
            }
            #self.variables.update({"Electrolyte concentration": c_e})
            self.variables.update({
                "Electrolyte concentration neg": c_e_n,
                "Electrolyte concentration sep": c_e_s,
                "Electrolyte concentration pos": c_e_p,
                                   })
            ic(self.variables)

    model = Model()
    variables = list(model.variables.keys())
    return model, variables


def reduced_c_phi():
    class Model(BaseModel):
        def __init__(self, name="reduced_c_phi"):
            super().__init__(name)
            param = self.param
            #c_e = pybamm.Variable(
                #"Electrolyte concentration",
                #["negative electrode", "separator", "positive electrode"],
            #)
            #c_e.short_name = "c_e"
            #c_e.print_name = "c_e"
            c_e_n = pybamm.Variable(
                "Electrolyte concentration neg", "negative electrode"
            )
            c_e_s = pybamm.Variable("Electrolyte concentration sep", "separator")
            c_e_p = pybamm.Variable(
                "Electrolyte concentration pos", "positive electrode"
            )
            c_e_n.short_name = "c_e_n"
            c_e_n.print_name = "c_e_n"
            c_e_s.short_name = "c_e_s"
            c_e_s.print_name = "c_e_s"
            c_e_p.short_name = "c_e_p"
            c_e_p.print_name = "c_e_p"
            c_e = pybamm.concatenation(c_e_n, c_e_s, c_e_p)
            c_e.short_name = "c_e"
            c_e.print_name = "c_e"

            #phi_e = pybamm.Variable(
                #"Electrolyte potential",
                #["negative electrode", "separator", "positive electrode"],
            #)
            #phi_e.short_name = "phi_e"
            #phi_e.print_name = "phi_e"
            phi_e_n = pybamm.Variable("Negative electrolyte potential", "negative electrode")
            phi_e_s = pybamm.Variable("Separator electrolyte potential", "separator")
            phi_e_p = pybamm.Variable("Positive electrolyte potential", "positive electrode")
            phi_e_n.short_name = "phi_e_n"
            phi_e_n.print_name = "phi_e_n"
            phi_e_s.short_name = "phi_e_s"
            phi_e_s.print_name = "phi_e_s"
            phi_e_p.short_name = "phi_e_p"
            phi_e_p.print_name = "phi_e_p"
            phi_e = pybamm.concatenation(phi_e_n, phi_e_s, phi_e_p)
            phi_e.short_name = "phi_e"
            phi_e.print_name = "phi_e"

            j = pybamm.concatenation(
                pybamm.PrimaryBroadcast(1 / param.l_n, "negative electrode"),
                pybamm.PrimaryBroadcast(0, "separator"),
                pybamm.PrimaryBroadcast(-1 / param.l_p, "positive electrode"),
            )
            self.rhs[c_e] = pybamm.div(pybamm.grad(c_e)) + j
            self.algebraic[phi_e] = (
                pybamm.div(pybamm.grad(c_e) / c_e - pybamm.grad(phi_e)) - j
            )

            self.initial_conditions[c_e] = pybamm.Scalar(1)
            self.initial_conditions[phi_e] = pybamm.Scalar(0)
            self.boundary_conditions[c_e] = {
                "left": (pybamm.Scalar(0), "Neumann"),
                "right": (pybamm.Scalar(0), "Neumann"),
            }
            self.boundary_conditions[phi_e] = {
                "left": (pybamm.Scalar(0), "Dirichlet"),
                "right": (pybamm.Scalar(0), "Neumann"),
            }
            #self.variables.update(
                #{"Electrolyte concentration": c_e, "Electrolyte potential": phi_e}
            #)
            self.variables.update({
                "Electrolyte concentration neg": c_e_n,
                "Electrolyte concentration sep": c_e_s,
                "Electrolyte concentration pos": c_e_p,
                "Negative electrolyte potential": phi_e_n,
                "Separator electrolyte potential": phi_e_s,
                "Positive electrolyte potential": phi_e_p,
                                   })
            ic(self.variables)

    model = Model()
    variables = list(model.variables.keys())
    return model, variables


def reduced_c_phi_j():
    class Model(BaseModel):
        def __init__(self, name="reduced_c_phi_j"):
            super().__init__(name)
            param = self.param
            c_e_n = pybamm.Variable(
                "Electrolyte concentration neg", "negative electrode"
            )
            c_e_s = pybamm.Variable("Electrolyte concentration sep", "separator")
            c_e_p = pybamm.Variable(
                "Electrolyte concentration pos", "positive electrode"
            )
            c_e_n.short_name = "c_e_n"
            c_e_n.print_name = "c_e_n"
            c_e_s.short_name = "c_e_s"
            c_e_s.print_name = "c_e_s"
            c_e_p.short_name = "c_e_p"
            c_e_p.print_name = "c_e_p"
            c_e = pybamm.concatenation(c_e_n, c_e_s, c_e_p)
            c_e.short_name = "c_e"
            c_e.print_name = "c_e"

            phi_e_n = pybamm.Variable("Negative electrolyte potential", "negative electrode")
            phi_e_s = pybamm.Variable("Separator electrolyte potential", "separator")
            phi_e_p = pybamm.Variable("Positive electrolyte potential", "positive electrode")
            phi_e_n.short_name = "phi_e_n"
            phi_e_n.print_name = "phi_e_n"
            phi_e_s.short_name = "phi_e_s"
            phi_e_s.print_name = "phi_e_s"
            phi_e_p.short_name = "phi_e_p"
            phi_e_p.print_name = "phi_e_p"
            phi_e = pybamm.concatenation(phi_e_n, phi_e_s, phi_e_p)
            phi_e.short_name = "phi_e"
            phi_e.print_name = "phi_e"

            phi_s_n = pybamm.Variable("Negative electrode potential", "negative electrode")
            phi_s_p = pybamm.Variable("Positive electrode potential", "positive electrode")
            phi_s_n.short_name = "phi_s_n"
            phi_s_n.print_name = "phi_s_n"
            phi_s_p.short_name = "phi_s_p"
            phi_s_p.print_name = "phi_s_p"

            j_n = c_e_n ** 0.5 * pybamm.sinh(phi_s_n - phi_e_n + 1)
            j_p = c_e_p ** 0.5 * pybamm.sinh(phi_s_p - phi_e_p - 4)
            j = pybamm.concatenation(j_n, pybamm.PrimaryBroadcast(0, "separator"), j_p)

            self.rhs[c_e] = pybamm.div(pybamm.grad(c_e)) + j
            self.algebraic[phi_e] = (
                pybamm.div(pybamm.grad(c_e) / c_e - pybamm.grad(phi_e)) - j
            )
            self.algebraic[phi_s_n] = pybamm.div(pybamm.grad(phi_s_n)) + j_n
            self.algebraic[phi_s_p] = pybamm.div(pybamm.grad(phi_s_p)) + j_p

            self.initial_conditions[c_e] = pybamm.Scalar(1)
            self.initial_conditions[phi_e] = pybamm.Scalar(0)
            self.initial_conditions[phi_s_n] = pybamm.Scalar(0)
            self.initial_conditions[phi_s_p] = pybamm.Scalar(2)

            self.boundary_conditions[c_e] = {
                "left": (pybamm.Scalar(0), "Neumann"),
                "right": (pybamm.Scalar(0), "Neumann"),
            }
            self.boundary_conditions[phi_e] = {
                "left": (pybamm.Scalar(0), "Neumann"),
                "right": (pybamm.Scalar(0), "Neumann"),
            }
            self.boundary_conditions[phi_s_n] = {
                "left": (pybamm.Scalar(0), "Dirichlet"),
                "right": (pybamm.Scalar(0), "Neumann"),
            }
            self.boundary_conditions[phi_s_p] = {
                "left": (pybamm.Scalar(0), "Neumann"),
                "right": (pybamm.Scalar(1), "Neumann"),
            }

            self.variables.update(
                {
                    #"Electrolyte concentration": c_e,
                    #"Electrolyte potential": phi_e,
                    "Electrolyte concentration neg": c_e_n,
                    "Electrolyte concentration sep": c_e_s,
                    "Electrolyte concentration pos": c_e_p,
                    "Negative electrolyte potential": phi_e_n,
                    "Separator electrolyte potential": phi_e_s,
                    "Positive electrolyte potential": phi_e_p,
                    "Negative electrode potential": phi_s_n,
                    "Positive electrode potential": phi_s_p,
                    "Interfacial current density": j,
                }
            )
            ic(self.variables)

    model = Model()
    variables = list(model.variables.keys())
    return model, variables


def solve_plot_generate(model, variables, current_input=False, include_q=False, num_pts=100, num_tsteps=1000):
    #ic(model.rhs.keys())
    #ic(model.rhs.items())
    #ic(model.algebraic)
    #ic(include_q)
    rhs_to_pop = []
    ic_to_pop = []
    if include_q:
        common_vars = ["Time", "x", "x_n", "x_p", "r_n", "r_p", "Discharge capacity [A.h]"]
        
    else:
        common_vars = ["Time", "x", "x_n", "x_p", "r_n", "r_p"]
        for i, var in enumerate(model.rhs.keys()):
            if var.name == "Discharge capacity [A.h]":
                rhs_to_pop.append((var, model.rhs[var]))
        for popee, val in rhs_to_pop:
            #ic(popee)
            #ic(val)
            model.rhs.pop(popee)
        #ic("ic's")
        for i, var in enumerate(model.initial_conditions):
            if var.name == "Discharge capacity [A.h]":
                ic_to_pop.append((var, model.initial_conditions[var]))
        for popee, val in ic_to_pop:
            #ic(popee)
            #ic(val)
            model.initial_conditions.pop(popee)
        #ic("bc's")
        #for i, var in enumerate(model.boundary_conditions):
            #ic(var.name)
            #ic(model.boundary_conditions[var])
    #ic(common_vars)
    variables = list(set(common_vars + variables))
    #ic(variables)
    dep_vars = [x for x in variables if x not in common_vars]
    #ic(dep_vars)

    parameter_values = model.default_parameter_values
    # parameter_values["Electrolyte diffusivity [m2.s-1]"] = 1e-10
    # parameter_values["Electrolyte conductivity [S.m-1]"] = 1
    # parameter_values["Negative electrode exchange-current density [A.m-2]"] = 1e-6
    # parameter_values["Positive electrode exchange-current density [A.m-2]"] = 1e-6
    # parameter_values["Negative electrode OCP [V]"] = 0.5
    # parameter_values["Positive electrode OCP [V]"] = 4
    parameter_values._replace_callable_function_parameters = False
    if current_input is True:
        inputs = {"I": parameter_values["Current function [A]"]}
        parameter_values["Current function [A]"] = pybamm.InputParameter("I")
    else:
        inputs = {}

    var = pybamm.standard_spatial_vars
    var_pts = {var.x_n: num_pts, var.x_s: num_pts, var.x_p: num_pts, var.r_n: num_pts, var.r_p: num_pts}

    sim = pybamm.Simulation(model, var_pts=var_pts, parameter_values=parameter_values)
    sim.set_parameters()

    #embed()

    # Print MTK (to be copy-pasted to correct files)
    mtk_str = pybamm.get_julia_mtk_model(
        sim.model, geometry=sim.geometry, tspan=(0, 3600)
    )
    # save mtk_str to file
    outfile = "./models/playground/model.jl"
    ic(outfile)
    with open(outfile, "w") as f:
        ic("writing outfile")
        f.write(mtk_str)
        ic("writing successful")
    #print(mtk_str)


    # Solve
    parameter_values._replace_callable_function_parameters = True
    for pushee, val in rhs_to_pop:
        model.rhs[pushee] = val
    for pushee, val in ic_to_pop:
        model.initial_conditions[pushee] = val
    sim = pybamm.Simulation(model, var_pts=var_pts, parameter_values=parameter_values)
    #sim.solve([0, 3600], inputs=inputs)
    #sim.mesh.add_ghost_meshes()
    pybamm.set_logging_level("INFO")
    time_steps = np.linspace(0, 3600, num=num_tsteps)
    sim.solve(time_steps, inputs=inputs)


    # Plot
    # sim.plot(dep_vars)

    # Save to pickle
    #sim.solution.save_data(
        #"pybamm/hardcoded_models/MTK_format/pybamm_solutions/" + model.name + ".pickle",
        #variables,
    #)
    # solve_plot_generate(*spm_no_r(), current_input=False)
    # solve_plot_generate(*reduced_c(), current_input=False)
    # solve_plot_generate(*spme(), current_input=False)
    # solve_plot_generate(*reduced_c_phi(), current_input=False)
    # solve_plot_generate(*reduced_c_phi_j(), current_input=False)
    # solve_plot_generate(*dfn_no_r(), current_input=False)
    # solve_plot_generate(*dfn(), current_input=False)
    # solve_plot_generate(*dfn(), current_input=True)

    return sim, mtk_str, variables


#solve_plot_generate(*spm(), current_input=False, include_q=False)
#solve_plot_generate(*spme(), current_input=False)
#solve_plot_generate(*reduced_c(), current_input=False)
#solve_plot_generate(*reduced_c_phi(), current_input=False)
#solve_plot_generate(*dfn(), current_input=False)
