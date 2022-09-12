begin
using IfElse
# ('negative particle',) -> r_n
# ('positive particle',) -> r_p
# ('negative electrode',) -> x_n
# ('separator',) -> x_s
# ('positive electrode',) -> x_p
@parameters t r_n r_p x_n x_s x_p
independent_variables_to_pybamm_names = Dict(
  :t => "Time",
  :r_n => "negative particle",
  :r_p => "positive particle",
  :x_n => "negative electrode",
  :x_s => "separator",
  :x_p => "positive electrode",
)
# 'X-averaged negative particle concentration' -> c_s_n_xav
# 'X-averaged positive particle concentration' -> c_s_p_xav
# 'Negative electrolyte concentration' -> c_e_n
# 'Separator electrolyte concentration' -> c_e_s
# 'Positive electrolyte concentration' -> c_e_p
@variables c_s_n_xav(..) c_s_p_xav(..) c_e_n(..) c_e_s(..) c_e_p(..)
dependent_variables_to_pybamm_names = Dict(
  :c_s_n_xav => "X-averaged negative particle concentration",
  :c_s_p_xav => "X-averaged positive particle concentration",
  :c_e_n => "Negative electrolyte concentration",
  :c_e_s => "Separator electrolyte concentration",
  :c_e_p => "Positive electrolyte concentration",
)
dependent_variables_to_dependencies = Dict(
  :c_s_n_xav => (:t, :r_n),
  :c_s_p_xav => (:t, :r_p),
  :c_e_n => (:t, :x_n),
  :c_e_s => (:t, :x_s),
  :c_e_p => (:t, :x_p),
)
Dt = Differential(t)
Dr_n = Differential(r_n)
Dr_p = Differential(r_p)
Dx_n = Differential(x_n)
Dx_s = Differential(x_s)
Dx_p = Differential(x_p)

# 'X-averaged negative particle concentration' equation
cache_6254820856749760283 = 8.813457647415216 * (1 / r_n^2 * Dr_n(r_n^2 * Dr_n(c_s_n_xav(t, r_n))))

# 'X-averaged positive particle concentration' equation
cache_8262395471222958549 = 22.598609352346717 * (1 / r_p^2 * Dr_p(r_p^2 * Dr_p(c_s_p_xav(t, r_p))))

# 'Electrolyte concentration' equation

function D_e(c_e, T)
   5.34e-10 * exp(-0.00065 * c_e) * exp(14.941767670181717 - (4454.8880308646785 * 1.0 / T))
end

cache_m8390082798472572148_n = -589429731.3893579
cache_m8390082798472572148_s = -3587155110.512914
cache_m8390082798472572148_p = -589429731.3893579
cache_m7799770279103888029_n = 
   0.1806862603261852 * x_n

cache_m7799770279103888029_s = 
   0.08030500458941565

cache_m7799770279103888029_p = 
   (0.007232292579356106 - (0.007232292579356106 * x_p)) / 0.04002679875215723

cache_9117425992584884052_n = 
   56.212339486148316

cache_9117425992584884052_s = 0.0
cache_9117425992584884052_p = 
   -56.212339486148316

cache_m7507122597998692872_n = 0.3
cache_m7507122597998692872_s = 1.0
cache_m7507122597998692872_p = 0.3
cache_m406707742878746567_n = (((Dx_n(((cache_m8390082798472572148_n * (D_e(c_e_n(t, x_n) * 1000.0, 298.15))) * Dx_n(c_e_n(t, x_n))) + cache_m7799770279103888029_n)) / -0.008035880643729006) + cache_9117425992584884052_n) / cache_m7507122597998692872_n
cache_m406707742878746567_s = (((Dx_s(((cache_m8390082798472572148_s * (D_e(c_e_s(t, x_s) * 1000.0, 298.15))) * Dx_s(c_e_s(t, x_s))) + cache_m7799770279103888029_s)) / -0.008035880643729006) + cache_9117425992584884052_s) / cache_m7507122597998692872_s
cache_m406707742878746567_p = (((Dx_p(((cache_m8390082798472572148_p * (D_e(c_e_p(t, x_p) * 1000.0, 298.15))) * Dx_p(c_e_p(t, x_p))) + cache_m7799770279103888029_p)) / -0.008035880643729006) + cache_9117425992584884052_p) / cache_m7507122597998692872_p
eqs = [
   Dt(c_s_n_xav(t, r_n)) ~ cache_6254820856749760283,
   Dt(c_s_p_xav(t, r_p)) ~ cache_8262395471222958549,
   Dt(c_e_n(t, x_n)) ~ cache_m406707742878746567_n,
   Dt(c_e_s(t, x_s)) ~ cache_m406707742878746567_s,
   Dt(c_e_p(t, x_p)) ~ cache_m406707742878746567_p,
]


ics_bcs = [
   # initial conditions
   c_s_n_xav(0, r_n) ~ 0.8000000000000016,
   c_s_p_xav(0, r_p) ~ 0.6000000000000001,
   c_e_n(0, x_n) ~ 1.0,
   c_e_s(0, x_s) ~ 1.0,
   c_e_p(0, x_p) ~ 1.0,
   # boundary conditions
   Dr_n(c_s_n_xav(t, 0.0)) ~ 0.0,
   Dr_n(c_s_n_xav(t, 1.0)) ~ -0.14182855923368468,
   Dr_p(c_s_p_xav(t, 0.0)) ~ 0.0,
   Dr_p(c_s_p_xav(t, 1.0)) ~ 0.03237700710041634,
   Dx_n(c_e_n(t, 0.0)) ~ 0.0,
   c_e_n(t, 0.4444444444444445) ~ c_e_s(t, 0.4444444444444445),
   c_e_p(t, 0.5555555555555556) ~ c_e_s(t, 0.5555555555555556),
   Dx_p(c_e_p(t, 1.0)) ~ 0.0,
]

t_domain = Interval(0.000,0.159)
r_n_domain = Interval(0.0, 1.0)
r_p_domain = Interval(0.0, 1.0)
x_n_domain = Interval(0.0, 0.4444444444444445)
x_s_domain = Interval(0.4444444444444445, 0.5555555555555556)
x_p_domain = Interval(0.5555555555555556, 1.0)

domains = [
   t in t_domain,
   r_n in r_n_domain,
   r_p in r_p_domain,
   x_n in x_n_domain,
   x_s in x_s_domain,
   x_p in x_p_domain,
]

ind_vars = [t, r_n, r_p, x_n, x_s, x_p]
dep_vars = [c_s_n_xav(t, r_n), c_s_p_xav(t, r_p), c_e_n(t, x_n), c_e_s(t, x_s), c_e_p(t, x_p)]

@named SPMe_pde_system = PDESystem(eqs, ics_bcs, domains, ind_vars, dep_vars)
pde_system = SPMe_pde_system

nothing
end
