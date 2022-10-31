begin
using IfElse
# ('negative electrode',) -> x_n
# ('separator',) -> x_s
# ('positive electrode',) -> x_p
@parameters t x_n x_s x_p
independent_variables_to_pybamm_names = Dict(
  :t => "Time",
  :x_n => "negative electrode",
  :x_s => "separator",
  :x_p => "positive electrode",
)
# 'R-averaged negative particle concentration' -> c_s_n_rav
# 'R-averaged positive particle concentration' -> c_s_p_rav
# 'Negative electrolyte concentration' -> c_e_n
# 'Separator electrolyte concentration' -> c_e_s
# 'Positive electrolyte concentration' -> c_e_p
# 'Negative electrode potential' -> phi_s_n
# 'Positive electrode potential' -> phi_s_p
# 'Negative electrolyte potential' -> phi_e_n
# 'Separator electrolyte potential' -> phi_e_s
# 'Positive electrolyte potential' -> phi_e_p
@variables c_s_n_rav(..) c_s_p_rav(..) c_e_n(..) c_e_s(..) c_e_p(..) phi_s_n(..) phi_s_p(..) phi_e_n(..) phi_e_s(..) phi_e_p(..)
dependent_variables_to_pybamm_names = Dict(
  :c_s_n_rav => "R-averaged negative particle concentration",
  :c_s_p_rav => "R-averaged positive particle concentration",
  :c_e_n => "Negative electrolyte concentration",
  :c_e_s => "Separator electrolyte concentration",
  :c_e_p => "Positive electrolyte concentration",
  :phi_s_n => "Negative electrode potential",
  :phi_s_p => "Positive electrode potential",
  :phi_e_n => "Negative electrolyte potential",
  :phi_e_s => "Separator electrolyte potential",
  :phi_e_p => "Positive electrolyte potential",
)
dependent_variables_to_dependencies = Dict(
  :c_s_n_rav => (:t, :x_n),
  :c_s_p_rav => (:t, :x_p),
  :c_e_n => (:t, :x_n),
  :c_e_s => (:t, :x_s),
  :c_e_p => (:t, :x_p),
  :phi_s_n => (:t, :x_n),
  :phi_s_p => (:t, :x_p),
  :phi_e_n => (:t, :x_n),
  :phi_e_s => (:t, :x_s),
  :phi_e_p => (:t, :x_p),
)
Dt = Differential(t)
Dx_n = Differential(x_n)
Dx_s = Differential(x_s)
Dx_p = Differential(x_p)

# 'R-averaged negative particle concentration' equation

function j0_n(c_e, c_s_surf, T)
   2e-05 * exp(15.119261670583445 - (4507.807867084453 * 1.0 / T)) * (c_e ^ 0.5) * (c_s_surf ^ 0.5) * ((24983.2619938437 - c_s_surf) ^ 0.5)
end

function U_n(sto)
   0.194 + 1.5 * exp(-120.0 * sto) + 0.0351 * tanh((sto - 0.286) / 0.083) - (0.0045 * tanh((sto - 0.849) / 0.119)) - (0.035 * tanh((sto - 0.9233) / 0.05)) - (0.0147 * tanh((sto - 0.5) / 0.034)) - (0.102 * tanh((sto - 0.194) / 0.142)) - (0.022 * tanh((sto - 0.9) / 0.0164)) - (0.011 * tanh((sto - 0.124) / 0.0226)) + 0.0155 * tanh((sto - 0.105) / 0.029)
end

cache_8471184973576659264 = -1.6666666666666667 * ((((2.0 * max(1.0,1e-08)) / 0.5925925925925926) * (j0_n(max(1e-08,c_e_n(t, x_n)) * 1000.0, (max(1e-08,min(c_s_n_rav(t, x_n),0.99999999)) * 24983.2619938437), 298.15))) * (sinh(0.5 * ((phi_s_n(t, x_n) - phi_e_n(t, x_n)) - (((U_n(c_s_n_rav(t, x_n)) + ((1e-06 * (1.0 / c_s_n_rav(t, x_n))) + (1e-06 * (1.0 / (c_s_n_rav(t, x_n) - 1.0))))) - 0.175189434028335) / 0.025692579121493725)))))

# 'R-averaged positive particle concentration' equation

function j0_p(c_e, c_s_surf, T)
   6e-07 * exp(15.962358172491644 - (4759.177089128383 * 1.0 / T)) * (c_e ^ 0.5) * (c_s_surf ^ 0.5) * ((51217.9257309275 - c_s_surf) ^ 0.5)
end

function U_p(sto)
   2.16216 + 0.07645 * tanh(30.834 - (57.858397200000006 * sto)) + 2.1581 * tanh(52.294 - (53.412228 * sto)) - (0.14169 * tanh(11.0923 - (21.0852666 * sto))) + 0.2051 * tanh(1.4684 - (5.829105600000001 * sto)) + 0.2531 * tanh((-1.062 * sto + 0.56478) / 0.1316) - (0.02167 * tanh(((1.062 * sto) - 0.525) / 0.006))
end

cache_m6738177367957728228 = -0.9755671139472863 * ((((2.0 * max(1.0,1e-08)) * 1.4062500000000002) * (j0_p(max(1e-08,c_e_p(t, x_p)) * 1000.0, (max(1e-08,min(c_s_p_rav(t, x_p),0.99999999)) * 51217.9257309275), 298.15))) * (sinh(0.5 * ((phi_s_p(t, x_p) - phi_e_p(t, x_p)) - (((U_p(c_s_p_rav(t, x_p)) + ((1e-06 * (1.0 / c_s_p_rav(t, x_p))) + (1e-06 * (1.0 / (c_s_p_rav(t, x_p) - 1.0))))) - 4.027013014008729) / 0.025692579121493725)))))

# 'Electrolyte concentration' equation

function D_e(c_e, T)
   5.34e-10 * exp(-0.00065 * c_e) * exp(14.941767670181717 - (4454.8880308646785 * 1.0 / T))
end

function kappa_e(c_e, T)
   (0.0911 + 0.0019100999999999999 * c_e - (1.052e-06 * (c_e ^ 2.0)) + 1.554e-10 * (c_e ^ 3.0)) * exp(13.9978223044089 - (4173.450720059513 * 1.0 / T))
end

cache_6689922603821238635_n = -589429731.3893579
cache_6689922603821238635_s = -3587155110.512914
cache_6689922603821238635_p = -589429731.3893579
cache_2850548738456284222_n = 0.7818002858515763
cache_2850548738456284222_s = 4.757885022498838
cache_2850548738456284222_p = 0.7818002858515763
cache_m2900951515657641564_n = 
   ((((2.0 * max(1.0,1e-08)) / 0.5925925925925926) * (j0_n(max(1e-08,c_e_n(t, x_n)) * 1000.0, (max(1e-08,min(c_s_n_rav(t, x_n),0.99999999)) * 24983.2619938437), 298.15))) * (sinh(0.5 * ((phi_s_n(t, x_n) - phi_e_n(t, x_n)) - (((U_n(c_s_n_rav(t, x_n)) + ((1e-06 * (1.0 / c_s_n_rav(t, x_n))) + (1e-06 * (1.0 / (c_s_n_rav(t, x_n) - 1.0))))) - 0.175189434028335) / 0.025692579121493725))))) / 0.04002679875215723

cache_m2900951515657641564_s = 0.0
cache_m2900951515657641564_p = 
   ((((2.0 * max(1.0,1e-08)) * 1.4062500000000002) * (j0_p(max(1e-08,c_e_p(t, x_p)) * 1000.0, (max(1e-08,min(c_s_p_rav(t, x_p),0.99999999)) * 51217.9257309275), 298.15))) * (sinh(0.5 * ((phi_s_p(t, x_p) - phi_e_p(t, x_p)) - (((U_p(c_s_p_rav(t, x_p)) + ((1e-06 * (1.0 / c_s_p_rav(t, x_p))) + (1e-06 * (1.0 / (c_s_p_rav(t, x_p) - 1.0))))) - 4.027013014008729) / 0.025692579121493725))))) / 0.04002679875215723

cache_m2387532842543811942_n = 0.3
cache_m2387532842543811942_s = 1.0
cache_m2387532842543811942_p = 0.3
cache_6342662750096404049_n = (((Dx_n(((cache_6689922603821238635_n * (D_e(c_e_n(t, x_n) * 1000.0, 298.15))) * Dx_n(c_e_n(t, x_n))) + (0.08030500458941565 * (((kappa_e(c_e_n(t, x_n) * 1000.0, 298.15)) * cache_2850548738456284222_n) * (((1.2 * Dx_n(c_e_n(t, x_n))) / c_e_n(t, x_n)) - Dx_n(phi_e_n(t, x_n))))))) / -0.008035880643729006) + cache_m2900951515657641564_n) / cache_m2387532842543811942_n
cache_6342662750096404049_s = (((Dx_s(((cache_6689922603821238635_s * (D_e(c_e_s(t, x_s) * 1000.0, 298.15))) * Dx_s(c_e_s(t, x_s))) + (0.08030500458941565 * (((kappa_e(c_e_s(t, x_s) * 1000.0, 298.15)) * cache_2850548738456284222_s) * (((1.2 * Dx_s(c_e_s(t, x_s))) / c_e_s(t, x_s)) - Dx_s(phi_e_s(t, x_s))))))) / -0.008035880643729006) + cache_m2900951515657641564_s) / cache_m2387532842543811942_s
cache_6342662750096404049_p = (((Dx_p(((cache_6689922603821238635_p * (D_e(c_e_p(t, x_p) * 1000.0, 298.15))) * Dx_p(c_e_p(t, x_p))) + (0.08030500458941565 * (((kappa_e(c_e_p(t, x_p) * 1000.0, 298.15)) * cache_2850548738456284222_p) * (((1.2 * Dx_p(c_e_p(t, x_p))) / c_e_p(t, x_p)) - Dx_p(phi_e_p(t, x_p))))))) / -0.008035880643729006) + cache_m2900951515657641564_p) / cache_m2387532842543811942_p

# 'Negative electrode potential' equation
cache_m8687180816311917355 = (-Dx_n(221.12651346369245 * Dx_n(phi_s_n(t, x_n)))) + ((((2.0 * max(1.0,1e-08)) / 0.5925925925925926) * (j0_n(max(1e-08,c_e_n(t, x_n)) * 1000.0, (max(1e-08,min(c_s_n_rav(t, x_n),0.99999999)) * 24983.2619938437), 298.15))) * (sinh(0.5 * ((phi_s_n(t, x_n) - phi_e_n(t, x_n)) - (((U_n(c_s_n_rav(t, x_n)) + ((1e-06 * (1.0 / c_s_n_rav(t, x_n))) + (1e-06 * (1.0 / (c_s_n_rav(t, x_n) - 1.0))))) - 0.175189434028335) / 0.025692579121493725)))))

# 'Positive electrode potential' equation
cache_8719332745503186028 = (-Dx_p(16.821663817574187 * Dx_p(phi_s_p(t, x_p)))) + ((((2.0 * max(1.0,1e-08)) * 1.4062500000000002) * (j0_p(max(1e-08,c_e_p(t, x_p)) * 1000.0, (max(1e-08,min(c_s_p_rav(t, x_p),0.99999999)) * 51217.9257309275), 298.15))) * (sinh(0.5 * ((phi_s_p(t, x_p) - phi_e_p(t, x_p)) - (((U_p(c_s_p_rav(t, x_p)) + ((1e-06 * (1.0 / c_s_p_rav(t, x_p))) + (1e-06 * (1.0 / (c_s_p_rav(t, x_p) - 1.0))))) - 4.027013014008729) / 0.025692579121493725)))))

# 'Electrolyte potential' equation
cache_2850548738456284222_n = 0.7818002858515763
cache_2850548738456284222_s = 4.757885022498838
cache_2850548738456284222_p = 0.7818002858515763
cache_m4198028591950717818_n = 
   (((2.0 * max(1.0,1e-08)) / 0.5925925925925926) * (j0_n(max(1e-08,c_e_n(t, x_n)) * 1000.0, (max(1e-08,min(c_s_n_rav(t, x_n),0.99999999)) * 24983.2619938437), 298.15))) * (sinh(0.5 * ((phi_s_n(t, x_n) - phi_e_n(t, x_n)) - (((U_n(c_s_n_rav(t, x_n)) + ((1e-06 * (1.0 / c_s_n_rav(t, x_n))) + (1e-06 * (1.0 / (c_s_n_rav(t, x_n) - 1.0))))) - 0.175189434028335) / 0.025692579121493725))))

cache_m4198028591950717818_s = 0.0
cache_m4198028591950717818_p = 
   (((2.0 * max(1.0,1e-08)) * 1.4062500000000002) * (j0_p(max(1e-08,c_e_p(t, x_p)) * 1000.0, (max(1e-08,min(c_s_p_rav(t, x_p),0.99999999)) * 51217.9257309275), 298.15))) * (sinh(0.5 * ((phi_s_p(t, x_p) - phi_e_p(t, x_p)) - (((U_p(c_s_p_rav(t, x_p)) + ((1e-06 * (1.0 / c_s_p_rav(t, x_p))) + (1e-06 * (1.0 / (c_s_p_rav(t, x_p) - 1.0))))) - 4.027013014008729) / 0.025692579121493725))))

cache_m7636102912382120110_n = (Dx_n(((kappa_e(c_e_n(t, x_n) * 1000.0, 298.15)) * cache_2850548738456284222_n) * (((1.2 * Dx_n(c_e_n(t, x_n))) / c_e_n(t, x_n)) - Dx_n(phi_e_n(t, x_n))))) - cache_m4198028591950717818_n
cache_m7636102912382120110_s = (Dx_s(((kappa_e(c_e_s(t, x_s) * 1000.0, 298.15)) * cache_2850548738456284222_s) * (((1.2 * Dx_s(c_e_s(t, x_s))) / c_e_s(t, x_s)) - Dx_s(phi_e_s(t, x_s))))) - cache_m4198028591950717818_s
cache_m7636102912382120110_p = (Dx_p(((kappa_e(c_e_p(t, x_p) * 1000.0, 298.15)) * cache_2850548738456284222_p) * (((1.2 * Dx_p(c_e_p(t, x_p))) / c_e_p(t, x_p)) - Dx_p(phi_e_p(t, x_p))))) - cache_m4198028591950717818_p
eqs = [
   Dt(c_s_n_rav(t, x_n)) ~ cache_8471184973576659264,
   Dt(c_s_p_rav(t, x_p)) ~ cache_m6738177367957728228,
   Dt(c_e_n(t, x_n)) ~ cache_6342662750096404049_n,
   Dt(c_e_s(t, x_s)) ~ cache_6342662750096404049_s,
   Dt(c_e_p(t, x_p)) ~ cache_6342662750096404049_p,
   0 ~ cache_m8687180816311917355,
   0 ~ cache_8719332745503186028,
   0 ~ cache_m7636102912382120110_n,
   0 ~ cache_m7636102912382120110_s,
   0 ~ cache_m7636102912382120110_p,
]


ics_bcs = [
   # initial conditions
   c_s_n_rav(0, x_n) ~ 0.8000000000000016,
   c_s_p_rav(0, x_p) ~ 0.6000000000000001,
   c_e_n(0, x_n) ~ 1.0,
   c_e_s(0, x_s) ~ 1.0,
   c_e_p(0, x_p) ~ 1.0,
   # boundary conditions
   phi_s_n(t, 0.0) ~ 0.0,
   Dx_n(phi_s_n(t, 0.4444444444444445)) ~ 0.0,
   Dx_p(phi_s_p(t, 0.5555555555555556)) ~ 0.0,
   Dx_p(phi_s_p(t, 1.0)) ~ -0.05944715165186363,
   Dx_n(c_e_n(t, 0.0)) ~ 0.0,
   c_e_n(t, 0.4444444444444445) ~ c_e_s(t, 0.4444444444444445),
   c_e_p(t, 0.5555555555555556) ~ c_e_s(t, 0.5555555555555556),
   Dx_p(c_e_p(t, 1.0)) ~ 0.0,
   Dx_n(phi_e_n(t, 0.0)) ~ 0.0,
   phi_e_n(t, 0.4444444444444445) ~ phi_e_s(t, 0.4444444444444445),
   phi_e_p(t, 0.5555555555555556) ~ phi_e_s(t, 0.5555555555555556),
   Dx_p(phi_e_p(t, 1.0)) ~ 0.0,
]

t_domain = Interval(0.000,0.159)
x_n_domain = Interval(0.0, 0.4444444444444445)
x_s_domain = Interval(0.4444444444444445, 0.5555555555555556)
x_p_domain = Interval(0.5555555555555556, 1.0)

domains = [
   t in t_domain,
   x_n in x_n_domain,
   x_s in x_s_domain,
   x_p in x_p_domain,
]

ind_vars = [t, x_n, x_s, x_p]
dep_vars = [c_s_n_rav(t, x_n), c_s_p_rav(t, x_p), c_e_n(t, x_n), c_e_s(t, x_s), c_e_p(t, x_p), phi_s_n(t, x_n), phi_s_p(t, x_p), phi_e_n(t, x_n), phi_e_s(t, x_s), phi_e_p(t, x_p)]

@named DFN_no_r_pde_system = PDESystem(eqs, ics_bcs, domains, ind_vars, dep_vars)
pde_system = DFN_no_r_pde_system

nothing
end
