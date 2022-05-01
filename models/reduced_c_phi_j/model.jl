begin
using IfElse
# ('negative electrode', 'separator', 'positive electrode') -> x
# ('negative electrode',) -> x_n
# ('positive electrode',) -> x_p
@parameters t x x_n x_p
independent_variables_to_pybamm_names = Dict(
  :t => "Time",
  :x => "negative electrode",
  :x_n => "negative electrode",
  :x_p => "positive electrode",
)
# 'Electrolyte concentration' -> u1
# 'Electrolyte potential' -> u2
# 'Negative electrode potential' -> u3
# 'Positive electrode potential' -> u4
@variables u1(..) u2(..) u3(..) u4(..)
dependent_variables_to_pybamm_names = Dict(
  :u1 => "Electrolyte concentration",
  :u2 => "Electrolyte potential",
  :u3 => "Negative electrode potential",
  :u4 => "Positive electrode potential",
)
dependent_variables_to_dependencies = Dict(
  :u1 => (:t, :x),
  :u2 => (:t, :x),
  :u3 => (:t, :x_n),
  :u4 => (:t, :x_p),
)
Dt = Differential(t)
Dx = Differential(x)
Dx_n = Differential(x_n)
Dx_p = Differential(x_p)

# 'Electrolyte concentration' equation

function concatenation(x, n, s, p)
   # A concatenation in the electrolyte domain
   IfElse.ifelse(
      x < 0.4444444444444445, n, IfElse.ifelse(
         x < 0.5555555555555556, s, p
      )
   )
end

#cache_8017838820490660438 = concatenation(x, (u1(t, x) ^ 0.5) * (sinh((u3(t, x_n) - u2(t, x)) + 1.0)), 0.0, 
   #(u1(t, x) ^ 0.5) * (sinh((u4(t, x_p) - u2(t, x)) - 4.0))
#)
cache_8017838820490660438 = concatenation(x, (u1(t, x) ^ 0.5) * (sinh((u3(t, x) - u2(t, x)) + 1.0)), 0.0, 
   (u1(t, x) ^ 0.5) * (sinh((u4(t, x) - u2(t, x)) - 4.0))
)
cache_4764846323173313089 = (Dx(Dx(u1(t, x)))) + cache_8017838820490660438

# 'Electrolyte potential' equation
#cache_8017838820490660438 = concatenation(x, (u1(t, x) ^ 0.5) * (sinh((u3(t, x_n) - u2(t, x)) + 1.0)), 0.0, 
   #(u1(t, x) ^ 0.5) * (sinh((u4(t, x_p) - u2(t, x)) - 4.0))
#)
cache_8017838820490660438 = concatenation(x, (u1(t, x) ^ 0.5) * (sinh((u3(t, x) - u2(t, x)) + 1.0)), 0.0, 
   (u1(t, x) ^ 0.5) * (sinh((u4(t, x) - u2(t, x)) - 4.0))
)
cache_3358561361070727352 = (Dx(Dx(u1(t, x)) / u1(t, x) - Dx(u2(t, x)))) - cache_8017838820490660438

# 'Negative electrode potential' equation
#cache_m1434052808539415122 = (Dx_n(Dx_n(u3(t, x_n)))) + ((u1(t, x) ^ 0.5) * (sinh((u3(t, x_n) - u2(t, x)) + 1.0)))
cache_m1434052808539415122 = (Dx_n(Dx_n(u3(t, x_n)))) + ((u1(t, x_n) ^ 0.5) * (sinh((u3(t, x_n) - u2(t, x_n)) + 1.0)))

# 'Positive electrode potential' equation
#cache_m3885295840502231583 = (Dx_p(Dx_p(u4(t, x_p)))) + ((u1(t, x) ^ 0.5) * (sinh((u4(t, x_p) - u2(t, x)) - 4.0)))
cache_m3885295840502231583 = (Dx_p(Dx_p(u4(t, x_p)))) + ((u1(t, x_p) ^ 0.5) * (sinh((u4(t, x_p) - u2(t, x_p)) - 4.0)))


eqs = [
   Dt(u1(t, x)) ~ cache_4764846323173313089,
   0 ~ cache_3358561361070727352,
   0 ~ cache_m1434052808539415122,
   0 ~ cache_m3885295840502231583,
]


ics_bcs = [
   # initial conditions
   u1(0, x) ~ 1.0,
   u2(0, x) ~ 0.0,
   u3(0, x_n) ~ 0.0,
   u4(0, x_p) ~ 2.0,
   # boundary conditions
   Dx(u1(t, 0.0)) ~ 0.0,
   Dx(u1(t, 1.0)) ~ 0.0,
   Dx(u2(t, 0.0)) ~ 0.0,
   Dx(u2(t, 1.0)) ~ 0.0,
   u3(t, 0.0) ~ 0.0,
   Dx_n(u3(t, 0.4444444444444445)) ~ 0.0,
   Dx_p(u4(t, 0.5555555555555556)) ~ 0.0,
   Dx_p(u4(t, 1.0)) ~ 1.0,
]

t_domain = Interval(0.000,1.000)
x_domain = Interval(0.0, 1.0)
x_n_domain = Interval(0.0, 0.4444444444444445)
x_p_domain = Interval(0.5555555555555556, 1.0)

domains = [
   t in t_domain,
   x in x_domain,
   x_n in x_n_domain,
   x_p in x_p_domain,
]

ind_vars = [t, x, x_n, x_p]
dep_vars = [u1(t, x), u2(t, x), u3(t, x_n), u4(t, x_p)]

@named reduced_c_phi_j_pde_system = PDESystem(eqs, ics_bcs, domains, ind_vars, dep_vars)
pde_system = reduced_c_phi_j_pde_system

nothing
end
