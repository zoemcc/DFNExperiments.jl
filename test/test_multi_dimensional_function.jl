begin 
    using DFNExperiments
    using Symbolics
    using Flux, DiffEqFlux
    using LinearAlgebra
end


begin
x1s = 0:0.1:1
x1i = 0.1
x2s = 0:0.2:1
x2i = 0.4
ys = 0:0.1:0.2
@show size(cartesian_product(x1s, x2i, ys; flat=Val{false}))
@show size(cartesian_product(x1s, x2i, ys; flat=Val{true}))
@show size(cartesian_product(ys, ys; flat=Val{true}))
end

begin
ivs = (:x1, :x2)
dvs = (:u1,)
dv_deps = ((:x1, :x2),)
chain = FastChain(FastDense(2,8,Flux.tanh),FastDense(8,1))
initθ = DiffEqFlux.initial_params.(chain)
mdf = MultiDimensionalFunction(chain, ParameterizedMatrixApplyFuncType(), ivs, dvs, dv_deps)

da = mdf(initθ, ys, ys; flat=Val{false})
end

begin
ivs = (:x1, :x2)
dvs = (:u1, :u2)
dv_deps = ((:x1, :x2), (:x1, :x2))
chain = [FastChain(FastDense(2,8,Flux.tanh),FastDense(8,1)), FastChain(FastDense(2,6,Flux.tanh), FastDense(6,1))]
initθ = DiffEqFlux.initial_params.(chain)
flat_initθ = reduce(vcat, initθ)
mdf = MultiDimensionalFunction(chain, VectorOfParameterizedMDFApplyFuncType(), ivs, dvs, dv_deps)
i = 1
mdf.f[i].ivs
mdf.f[i].dvs
mdf.f[i].dv_deps

da = mdf(flat_initθ, ys, ys; flat=Val{false})
end

begin
ivs = (:x1, :x2)
dvs = (:u1, :u2, :u3)
dv_deps = ((:x2,), (:x1,), (:x1, :x2))
chain = [FastChain(FastDense(1,8,Flux.tanh),FastDense(8,1)), FastChain(FastDense(1,6,Flux.tanh), FastDense(6,1)), FastChain(FastDense(2,4,Flux.tanh), FastDense(4,1))]
initθ = DiffEqFlux.initial_params.(chain)
flat_initθ = reduce(vcat, initθ)
mdf = MultiDimensionalFunction(chain, VectorOfParameterizedMDFApplyFuncType(), ivs, dvs, dv_deps)

i = 3
mdf.f[i].ivs
mdf.f[i].dvs
mdf.f[i].dv_deps

da = mdf(flat_initθ, ys, ys; flat=Val{false})
end


"""
things I want to plot:

pde values and slices of them (heterogenous)



"""


"""
desired interface:

AbstractPDESolution <: AbstractNoTimeSolution
FunctionPDESolution <: AbstractPDESolution
PointwisePDESolution <: AbstractPDESolution

would like to treat both as a function and both as a grid

res::FunctionPDESolution = solve(prob, opt)
sample_slice()


"""



"""
const TYPE_LENGHT_LIMIT = Ref(20)
function Base.print_type_stacktrace(io, type; color=:normal)
str = first(sprint(show, type, context=io), TYPE_LENGHT_LIMIT[])
i = findfirst('{', str)
if isnothing(i) || !get(io, :backtrace, false)::Bool
printstyled(io, str; color=color)
else
printstyled(io, str[1:prevind(str,i)]; color=color)
printstyled(io, str[i:end]; color=:light_black)
end
end """
