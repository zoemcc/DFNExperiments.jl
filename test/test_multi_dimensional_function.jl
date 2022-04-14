begin 
    #using DFNExperiments
    #using ModelingToolkit, DomainSets
    using Symbolics
    #import ModelingToolkit: Interval, infimum, supremum
    #using LabelledArrays
    #using NeuralPDE
    #using JSON
    #using Random
    using Flux, DiffEqFlux
    using LinearAlgebra
    #using Plots
    #using GalacticOptim
end

abstract type AbstractApplyFuncType end

struct ParameterizedMatrixApplyFuncType <: AbstractApplyFuncType end
struct VectorOfParameterizedMDFApplyFuncType <: AbstractApplyFuncType end

struct MultiDimensionalFunction{Func, FType <: AbstractApplyFuncType, Extra, IVsNum, DVsNum, DVDepsTuple} 
    f::Func
    ftype::FType
    extra_data::Extra
    ivs::NTuple{IVsNum, Symbol}
    dvs::NTuple{DVsNum, Symbol}
    dv_deps::DVDepsTuple 

    function MultiDimensionalFunction(f, ftype::FType, ivs::NTuple{IVsNum, Symbol}, dvs::NTuple{DVsNum, Symbol}, dv_deps::Tuple) where
        {FType <: AbstractApplyFuncType, IVsNum, DVsNum}

        @assert length(dv_deps) == DVsNum
        for dv_dep in dv_deps
            @assert dv_dep isa Tuple
            @assert length(dv_dep) <= IVsNum
            for iv in dv_dep
                @assert iv isa Symbol
            end
        end

        if FType == VectorOfParameterizedMDFApplyFuncType
            f_types = typeof.(f)
            @assert length(f) == DVsNum
            all_fastchains = all(map(ft_i->ft_i <: FastChain, f_types))
            if all_fastchains
                initθ = DiffEqFlux.initial_params.(chain)
                param_lengths = length.(initθ)
                param_lengths_cumulative = cumsum(vcat([0], param_lengths...)) 
                param_indices = [param_lengths_cumulative[i] + 1 : param_lengths_cumulative[i+1] for i in 1:length(param_lengths)]
                extra = Dict(:param_indices => param_indices)
                fs_mdf = map(1:DVsNum) do i
                    MultiDimensionalFunction(f[i], ParameterizedMatrixApplyFuncType(), dv_deps[i], dvs[i:i], dv_deps[i:i])
                end
                f = fs_mdf
            else
                throw("currently only supports vector of FastChains")
            end

        else
            extra = Nothing
        end

        new{typeof(f), FType, typeof(extra), IVsNum, DVsNum, typeof(dv_deps)}(f, ftype, extra, ivs, dvs, dv_deps)
    end
end

# Neural network

@generated function (f::MultiDimensionalFunction{Func, VectorOfParameterizedMDFApplyFuncType})(θ, x::Union{<:Real, AbstractVector{<:Real}}...; flat=Val{false}) where {Func}
    quote
        # slice out the portion of θ that is appropriate and then call the sub MDF with that slice on the appropriate data. currently only supports homogenous
        param_indices = f.extra_data[:param_indices]
        f_evals = Vector{Array{eltype(θ)}}(undef, length(f.dvs))
        for i in 1:length(f.dvs)
            param_indices_i = param_indices[i]
            @show param_indices_i
            f_evals[i] = f.f[i]((@view θ[param_indices_i]), x...; flat=flat)
        end
        return f_evals
    end
end

@generated function (f::MultiDimensionalFunction{Func, ParameterizedMatrixApplyFuncType})(θ, x::Union{<:Real, AbstractVector{<:Real}}...; flat=Val{false}) where {Func}
    quote
        # need to make a 2D array with the cartesian product of all the broadcasted 
        # and then apply the function to the point array
        # and reshape into the correct shape if flat is nothing
        if flat == Val{false}
            point_array, resize_info = cartesian_product(x...; flat=Val{true}, resize_info=Val{true}, debug=Val{false}) 
            transformed_array = f.f(point_array, θ)
            reshaped_transformed_array = reshape(transformed_array, :, resize_info...)
            return reshaped_transformed_array
        elseif flat == Val{true}
            point_array = cartesian_product(x...; flat=Val{true}, resize_info=Val{false}, debug=Val{false}) 
            transformed_array = f.f(point_array, θ)
            return transformed_array
        else
            throw("Flat must be either Val{false} or Val{true}")
        end
    end
end


@generated function cartesian_product(x::Union{<:Real, AbstractVector{<:Real}}...; flat=Val{false}, debug=Val{false}, resize_info=Val{false}) 
    broadcast_dims = map(x_i->x_i <: AbstractVector, x)
    broadcast_indices = Set(map(first, filter(i_b->i_b[2], collect(enumerate(broadcast_dims)))))
    N = length(x)
    eltypes = map(x_i -> x_i <: AbstractVector ? eltype(x_i) : x_i, x)
    promotion_type = promote_type(eltypes...)
    quote
        # need to make a 2D array with the cartesian product of all the broadcasted 
        broadcast_lengths = length.(x)
        num_points = prod(broadcast_lengths)
        iter_indices = CartesianIndices(tuple(broadcast_lengths...))

        if flat == Val{false} 
            point_array = Array{$promotion_type, ($N + 1)}(undef, $N, broadcast_lengths...) # TODO: make the AbstractArray type be more general
            for (i, index) in enumerate(iter_indices)
                for j in 1:$N
                    if j in $(broadcast_indices)
                        point_array[j, index] = x[j][index.I[j]]
                    else
                        point_array[j, index] = x[j]
                    end
                end
            end
        elseif flat == Val{true}
            point_array = Array{$promotion_type, 2}(undef, $N, num_points) # TODO: make the AbstractArray type be more general
            for (i, index) in enumerate(iter_indices)
                for j in 1:$N
                    if j in $(broadcast_indices)
                        point_array[j, i] = x[j][index.I[j]]
                    else
                        point_array[j, i] = x[j]
                    end
                end
            end
        else
            throw("Flat must be either Val{false} or Val{true}")
        end

        if debug == Val{false}
            if resize_info == Val{false}
                return point_array
            elseif resize_info == Val{true}
                return (point_array=point_array, resize_shape=broadcast_lengths)
            else
                throw("resize_info must be either Val{false} or Val{true}")
            end
        elseif debug == Val{true}
            return (N=$N, broadcast_dims=$(broadcast_dims), broadcast_indices=$(broadcast_indices), broadcast_lengths=broadcast_lengths, 
                num_points=num_points, point_array=point_array, shape=size(point_array), iter_indices=iter_indices)
        else
            throw("Debug must be either Val{false} or Val{true}")
        end
    end

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

da = mdf(flat_initθ, ys, ys; flat=Val{false})
end
nothing


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
