
function recursive_array_size(data) 
    array_size = []
    cur_data = data
    while true 
        if cur_data isa Vector
            push!(array_size, size(cur_data, 1))
            cur_data = cur_data[1]
        else
            break
        end
    end
    array_size
end

function recursive_get_index(recursive_data, index)
    if length(index) == 1
        return recursive_data[index[1]]
    else
        return recursive_get_index(recursive_data[index[1]], index[2:end])
    end
end

function recursive_array_to_array(recursive_data; do_reverse=true)
    array_size = recursive_array_size(recursive_data)
    array_size_possibly_reversed = do_reverse ? reverse(array_size) : array_size
    data_array = Array{Float64, length(array_size)}(undef, array_size_possibly_reversed...)
    indices = CartesianIndices(axes(data_array))
    for index in indices
        rec_index = do_reverse ? reverse(index.I) : index.I
        data_array[index] = recursive_get_index(recursive_data, rec_index)
    end
    data_array
end
