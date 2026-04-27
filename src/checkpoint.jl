using HDF5
"""
    count_existing_steps(filename::String, tag::String) -> Int

Returns the number of completed METTS steps already stored in the HDF5 file
by counting the entries in the `product_state` dataset. Returns 0 if the file
does not exist. Throws an error if it is corrupt or has no
`product_state` dataset.
"""
function count_existing_steps(filename::String, product_state_name::String)::Int
    isfile(filename) || return 0 # in case the file does not exist
    try
        result = h5open(filename, "r") do f
            if !haskey(f, product_state_name)
                existing_datasets = keys(f)

                if isempty(existing_datasets)
                    @warn "Checkpoint file exists but is empty" filename
                    return 0
                else
                    found_keys_str = join(existing_datasets, ", ")
                    error("Checkpoint file '$filename' contains data ($found_keys_str) but is missing the mandatory '$product_state_name' dataset.")
                end
            end

            # Get dimensions
            dataset_dims = size(f[product_state_name])

            if length(dataset_dims) ∉ (1, 2)
                error("Data shape mismatch in '$filename': expected a 1D or 2D array for '$product_state_name', but got dimensions $(dataset_dims).")
            elseif length(dataset_dims) == 2
                return dataset_dims[2]
            else
                return 1
            end
            return result
        end
    catch e
        @error "Failed to count existing steps due to an error:" filename exception = (e, catch_backtrace())
        rethrow(e)
    end
end


"""
    read_last_product_state(filename::String)
Reads the last collapsed product state from the HDF5 file. Errors if
the file cannot be read or the dataset is empty.
"""
function read_last_product_state(filename::String, product_state_name::String)::Vector{Int}
    isfile(filename) || error("File '$filename' not found.")
    try
        product_state = h5open(filename, "r") do f
            if !haskey(f, product_state_name)
                error("No '$product_state_name' dataset found in '$filename'.")
            end

            dataset = f[product_state_name]
            dataset_dims = size(dataset)

            if 0 in dataset_dims || isempty(dataset_dims)
                error("Empty '$product_state_name' dataset in '$filename'.")
            end

            if length(dataset_dims) ∉ (1, 2)
                error("Data shape mismatch in '$filename': expected a 1D or 2D array for '$product_state_name', but got dimensions $(dataset_dims).")
            elseif length(dataset_dims) == 2
                return dataset[:, end]
            else
                return dataset[:]
            end
        end
        return product_state
    catch e
        @error "Failed to read product state due to an error:" filename exception = (e, catch_backtrace())
        rethrow(e)
    end
end