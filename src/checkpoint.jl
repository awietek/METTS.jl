using HDF5
using HDF5

"""
    count_existing_steps(filename::String) -> Int

Returns the number of completed METTS steps already stored in the HDF5 file
by counting the entries in the `product_state` dataset. Returns 0 if the file
does not exist. Deletes the file and warns if it is corrupt or has no
`product_state` dataset.
"""
function count_existing_steps(filename::String)::Int
    isfile(filename) || return 0
    try
        result = h5open(filename, "r") do f
            if !haskey(f, "product_state")
                @warn "Checkpoint file exists but contains no `product_state` dataset, deleting and starting fresh." filename
                rm(filename)
                return 0
            end
            length(read(f["product_state"]))
        end
        return result
    catch e
        if isa(e, HDF5.Exception)
            @warn "Checkpoint file appears corrupted, please inspect it manually before continuing." filename exception=e
        else
            @warn "Checkpoint file is unavailable." filename exception=e
        end
        return 0
    end
end


"""
    read_last_product_state(filename::String)

Reads the last collapsed product state from the HDF5 file. Errors if
the file cannot be read or the dataset is empty.
"""
function read_last_product_state(filename::String)
    isfile(filename) || error("Checkpoint file '$filename' not found.")
    h5open(filename, "r") do f
        haskey(f, "product_state") || error("No `product_state` dataset found in '$filename'.")
        states = read(f["product_state"])
        isempty(states) && error("Empty `product_state` dataset in '$filename'.")
        states[end]
    end
end