"""
    local_state_index(s::Index, name::AbstractString) -> Int

Return the integer state index for the state named `name` on site index `s`,
as defined by ITensors for the site's type.
"""
local_state_index(s::Index, name::AbstractString) = ITensors.val(s, name)

"""
    local_state_string(s::Index, n::Int) -> String

Return the state name for integer index `n` on site index `s`.
This is the inverse of `local_state_index`.
"""
local_state_string(s::Index, n::Int) = local_state_strings(s)[n]

"""
    local_state_strings(s::Index) -> Vector{String}

Return all allowed state names for site index `s`, ordered by their ITensors
integer index. Supported site types: `"S=1/2"`, `"tJ"`, `"Electron"`.
"""
function local_state_strings(s::Index)
    if hastags(s, "S=1/2") || hastags(s, "Spinhalf")
        return ["Up", "Dn"]
    elseif hastags(s, "tJ")
        return ["Emp", "Up", "Dn"]
    elseif hastags(s, "Electron")
        return ["Emp", "Up", "Dn", "UpDn"]
    else
        error("Unsupported site type: $(tags(s))")
    end
end

"""
    local_state_integers(s::Index) -> Vector{Int}

Return the integer indices of all allowed states for site index `s`,
ordered consistently with `local_state_strings`.
"""
local_state_integers(s::Index) = [local_state_index(s, name) for name in local_state_strings(s)]
