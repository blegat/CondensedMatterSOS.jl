module CondensedMatterSOS

using DataStructures
using Combinatorics

import MutableArithmetics
const MA = MutableArithmetics

import MultivariatePolynomials
const MP = MultivariatePolynomials

include("types.jl")

# Used to print variable
function _name_splits(var::SpinVariable, suffixes)
    splits = split(NAMES[var.id], r"[\[,\]]\s*", keepempty=false)
    suffix = suffixes[var.index + 1]
    name = splits[1] * suffix
    return name, splits
end
function MP.name(var::SpinVariable)
    name, splits = _name_splits(var, ["x", "y", "z"])
    if length(splits) == 1
        return name
    else
        return _spin_name(name, splits[2:end])
    end
end
function MP.name_base_indices(var::SpinVariable) # Used to print variable
    name, splits = _name_splits(var, ["ˣ", "ʸ", "ᶻ"])
    if length(splits) == 1
        return name, Int[]
    else
        return name, parse.(Int, splits[2:end])
    end
end

include("operators.jl")
include("monom-set.jl")

include("ising.jl")

include("symmetry.jl")
include("action.jl")
include("sos.jl")

export @spin
end # module
