module CondensedMatterSOS

using MultivariatePolynomials
using DataStructures

const MP = MultivariatePolynomials
include("types.jl")

# TODO remove when MP v0.3.9 is out
function Base.convert(::Type{TT}, t::TT) where {T, TT <: AbstractTerm{T}}
    return t
end
function Base.promote_rule(::Type{TT}, t::TT) where {T, TT <: AbstractTerm{T}}
    return t
end
# Needed in Julia v1.0 to avoid clashin with promotion.jl:236.
# See https://travis-ci.com/github/blegat/CondensedMatterSOS.jl/jobs/333686069
Base.promote_rule(::Type{PT}, ::Type{Any}) where {PT<:MP.APL} = promote_rule_constant(Any, PT)

function MP.name_base_indices(var::SpinVariable) # Used to print variable
    splits = split(NAMES[var.id], r"[\[,\]]\s*", keepempty=false)
    suffix = ("ˣ", "ʸ", "ᶻ")[var.index + 1]
    name = splits[1] * suffix
    if length(splits) == 1
        return name, Int[]
    else
        return name, parse.(Int, splits[2:end])
    end
end


include("operators.jl")

export @spin
end # module
