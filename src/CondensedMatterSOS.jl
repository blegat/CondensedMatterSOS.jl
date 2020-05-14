module CondensedMatterSOS

using MultivariatePolynomials
using DataStructures

const MP = MultivariatePolynomials
include("types.jl")

# TODO remove when MP v0.3.9 is out
function Base.convert(::Type{TT}, t::TT) where {T, TT <: AbstractTerm{T}}
    return t
end
# Needed in Julia v1.0
function Base.promote_rule(TT::Type{<:AbstractTerm{T}}, ST::Type{<:AbstractTermLike{S}}) where {S, T}
    U = promote_type(S, T)
    UT = termtype(ST, U)
    if UT != termtype(TT, U)
        error("Cannot promote `$ST` and `$TT` to the same type.")
    end
    return UT
end

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
