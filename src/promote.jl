const SpinMonomialLike = Union{SpinVariable, SpinMonomial}

# Promotion with Term
Base.promote_rule(::Type{SpinTerm{S}}, ::Type{SpinTerm{T}}) where {S, T} = SpinTerm{promote_type(S, T)}
Base.promote_rule(::Type{<:SpinMonomialLike}, ::Type{SpinTerm{T}}) where {T} = SpinTerm{promote_type(T, Int)}
