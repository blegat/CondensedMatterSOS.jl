struct SpinPolynomial{T} <: MP.AbstractPolynomial{T}
    terms::Vector{SpinTerm{T}}
end
MP.terms(p::SpinPolynomial) = p.terms
monomialtype(::Type{<:Union{SpinVariable, SpinMonomial, SpinTerm, SpinPolynomial}}) = SpinMonomial

#With this I solve 2*sx[1]<sx[1]
function SpinTerm{T}(spin::Union{SpinVariable, SpinMonomial}) where T
    return SpinTerm(one(T),monomial(spin))
end

#With this I solve sx[1]+sx[2]
function MP.polynomial(vterm::Array{SpinTerm{T},1}, s::MP.SortedUniqState) where T
    return SpinPolynomial(vterm);
end