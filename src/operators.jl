function _coef_variable(a::SpinVariable, b::SpinVariable)
    i = a.index
    j = b.index
    pos = mod(i - j, 3)
    if pos == 2
        return 1im, SpinVariable(a.id, mod(i+pos,3))
    elseif pos == 1
        return -im, SpinVariable(a.id, mod(i+pos,3))
    else
        error("Invalid `index` field for variables, we should have 0 <= `a.index`, `b.index` < 3.")
    end
end
function Base.:*(a::SpinVariable, b::SpinVariable)
    if a.id == b.id
        if a.index == b.index
            return (1 + 0im) * constantmonomial(a)
            # We want to return `1` but in which type ?
            # We could use `Bool` type as it the type compatible with the most other types in Julia
            # but currently the convention is `Int` in MP, i.e. variables and monomials are `AbstractTerm{Int}`.
        else
            return prod(_coef_variable(a, b))
        end
    else
        return (1 + 0im) * SpinMonomial([a, b])
    end
end

function var_op!(op::Function, m::SpinMonomial, variable::SpinVariable)
    site = variable.id
    m_variable = get(m.variables, site, nothing)
    if m_variable === nothing
        coef = 1 + 0im
        m.variables[site] = variable
    else
        if m_variable.index == variable.index
            coef = 1 + 0im
            delete!(m.variables, site)
        else
            coef, m.variables[site] = op(m_variable, variable)
        end
    end
    return coef
end

function Base.:*(a::SpinMonomial, b::SpinVariable)
    c = deepcopy(a)
    coef = var_op!(_coef_variable, c, b)
    return coef * c
end
function Base.:*(a::SpinVariable, b::SpinMonomial)
    c = deepcopy(b)
    coef = var_op!((x, y) -> _coef_variable(y, x), c, a)
    return coef * c
end

function Base.:*(a::SpinMonomial, b::SpinMonomial)
    c = deepcopy(a)
    coef = 1 + 0im
    for variable in values(b.variables)
        coef *= var_op!(_coef_variable, c, variable)
    end
    return coef * c
end
MP.multconstant(α, m::SpinMonomial) = SpinTerm(α, m)
MP.multconstant(m::SpinMonomial, α) = SpinTerm(α, m)
function Base.:(==)(a::SpinVariable, b::SpinVariable)
    return (a.id==b.id) && (a.index==b.index)
end
function Base.:(==)(a::SpinMonomial, b::SpinMonomial)
    return a.variables == b.variables
end

# `σx * σx` is a `Bool` so we need to implement comparison with `Bool`s too.
Base.isless(a::Bool, b::SpinVariable) = true
Base.isless(a::SpinVariable, b::Bool) = false
function Base.isless(a::SpinVariable, b::SpinVariable)
    return a.id > b.id || (a.id == b.id && a.index > b.index)
end
function Base.isless(a::SpinVariable, b::SpinMonomial)
    return length(b.variables) > 1 ||
        (isone(length(b.variables)) && a < first(values(b.variables)))
end
function Base.isless(a::SpinMonomial, b::SpinVariable)
    return isconstant(a) ||
        (isone(length(a.variables)) && first(values(a.variables)) < b)
end
function Base.isless(a::SpinMonomial, b::SpinMonomial)
    dica = a.variables
    dicb = b.variables
    la = length(a.variables)
    lb = length(b.variables)
    if la != lb
        return la < lb
    end
    pa = startof(dica)
    pb = startof(dicb)
    # `la == lb` so we only don't need to check `pb != pastendsemitoken(dicb)`.
    while pa != pastendsemitoken(dica)
        ka, va = deref((dica, pa))
        kb, vb = deref((dicb, pb))
        if va < vb
            return true
        elseif vb < va
            return false
        end
        pa = advance((dica, pa))
        pb = advance((dicb, pb))
    end
    return false
end

combine_plus(t1::SpinTerm, t2::SpinTerm) = (coefficient(t1) + coefficient(t2)) * monomial(t1)
combine_plus(t::SpinTerm, ::MA.Zero) = MA.mutable_copy(t)
combine_plus(::MA.Zero, t::SpinTerm) = MA.mutable_copy(t)
combine_minus(t1::SpinTerm, t2::SpinTerm) = (coefficient(t1) - coefficient(t2)) * monomial(t1)
combine_minus(t::SpinTerm, ::MA.Zero) = MA.mutable_copy(t)
combine_minus(::MA.Zero, t::SpinTerm) = -t
function MA.promote_operation(::typeof(combine_plus), ::Type{SpinTerm{S}}, ::Type{SpinTerm{T}}) where {S, T}
    return SpinTerm{MA.promote_operation(+, S, T)}
end
function MA.promote_operation(::typeof(combine_minus), ::Type{SpinTerm{S}}, ::Type{SpinTerm{T}}) where {S, T}
    return SpinTerm{MA.promote_operation(-, S, T)}
end
compare(t1::SpinTerm, t2::SpinTerm) = monomial(t1) > monomial(t2)

# TODO taken from TypedPolynomials
join_terms_plus(terms1::AbstractArray{<:SpinTerm}, terms2::AbstractArray{<:SpinTerm}) = Sequences.mergesorted(terms1, terms2, compare, combine_plus)
join_terms_minus(terms1::AbstractArray{<:SpinTerm}, terms2::AbstractArray{<:SpinTerm}) = Sequences.mergesorted(terms1, terms2, compare, combine_minus)

Base.:+(p1::SpinPolynomial, p2::SpinPolynomial) = SpinPolynomial(join_terms_plus(terms(p1), terms(p2)))
Base.:+(p::SpinPolynomial, t::SpinTerm) = SpinPolynomial(join_terms_plus(terms(p), [t]))
Base.:+(t::SpinTerm, p::SpinPolynomial) = p + t
Base.:-(p1::SpinPolynomial, p2::SpinPolynomial) = SpinPolynomial(join_terms_minus(terms(p1), terms(p2)))
Base.:-(p::SpinPolynomial, t::SpinTerm) = SpinPolynomial(join_terms_minus(terms(p), [t]))
Base.:-(t::SpinTerm, p::SpinPolynomial) = SpinPolynomial(join_terms_minus([t], terms(p)))
