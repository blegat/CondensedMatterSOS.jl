function Base.:*(a::SpinVariable, b::SpinVariable)
    if a.id == b.id
        if a.index == b.index
            return true # We want to return `1` but in which type ?
                        # We use `Bool` type as it the type compatible with the most other types in Julia.
        else
            i = a.index
            j = b.index
            pos = mod(i-j,3);
            if pos==2
                return  im*SpinVariable(a.id, mod(i+pos,3))
            elseif pos==1
                return  -im*SpinVariable(a.id, mod(i+pos,3))
            else
                error("Invalid `index` field for variables, we should have 0 <= `a.index`, `b.index` < 3.")
            end
        end
    else
        SpinMonomial([a,b]);
    end
end

function var_op!(op::Function, m::SpinMonomial, variable::SpinVariable)
    site = variable.id
    m_variable = get(m.variables, site, nothing)
    if m_variable === nothing
        coef = true
        m.variables[site] = variable
    else
        term = op(m_variable, variable) # It is either a `Term` or a `Bool`
        if term isa SpinTerm
            coef = coefficient(term)
            m.variables[site] = first(variables(term))
        else
            coef = term
            delete!(m.variables, site)
        end
    end
    return coef
end

function Base.:*(a::SpinMonomial, b::SpinVariable)
    c = deepcopy(a)
    coef = var_op!(*, c, b)
    return coef * c
end
function Base.:*(a::SpinVariable, b::SpinMonomial)
    c = deepcopy(b)
    coef = var_op!((x, y) -> y * x, c, a)
    return coef * c
end

function Base.:*(a::SpinMonomial, b::SpinMonomial)
    c = deepcopy(a)
    coef = true
    for variable in values(b.variables)
        coef *= var_op!(*, c, variable)
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
