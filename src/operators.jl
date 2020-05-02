
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
                error("not valid index")
            end
        end
    else
        SpinMonomial([a,b]);
    end
end


function Base.:*(n::Number,spin::SpinVariable)
    return SpinTerm(n,SpinMonomial([spin]))
end



function Base.:*(b::SpinTerm,a::SpinVariable)
    return a*b;
end
function Base.:*(a::SpinVariable, b::SpinTerm)
    return coefficient(b)*(a*monomial(b));
end
function Base.:*(a::SpinMonomial, b::SpinVariable)
    return b*a;
end
function Base.:*(a::SpinVariable, b::SpinMonomial)
    c = deepcopy(b);
    if in(a.id,keys(b.variables))
        delete!(c.variables,a.id);
        site_mult = b.variables[a.id]*a; #var*var=term or Bool;
        if typeof(site_mult) == Bool
            return c;
        end
        push!(c.variable, a.id=>monomial(site_mult)[a.id]);
        return SpinTerm(coefficient(site_mult), c);
    else
        push!(c.variables, a.id=>a);
        return c;
    end
end


function Base.:*(a::SpinMonomial, b::SpinMonomial)
    c = deepcopy(a)
    coef = true
    for (site, variable) in b.variables
        c_variable = get(c.variables, site, nothing)
        if c_variable === nothing
            c.variables[site] = variable
        else
            term = c_variable * variable # It is either a `Term` or a `Bool`
            if term isa SpinTerm
                coef *= coefficient(term)
                c.variables[site] = first(variables(term))
            else
                coef *= term
                delete!(c.variables, site)
            end
        end
    end
    return coef * c
end
function Base.:*(b::SpinTerm,a::Number)
    return a*b;
end
function Base.:*(a::Number,b::SpinTerm)
    return SpinTerm(a*coefficient(b), monomial(b));
end
function Base.:*(b::SpinMonomial,a::Number)
    return a*b;
end
function Base.:*(a::Number,b::SpinMonomial)
    return SpinTerm(a, b);
end
function Base.:*(b::SpinTerm,a::SpinMonomial)
    return a*b;
end
function Base.:*(a::SpinMonomial, b::SpinTerm)
    return coefficient(b)*(a*monomial(b))
end
function Base.:*(a::SpinTerm, b::SpinTerm)
    return (coefficient(a)*coefficient(b))*(monomial(a)*monomial(b));
end
function Base.:(==)(a::SpinVariable, b::SpinVariable)
    return (a.id==b.id)&&(a.index==b.index)
end
function Base.:(==)(a::SpinMonomial, b::SpinMonomial)
    return a.variables==b.variables;
end
function Base.:(==)(a::SpinTerm, b::SpinTerm)
    return (a.monomial==b.monomial)&&(a.coefficient==b.coefficient);
end
function Base.:(==)(b::SpinTerm, a::SpinMonomial)
    return a==b;
end
function Base.:(==)(a::SpinMonomial,b::SpinTerm)
    c = isone(coefficient(b));
    return c&&(a==monomial(b));
end
