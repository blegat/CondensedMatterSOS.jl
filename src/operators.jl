
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
    println("term*variables");
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
        println("DENTRO IF")
        delete!(c.variables,a.id);
        site_mult = b.variables[a.id]*a; #var*var=term or Bool;
        if typeof(site_mult) == Bool
            return c;
        end
        push!(c.variable, a.id=>monomial(site_mult)[a.id]);
        return SpinTerm(coefficient(site_mult), c);
    else
        println("DENTRO ELSE")
        push!(c.variables, a.id=>a);
        return c;
    end
end


function Base.:*(a::SpinMonomial, b::SpinMonomial)
    c = deepcopy(a);
    sites_a = keys(a.variables);
    sites_b = keys(b.variables);
    same_sites = intersect(sites_a,sites_b);
    sda = setdiff(sites_a,sites_b);
    sdb = setdiff(sites_b,sites_a);
    for (site) in sdb
        push!(c.variables, site=>b.variables[site]);
    end
    if isempty(same_sites)
        return c;
    else
    coeff = 1;
        for (site) in (same_sites)
            delete!(c.variables, site)
            temp_term = a.variables[site]*b.variables[site]; #term or Bool
            if typeof(temp_term) != Bool
                coeff = coeff*coefficient(temp_term);
                push!(c.variables, site=>variables(monomial(temp_term))[1])
            end
        end
        return coeff*c;
    end
end
function Base.:*(b::SpinTerm,a::Number)
    return a*b;
end
function Base.:*(a::Number,b::SpinTerm)
    println("int*term")
    return SpinTerm(a*coefficient(b), monomial(b));
end
function Base.:*(b::SpinMonomial,a::Number)
    return a*b;
end
function Base.:*(a::Number,b::SpinMonomial)
    println("int*term")
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
