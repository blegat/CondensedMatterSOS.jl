module CondensedMatterSOS

using MultivariatePolynomials
using DataStructures

const MP = MultivariatePolynomials
# The names ares stored outside so that `isbits(::SpinVariable)` is `true`.
# const NAMES = String[]
#
# struct SpinVariable <: MP.AbstractVariable
#     id::Int # Spin id, spins with different id commute
#     index::Int # 0 means x, 1 means y and 2 means z
# end
#
# function spin(name::String)
#     push!(NAMES, name)
#     id = length(NAMES)
#     return SpinVariable(id, 0), SpinVariable(id, 1), SpinVariable(id, 2)
# end
#
#
#
# struct SpinMonomial <: MP.AbstractMonomial
#     variables::SortedDict{Int,SpinVariable};
#     function SpinMonomial(vec::Vector{SpinVariable})
#         variables   = SortedDict{Int, SpinVariable}()
#         for spin in vec
#             if in(keys(variables),spin.id)
#                 error("Monomial with repeated variable")
#             end
#             push!(variables, spin.id => spin);
#         end
#         new(variables)
#     end
# end
#
# function MP.exponents(spin::SpinMonomial)
#        #  #Exponents is a vector long as the the biggest site index of the spin
#        #  exponents = zeros(Int, last(spin.variables)[1])
#        #  for (key,value) in spin.variables
#        #      #if the spin is actually present, exponent=1 that is the maximum
#        #      exponents[key] = 1;
#        #  end
#        # return exponents;
#
#        #There must be a 1to1 corresponendence between variables and exponents
#        return ones(Int, length(spin.variables))
# end
#
# function MP.variables(spin::SpinMonomial)
#     var = [];
#     for (key,value) in spin.variables
#         push!(var,value)
#     end
#     return var;
# end
#
# struct SpinTerm{T} <: MP.AbstractTerm{T}
#     coefficient::T
#     monomial::SpinMonomial
# end
#
# function MP.monomial(term::SpinTerm)
#     return term.monomial
# end
#
# function MP.coefficient(term::SpinTerm)
#     return term.coefficient
# end
#
include("types.jl")
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

# function Base.:*(a::SpinVariable, b::SpinVariable)
#     if a.id == b.id
#         if a.index == b.index
#             return true # We want to return `1` but in which type ?
#                         # We use `Bool` type as it the type compatible with the most other types in Julia.
#         else
#             i = a.index
#             j = b.index
#             pos = mod(i-j,3);
#             if pos==2
#                 return  im*SpinVariable(a.id, mod(i+pos,3))
#             elseif pos==1
#                 return  -im*SpinVariable(a.id, mod(i+pos,3))
#             else
#                 error("not valid index")
#             end
#         end
#     else
#         SpinMonomial([a,b]);
#     end
# end
#
#
# function Base.:*(n::Number,spin::SpinVariable)
#     return SpinTerm(n,SpinMonomial([spin]))
# end

function spin_index(prefix::String, indices)
    return spin(prefix * "[" * join(indices, ",") * "]")
end

function array_spin(prefix, indices...)
    σs = map(i -> spin_index(prefix, i), Iterators.product(indices...))
    return [σ[1] for σ in σs], [σ[2] for σ in σs], [σ[3] for σ in σs]
end

# function build_spin(var)
#     if isa(var, Symbol)
#         σx = Symbol(string(var) * "x")
#         σy = Symbol(string(var) * "y")
#         σz = Symbol(string(var) * "z")
#         return [σx, σy, σz], :(($(esc(σx)), $(esc(σy)), $(esc(σz))) = spin($"$var"))
#     else
#         isa(var, Expr) || error("Expected $var to be a variable name")
#         Base.Meta.isexpr(var, :ref) || error("Expected $var to be of the form varname[idxset]")
#         (2 ≤ length(var.args)) || error("Expected $var to have at least one index set")
#         varname = var.args[1]
#         prefix = string(varname)
#         σx = Symbol(prefix * "x")
#         σy = Symbol(prefix * "y")
#         σz = Symbol(prefix * "z")
#         return [σx, σy, σz], :(($(esc(σx)), $(esc(σy)), $(esc(σz))) = array_spin($prefix, $(esc.(var.args[2:end])...)))
#     end
# end
#
# function build_spins(args)
#     vars = Symbol[]
#     exprs = []
#     for arg in args
#         var, expr = build_spin(arg)
#         append!(vars, var)
#         push!(exprs, expr)
#     end
#     return vars, exprs
# end
#
# # Variable vector x returned garanteed to be sorted so that if p is built with x then vars(p) == x
# macro spin(args...)
#     vars, exprs = build_spins(args)
#     :($(foldl((x,y) -> :($x; $y), exprs, init=:())); $(Expr(:tuple, esc.(vars)...)))
# end


#
#
# function Base.:*(b::SpinTerm,a::SpinVariable)
#     println("term*variables");
#     return a*b;
# end
# function Base.:*(a::SpinVariable, b::SpinTerm)
#     return coefficient(b)*(a*monomial(b));
# end
# function Base.:*(a::SpinMonomial, b::SpinVariable)
#     return b*a;
# end
# function Base.:*(a::SpinVariable, b::SpinMonomial)
#     c = deepcopy(b);
#     if in(a.id,keys(b.variables))
#         println("DENTRO IF")
#         delete!(c.variables,a.id);
#         site_mult = b.variables[a.id]*a; #var*var=term or Bool;
#         if typeof(site_mult) == Bool
#             return c;
#         end
#         push!(c.variable, a.id=>monomial(site_mult)[a.id]);
#         return SpinTerm(coefficient(site_mult), c);
#     else
#         println("DENTRO ELSE")
#         push!(c.variables, a.id=>a);
#         return c;
#     end
# end
#
#
# function Base.:*(a::SpinMonomial, b::SpinMonomial)
#     c = deepcopy(a);
#     sites_a = keys(a.variables);
#     sites_b = keys(b.variables);
#     same_sites = intersect(sites_a,sites_b);
#     sda = setdiff(sites_a,sites_b);
#     sdb = setdiff(sites_b,sites_a);
#     for (site) in sdb
#         push!(c.variables, site=>b.variables[site]);
#     end
#     if isempty(same_sites)
#         return c;
#     else
#     coeff = 1;
#         for (site) in (same_sites)
#             delete!(c.variables, site)
#             temp_term = a.variables[site]*b.variables[site]; #term or Bool
#             if typeof(temp_term) != Bool
#                 coeff = coeff*coefficient(temp_term);
#                 push!(c.variables, site=>variables(monomial(temp_term))[1])
#             end
#         end
#         return coeff*c;
#     end
# end
# function Base.:*(b::SpinTerm,a::Number)
#     return a*b;
# end
# function Base.:*(a::Number,b::SpinTerm)
#     println("int*term")
#     return SpinTerm(a*coefficient(b), monomial(b));
# end
# function Base.:*(b::SpinMonomial,a::Number)
#     return a*b;
# end
# function Base.:*(a::Number,b::SpinMonomial)
#     println("int*term")
#     return SpinTerm(a, b);
# end
# function Base.:*(b::SpinTerm,a::SpinMonomial)
#     return a*b;
# end
# function Base.:*(a::SpinMonomial, b::SpinTerm)
#     return coefficient(b)*(a*monomial(b))
# end
# function Base.:*(a::SpinTerm, b::SpinTerm)
#     return (coefficient(a)*coefficient(b))*(monomial(a)*monomial(b));
# end
# function Base.:(==)(a::SpinVariable, b::SpinVariable)
#     return (a.id==b.id)&&(a.index==b.index)
# end
# function Base.:(==)(a::SpinMonomial, b::SpinMonomial)
#     return a.variables==b.variables;
# end
# function Base.:(==)(a::SpinTerm, b::SpinTerm)
#     return (a.monomial==b.monomial)&&(a.coefficient==b.coefficient);
# end
# function Base.:(==)(b::SpinTerm, a::SpinMonomial)
#     return a==b;
# end
# function Base.:(==)(a::SpinMonomial,b::SpinTerm)
#     c = isone(coefficient(b));
#     return c&&(a==monomial(b));
# end


#
# @spin sigma[1:4];
# A = sigmax[1];
# B = sigmax[1]*sigmax[2];
# C = sigmax[3]*sigmax[4];
# D = 3*sigmax[1]*sigmax[2];
# E = 4*sigmax[3]*sigmax[4];
# B*C
# D*E
# A==A
# C==C
# println();
#
# A = -im*7*sigmax[2]*sigmaz[3];
# B = sigmax[1]*sigmax[2];
# C = sigmax[3]*sigmay[4];
# D = 3*sigmax[1]*sigmax[2];
#
# sigmax[1]*sigmax[1] == true
# sigmax[2]*sigmax[2] == true
# typeof(sigmax[1]*sigmax[2]) == SpinMonomial
# variables((sigmax[1]*sigmax[2])) == [sigmax[1], sigmax[2]]
# sigmax[1]*sigmay[1] == im*sigmaz[1]
# sigmay[1]*sigmax[1] == -im*sigmaz[1]
# C*D==3*sigmax[1]*sigmax[2]*sigmax[3]*sigmay[4]
# C == sigmax[3]*sigmay[4]
# D == 3*sigmax[1]*sigmax[2]
# B*D == 3
# B == sigmax[1]*sigmax[2]
# D == 3*sigmax[1]*sigmax[2]
# B*A*C == -7*sigmax[1]*sigmay[3]*sigmay[4]
# A == -im*7*sigmax[2]*sigmaz[3]
# B == sigmax[1]*sigmax[2]
# C == sigmax[3]*sigmay[4]
# (3*B)*A*C == -21*sigmax[1]*sigmay[3]*sigmay[4]
# A == -im*7*sigmax[2]*sigmaz[3]
# B == sigmax[1]*sigmax[2]
# C == sigmax[3]*sigmay[4]

export @spin
end # module
