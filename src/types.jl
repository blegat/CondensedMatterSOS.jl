const NAMES = String[]

struct SpinVariable <: MP.AbstractVariable
    id::Int # Spin id, spins with different id commute
    index::Int # 0 means x, 1 means y and 2 means z
end

function spin(name::String)
    push!(NAMES, name)
    id = length(NAMES)
    return SpinVariable(id, 0), SpinVariable(id, 1), SpinVariable(id, 2)
end



struct SpinMonomial <: MP.AbstractMonomial
    variables::SortedDict{Int, SpinVariable}
end
function SpinMonomial(vec::Vector{SpinVariable})
    variables   = SortedDict{Int, SpinVariable}()
    for spin in vec
        if in(keys(variables),spin.id)
            error("Monomial with repeated variable")
        end
        push!(variables, spin.id => spin);
    end
    return SpinMonomial(variables)
end

Base.copy(m::SpinMonomial) = SpinMonomial(copy(m.variables))
MP.isconstant(mono::SpinMonomial) = iszero(length(mono.variables))
MP.monomial(spin::SpinVariable) = SpinMonomial([spin])
function MP.exponents(spin::SpinMonomial)
   #There must be a 1to1 corresponendence between variables and exponents
   return ones(Int, length(spin.variables))
   # TODO It would be cheaper to return `FillArrays.Ones{Int}(length(spin.variables))`.
   #      but this is currently blocked by https://github.com/JuliaArrays/FillArrays.jl/issues/96
end
function MP.degree(spin::SpinMonomial, variable::SpinVariable)
    var = get(spin.variables, variable.id, nothing)
    return (var === nothing || var.index != variable.index) ? 0 : 1
end
function MP.powers(m::CondensedMatterSOS.SpinMonomial)
    # TODO maybe use MappedArrays.jl here so that ẁe can hardcode the `eltype`
    #      as `eltype(::Base.Generator)` is `Any`.
    return Base.Generator(v -> (v, 1), variables(m))
end

MP.variables(spin::SpinMonomial) = collect(values(spin.variables))

struct SpinTerm{T} <: MP.AbstractTerm{T}
    coefficient::T
    monomial::SpinMonomial
end

# TODO move to MP
MP.convertconstant(::Type{SpinTerm{T}}, α) where {T} = convert(T, α) * constantmonomial(SpinTerm{T})
Base.copy(t::SpinTerm) = SpinTerm(coefficient(t), copy(monomial(t)))
MA.mutable_copy(t::SpinTerm) = SpinTerm(MA.copy_if_mutable(coefficient(t)), copy(monomial(t)))

_spin_name(prefix::String, indices) = prefix * "[" * join(indices, ",") * "]"
function spin_index(prefix::String, indices)
    return spin(_spin_name(prefix, indices))
end

function array_spin(prefix, indices...)
    σs = map(i -> spin_index(prefix, i), Iterators.product(indices...))
    return [σ[1] for σ in σs], [σ[2] for σ in σs], [σ[3] for σ in σs]
end

function build_spin(var)
    if isa(var, Symbol)
        σx = Symbol(string(var) * "x")
        σy = Symbol(string(var) * "y")
        σz = Symbol(string(var) * "z")
        return var, :($(esc(var)) = ($(esc(σx)), $(esc(σy)), $(esc(σz))) = spin($"$var"))
    else
        isa(var, Expr) || error("Expected $var to be a variable name")
        Base.Meta.isexpr(var, :ref) || error("Expected $var to be of the form varname[idxset]")
        (2 ≤ length(var.args)) || error("Expected $var to have at least one index set")
        varname = var.args[1]
        prefix = string(varname)
        σx = Symbol(prefix * "x")
        σy = Symbol(prefix * "y")
        σz = Symbol(prefix * "z")
        return varname, :($(esc(varname)) = ($(esc(σx)), $(esc(σy)), $(esc(σz))) = array_spin($prefix, $(esc.(var.args[2:end])...)))
    end
end

function build_spins(args)
    vars = Symbol[]
    exprs = []
    for arg in args
        var, expr = build_spin(arg)
        push!(vars, var)
        push!(exprs, expr)
    end
    return vars, exprs
end

# Variable vector x returned garanteed to be sorted so that if p is built with x then vars(p) == x
macro spin(args...)
    vars, exprs = build_spins(args)
    :($(foldl((x,y) -> :($x; $y), exprs, init=:())); $(Expr(:tuple, esc.(vars)...)))
end



struct SpinPolynomial{T} <: MP.AbstractPolynomial{T}
    terms::Vector{SpinTerm{T}}
end
MP.terms(p::SpinPolynomial) = p.terms
MP.zero(::Type{SpinPolynomial{T}}) where {T} = SpinPolynomial(SpinTerm{T}[])

function MP.variables(monos::Vector{SpinMonomial})
    vars = Set{SpinVariable}()
    for mono in monos
        union!(vars, variables(mono))
    end
    return sort(collect(vars), rev=true)
end
MP.variables(p::SpinPolynomial) = MP.variables(MP.monomials(p))

# TODO move to MP
function MP.polynomial(m::Union{SpinMonomial, SpinVariable}, T::Type)
    return MP.polynomial(one(T) * m, T)
end
function MP.polynomial(t::SpinTerm{T}, ::Type{T}) where T
    return SpinPolynomial([t])
end
function MP.polynomial(t::SpinTerm, T::Type)
    return MP.polynomial(MP.changecoefficienttype(t, T), T)
end
function MP.polynomial(p::SpinPolynomial, T::Type)
    return SpinPolynomial(MP.changecoefficienttype.(terms(p), T))
end

const SpinLike = Union{SpinVariable, SpinMonomial, SpinTerm, SpinPolynomial}
MP.variable_union_type(::Union{SpinLike, Type{<:SpinLike}}) = SpinVariable
MP.monomialtype(::Type{<:SpinLike}) = SpinMonomial
function MP.constantmonomial(::Union{SpinLike, Type{<:SpinLike}})
    return SpinMonomial(SpinVariable[])
end
MP.termtype(::Union{SpinLike, Type{<:SpinLike}}, T::Type) = SpinTerm{T}
MP.polynomialtype(::Union{SpinLike, Type{<:SpinLike}}, T::Type) = SpinPolynomial{T}

# #With this I solve 2*sx[1]<sx[1]
# function SpinTerm{T}(spin::Union{SpinVariable, SpinMonomial}) where T
#     return SpinTerm(one(T),monomial(spin))
# end
MP.polynomial(terms::Vector{SpinTerm{T}}, ::MP.SortedUniqState) where {T} = SpinPolynomial{T}(terms)
MP.polynomial!(terms::Vector{SpinTerm{T}}, ::MP.SortedUniqState) where {T} = SpinPolynomial{T}(terms)
