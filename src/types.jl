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
Base.convert(::Type{SpinMonomial}, spin::SpinVariable) = MP.monomial(spin)
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
function MP.powers(m::SpinMonomial)
    return MP.LazyMap{Tuple{SpinVariable,Int}}(Base.Fix2(tuple, 1), MP.variables(m))
end

MP.variables(spin::SpinMonomial) = collect(values(spin.variables))

# Use default `MP.Term`
MP.term_type(::Union{Type{M}}, ::Type{T}) where {M<:SpinMonomial,T} = MP.Term{T,M}

# Use default `MP.Polynomial` with `Vector`
MP.polynomial_type(::Type{MP.Term{C,M}}) where {C,M<:SpinMonomial} = MP.Polynomial{C,MP.Term{C,M},Vector{MP.Term{C, M}}}
MP.polynomial_type(::Type{MP.Term{C,M} where C}) where {M<:SpinMonomial} = MP.Polynomial

function MP.variables(monos::AbstractVector{SpinMonomial})
    vars = Set{SpinVariable}()
    for mono in monos
        union!(vars, MP.variables(mono))
    end
    return sort(collect(vars), rev=true)
end
MP.variables(p::MP.Polynomial{C,MP.Term{C,SpinMonomial}}) where {C} = MP.variables(MP.monomials(p))

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

"""
    @spin(σ[1:N1, 1:N2, ...], ...)

Return a tuple of 3-tuples `σ = (σx, σy, σz)` where the three elements of this tuple
are arrays of size `(N1, N2, ...)`. Moreover, the product `σ[i][I...] * σ[j][J...]`
* commutes if `I != J`,
* if `I == J`, it
  - commutes if `i == j`, moreover it satisfies the identity: `σx[I...] * σx[I...] == σy[I...] * σy[I...] == σz[I...] * σz[I...] = 1`.
  - anticommutes if `i != j`, moreover it satisifes the identities:
    * `σx * σy == -σy * σx == im * σz`;
    * `σy * σz == -σz * σy == im * σx`;
    * `σz * σx == -σx * σz == im * σy`.

It also sets the 3-tuple to the local Julia variable with name `σ`.

The macro can either be used to create several groups that commute with each other
```jldoctest
julia> using CondensedMatterSOS

julia> (σ1x, σ1y, σ1z), (σ2x, σ2y, σ2z) = @spin(σ1, σ2)
((σ1ˣ, σ1ʸ, σ1ᶻ), (σ2ˣ, σ2ʸ, σ2ᶻ))

julia> σ1x * σ2y
σ1ˣσ2ʸ

julia> σ1x * σ1y
(0 + 1im)σ1ᶻ

julia> σ2z * σ2y
(0 - 1im)σ2ˣ

julia> (σ1x + σ2y) * (σ2x + im * σ1z)
σ1ˣσ2ˣ + (0 + 1im)σ1ᶻσ2ʸ + σ1ʸ + (0 - 1im)σ2ᶻ
```

Or also a vector of groups commuting with each other. Note that it returns a 1-tuple
containing a 3-tuple hence the needed comma after `(σx, σy, σz)`

```jldoctest
julia> using CondensedMatterSOS

julia> (σx, σy, σz), = @spin(σ[1:2])
((CondensedMatterSOS.SpinVariable[σˣ₁, σˣ₂], CondensedMatterSOS.SpinVariable[σʸ₁, σʸ₂], CondensedMatterSOS.SpinVariable[σᶻ₁, σᶻ₂]),)

julia> σx[1] * σx[2]
σˣ₁σˣ₂

julia> σx[2] * σx[1]
σˣ₁σˣ₂

julia> σx[1] * σx[1]
(1 + 0im)

julia> σy[1] * σy[1]
(1 + 0im)

julia> σx[1] * σy[1]
(0 + 1im)σᶻ₁

julia> σy[1] * σx[1]
(0 - 1im)σᶻ₁

julia> σz[1] * σx[1]
(0 + 1im)σʸ₁
```
"""
macro spin(args...)
    # Variable vector x returned garanteed to be sorted so that if p is built with x then vars(p) == x
    vars, exprs = build_spins(args)
    :($(foldl((x,y) -> :($x; $y), exprs, init=:())); $(Expr(:tuple, esc.(vars)...)))
end

const SpinLike = Union{SpinVariable, SpinMonomial}
MP.variable_union_type(::Union{SpinLike, Type{<:SpinLike}}) = SpinVariable
MP.monomial_type(::Type{<:SpinLike}) = SpinMonomial
function MP.constant_monomial(::Union{SpinLike, Type{<:SpinLike}})
    return SpinMonomial(SpinVariable[])
end
MP.multiplication_preserves_monomial_order(::Type{SpinMonomial}) = false
