using SumOfSquares
import Random
using GroupsCore

using PermutationGroups
function __symmetric_grp_gens(n::Integer)
    return (Perm(UInt16[2, 1], false), Perm(UInt16[2:n; 1], false))
end
symmetric_group(n::Integer) = PermGroup(__symmetric_grp_gens(n)...)

export Lattice1Group

## Klein C₂⊕C₂ group
struct KleinGroup <: Group end
"""
    struct KleinElement <: GroupElement
        id::Int
    end

The transformation:
* (σx, σy) → (-σx, -σy)
* (σy, σz) → (-σy, -σz)
* (σz, σx) → (-σz, -σx)
form a the Klein group.
"""
struct KleinElement <: GroupElement
    id::Int
end

GroupsCore.gens(::KleinGroup) = [KleinElement(1), KleinElement(2)]
GroupsCore.order(::Type{T}, G::KleinGroup) where {T} = convert(T, 4)
Base.IteratorSize(::Type{KleinGroup}) = Base.HasLength()
Base.eltype(::Type{KleinGroup}) = KleinElement
Base.iterate(::KleinGroup, st=0) = ifelse(st ≥ 4, nothing, (KleinElement(st), st + 1))
Base.one(::KleinGroup) = KleinElement(0)

Base.parent(::KleinElement) = KleinGroup()
GroupsCore.order(::Type{T}, el::KleinElement) where {T} = convert(T, ifelse(iszero(el.id), 1, 2))
Base.:(==)(a::KleinElement, b::KleinElement) = a.id == b.id
Base.inv(el::KleinElement) = el

function Base.:*(a::KleinElement, b::KleinElement)
    a.id > b.id && return b * a
    a.id == 0 && return b
    a.id == b.id && return one(a)
    a.id == 2 && return KleinElement(1)
    b.id == 2 && return KleinElement(3)
    return KleinElement(2)
end

function Base.rand(rng::Random.AbstractRNG, st::Random.SamplerTrivial{KleinGroup})
    klein = st[]
    id = rand(rng, 1:order(Int, klein))
    return KleinElement(id - 1)
end

## PermKleinGroup, aka the semi-direct product of Klein and S₃ over ϕ,
# the homomorphism ϕ: S₃ → Aut(Klein) given by the permutation
# of the three non-trivial elements in Klein.
struct KleinPermGroup <: Group end

"""
Group element `k * p = p * k^inv(p)`.
"""
struct KleinPermElement <: GroupElement
    p::PermutationGroups.Perm{UInt16}
    k::KleinElement
end

function GroupsCore.gens(::KleinPermGroup)
    Kid = one(KleinGroup())
    return [KleinPermElement(g, Kid) for g in __symmetric_grp_gens(3)]
end
GroupsCore.order(::Type{T}, G::KleinPermGroup) where {T} = convert(T, 6 * 4)
Base.IteratorSize(::Type{KleinPermGroup}) = Base.HasShape{2}()
Base.size(::KleinPermGroup) = (6, 4)
Base.eltype(::Type{KleinPermGroup}) = KleinPermElement

function Base.iterate(::KleinPermGroup)
    P3 = map(PermutationGroups.Perms.perm, symmetric_group(3))
    IT = Iterators.product(P3, KleinGroup())
    (p, k), st = iterate(IT)
    return KleinPermElement(p, k), (IT, st)
end

function Base.iterate(::KleinPermGroup, it_st)
    IT, st = it_st
    el_st = iterate(IT, st)
    isnothing(el_st) && return nothing
    (p, k), st = el_st
    return KleinPermElement(p, k), (IT, st)
end

Base.one(::KleinPermGroup) = KleinPermElement(one(Perm{UInt16}), one(KleinGroup()))

Base.parent(::KleinPermElement) = KleinPermGroup()
Base.isone(el::KleinPermElement) = isone(el.k) && isone(el.p)

function Base.hash(el::KleinPermElement, u::UInt64)
    return hash(el.k, hash(el.p, u))
end
function Base.:(==)(a::KleinPermElement, b::KleinPermElement)
    return a.p == b.p && a.k == b.k
end

function __klein_ϕ(k::KleinElement, p::PermutationGroups.AbstractPermutation)
    return ifelse(k.id == 0, k, KleinElement(k.id^p))
end

# k^(inv(p)) * inv(p) * k * p = k^(inv(p)) * k * inv(p) * p = 1
function Base.inv(el::KleinPermElement)
    inv_p = inv(el.p)
    KleinPermElement(inv_p, __klein_ϕ(el.k, inv_p))
end

# p * k * q * k' = p * q * k^p * k'
# k * p * k' * q = k * k'^p * p * q
function Base.:*(a::KleinPermElement, b::KleinPermElement)
    return KleinPermElement(a.p * b.p, __klein_ϕ(a.k, b.p) * b.k)
end

function Base.rand(rng::Random.AbstractRNG, ::Random.SamplerTrivial{KleinPermGroup})
    p = Perm{UInt16}(Random.randperm(rng, 3), false)
    k = rand(rng, KleinGroup())
    return KleinPermElement(p, k)
end

function Base.deepcopy_internal(g::KleinPermElement, dict::IdDict)
    p = if haskey(dict, objectid(g.p))
        dict[objectid(g.p)]
    else
        Base.deepcopy_internal(g.p, dict)
    end
    return KleinPermElement(p, g.k)
end

# CyclicGroup
struct CyclicGroup <: Group
    n::Int
end
struct CyclicElem <: GroupElement
    n::Int
    id::Int
end

GroupsCore.gens(c::CyclicGroup) = [CyclicElem(c.n, 1)]
GroupsCore.order(::Type{T}, c::CyclicGroup) where {T} = convert(T, c.n)
Base.IteratorSize(::Type{CyclicGroup}) = Base.HasLength()
Base.eltype(::Type{CyclicGroup}) = CyclicElem
Base.iterate(c::CyclicGroup, st=0) = ifelse(st ≥ c.n, nothing, (CyclicElem(c.n, st), st + 1))
Base.one(C::CyclicGroup) = CyclicElem(C.n, 0)

Base.parent(c::CyclicElem) = CyclicGroup(c.n)
GroupsCore.order(el::CyclicElem) = div(el.n, gcd(el.n, el.id))
Base.:(==)(a::CyclicElem, b::CyclicElem) = a.n == b.n && a.id == b.id
Base.inv(el::CyclicElem) = CyclicElem(el.n, (el.n - el.id) % el.n)
Base.:*(a::CyclicElem, b::CyclicElem) = CyclicElem(a.n, (a.id + b.id) % a.n)
Base.:^(el::CyclicElem, k::Integer) = CyclicElem(el.n, (el.id * k) % el.n)

function Base.rand(rng::Random.AbstractRNG, st::Random.SamplerTrivial{CyclicGroup})
    C = st[]
    id = rand(rng, 1:order(Int, C))
    return CyclicElem(C.n, id - 1)
end

# General direct sum of two groups
struct DirectSumGroup{G1<:Group,G2<:Group} <: Group
    H::G1
    K::G2
end
struct DirectSumElem{G1El<:GroupElement,G2El<:GroupElement} <: GroupElement
    h::G1El
    k::G2El
end
function GroupsCore.gens(G::DirectSumGroup)
    elts = [DirectSumElem(h, one(G.K)) for h in gens(G.H)]
    append!(elts, [DirectSumElem(one(G.H), k) for k in gens(G.K)])
    return elts
end
GroupsCore.order(::Type{T}, G::DirectSumGroup) where {T} =
    convert(T, order(T, G.H) * order(T, G.K))

function Base.IteratorSize(::Type{DirectSumGroup{H,K}}) where {H,K}
    if Base.IteratorSize(H) isa Base.IsInfinite || Base.IteratorSize(K) isa Base.IsInfinite
        return Base.IsInfinite()
    elseif Base.IteratorSize(H) isa Base.SizeUnknown || Base.IteratorSize(K) isa Base.SizeUnknown
        return Base.SizeUnknown()
    else
        return Base.HasShape{2}()
    end
end

Base.size(G::DirectSumGroup) = (order(Int, G.H), order(Int, G.K))

function Base.eltype(::Type{DirectSumGroup{H,K}}) where {H,K}
    return DirectSumElem{eltype(H),eltype(K)}
end

function Base.iterate(G::DirectSumGroup)
    IT = Iterators.product(G.H, G.K)
    el, st = iterate(IT)
    return DirectSumElem(el...), (IT, st)
end

function Base.iterate(::DirectSumGroup, it_st)
    el_st = iterate(it_st...)
    isnothing(el_st) && return nothing
    el, st = el_st
    return DirectSumElem(el...), (first(it_st), st)
end

Base.one(G::DirectSumGroup) = DirectSumElem(one(G.H), one(G.K))

Base.parent(g::DirectSumElem) = DirectSumGroup(parent(g.h), parent(g.k))
function GroupsCore.order(::Type{T}, g::DirectSumElem) where {T}
    return convert(T, lcm(order(T, g.h), order(T, g.k)))
end
Base.:(==)(g1::DirectSumElem, g2::DirectSumElem) = g1.h == g2.h && g1.k == g2.k
Base.hash(el::DirectSumElem, h::UInt) = hash(el.h, hash(el.k, hash(eltype(el), h)))
Base.inv(g::DirectSumElem) = DirectSumElem(inv(g.h), inv(g.k))
Base.:*(g1::DirectSumElem, g2::DirectSumElem) = DirectSumElem(g1.h * g2.h, g1.k * g2.k)
Base.:^(g::DirectSumElem, n::Integer) = DirectSumElem(g.h^n, g.k^n)

function Base.rand(rng::Random.AbstractRNG, st::Random.SamplerTrivial{<:DirectSumGroup})
    G = st[]
    h = rand(rng, G.H)
    k = rand(rng, G.K)
    return DirectSumElem(h, k)
end

function Base.deepcopy_internal(g::DirectSumElem, dict::IdDict)
    h = if isbits(g.h)
        g.h
    elseif haskey(dict, objectid(g.h))
        dict[objectid(g.h)]
    else
        Base.deepcopy_internal(g.h, dict)
    end

    k = if isbits(g.k)
        g.k
    elseif haskey(dict, objectid(g.k))
        dict[objectid(g.k)]
    else
        Base.deepcopy_internal(g.k, dict)
    end
    return DirectSumElem(h, k)
end

# two aliases
const Lattice1Group = DirectSumGroup{CyclicGroup,KleinPermGroup}
const DirectSum = DirectSumElem{CyclicElem,KleinPermElement}

"""
    Lattice1Group(n::Integer)

The direct sum of the cyclic group of order `n` and `KleinPermGroup()`.
"""
(::Type{Lattice1Group})(n::Integer) = DirectSumGroup(CyclicGroup(n), KleinPermGroup())
