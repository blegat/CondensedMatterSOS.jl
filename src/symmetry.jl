using SumOfSquares
using GroupsCore
using PermutationGroups

export Lattice1Group

struct KleinElement <: GroupElement
    id::Int
end
Base.:(==)(a::KleinElement, b::KleinElement) = a.id == b.id

PermutationGroups.order(el::KleinElement) = iszero(el.id) ? 1 : 2
Base.inv(el::KleinElement) = el

function Base.:*(a::KleinElement, b::KleinElement)
    a.id > b.id && return b * a
    a.id == 0 && return b
    a.id == b.id && return one(a)
    a.id == 2 && return KleinElement(1)
    b.id == 2 && return KleinElement(3)
    return KleinElement(2)
end

Base.conj(a::KleinElement, b::KleinElement) = inv(b) * a * b
Base.:^(a::KleinElement, b::KleinElement) = conj(a, b)

struct KleinGroup <: Group end
Base.one(::Union{KleinGroup, KleinElement}) = KleinElement(0)
PermutationGroups.gens(::KleinGroup) = [KleinElement(1), KleinElement(2)]
PermutationGroups.order(::Type{T}, G::KleinGroup) where {T} = convert(T, 4)
function Base.iterate(::KleinGroup, prev::KleinElement=KleinElement(-1))
    id = prev.id + 1
    if id > 4
        return nothing
    else
        next = KleinElement(id)
        return next, next
    end
end

#SymbolicWedderburn.conjugacy_classes_orbit(KleinGroup())

function perm_klein(k::KleinElement, p::Perm)
    if k.id == 0
        return k
    else
        return KleinElement(k.id^p)
    end
end

"""
Group element `k * p = p * k^inv(p)`.
"""
struct KleinPermElement <: GroupElement
    p::Perm{Int}
    k::KleinElement
end
Base.isone(el::KleinPermElement) = isone(el.k) && isone(el.p)

function Base.hash(el::KleinPermElement, u::UInt64)
    return hash(el.k, hash(el.p, u))
end
function Base.:(==)(a::KleinPermElement, b::KleinPermElement)
    return a.p == b.p && a.k == b.k
end

# k^(inv(p)) * inv(p) * k * p = k^(inv(p)) * k * inv(p) * p = 1
function Base.inv(el::KleinPermElement)
    inv_p = inv(el.p)
    KleinPermElement(inv_p, perm_klein(el.k, inv_p))
end

# p * k * q * k' = p * q * k^p * k'
# k * p * k' * q = k * k'^p * p * q
function Base.:*(a::KleinPermElement, b::KleinPermElement)
    return KleinPermElement(a.p * b.p, perm_klein(a.k, b.p) * b.k)
end
function Base.:^(el::KleinPermElement, k::Integer)
    return Base.power_by_squaring(el, k)
end

Base.conj(a::KleinPermElement, b::KleinPermElement) = inv(b) * a * b
Base.:^(a::KleinPermElement, b::KleinPermElement) = conj(a, b)

function PermutationGroups.order(el::KleinPermElement)
    cur = el
    i = 1
    while !isone(cur)
        i += 1
        cur *= el
    end
    return i
end

struct KleinPermGroup <: Group end
Base.one(::Union{KleinPermGroup, KleinPermElement}) = KleinPermElement(Perm(3), one(KleinGroup()))
# See https://github.com/blegat/CondensedMatterSOS.jl/pull/31#issuecomment-1717665659
symmetric_group(n) = PermutationGroups.PermGroup(PermutationGroups.Perm([2,1]), PermutationGroups.Perm([2:n;1]))
function PermutationGroups.gens(::KleinPermGroup)
    els = [KleinPermElement(Perm(3), KleinElement(1))]
    for g in gens(symmetric_group(3))
        push!(els, KleinPermElement(g, one(KleinGroup())))
    end
    return els
end
PermutationGroups.order(::Type{T}, G::KleinPermGroup) where {T} = convert(T, 6 * 4)
const P3 = collect(symmetric_group(3))
const IT = Iterators.product(1:6, 0:3)
function Base.iterate(::KleinPermGroup, args...)
    el_st = iterate(IT, args...)
    if el_st === nothing
        return nothing
    else
        el, st = el_st
        return KleinPermElement(P3[el[1]], KleinElement(el[2])), st
    end
end

struct CyclicElem <: GroupElement
    n::Int
    id::Int
end
Base.:(==)(a::CyclicElem, b::CyclicElem) = a.n == b.n && a.id == b.id
Base.inv(el::CyclicElem) = CyclicElem(el.n, (el.n - el.id) % el.n)

function Base.:*(a::CyclicElem, b::CyclicElem)
    return CyclicElem(a.n, (a.id + b.id) % a.n)
end
Base.:^(el::CyclicElem, k::Integer) = CyclicElem(el.n, (el.id * k) % el.n)

Base.conj(a::CyclicElem, b::CyclicElem) = inv(b) * a * b
Base.:^(a::CyclicElem, b::CyclicElem) = conj(a, b)

function PermutationGroups.order(el::CyclicElem)
    return div(el.n, gcd(el.n, el.id))
end

struct CyclicGroup <: Group
    n::Int
end
Base.one(c::Union{CyclicGroup, CyclicElem}) = CyclicElem(c.n, 0)
PermutationGroups.gens(c::CyclicGroup) = [CyclicElem(c.n, 1)]
PermutationGroups.order(::Type{T}, c::CyclicGroup) where {T} = convert(T, c.n)
function Base.iterate(c::CyclicGroup, prev::CyclicElem=CyclicElem(c.n, -1))
    id = prev.id + 1
    if id >= c.n
        return nothing
    else
        next = CyclicElem(c.n, id)
        return next, next
    end
end

#SymbolicWedderburn.conjugacy_classes_orbit(KleinPermGroup())

struct DirectSum <: GroupElement
    c::CyclicElem
    k::KleinPermElement
end
function PermutationGroups.order(el::DirectSum)
    return lcm(PermutationGroups.order(el.c), PermutationGroups.order(el.k))
end

function Base.hash(el::DirectSum, u::UInt64)
    return hash(el.k, hash(el.c, u))
end
function Base.:(==)(a::DirectSum, b::DirectSum)
    return a.c == b.c && a.k == b.k
end

function Base.inv(el::DirectSum)
    DirectSum(inv(el.c), inv(el.k))
end

# k * p * k' * q = k * k'^p * p * q
function Base.:*(a::DirectSum, b::DirectSum)
    return DirectSum(a.c * b.c, a.k * b.k)
end
Base.:^(a::DirectSum, k::Integer) = DirectSum(a.c^k, a.k^k)

Base.conj(a::DirectSum, b::DirectSum) = inv(b) * a * b
Base.:^(a::DirectSum, b::DirectSum) = conj(a, b)

struct Lattice1Group <: Group
    n::Int
end
Base.eltype(::Lattice1Group) = DirectSum
Base.one(el::DirectSum) = DirectSum(one(el.c), one(el.k))
Base.one(L::Lattice1Group) = DirectSum(CyclicElem(L.n, 0), one(KleinPermGroup()))
function PermutationGroups.gens(L::Lattice1Group)
    els = DirectSum[]
    for g in gens(CyclicGroup(L.n))
        push!(els, DirectSum(g, one(KleinPermGroup())))
    end
    for g in gens(KleinPermGroup())
        push!(els, DirectSum(one(CyclicGroup(L.n)), g))
    end
    return els
end
PermutationGroups.order(L::Lattice1Group) = PermutationGroups.order(Int, L)
function PermutationGroups.order(::Type{T}, L::Lattice1Group) where {T}
    return order(T, CyclicGroup(L.n)) * order(T, KleinPermGroup())
end
Base.length(L::Lattice1Group) = PermutationGroups.order(L)
function Base.iterate(L::Lattice1Group)
    el_p, st_p = iterate(CyclicGroup(L.n))
    el_k, st_k = iterate(KleinPermGroup())
    return DirectSum(el_p, el_k), ((el_p, st_p), st_k)
end
function Base.iterate(L::Lattice1Group, st)
    el_st_k = iterate(KleinPermGroup(), st[2])
    if el_st_k === nothing
        el_st_p = iterate(CyclicGroup(L.n), st[1][2])
        el_st_p === nothing && return nothing
        el_k, st_k = iterate(KleinPermGroup())
        return DirectSum(el_st_p[1], el_k), (el_st_p, st_k)
    end
    return DirectSum(st[1][1], el_st_k[1]), (st[1], el_st_k[2])
end
