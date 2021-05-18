"""
    monomials(vars::Vector{CondensedMatterSOS.SpinVariable}, deg::Int)

This version returns all the monomials of degree deg with all the possible orders
```
E1. (sx[1],sx[2])->[sx[1]*sx[2],sx[2]*sx[1]];
```
```
E2. (sx[1],sy[1])->[sx[1]*sy[1],sy[1]*sx[1]];
```
"""
function monomials(vars::Vector{CondensedMatterSOS.SpinVariable}, deg::Int)
    CC = combinations(vars,deg)
    BB = permutations.(collect(Iterators.take(CC,length(CC))))
    AA = prod.(vcat(map(x->collect(Iterators.take(x,length(x))), vcat(BB...))...))
    return unique(AA)
end

function _add_monomials!(monos, sites, current_mono, i, filter)
    if i > length(sites)
        if filter(current_mono)
            push!(monos, current_mono)
        end
    else
        for xyz in 0:2
            var = SpinVariable(sites[i], xyz)
            _add_monomials!(monos, sites, current_mono * var, i + 1, filter)
        end
    end
end
function _monomials(sites::Vector{Int}, deg::Int, filter::Function, monos, consecutive)
    if consecutive
        if deg == 0
            _add_monomials!(monos, Int[], MP.constantmonomial(SpinMonomial), 1, filter)
            return monos
        end
        if deg >= length(sites)
            if deg == length(sites)
                _add_monomials!(monos, sites, MP.constantmonomial(SpinMonomial), 1, filter)
            end
            return monos
        end
        for i in eachindex(sites)
            n = i+deg-1
            m = length(sites)
            active_sites = sites[i:min(n, m)]
            if n > m
                append!(active_sites, sites[1:(n-m)])
            end
            @assert length(active_sites) == deg
            _add_monomials!(monos, active_sites, MP.constantmonomial(SpinMonomial), 1, filter)
        end
    else
        for active_sites in combinations(sites, deg)
            _add_monomials!(monos, active_sites, MP.constantmonomial(SpinMonomial), 1, filter)
        end
    end
    return monos
end
function _monomials(vars::Vector{SpinVariable}, deg::Int, filter::Function, monos, consecutive)
    sites = unique!(sort!(map(var -> var.id, vars)))
    return _monomials(sites, deg, filter, monos, consecutive)
end

function MP.monomials(vars::Vector{SpinVariable}, deg::Int, filter::Function = x -> true; consecutive=false)
    return _monomials(vars, deg, filter, SpinMonomial[], consecutive)
end

function MP.monomials(vars::Vector{SpinVariable}, degs::AbstractVector{Int}, filter::Function = x -> true; consecutive=false)
    monos = SpinMonomial[]
    for deg in sort(degs)
        _monomials(vars, deg, filter, monos, consecutive)
    end
    return monos
end

# # using CondensedMatterSOS
# # function monomials(vars::Vector{CondensedMatterSOS.SpinVariable}, deg::Int)
# #
# #
# # using CondensedMatterSOS
# # using Combinatorics
# # using Flux: batch
# @spin s[1:4]
# # CC = combinations([sx[1],sy[1],sz[1]],2)
# #
# # prod.(collect(Iterators.take(CC,length(CC))))
# # collect(Iterators.take(CC,length(CC))) #Un vettore di vettori in cui ogni vettore è una combinazione
# # #ora devo permutare ognuna di queste combinazioni
# # AA = permutations.(collect(Iterators.take(CC,length(CC))))
# # for i in permutations.([[1,2,3],[4,5]])[1]
# #     println(i)
# # end
# #
# # collect(Iterators.take.(batch(AA)))
# # println("--")
# # println.(batch(AA))
# # length(AA)
# # AA
# #
# # prod.(batch(map(x->collect(Iterators.take(x,length(x))), batch(AA))))
# #
# monomials([sx[1],sy[1]],2)
# #
# #
# #
# # hash(sx[1]*sx[2])
# #
# # variables(sx[1]*sx[2])
# #
# # collect(values((sx[1]*sx[2]).variables))
# #
# # (sx[1]*sx[2]).variables
# unique([2*sx[1],sx[1]])
# hash(2*sx[1])
# hash(monomial(2*sx[1]))
#
# hash(sx[1]*sy[3])
# hash(collect(values(monomial(2*sx[1]).variables)))
#
#
#
# typeof(monomial(2*sx[1]))
#
# variable(monomial(sx[1]*sy[2]))
#  variable(monomial(sx[1]))
#
#
# hash(monomial(2*sx[1]*sy[2]))
# hash(monomial(2*sx[1]))
#
# degree(monomial(2*sx[1]*sy[2]))
# degree(monomial(2*sx[1]))
#
# powers(2*sx[1]*sy[2])
#
# CC = monomials([sx[1],sy[1]],2)
# CC[1] == im*sz[1]
# CC[2] == -im*sz[1]
# monomials([sx[1],sy[1]],2)==[im*sz[1],-im*sz[1]]
#
# map((x,y)->x==y, (CC,[im*sz[1],-im*sz[1]])))
#
# [im*sz[1],-im*sz[1]]
# typeof(im*sz[1])
# typeof(-im*sz[1])
