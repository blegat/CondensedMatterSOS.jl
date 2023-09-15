import SumOfSquares.Symmetry
export Action

struct Action <: Symmetry.OnMonomials
    σ
end
Symmetry.SymbolicWedderburn.coeff_type(::Action) = Float64
function Symmetry.SymbolicWedderburn.action(a::Action, el::DirectSum, mono::CondensedMatterSOS.SpinMonomial)
    isempty(mono.variables) && return 1 * mono
    sign = 1
    vars = map(values(mono.variables)) do var
        rel_id = var.id - a.σ[1][1].id
        rel_index = var.index + 1
        @assert a.σ[rel_index][rel_id + 1] == var
        id = ((rel_id + el.h.id) % el.h.n) + a.σ[1][1].id
        index = (rel_index^el.k.p) - 1
        new_var = CondensedMatterSOS.SpinVariable(id, index)
        if el.k.k.id != 0 && el.k.k.id != index + 1
            sign *= -1
        end
        return new_var
    end
    return sign * CondensedMatterSOS.SpinMonomial(vars)
end
