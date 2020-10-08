using SumOfSquares
using MultivariatePolynomials
const MP = MultivariatePolynomials
using ComplexOptInterface
const COI = ComplexOptInterface
using CondensedMatterSOS

Base.promote_rule(::Type{T}, ::Type{MOI.SingleVariable}) where {T<:Number} = MOI.ScalarAffineFunction{T}
PolyJuMP.non_constant(a::Vector{<:MOI.AbstractFunction}) = a
function Base.convert(::Type{MOI.ScalarAffineFunction{T}}, f::MOI.ScalarAffineFunction) where T
    return MOI.ScalarAffineFunction(
        [MOI.ScalarAffineTerm(convert(T, t.coefficient), t.variable_index) for t in f.terms],
        convert(T, f.constant)
    )
end
function MOI.Utilities.operate(::typeof(*), ::Type{Complex{Float64}}, a::Complex{Float64}, f::MOI.ScalarAffineFunction{Float64})
    MOI.Utilities.operate(*, Complex{Float64}, a, convert(MOI.ScalarAffineFunction{Complex{Float64}}, f))
end

Base.:*(f::MOI.ScalarAffineFunction{Complex{Float64}}, a::Number) = f * convert(Complex{Float64}, a)
Base.:*(f::MOI.ScalarAffineFunction{Complex{Float64}}, a::Complex{Float64}) = MOI.Utilities.operate(*, Complex{Float64}, a, f)
Base.:+(a::Int, f::MOI.ScalarAffineFunction{Complex{Float64}}) = convert(Complex{Float64}, a) + f
Base.:+(a::Complex{Float64}, f::MOI.ScalarAffineFunction{Complex{Float64}}) = MOI.Utilities.operate(+, Complex{Float64}, a, f)

function SumOfSquares.matrix_cone(::Type{COI.HermitianPositiveSemidefiniteConeTriangle}, d)
    return COI.HermitianPositiveSemidefiniteConeTriangle(d)
end
using MutableArithmetics
const MA = MutableArithmetics
function MOI.Utilities._add_sub_affine_terms(
    op::Union{typeof(+), typeof(-)}, terms::Vector{MOI.ScalarAffineTerm{T}},
    α::T, f::MOI.SingleVariable) where T
    push!(terms, MOI.ScalarAffineTerm(op(α), f.variable))
end
MOI.Utilities._constant(::Type{T}, ::MOI.SingleVariable) where T = zero(T)
#function SumOfSquares.add_gram_matrix(
#        model::MOI.ModelLike,
#        matrix_cone_type::Type{COI.HermitianPositiveSemidefiniteConeTriangle},
#        basis::AbstractPolynomialBasis, T::Type{<:Complex})
#    n = length(basis)
#    Q, cQ = MOI.add_constrained_variables(model, SumOfSquares.matrix_cone(matrix_cone_type, n))
#    N = MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(n))
#    C = [zero(MOI.ScalarAffineFunction{T}) for i in 1:N]
#    k_real = 0
#    k_imag = 0
#    for j in 1:n
#        for i in 1:j
#            k_real += 1
#            MA.mutable_operate!(MA.add_mul, C[k_real], one(T), MOI.SingleVariable(Q[k_real]))
#            if i != j
#                k_imag += 1
#                MA.mutable_operate!(MA.add_mul, C[k_real], one(T) * im, MOI.SingleVariable(Q[N + k_imag]))
#            end
#        end
#    end
#    q = SumOfSquares.build_gram_matrix(C, basis, T)
#    return q, Q, cQ
#end
function Base.getindex(Q::SymMatrix, i::Integer, j::Integer)
    if i <= j
        return Q.Q[MultivariateMoments.trimap(j, i)]
    else
        return conj(Q[j, i])
    end
end

function COI.Bridges.Constraint.operate_coefficients(f, T::Type, func::MOI.ScalarAffineFunction)
    return MOI.ScalarAffineFunction(
        [COI.Bridges.Constraint.operate_coefficient(f, T, term) for term in func.terms],
        f(func.constant)
    )
end
# We assume MOI variables are real
Base.conj(func::MOI.Utilities.TypedLike{T}) where T = COI.Bridges.Constraint.operate_coefficients(conj, T, func)

default_name(vi) = string("x[", vi.value, "]")
function function_string(::Type{JuMP.REPLMode}, v::MOI.VariableIndex, name = default_name)
    var_name = name(v)
    if !isempty(var_name)
        return var_name
    else
        return "noname"
    end
end
function function_string(::Type{JuMP.IJuliaMode}, v::MOI.VariableIndex, name = default_name)
    var_name = name(v)
    if !isempty(var_name)
        # TODO: This is wrong if variable name constains extra "]"
        return replace(replace(var_name, "[" => "_{", count = 1), "]" => "}")
    else
        return "noname"
    end
end
function_string(mode, func::MOI.SingleVariable) = function_string(mode, func.variable)
# Whether something is zero or not for the purposes of printing it
# oneunit is useful e.g. if coef is a Unitful quantity. The second `abs` is import if it is complex.
_is_zero_for_printing(coef) = abs(coef) < 1e-10 * abs(oneunit(coef))
# Whether something is one or not for the purposes of printing it.
_is_one_for_printing(coef) = _is_zero_for_printing(abs(coef) - oneunit(coef))
_is_one_for_printing(coef::Complex) = _is_one_for_printing(real(coef)) && _is_zero_for_printing(imag(coef))
_sign_string(coef) = coef < zero(coef) ? " - " : " + "
_sign_string(coef) = coef < zero(coef) ? " - " : " + "
function function_string(mode, func::MOI.ScalarAffineFunction, show_constant=true)
    # If the expression is empty, return the constant (or 0)
    if isempty(func.terms)
        return show_constant ? JuMP._string_round(MOI.constant(func)) : "0"
    end

    term_str = Vector{String}(undef, 2length(func.terms))
    elm = 1

    for term in func.terms
        pre = _is_one_for_printing(term.coefficient) ? "" : JuMP._string_round(term.coefficient) * " "

        term_str[2 * elm - 1] = "+"
        term_str[2 * elm] = string(pre, function_string(mode, term.variable_index))
        elm += 1
    end

    if elm == 1
        # Will happen with cancellation of all terms
        # We should just return the constant, if its desired
        return show_constant ? JuMP._string_round(a.constant) : "0"
    else
        # Correction for very first term - don't want a " + "/" - "
        term_str[1] = (term_str[1] == " - ") ? "-" : ""
        ret = join(term_str[1 : 2 * (elm - 1)])
        if !_is_zero_for_printing(MOI.constant(func)) && show_constant
            ret = string(ret, "+",
                         JuMP._string_round(MOI.constant(constant)))
        end
        return ret
    end
end
function Base.show(io::IO, f::MOI.AbstractScalarFunction)
    print(io, function_string(JuMP.REPLMode, f))
end
function Base.show(io::IO, ::MIME"text/latex", f::MOI.AbstractScalarFunction)
    print(io, JuMP._wrap_in_math_mode(function_string(JuMP.IJuliaMode, f)))
end

MOI.Bridges.is_bridged(b::MOI.Bridges.AbstractBridgeOptimizer, vi::MOI.VariableIndex) = vi.value < 0 && MOI.Bridges.Variable.has_bridges(MOI.Bridges.Variable.bridges(b))
