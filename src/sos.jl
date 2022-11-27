using SumOfSquares
export optimizer_with_attributes, MOI, MOIU, NoSparsity, MonomialSparsity, ChordalCompletion

export energy

function all_odd(mono::SpinMonomial)
    return all(0:2) do index
        isodd(count(v -> v.index == index, values(mono.variable)))
    end
end

"""
    energy(H, maxdegree, solver;
        cone=NonnegPolyInnerCone{SumOfSquares.COI.HermitianPositiveSemidefiniteConeTriangle}(),
        sparsity=MonomialSparsity(),
        non_sparse=SumOfSquares.Certificate.MaxDegree(cone, MonomialBasis, maxdegree),
        certificate=sparsity isa NoSparsity ? non_sparse : SumOfSquares.Certificate.SparseIdeal(sparsity, non_sparse),
        kws...
    )

Compute a lower bound to the Ground State Energe (GSE) of the Hamiltonian `H`
using a Sum-of-Squares certificate with monomials of maximum degree `maxdegree`
and `solver` as a solver. The cone used is `cone` which is Hermitian PSD
matrices by default. The `sparsity` is exploited to reduce the certificate.
The certificate for each reduced part is `non_sparse`. The overall certificate
is `certificate`. The rest of the keywords are passed to the `@constraint` macro
of the SumOfSquares package.
"""
function energy(H, maxdegree, solver;
    cone=NonnegPolyInnerCone{SumOfSquares.COI.HermitianPositiveSemidefiniteConeTriangle}(),
    sparsity=MonomialSparsity(),
    non_sparse=SumOfSquares.Certificate.MaxDegree(cone, MonomialBasis, maxdegree),
    certificate=sparsity isa NoSparsity ? non_sparse : SumOfSquares.Certificate.SparseIdeal(sparsity, non_sparse),
    odd_sym=false,
    kws...
)
    model = MOI.instantiate(solver, with_bridge_type=Float64)
    MOI.Bridges.add_bridge(model, SumOfSquares.Bridges.Constraint.SOSPolynomialBridge{Complex{Float64}})
    MOI.Bridges.add_bridge(model, SumOfSquares.Bridges.Constraint.EmptyBridge{Float64})
    MOI.Bridges.add_bridge(model, SumOfSquares.Bridges.Constraint.PositiveSemidefinite2x2Bridge{Float64})
    MOI.Bridges.add_bridge(model, PolyJuMP.ZeroPolynomialBridge{Complex{Float64}})
    MOI.Bridges.add_bridge(model, SumOfSquares.COI.Bridges.Variable.HermitianToSymmetricPSDBridge{Float64})
    MOI.Bridges.add_bridge(model, SumOfSquares.COI.Bridges.Constraint.SplitZeroBridge{Float64})
    γ = MOI.SingleVariable(MOI.add_variable(model))
    poly = convert(SpinPolynomial{Complex{Float64}}, H) - (1.0 + 0.0im) * γ
    if odd_sym
        gram_basis = SumOfSquares.Certificate.get(certificate, SumOfSquares.Certificate.GramBasis(), poly)
        MCT = SumOfSquares.COI.HermitianPositiveSemidefiniteConeTriangle
        g, Q, cQ = SumOfSquares.add_gram_matrix(model, MCT, gram_basis, Complex{Float64})
        q = MA.operate!(-, poly, g)
        for t in MP.terms(q)
            MOI.add_constraint(model, MOIU.operate(vcat, Complex{Float64}, MP.coefficient(t)), MOI.Zeros(1))
        end
    else
        c = SumOfSquares.add_constraint(model, poly, cone; ideal_certificate=certificate, kws...)
    end
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{MOI.SingleVariable}(), γ)
    MOI.optimize!(model)
    if MOI.get(model, MOI.TerminationStatus()) != MOI.OPTIMAL
        @warn("Termination status: $(MOI.get(model, MOI.TerminationStatus())), $(MOI.get(model, MOI.RawStatusString()))")
    end
    if odd_sym

        gram = SumOfSquares.Bridges.Constraint._gram(Q -> MOI.get(model, MOI.VariablePrimal(), Q),
                     Q, gram_basis, Complex{Float64}, MCT)
        ν = if cQ isa Vector{<:MOI.ConstraintIndex}
            SumOfSquares.build_moment_matrix(gram_basis) do i
                MOI.get(model, MOI.ConstraintDual(), cQ[i])
            end
        else
            SOS.build_moment_matrix(
                MOI.get(model, MOI.ConstraintDual(), cQ),
                gram_basis,
            )
        end
    else
        gram = MOI.get(model, SumOfSquares.GramMatrixAttribute(), c)
        ν = MOI.get(model, SumOfSquares.MomentMatrixAttribute(), c)
    end
    MOI.get(model, MOI.ObjectiveValue()), gram, ν
end
