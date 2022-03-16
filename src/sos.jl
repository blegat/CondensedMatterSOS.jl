using SumOfSquares
export optimizer_with_attributes, MOI, MOIU, NoSparsity, MonomialSparsity, ChordalCompletion

export energy

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
    sparsity=Sparsity.Monomial(),
    non_sparse=SumOfSquares.Certificate.MaxDegree(cone, MonomialBasis, maxdegree),
    certificate=sparsity isa Sparsity.NoPattern ? non_sparse : SumOfSquares.Certificate.Sparsity.Ideal(sparsity, non_sparse),
    kws...
)
    model = MOI.instantiate(solver, with_bridge_type=Float64)
    MOI.Bridges.add_bridge(model, SumOfSquares.Bridges.Constraint.SOSPolynomialBridge{Complex{Float64}})
    MOI.Bridges.add_bridge(model, SumOfSquares.Bridges.Constraint.EmptyBridge{Float64})
    MOI.Bridges.add_bridge(model, SumOfSquares.Bridges.Constraint.PositiveSemidefinite2x2Bridge{Float64})
    MOI.Bridges.add_bridge(model, PolyJuMP.ZeroPolynomialBridge{Complex{Float64}})
    MOI.Bridges.add_bridge(model, SumOfSquares.COI.Bridges.Variable.HermitianToSymmetricPSDBridge{Float64})
    MOI.Bridges.add_bridge(model, SumOfSquares.COI.Bridges.Constraint.SplitZeroBridge{Float64})
    γ = MOI.add_variable(model)
    poly = convert(SpinPolynomial{Complex{Float64}}, H) - (1.0 + 0.0im) * γ
    c = SumOfSquares.add_constraint(model, poly, SOSCone(); ideal_certificate=certificate, kws...)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{MOI.VariableIndex}(), γ)
    MOI.optimize!(model)
    if MOI.get(model, MOI.TerminationStatus()) != MOI.OPTIMAL
        @warn("Termination status: $(MOI.get(model, MOI.TerminationStatus())), $(MOI.get(model, MOI.RawStatusString()))")
    end
    gram = MOI.get(model, SumOfSquares.GramMatrixAttribute(), c)
    ν = MOI.get(model, SumOfSquares.MomentMatrixAttribute(), c)
    MOI.get(model, MOI.ObjectiveValue()), gram, ν
end
