using SumOfSquares
export optimizer_with_attributes, MOI, MOIU

export energy

function energy(H, maxdegree, solver;
    cone=NonnegPolyInnerCone{SumOfSquares.COI.HermitianPositiveSemidefiniteConeTriangle}(),
    sparsity=MonomialSparsity(),
    non_sparse=SumOfSquares.Certificate.MaxDegree(cone, MonomialBasis, maxdegree),
    certificate=sparsity isa NoSparsity ? non_sparse : SumOfSquares.Certificate.SparseIdeal(sparsity, non_sparse),
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
    c = SumOfSquares.add_constraint(model, H - (1.0 + 0.0im) * γ, SOSCone(); ideal_certificate=certificate, kws...)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{MOI.SingleVariable}(), γ)
    MOI.optimize!(model)
    if MOI.get(model, MOI.TerminationStatus()) != MOI.OPTIMAL
        @warn("Termination status: $(MOI.get(model, MOI.TerminationStatus())), $(MOI.get(model, MOI.RawStatusString()))")
    end
    gram = MOI.get(model, SumOfSquares.GramMatrixAttribute(), c)
    ν = MOI.get(model, SumOfSquares.MomentMatrixAttribute(), c)
    MOI.get(model, MOI.ObjectiveValue()), gram, ν
end
