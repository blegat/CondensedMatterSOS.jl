using Reexport
@reexport using SumOfSquares
export optimizer_with_attributes, MOI, MOIU, NoSparsity, MonomialSparsity, ChordalCompletion

export energy

"""
    energy(H, maxdegree, solver;
        cone=NonnegPolyInnerCone{MOI.HermitianPositiveSemidefiniteConeTriangle}(),
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
    cone=NonnegPolyInnerCone{MOI.HermitianPositiveSemidefiniteConeTriangle}(),
    sparsity=Sparsity.Monomial(),
    non_sparse=SumOfSquares.Certificate.MaxDegree(cone, MonomialBasis, maxdegree),
    certificate=sparsity isa Sparsity.NoPattern ? non_sparse : SumOfSquares.Certificate.Sparsity.Ideal(sparsity, non_sparse),
    kws...
)
    model = Model(solver)
    @variable(model, γ)
    poly = convert(MP.Polynomial{Complex{Float64}}, H) - (1.0 + 0.0im) * γ
    cone = NonnegPolyInnerCone{MOI.HermitianPositiveSemidefiniteConeTriangle}()
    c = @constraint(model, poly in cone, ideal_certificate = certificate, kws...)
    @objective(model, Max, γ)
    optimize!(model)
    if termination_status(model) != MOI.OPTIMAL
        @warn("Termination status: $(termination_status(model)), $(raw_status(model))")
    end
    gram = gram_matrix(c)
    ν = MOI.get(model, SumOfSquares.MomentMatrixAttribute(), c)
    objective_value(model), gram, ν
end
