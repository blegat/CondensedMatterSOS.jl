# # Benchmark 1

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Benchmark_1.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Benchmark_1.ipynb)

# We study the Hamiltonian of the Heisenberg model with periodic boundary conditions.

using Test #src
using CondensedMatterSOS
import MultivariatePolynomials
const MP = MultivariatePolynomials
@spin σ[1:3]
heisenberg_hamiltonian(σ, true)

## Let's pick a solver from [this list](https://jump.dev/JuMP.jl/dev/installation/#Getting-Solvers).

using CSDP
solver = optimizer_with_attributes(
    () -> MOIU.CachingOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()), CSDP.Optimizer()),
    MOI.Silent() => false,
);

# We can compute a lower bound `-2√2` to the ground state energy as follow:

function hamiltonian_energy(N, maxdegree, solver; symmetry=true, consecutive=false, kws...)
    @spin σ[1:N]
    H = heisenberg_hamiltonian(σ, true)
    G = Lattice1Group(N)
    cone = NonnegPolyInnerCone{MOI.HermitianPositiveSemidefiniteConeTriangle}()
    @assert iseven(maxdegree)
    cert = SumOfSquares.Certificate.FixedBasis(
        cone,
        MonomialBasis(MP.monomials(σ[1], 0:div(maxdegree, 2), consecutive=consecutive)),
    )
    certificate = Symmetry.Ideal(
        Symmetry.Pattern(G, Action(σ)),
        cert,
    )
    if symmetry
        energy(H, maxdegree, solver; certificate = certificate, kws...)
    else
        energy(H, maxdegree, solver; kws...)
    end
end
bound, gram, ν = hamiltonian_energy(
    2,
    2,
    solver,
    symmetry = false,
    sparsity = SumOfSquares.Sparsity.NoPattern(),
)
@test bound ≈ -6 rtol=1e-6 #src
bound

# We can see that the moment matrix uses all monomials:

@test length(ν.basis.monomials) == 7 #src
ν.basis.monomials

# # Symmetry reduction
#
# We can reduce the computation using symmetry reduction as follows.

using CondensedMatterSOS

bound, gram, ν = hamiltonian_energy(
    2,
    2,
    solver,
)
@test bound ≈ -6 rtol=1e-6 #src
bound

# The reduction is obtained by block diagonalizing with a change of polynomial
# basis to the isotypical basis.

display([M.basis.polynomials for M in ν.blocks])

@test length(ν.blocks) == 7 #src
[M.basis.polynomials for M in ν.blocks]

# Let's try this for 3 sites. First without symmetry.

bound, gram, ν = hamiltonian_energy(
    3,
    2,
    solver,
    symmetry = false,
)
@show bound
@test bound ≈ -4.5 rtol=1e-6 #src

# Now with symmetry.

bound, gram, ν = hamiltonian_energy(
    3,
    2,
    solver,
)
@show bound
@test bound ≈ -4.5 rtol=1e-6 #src

# Let's look at the isotypical basis.

display([M.basis.polynomials for M in ν.blocks])

# Now let's define a function for our common use case.

function f(L, d=1, consecutive=false; symmetry=true)
    @show L
    println("***")
    @show d
    bound, gram, ν = @time hamiltonian_energy(
        L,
        2d,
        solver,
        consecutive=consecutive,
        symmetry=symmetry,
    )
    @show bound
    for M in ν.blocks
        display(round.(M.basis.polynomials, digits=6))
    end
    println("E/N = ", bound / L)
    println("------------------------------------")
end

f(6, 1, true)

# Now with `d = 2`.

f(6, 2, true)

# | id     | irep 1 | irep 2 | irep 3 | irep 4 |
# |--------|--------|--------|--------|--------|
# | degree | 2      | 3      | 1      | 3      |
# | mult 2 | 1      | 3      | 2      | 1      |
# | mult 3 | 3      | 6      | 4      | 3      |
# | mult 4 | 6      | 10     | 7      | 6      |
# | mult 5 | 10     | 15     | 11     | 10     |
# | mult 6 | 15     | 21     | 16     | 15     |
# | mult 7 | 21     | 28     | 22     | 21     |
