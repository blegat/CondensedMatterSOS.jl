# # Benchmark 1

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Benchmark_1.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Benchmark_1.ipynb)

# We study the Hamiltonian of the Heisenberg model with periodic boundary conditions.

using Test #src
using CondensedMatterSOS
import MultivariatePolynomials as MP
@spin σ[1:3]
MP.monomials(σ[1], 0:2, consecutive=true)
heisenberg_hamiltonian(σ, true)

# Let's pick a solver from [this list](https://jump.dev/JuMP.jl/dev/installation/#Getting-Solvers).

import Clarabel
solver = Clarabel.Optimizer

# We can compute a lower bound `-2√2` to the ground state energy as follow:

function hamiltonian_energy(N, maxdegree, solver; symmetry=true, consecutive=false, kws...)
    @spin σ[1:N]
    G = Lattice1Group(N)
    @assert iseven(maxdegree)
    H = heisenberg_hamiltonian(σ, true)
    cone = NonnegPolyInnerCone{MOI.HermitianPositiveSemidefiniteConeTriangle}()
    certificate = SumOfSquares.Certificate.FixedBasis(
        cone,
        MonomialBasis(MP.monomials(σ[1], 0:div(maxdegree, 2), consecutive=consecutive)),
    )
    if symmetry
        certificate = Symmetry.Ideal(
            Symmetry.Pattern(G, Action(σ)),
            certificate,
        )
    end
    energy(H, maxdegree, solver; certificate, kws...)
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

@test length(ν.blocks) == 10 #src
[M.basis.polynomials for M in ν.blocks]

# Now let's define a function for our common use case.

function f(N, d=1; verbose = 1, kws...)
    @show N
    println("***")
    @show d
    bound, _, ν = @time hamiltonian_energy(N, 2d, solver; kws...)
    @show bound
    block_sizes = map(ν.blocks) do M
        return length(M.basis.polynomials)
    end
    @show block_sizes
    if verbose >= 2
        for M in ν.blocks
            display(round.(M.basis.polynomials, digits=6))
        end
    end
    println("E/N = ", bound / N)
    println("------------------------------------")
    return bound / N
end

# With `d = 1`, we find a lower bound of `-3`:

lb = f(6, 1, consecutive = true, symmetry = true)
@test lb ≈ -3 rtol=1e-6 #src

# Now with `d = 2`, we find `-2`:

lb = f(6, 2, consecutive = true, symmetry = true)
@test lb ≈ -2 rtol=1e-6 #src

# Now with `d = 3`, we find `-1.8685`:

lb = f(6, 3, consecutive = true, symmetry = true)
@test lb ≈ -1.8685 rtol=1e-6 rtol = 1e-3 #src

# | id     | irep 1 | irep 2 | irep 3 | irep 4 |
# |--------|--------|--------|--------|--------|
# | degree | 2      | 3      | 1      | 3      |
# | mult 2 | 1      | 3      | 2      | 1      |
# | mult 3 | 3      | 6      | 4      | 3      |
# | mult 4 | 6      | 10     | 7      | 6      |
# | mult 5 | 10     | 15     | 11     | 10     |
# | mult 6 | 15     | 21     | 16     | 15     |
# | mult 7 | 21     | 28     | 22     | 21     |
