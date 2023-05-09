# # Ising spin glass model

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Ising.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Ising.ipynb)

# ## 1D

# We study the Hamiltonian of the transverse-field Ising spin glass model without periodic boundary conditions and transverse field `g` set to 0.4.

using Test #src
import Random
using CondensedMatterSOS
@spin σ[1:2]
function hamiltonian(σ)
    Random.seed!(1) # We set the seed so that the results of this example are reproducible.
    return ising_glass_hamiltonian(σ, 0.4, false)
end
hamiltonian(σ)

# Let's pick a solver from [this list](https://jump.dev/JuMP.jl/dev/installation/#Getting-Solvers).

using CSDP
solver = optimizer_with_attributes(
    () -> MOIU.CachingOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()), CSDP.Optimizer()),
    MOI.Silent() => true
)

# We can compute a lower bound `-0.8031` to the ground state energy as follow:

function hamiltonian_energy(N, maxdegree, solver; kws...)
    @spin σ[1:N]
    H = hamiltonian(σ)
    energy(H, maxdegree, solver; kws...)
end
bound, gram, ν = hamiltonian_energy(2, 2, solver, sparsity = NoSparsity())
@test bound ≈ -0.8031 rtol=1e-3 #src
bound

# We can see that the moment matrix uses all monomials:

@test length(ν.basis.monomials) == 7 #src
ν.basis.monomials

# Using term sparsity with block/cluster completion, we get the same bound:

bound, gram, ν = hamiltonian_energy(2, 2, solver)
@test bound ≈ -0.8031 rtol=1e-3 #src
bound

# But with a smaller basis:

@test length(ν.sub_moment_matrices) == 2 #src
[M.basis.monomials for M in ν.sub_moment_matrices]

# Using term sparsity with chordal completion, we get a smaller bound:

bound, gram, ν = hamiltonian_energy(2, 2, solver, sparsity = MonomialSparsity(ChordalCompletion()))
@test bound ≈ -0.87058 rtol=1e-3 #src
bound

# But with an even smaller basis:

@test length(ν.sub_moment_matrices) == 5 #src
[M.basis.monomials for M in ν.sub_moment_matrices]

# ## 2D

# We now study the same model but on a square lattice.

@spin σ[1:2, 1:2]
hamiltonian(σ)

# We can compute a lower bound `-4.15244` to the ground state energy as follow:

function hamiltonian_energy(N, M, maxdegree, solver; kws...)
    @spin σ[1:N, 1:M]
    H = hamiltonian(σ)
    energy(H, maxdegree, solver; kws...)
end
bound, gram, ν = hamiltonian_energy(2, 2, 2, solver, sparsity = NoSparsity())
@test bound ≈ -4.15244 rtol=1e-3 #src
bound

# We can see that the moment matrix uses all monomials:

@test length(ν.basis.monomials) == 13 #src
ν.basis.monomials

# Using term sparsity with block/cluster completion, we get the same bound:

bound, gram, ν = hamiltonian_energy(2, 2, 2, solver)
@test bound ≈ -4.15244 rtol=1e-3 #src
bound

# But with a smaller basis:

@test length(ν.sub_moment_matrices) == 2 #src
[M.basis.monomials for M in ν.sub_moment_matrices]

# Using term sparsity with chordal completion, we get a smaller bound:

bound, gram, ν = hamiltonian_energy(2, 2, 2, solver, sparsity = MonomialSparsity(ChordalCompletion()))
@test bound ≈ -5.1878 rtol=1e-3 #src
bound

# But with an even smaller basis:

@test length(ν.sub_moment_matrices) == 9 #src
[M.basis.monomials for M in ν.sub_moment_matrices]
