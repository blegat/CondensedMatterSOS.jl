# # Ising model

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Ising.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Ising.ipynb)

# We study the Hamiltonian of the transverse-field Ising model with periodic boundary conditions and transverse field `g` set to 1.

using Test #src
using CondensedMatterSOS
@spin σ[1:2]
ising_hamiltonian(σ, 1, true)

# Let's pick a solver from [this list](https://jump.dev/JuMP.jl/dev/installation/#Getting-Solvers).

using CSDP
solver = optimizer_with_attributes(
    () -> MOIU.CachingOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()), CSDP.Optimizer()),
    MOI.Silent() => true
)

# We can compute a lower bound `-2√2` to the ground state energy as follow:

function hamiltonian_energy(N, maxdegree, solver; kws...)
    @spin σ[1:N]
    H = ising_hamiltonian(σ, 1, true)
    energy(H, maxdegree, solver; kws...)
end
bound, gram, ν = hamiltonian_energy(2, 2, solver, sparsity = NoSparsity())
@test bound ≈ -2 * √2 rtol=1e-6 #src
bound

# We can see that the moment matrix uses all monomials:

@test length(ν.basis.monomials) == 7 #src
ν.basis.monomials

# Using term sparsity with block/cluster completion, we get the same bound:

bound, gram, ν = hamiltonian_energy(2, 2, solver)
@test bound ≈ -2 * √2 rtol=1e-6 #src
bound

# But with a smaller basis:

@test length(ν.sub_moment_matrices) == 2 #src
[M.basis.monomials for M in ν.sub_moment_matrices]

# Using term sparsity with chordal completion, we get a smaller bound:

bound, gram, ν = hamiltonian_energy(2, 2, solver, sparsity = MonomialSparsity(ChordalCompletion()))
@test bound ≈ -4 rtol=1e-6 #src
bound

# But with an even smaller basis:

@test length(ν.sub_moment_matrices) == 5 #src
[M.basis.monomials for M in ν.sub_moment_matrices]
