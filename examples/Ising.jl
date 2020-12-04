# # Ising model

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Ising.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Ising.ipynb)

using Test #src
using CondensedMatterSOS
@spin σ[1:2]
function hamiltonian(σ)
    σx, σy, σz = σ
    N = length(σx)
    return -sum(σx[n] * σx[n+1] for n in 1:(N-1)) - σx[N] * σx[1] - sum(σz)
end
hamiltonian(σ)

function hamiltonian_energy(N, maxdegree, solver; kws...)
    @spin σ[1:N]
    H = 1.0 * hamiltonian(σ)
    energy(H, maxdegree, solver; kws...)
end

using CSDP
solver = optimizer_with_attributes(
    () -> MOIU.CachingOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()), CSDP.Optimizer()),
    MOI.Silent() => true
)
bound, ν = hamiltonian_energy(2, 2, solver)
@test bound ≈ -2 * √2 rtol=1e-6 #src
