export ising_hamiltonian, ising_glass_hamiltonian

function _wrap(I::CartesianIndex, sz::Tuple, periodic::Bool)
    if periodic
        return CartesianIndex(mod1.(I.I, sz))
    else
        return I
    end
end
function interactions(σ::Array{T,N}, periodic::Bool, coef::Function) where {T,N}
    # Δ = [CartesianIndex(1, 0, ..., 0), ..., CartesianIndex(0, ..., 0, 1)]
    Δs = [CartesianIndex(ntuple(j -> Int(i == j), Val(N))) for i in 1:N]
    term(I, J) = coef(I, J) * σ[I] * σ[J]
    return sum([
        term(I, _wrap(I + Δ, size(σ), periodic))
        for I in CartesianIndices(σ), Δ in Δs
        if _wrap(I + Δ, size(σ), periodic) in CartesianIndices(σ)
    ])
end

"""
    ising_hamiltonian(σ, g, periodic::Bool, coef::Function = (I, J) -> 1)

Return the Hamiltonian of the transverse-field Ising model with periodic boundary conditions if `periodic` is `true` and transverse field `g`.
The input `σ` should be a triple `σx, σy, σz` of arrays as created with [`@spin`](@ref).

```jldoctest
julia> using CondensedMatterSOS

julia> @spin(σ[1:4])
((CondensedMatterSOS.SpinVariable[σˣ₁, σˣ₂, σˣ₃, σˣ₄], CondensedMatterSOS.SpinVariable[σʸ₁, σʸ₂, σʸ₃, σʸ₄], CondensedMatterSOS.SpinVariable[σᶻ₁, σᶻ₂, σᶻ₃, σᶻ₄]),)

julia> ising_hamiltonian(σ, 0.5, false)
σˣ₁σˣ₂ + σˣ₂σˣ₃ + σˣ₃σˣ₄ + (-0.5 + 0.0im)σᶻ₁ + (-0.5 + 0.0im)σᶻ₂ + (-0.5 + 0.0im)σᶻ₃ + (-0.5 + 0.0im)σᶻ₄
```
"""
function ising_hamiltonian(σ, g, periodic::Bool, coef::Function = (I, J) -> 1)
    σx, σy, σz = σ
    return interactions(σx, periodic, coef) - g * sum(σz)
end

"""
    ising_glass_hamiltonian(σ, g, periodic::Bool)

Return the Hamiltonian of the transverse-field Ising spin glass model with periodic boundary conditions if `periodic` is `true` and transverse field `g`.
The input `σ` should be a triple `σx, σy, σz` of arrays as created with [`@spin`](@ref).

```jldoctest
julia> using CondensedMatterSOS

julia> @spin(σ[1:3, 1:2])
((CondensedMatterSOS.SpinVariable[σˣ₁₋₁ σˣ₁₋₂; σˣ₂₋₁ σˣ₂₋₂; σˣ₃₋₁ σˣ₃₋₂], CondensedMatterSOS.SpinVariable[σʸ₁₋₁ σʸ₁₋₂; σʸ₂₋₁ σʸ₂₋₂; σʸ₃₋₁ σʸ₃₋₂], CondensedMatterSOS.SpinVariable[σᶻ₁₋₁ σᶻ₁₋₂; σᶻ₂₋₁ σᶻ₂₋₂; σᶻ₃₋₁ σᶻ₃₋₂]),)

julia> import Random; Random.seed!(0); # We set the seed so that the results of this example are reproducible.

julia> ising_glass_hamiltonian(σ, 1.5, false)
(0.6791074260357777 + 0.0im)σˣ₁₋₁σˣ₂₋₁ + (0.5866170746331097 + 0.0im)σˣ₁₋₁σˣ₁₋₂ + (0.8284134829000359 + 0.0im)σˣ₂₋₁σˣ₃₋₁ + (0.29733585084941616 + 0.0im)σˣ₂₋₁σˣ₂₋₂ + (0.06494754854834232 + 0.0im)σˣ₃₋₁σˣ₃₋₂ + (-0.3530074003005963 - 0.0im)σˣ₁₋₂σˣ₂₋₂ + (-0.13485387193052173 - 0.0im)σˣ₂₋₂σˣ₃₋₂ + (-1.5 + 0.0im)σᶻ₁₋₁ + (-1.5 + 0.0im)σᶻ₂₋₁ + (-1.5 + 0.0im)σᶻ₃₋₁ + (-1.5 + 0.0im)σᶻ₁₋₂ + (-1.5 + 0.0im)σᶻ₂₋₂ + (-1.5 + 0.0im)σᶻ₃₋₂

julia> ising_glass_hamiltonian((σ[1][1:2, 1:2], σ[2][1:2, 1:2], σ[3][1:2, 1:2]), 1, true)
(-0.6232277759150394 - 0.0im)σˣ₁₋₁σˣ₂₋₁ + (0.048825806096416735 + 0.0im)σˣ₁₋₁σˣ₁₋₂ + (0.05112780493936531 + 0.0im)σˣ₂₋₁σˣ₂₋₂ + (0.8854230743112911 + 0.0im)σˣ₁₋₂σˣ₂₋₂ + -σᶻ₁₋₁ + -σᶻ₂₋₁ + -σᶻ₁₋₂ + -σᶻ₂₋₂
```
"""
function ising_glass_hamiltonian(σ, g, periodic::Bool)
    return ising_hamiltonian(σ, g, periodic, (I, J) -> randn())
end
