var documenterSearchIndex = {"docs":
[{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"EditURL = \"../../../examples/Ising.jl\"","category":"page"},{"location":"generated/Ising/#Ising-model","page":"Ising","title":"Ising model","text":"","category":"section"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"(Image: ) (Image: )","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"We study the Hamiltonian of the transverse-field Ising model with periodic boundary conditions and transverse field g set to 1.","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"using CondensedMatterSOS\n@spin σ[1:2]\nising_hamiltonian(σ, 1, true)","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"Let's pick a solver from this list.","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"using CSDP\nsolver = optimizer_with_attributes(\n    () -> MOIU.CachingOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()), CSDP.Optimizer()),\n    MOI.Silent() => true\n)","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"We can compute a lower bound -2√2 to the ground state energy as follow:","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"function hamiltonian_energy(N, maxdegree, solver; kws...)\n    @spin σ[1:N]\n    H = ising_hamiltonian(σ, 1, true)\n    energy(H, maxdegree, solver; kws...)\nend\nbound, gram, ν = hamiltonian_energy(2, 2, solver, sparsity = NoSparsity())\nbound","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"We can see that the moment matrix uses all monomials:","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"ν.basis.monomials","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"Using term sparsity with block/cluster completion, we get the same bound:","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"bound, gram, ν = hamiltonian_energy(2, 2, solver)\nbound","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"But with a smaller basis:","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"[M.basis.monomials for M in ν.blocks]","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"Using term sparsity with chordal completion, we get a smaller bound:","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"bound, gram, ν = hamiltonian_energy(2, 2, solver, sparsity = MonomialSparsity(ChordalCompletion()))\nbound","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"But with an even smaller basis:","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"[M.basis.monomials for M in ν.blocks]","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"This page was generated using Literate.jl.","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"EditURL = \"../../../examples/Benchmark_1.jl\"","category":"page"},{"location":"generated/Benchmark_1/#Benchmark-1","page":"Benchmark_1","title":"Benchmark 1","text":"","category":"section"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"(Image: ) (Image: )","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"We study the Hamiltonian of the Heisenberg model with periodic boundary conditions.","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"using CondensedMatterSOS\nimport MultivariatePolynomials as MP\n@spin σ[1:3]\nMP.monomials(σ[1], 0:2, consecutive=true)\nheisenberg_hamiltonian(σ, true)","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"Let's pick a solver from this list.","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"import Clarabel\nsolver = Clarabel.Optimizer","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"We can compute a lower bound -2√2 to the ground state energy as follow:","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"function hamiltonian_energy(N, maxdegree, solver; symmetry=true, consecutive=false, kws...)\n    @spin σ[1:N]\n    G = Lattice1Group(N)\n    @assert iseven(maxdegree)\n    H = heisenberg_hamiltonian(σ, true)\n    cone = NonnegPolyInnerCone{MOI.HermitianPositiveSemidefiniteConeTriangle}()\n    certificate = SumOfSquares.Certificate.FixedBasis(\n        cone,\n        MonomialBasis(MP.monomials(σ[1], 0:div(maxdegree, 2), consecutive=consecutive)),\n    )\n    if symmetry\n        certificate = Symmetry.Ideal(\n            Symmetry.Pattern(G, Action(σ)),\n            certificate,\n        )\n    end\n    energy(H, maxdegree, solver; certificate, kws...)\nend\nbound, gram, ν = hamiltonian_energy(\n    2,\n    2,\n    solver,\n    symmetry = false,\n    sparsity = SumOfSquares.Sparsity.NoPattern(),\n)\nbound","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"We can see that the moment matrix uses all monomials:","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"ν.basis.monomials","category":"page"},{"location":"generated/Benchmark_1/#Symmetry-reduction","page":"Benchmark_1","title":"Symmetry reduction","text":"","category":"section"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"We can reduce the computation using symmetry reduction as follows.","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"using CondensedMatterSOS\n\nbound, gram, ν = hamiltonian_energy(\n    2,\n    2,\n    solver,\n)\nbound","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"The reduction is obtained by block diagonalizing with a change of polynomial basis to the isotypical basis.","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"[M.basis.polynomials for M in ν.blocks]","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"Let's try this for 3 sites. First without symmetry.","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"bound, gram, ν = hamiltonian_energy(\n    3,\n    2,\n    solver,\n    symmetry = false,\n)\n@show bound","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"Now with symmetry.","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"bound, gram, ν = hamiltonian_energy(\n    3,\n    2,\n    solver,\n)\n@show bound","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"Let's look at the isotypical basis.","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"[M.basis.polynomials for M in ν.blocks]","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"Now let's define a function for our common use case.","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"function f(N, d=1; verbose = 1, kws...)\n    @show N\n    println(\"***\")\n    @show d\n    bound, _, ν = @time hamiltonian_energy(N, 2d, solver; kws...)\n    @show bound\n    block_sizes = map(ν.blocks) do M\n        return length(M.basis.polynomials)\n    end\n    @show block_sizes\n    if verbose >= 2\n        for M in ν.blocks\n            display(round.(M.basis.polynomials, digits=6))\n        end\n    end\n    println(\"E/N = \", bound / N)\n    println(\"------------------------------------\")\n    return bound / N\nend","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"With d = 1, we find a lower bound of -3:","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"lb = f(6, 1, consecutive = true, symmetry = true)","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"Now with d = 2, we find -2:","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"lb = f(6, 2, consecutive = true, symmetry = true)","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"Now with d = 3, we find -1.8685:","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"lb = f(6, 3, consecutive = true, symmetry = true)","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"id irep 1 irep 2 irep 3 irep 4\ndegree 2 3 1 3\nmult 2 1 3 2 1\nmult 3 3 6 4 3\nmult 4 6 10 7 6\nmult 5 10 15 11 10\nmult 6 15 21 16 15\nmult 7 21 28 22 21","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"This page was generated using Literate.jl.","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"EditURL = \"../../../examples/Glass.jl\"","category":"page"},{"location":"generated/Glass/#Ising-spin-glass-model","page":"Glass","title":"Ising spin glass model","text":"","category":"section"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"(Image: ) (Image: )","category":"page"},{"location":"generated/Glass/#D","page":"Glass","title":"1D","text":"","category":"section"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"We study the Hamiltonian of the transverse-field Ising spin glass model without periodic boundary conditions and transverse field g set to 0.4.","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"import Random\nusing CondensedMatterSOS\n@spin σ[1:2]\nfunction hamiltonian(σ)\n    Random.seed!(1) # We set the seed so that the results of this example are reproducible.\n    return ising_glass_hamiltonian(σ, 0.4, false)\nend\nhamiltonian(σ)","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"Let's pick a solver from this list.","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"import Clarabel\nsolver = Clarabel.Optimizer","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"We can compute a lower bound -0.8031 to the ground state energy as follow:","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"function hamiltonian_energy(N, maxdegree, solver; kws...)\n    @spin σ[1:N]\n    H = hamiltonian(σ)\n    energy(H, maxdegree, solver; kws...)\nend\nbound, gram, ν = hamiltonian_energy(2, 2, solver, sparsity = NoSparsity())\nbound","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"We can see that the moment matrix uses all monomials:","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"ν.basis.monomials","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"Using term sparsity with block/cluster completion, we get the same bound:","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"bound, gram, ν = hamiltonian_energy(2, 2, solver)\nbound","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"But with a smaller basis:","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"[M.basis.monomials for M in ν.blocks]","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"Using term sparsity with chordal completion, we get a smaller bound:","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"bound, gram, ν = hamiltonian_energy(2, 2, solver, sparsity = MonomialSparsity(ChordalCompletion()))\nbound","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"But with an even smaller basis:","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"[M.basis.monomials for M in ν.blocks]","category":"page"},{"location":"generated/Glass/#D-2","page":"Glass","title":"2D","text":"","category":"section"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"We now study the same model but on a square lattice.","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"@spin σ[1:2, 1:2]\nhamiltonian(σ)","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"We can compute a lower bound -4.15244 to the ground state energy as follow:","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"function hamiltonian_energy(N, M, maxdegree, solver; kws...)\n    @spin σ[1:N, 1:M]\n    H = hamiltonian(σ)\n    energy(H, maxdegree, solver; kws...)\nend\nbound, gram, ν = hamiltonian_energy(2, 2, 2, solver, sparsity = NoSparsity())\nbound","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"We can see that the moment matrix uses all monomials:","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"ν.basis.monomials","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"Using term sparsity with block/cluster completion, we get the same bound:","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"bound, gram, ν = hamiltonian_energy(2, 2, 2, solver)\nbound","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"But with a smaller basis:","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"[M.basis.monomials for M in ν.blocks]","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"Using term sparsity with chordal completion, we get a smaller bound:","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"bound, gram, ν = hamiltonian_energy(2, 2, 2, solver, sparsity = MonomialSparsity(ChordalCompletion()))\nbound","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"But with an even smaller basis:","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"[M.basis.monomials for M in ν.blocks]","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"This page was generated using Literate.jl.","category":"page"},{"location":"#Condensed-Matter-with-Sum-of-Squares","page":"Introduction","title":"Condensed Matter with Sum-of-Squares","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"@spin\nising_hamiltonian\nising_glass_hamiltonian\nenergy","category":"page"},{"location":"#CondensedMatterSOS.@spin","page":"Introduction","title":"CondensedMatterSOS.@spin","text":"@spin(σ[1:N1, 1:N2, ...], ...)\n\nReturn a tuple of 3-tuples σ = (σx, σy, σz) where the three elements of this tuple are arrays of size (N1, N2, ...). Moreover, the product σ[i][I...] * σ[j][J...]\n\ncommutes if I != J,\nif I == J, it\ncommutes if i == j, moreover it satisfies the identity: σx[I...] * σx[I...] == σy[I...] * σy[I...] == σz[I...] * σz[I...] = 1.\nanticommutes if i != j, moreover it satisifes the identities:\nσx * σy == -σy * σx == im * σz;\nσy * σz == -σz * σy == im * σx;\nσz * σx == -σx * σz == im * σy.\n\nIt also sets the 3-tuple to the local Julia variable with name σ.\n\nThe macro can either be used to create several groups that commute with each other\n\njulia> using CondensedMatterSOS\n\njulia> (σ1x, σ1y, σ1z), (σ2x, σ2y, σ2z) = @spin(σ1, σ2)\n((σ1ˣ, σ1ʸ, σ1ᶻ), (σ2ˣ, σ2ʸ, σ2ᶻ))\n\njulia> σ1x * σ2y\nσ1ˣσ2ʸ\n\njulia> σ1x * σ1y\n(0 + 1im)σ1ᶻ\n\njulia> σ2z * σ2y\n(0 - 1im)σ2ˣ\n\njulia> (σ1x + σ2y) * (σ2x + im * σ1z)\nσ1ˣσ2ˣ + (0 + 1im)σ1ᶻσ2ʸ + σ1ʸ + (0 - 1im)σ2ᶻ\n\nOr also a vector of groups commuting with each other. Note that it returns a 1-tuple containing a 3-tuple hence the needed comma after (σx, σy, σz)\n\njulia> using CondensedMatterSOS\n\njulia> (σx, σy, σz), = @spin(σ[1:2])\n((CondensedMatterSOS.SpinVariable[σˣ₁, σˣ₂], CondensedMatterSOS.SpinVariable[σʸ₁, σʸ₂], CondensedMatterSOS.SpinVariable[σᶻ₁, σᶻ₂]),)\n\njulia> σx[1] * σx[2]\nσˣ₁σˣ₂\n\njulia> σx[2] * σx[1]\nσˣ₁σˣ₂\n\njulia> σx[1] * σx[1]\n(1 + 0im)\n\njulia> σy[1] * σy[1]\n(1 + 0im)\n\njulia> σx[1] * σy[1]\n(0 + 1im)σᶻ₁\n\njulia> σy[1] * σx[1]\n(0 - 1im)σᶻ₁\n\njulia> σz[1] * σx[1]\n(0 + 1im)σʸ₁\n\n\n\n\n\n","category":"macro"},{"location":"#CondensedMatterSOS.ising_hamiltonian","page":"Introduction","title":"CondensedMatterSOS.ising_hamiltonian","text":"ising_hamiltonian(σ, g, periodic::Bool, coef::Function = (I, J) -> 1)\n\nReturn the Hamiltonian of the transverse-field Ising model with periodic boundary conditions if periodic is true and transverse field g. The input σ should be a triple σx, σy, σz of arrays as created with @spin.\n\njulia> using CondensedMatterSOS\n\njulia> @spin(σ[1:4])\n((CondensedMatterSOS.SpinVariable[σˣ₁, σˣ₂, σˣ₃, σˣ₄], CondensedMatterSOS.SpinVariable[σʸ₁, σʸ₂, σʸ₃, σʸ₄], CondensedMatterSOS.SpinVariable[σᶻ₁, σᶻ₂, σᶻ₃, σᶻ₄]),)\n\njulia> ising_hamiltonian(σ, 0.5, false)\nσˣ₁σˣ₂ + σˣ₂σˣ₃ + σˣ₃σˣ₄ + (-0.5 + 0.0im)σᶻ₁ + (-0.5 + 0.0im)σᶻ₂ + (-0.5 + 0.0im)σᶻ₃ + (-0.5 + 0.0im)σᶻ₄\n\n\n\n\n\n","category":"function"},{"location":"#CondensedMatterSOS.ising_glass_hamiltonian","page":"Introduction","title":"CondensedMatterSOS.ising_glass_hamiltonian","text":"ising_glass_hamiltonian(σ, g, periodic::Bool)\n\nReturn the Hamiltonian of the transverse-field Ising spin glass model with periodic boundary conditions if periodic is true and transverse field g. The input σ should be a triple σx, σy, σz of arrays as created with @spin.\n\njulia> using CondensedMatterSOS\n\njulia> @spin(σ[1:3, 1:2])\n((CondensedMatterSOS.SpinVariable[σˣ₁₋₁ σˣ₁₋₂; σˣ₂₋₁ σˣ₂₋₂; σˣ₃₋₁ σˣ₃₋₂], CondensedMatterSOS.SpinVariable[σʸ₁₋₁ σʸ₁₋₂; σʸ₂₋₁ σʸ₂₋₂; σʸ₃₋₁ σʸ₃₋₂], CondensedMatterSOS.SpinVariable[σᶻ₁₋₁ σᶻ₁₋₂; σᶻ₂₋₁ σᶻ₂₋₂; σᶻ₃₋₁ σᶻ₃₋₂]),)\n\njulia> import Random; Random.seed!(0); # We set the seed so that the results of this example are reproducible.\n\njulia> ising_glass_hamiltonian(σ, 1.5, false)\n(0.6791074260357777 + 0.0im)σˣ₁₋₁σˣ₂₋₁ + (0.5866170746331097 + 0.0im)σˣ₁₋₁σˣ₁₋₂ + (0.8284134829000359 + 0.0im)σˣ₂₋₁σˣ₃₋₁ + (0.29733585084941616 + 0.0im)σˣ₂₋₁σˣ₂₋₂ + (0.06494754854834232 + 0.0im)σˣ₃₋₁σˣ₃₋₂ + (-0.3530074003005963 - 0.0im)σˣ₁₋₂σˣ₂₋₂ + (-0.13485387193052173 - 0.0im)σˣ₂₋₂σˣ₃₋₂ + (-1.5 + 0.0im)σᶻ₁₋₁ + (-1.5 + 0.0im)σᶻ₂₋₁ + (-1.5 + 0.0im)σᶻ₃₋₁ + (-1.5 + 0.0im)σᶻ₁₋₂ + (-1.5 + 0.0im)σᶻ₂₋₂ + (-1.5 + 0.0im)σᶻ₃₋₂\n\njulia> ising_glass_hamiltonian((σ[1][1:2, 1:2], σ[2][1:2, 1:2], σ[3][1:2, 1:2]), 1, true)\n(-0.6232277759150394 - 0.0im)σˣ₁₋₁σˣ₂₋₁ + (0.048825806096416735 + 0.0im)σˣ₁₋₁σˣ₁₋₂ + (0.05112780493936531 + 0.0im)σˣ₂₋₁σˣ₂₋₂ + (0.8854230743112911 + 0.0im)σˣ₁₋₂σˣ₂₋₂ + -σᶻ₁₋₁ + -σᶻ₂₋₁ + -σᶻ₁₋₂ + -σᶻ₂₋₂\n\n\n\n\n\n","category":"function"},{"location":"#CondensedMatterSOS.energy","page":"Introduction","title":"CondensedMatterSOS.energy","text":"energy(H, maxdegree, solver;\n    cone=NonnegPolyInnerCone{MOI.HermitianPositiveSemidefiniteConeTriangle}(),\n    sparsity=MonomialSparsity(),\n    non_sparse=SumOfSquares.Certificate.MaxDegree(cone, MonomialBasis, maxdegree),\n    certificate=sparsity isa NoSparsity ? non_sparse : SumOfSquares.Certificate.SparseIdeal(sparsity, non_sparse),\n    kws...\n)\n\nCompute a lower bound to the Ground State Energe (GSE) of the Hamiltonian H using a Sum-of-Squares certificate with monomials of maximum degree maxdegree and solver as a solver. The cone used is cone which is Hermitian PSD matrices by default. The sparsity is exploited to reduce the certificate. The certificate for each reduced part is non_sparse. The overall certificate is certificate. The rest of the keywords are passed to the @constraint macro of the SumOfSquares package.\n\n\n\n\n\n","category":"function"}]
}
