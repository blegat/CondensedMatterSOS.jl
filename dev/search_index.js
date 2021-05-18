var documenterSearchIndex = {"docs":
[{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"EditURL = \"https://github.com/blegat/CondensedMatterSOS.jl/blob/master/examples/Ising.jl\"","category":"page"},{"location":"generated/Ising/#Ising-model","page":"Ising","title":"Ising model","text":"","category":"section"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"(Image: ) (Image: )","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"We study the Hamiltonian of the transverse-field Ising model with periodic boundary conditions and transverse field g set to 1.","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"using CondensedMatterSOS\n@spin σ[1:2]\nising_hamiltonian(σ, 1, true)","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"Let's pick a solver from this list.","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"using CSDP\nsolver = optimizer_with_attributes(\n    () -> MOIU.CachingOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()), CSDP.Optimizer()),\n    MOI.Silent() => true\n)","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"We can compute a lower bound -2√2 to the ground state energy as follow:","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"function hamiltonian_energy(N, maxdegree, solver; kws...)\n    @spin σ[1:N]\n    H = ising_hamiltonian(σ, 1, true)\n    energy(H, maxdegree, solver; kws...)\nend\nbound, gram, ν = hamiltonian_energy(2, 2, solver, sparsity = NoSparsity())\nbound","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"We can see that the moment matrix uses all monomials:","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"ν.basis.monomials","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"Using term sparsity with block/cluster completion, we get the same bound:","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"bound, gram, ν = hamiltonian_energy(2, 2, solver)\nbound","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"But with a smaller basis:","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"[M.basis.monomials for M in ν.sub_moment_matrices]","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"Using term sparsity with chordal completion, we get a smaller bound:","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"bound, gram, ν = hamiltonian_energy(2, 2, solver, sparsity = MonomialSparsity(ChordalCompletion()))\nbound","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"But with an even smaller basis:","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"[M.basis.monomials for M in ν.sub_moment_matrices]","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"","category":"page"},{"location":"generated/Ising/","page":"Ising","title":"Ising","text":"This page was generated using Literate.jl.","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"EditURL = \"https://github.com/blegat/CondensedMatterSOS.jl/blob/master/examples/Benchmark_1.jl\"","category":"page"},{"location":"generated/Benchmark_1/#Benchmark-1","page":"Benchmark_1","title":"Benchmark 1","text":"","category":"section"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"(Image: ) (Image: )","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"We study the Hamiltonian of the Heisenberg model with periodic boundary conditions.","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"using CondensedMatterSOS\n@spin σ[1:3]\nheisenberg_hamiltonian(σ, true)\n\n# Let's pick a solver from [this list](https://jump.dev/JuMP.jl/dev/installation/#Getting-Solvers).\n\nusing CSDP\nsolver = optimizer_with_attributes(\n    () -> MOIU.CachingOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()), CSDP.Optimizer()),\n    MOI.Silent() => false,\n)","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"We can compute a lower bound -2√2 to the ground state energy as follow:","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"include(joinpath(dirname(dirname(pathof(CondensedMatterSOS))), \"examples\", \"symmetry.jl\"))\nfunction hamiltonian_energy(N, maxdegree, solver; symmetry=true, consecutive=false, kws...)\n    @spin σ[1:N]\n    function action(mono::CondensedMatterSOS.SpinMonomial, el::DirectSum)\n        isempty(mono.variables) && return 1 * mono\n        sign = 1\n        vars = map(values(mono.variables)) do var\n            rel_id = var.id - σ[1][1].id\n            rel_index = var.index + 1\n            @assert σ[rel_index][rel_id + 1] == var\n            id = ((rel_id + el.c.id) % el.c.n) + σ[1][1].id\n            index = (rel_index^el.k.p) - 1\n            new_var = CondensedMatterSOS.SpinVariable(id, index)\n            if el.k.k.id != 0 && el.k.k.id != index + 1\n                sign *= -1\n            end\n            return new_var\n        end\n        return sign * CondensedMatterSOS.SpinMonomial(vars)\n    end\n    function action(term::CondensedMatterSOS.SpinTerm, el::DirectSum)\n        return MP.coefficient(term) * action(MP.monomial(term), el)\n    end\n    function action(poly::CondensedMatterSOS.SpinPolynomial, el::DirectSum)\n        return MP.polynomial([action(term, el) for term in MP.terms(poly)])\n    end\n    H = heisenberg_hamiltonian(σ, true)\n    G = Lattice1Group(N)\n    cone = NonnegPolyInnerCone{SumOfSquares.COI.HermitianPositiveSemidefiniteConeTriangle}()\n    @assert iseven(maxdegree)\n    cert = SumOfSquares.Certificate.FixedBasis(\n        cone,\n        MonomialBasis(MP.monomials(σ[1], 0:div(maxdegree, 2), consecutive=consecutive)),\n    )\n    display(cert.basis.monomials)\n    certificate = SymmetricIdeal(\n        cert,\n        G,\n        action,\n    )\n    if symmetry\n        energy(H, maxdegree, solver; certificate = certificate, kws...)\n    else\n        energy(H, maxdegree, solver; kws...)\n    end\nend\nbound, gram, ν = hamiltonian_energy(\n    2,\n    2,\n    solver,\n    symmetry = false,\n    sparsity = NoSparsity(),\n)\nbound","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"We can see that the moment matrix uses all monomials:","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"ν.basis.monomials","category":"page"},{"location":"generated/Benchmark_1/#Symmetry-reduction","page":"Benchmark_1","title":"Symmetry reduction","text":"","category":"section"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"We can reduce the computation using symmetry reduction as follows.","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"using CondensedMatterSOS\n\nbound, gram, ν = hamiltonian_energy(\n    2,\n    2,\n    solver,\n)\nbound","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"The reduction is obtained by block diagonalizing with a change of polynomial basis to the isotypical basis.","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"display([M.basis.polynomials for M in ν.sub_moment_matrices])\n\n[M.basis.polynomials for M in ν.sub_moment_matrices]","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"Let's try this for 3 sites. First without symmetry.","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"bound, gram, ν = hamiltonian_energy(\n    3,\n    2,\n    solver,\n    symmetry = false,\n)\n@show bound","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"Now with symmetry.","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"bound, gram, ν = hamiltonian_energy(\n    3,\n    2,\n    solver,\n)\n@show bound","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"Let's look at the isotypical basis.","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"display([M.basis.polynomials for M in ν.sub_moment_matrices])","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"Now let's define a function for our common use case.","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"function f(L, d=1, consecutive=false; symmetry=true)\n    @show L\n    println(\"***\")\n    @show d\n    bound, gram, ν = @time hamiltonian_energy(\n        L,\n        2d,\n        solver,\n        consecutive=consecutive,\n        symmetry=symmetry,\n    )\n    @show bound\n    for M in ν.sub_moment_matrices\n        display(round.(M.basis.polynomials, digits=6))\n    end\n    println(\"E/N = \", bound / L)\n    println(\"------------------------------------\")\nend\n\nf(6, 1, true)","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"Now with d = 2.","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"f(6, 2, true)","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"id irep 1 irep 2 irep 3 irep 4\ndegree 2 3 1 3\nmult 2 1 3 2 1\nmult 3 3 6 4 3\nmult 4 6 10 7 6\nmult 5 10 15 11 10\nmult 6 15 21 16 15\nmult 7 21 28 22 21","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"","category":"page"},{"location":"generated/Benchmark_1/","page":"Benchmark_1","title":"Benchmark_1","text":"This page was generated using Literate.jl.","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"EditURL = \"https://github.com/blegat/CondensedMatterSOS.jl/blob/master/examples/Glass.jl\"","category":"page"},{"location":"generated/Glass/#Ising-spin-glass-model","page":"Glass","title":"Ising spin glass model","text":"","category":"section"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"(Image: ) (Image: )","category":"page"},{"location":"generated/Glass/#D","page":"Glass","title":"1D","text":"","category":"section"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"We study the Hamiltonian of the transverse-field Ising spin glass model without periodic boundary conditions and transverse field g set to 0.4.","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"import Random\nusing CondensedMatterSOS\n@spin σ[1:2]\nfunction hamiltonian(σ)\n    Random.seed!(1) # We set the seed so that the results of this example are reproducible.\n    return ising_glass_hamiltonian(σ, 0.4, false)\nend\nhamiltonian(σ)","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"Let's pick a solver from this list.","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"using CSDP\nsolver = optimizer_with_attributes(\n    () -> MOIU.CachingOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()), CSDP.Optimizer()),\n    MOI.Silent() => true\n)","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"We can compute a lower bound -0.853452 to the ground state energy as follow:","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"function hamiltonian_energy(N, maxdegree, solver; kws...)\n    @spin σ[1:N]\n    H = hamiltonian(σ)\n    energy(H, maxdegree, solver; kws...)\nend\nbound, gram, ν = hamiltonian_energy(2, 2, solver, sparsity = NoSparsity())\nbound","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"We can see that the moment matrix uses all monomials:","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"ν.basis.monomials","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"Using term sparsity with block/cluster completion, we get the same bound:","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"bound, gram, ν = hamiltonian_energy(2, 2, solver)\nbound","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"But with a smaller basis:","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"[M.basis.monomials for M in ν.sub_moment_matrices]","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"Using term sparsity with chordal completion, we get a smaller bound:","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"bound, gram, ν = hamiltonian_energy(2, 2, solver, sparsity = MonomialSparsity(ChordalCompletion()))\nbound","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"But with an even smaller basis:","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"[M.basis.monomials for M in ν.sub_moment_matrices]","category":"page"},{"location":"generated/Glass/#D-2","page":"Glass","title":"2D","text":"","category":"section"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"We now study the same model but on a square lattice.","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"@spin σ[1:2, 1:2]\nhamiltonian(σ)","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"We can compute a lower bound -1.990513 to the ground state energy as follow:","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"function hamiltonian_energy(N, M, maxdegree, solver; kws...)\n    @spin σ[1:N, 1:M]\n    H = hamiltonian(σ)\n    energy(H, maxdegree, solver; kws...)\nend\nbound, gram, ν = hamiltonian_energy(2, 2, 2, solver, sparsity = NoSparsity())\nbound","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"We can see that the moment matrix uses all monomials:","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"ν.basis.monomials","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"Using term sparsity with block/cluster completion, we get the same bound:","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"bound, gram, ν = hamiltonian_energy(2, 2, 2, solver)\nbound","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"But with a smaller basis:","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"[M.basis.monomials for M in ν.sub_moment_matrices]","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"Using term sparsity with chordal completion, we get a smaller bound:","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"bound, gram, ν = hamiltonian_energy(2, 2, 2, solver, sparsity = MonomialSparsity(ChordalCompletion()))\nbound","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"But with an even smaller basis:","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"[M.basis.monomials for M in ν.sub_moment_matrices]","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"","category":"page"},{"location":"generated/Glass/","page":"Glass","title":"Glass","text":"This page was generated using Literate.jl.","category":"page"},{"location":"#Condensed-Matter-with-Sum-of-Squares","page":"Introduction","title":"Condensed Matter with Sum-of-Squares","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"@spin\nising_hamiltonian\nising_glass_hamiltonian\nenergy","category":"page"},{"location":"#CondensedMatterSOS.@spin","page":"Introduction","title":"CondensedMatterSOS.@spin","text":"@spin(σ[1:N1, 1:N2, ...], ...)\n\nReturn a tuple of 3-tuples σ = (σx, σy, σz) where the three elements of this tuple are arrays of size (N1, N2, ...). Moreover, the product σ[i][I...] * σ[j][J...]\n\ncommutes if I != J,\nif I == J, it\ncommutes if i == j, moreover it satisfies the identity: σx[I...] * σx[I...] == σy[I...] * σy[I...] == σz[I...] * σz[I...] = 1.\nanticommutes if i != j, moreover it satisifes the identities:\nσx * σy == -σy * σx == im * σz;\nσy * σz == -σz * σy == im * σx;\nσz * σx == -σx * σz == im * σy.\n\nIt also sets the 3-tuple to the local Julia variable with name σ.\n\nThe macro can either be used to create several groups that commute with each other\n\njulia> using CondensedMatterSOS\n\njulia> (σ1x, σ1y, σ1z), (σ2x, σ2y, σ2z) = @spin(σ1, σ2)\n((σ1ˣ, σ1ʸ, σ1ᶻ), (σ2ˣ, σ2ʸ, σ2ᶻ))\n\njulia> σ1x * σ2y\nσ1ˣσ2ʸ\n\njulia> σ1x * σ1y\n(0 + 1im)σ1ᶻ\n\njulia> σ2z * σ2y\n(0 - 1im)σ2ˣ\n\njulia> (σ1x + σ2y) * (σ2x + im * σ1z)\nσ1ˣσ2ˣ + (0 + 1im)σ1ᶻσ2ʸ + σ1ʸ + (0 - 1im)σ2ᶻ\n\nOr also a vector of groups commuting with each other. Note that it returns a 1-tuple containing a 3-tuple hence the needed comma after (σx, σy, σz)\n\njulia> using CondensedMatterSOS\n\njulia> (σx, σy, σz), = @spin(σ[1:2])\n((CondensedMatterSOS.SpinVariable[σˣ₁, σˣ₂], CondensedMatterSOS.SpinVariable[σʸ₁, σʸ₂], CondensedMatterSOS.SpinVariable[σᶻ₁, σᶻ₂]),)\n\njulia> σx[1] * σx[2]\nσˣ₁σˣ₂\n\njulia> σx[2] * σx[1]\nσˣ₁σˣ₂\n\njulia> σx[1] * σx[1]\n(1 + 0im)\n\njulia> σy[1] * σy[1]\n(1 + 0im)\n\njulia> σx[1] * σy[1]\n(0 + 1im)σᶻ₁\n\njulia> σy[1] * σx[1]\n(0 - 1im)σᶻ₁\n\njulia> σz[1] * σx[1]\n(0 + 1im)σʸ₁\n\n\n\n\n\n","category":"macro"},{"location":"#CondensedMatterSOS.ising_hamiltonian","page":"Introduction","title":"CondensedMatterSOS.ising_hamiltonian","text":"ising_hamiltonian(σ, g, periodic::Bool, coef::Function = (I, J) -> 1)\n\nReturn the Hamiltonian of the transverse-field Ising model with periodic boundary conditions if periodic is true and transverse field g. The input σ should be a triple σx, σy, σz of arrays as created with @spin.\n\njulia> using CondensedMatterSOS\n\njulia> @spin(σ[1:4])\n((CondensedMatterSOS.SpinVariable[σˣ₁, σˣ₂, σˣ₃, σˣ₄], CondensedMatterSOS.SpinVariable[σʸ₁, σʸ₂, σʸ₃, σʸ₄], CondensedMatterSOS.SpinVariable[σᶻ₁, σᶻ₂, σᶻ₃, σᶻ₄]),)\n\njulia> ising_hamiltonian(σ, 0.5, false)\nσˣ₁σˣ₂ + σˣ₂σˣ₃ + σˣ₃σˣ₄ + (-0.5 + 0.0im)σᶻ₁ + (-0.5 + 0.0im)σᶻ₂ + (-0.5 + 0.0im)σᶻ₃ + (-0.5 + 0.0im)σᶻ₄\n\n\n\n\n\n","category":"function"},{"location":"#CondensedMatterSOS.ising_glass_hamiltonian","page":"Introduction","title":"CondensedMatterSOS.ising_glass_hamiltonian","text":"ising_glass_hamiltonian(σ, g, periodic::Bool)\n\nReturn the Hamiltonian of the transverse-field Ising spin glass model with periodic boundary conditions if periodic is true and transverse field g. The input σ should be a triple σx, σy, σz of arrays as created with @spin.\n\njulia> using CondensedMatterSOS\n\njulia> @spin(σ[1:3, 1:2])\n((CondensedMatterSOS.SpinVariable[σˣ₁₋₁ σˣ₁₋₂; σˣ₂₋₁ σˣ₂₋₂; σˣ₃₋₁ σˣ₃₋₂], CondensedMatterSOS.SpinVariable[σʸ₁₋₁ σʸ₁₋₂; σʸ₂₋₁ σʸ₂₋₂; σʸ₃₋₁ σʸ₃₋₂], CondensedMatterSOS.SpinVariable[σᶻ₁₋₁ σᶻ₁₋₂; σᶻ₂₋₁ σᶻ₂₋₂; σᶻ₃₋₁ σᶻ₃₋₂]),)\n\njulia> import Random; Random.seed!(0); # We set the seed so that the results of this example are reproducible.\n\njulia> ising_glass_hamiltonian(σ, 1.5, false)\n(0.6791074260357777 + 0.0im)σˣ₁₋₁σˣ₂₋₁ + (0.5866170746331097 + 0.0im)σˣ₁₋₁σˣ₁₋₂ + (0.8284134829000359 + 0.0im)σˣ₂₋₁σˣ₃₋₁ + (0.29733585084941616 + 0.0im)σˣ₂₋₁σˣ₂₋₂ + (0.06494754854834232 + 0.0im)σˣ₃₋₁σˣ₃₋₂ + (-0.3530074003005963 - 0.0im)σˣ₁₋₂σˣ₂₋₂ + (-0.13485387193052173 - 0.0im)σˣ₂₋₂σˣ₃₋₂ + (-1.5 + 0.0im)σᶻ₁₋₁ + (-1.5 + 0.0im)σᶻ₂₋₁ + (-1.5 + 0.0im)σᶻ₃₋₁ + (-1.5 + 0.0im)σᶻ₁₋₂ + (-1.5 + 0.0im)σᶻ₂₋₂ + (-1.5 + 0.0im)σᶻ₃₋₂\n\njulia> ising_glass_hamiltonian((σ[1][1:2, 1:2], σ[2][1:2, 1:2], σ[3][1:2, 1:2]), 1, true)\n(-0.6232277759150394 - 0.0im)σˣ₁₋₁σˣ₂₋₁ + (0.048825806096416735 + 0.0im)σˣ₁₋₁σˣ₁₋₂ + (0.05112780493936531 + 0.0im)σˣ₂₋₁σˣ₂₋₂ + (0.8854230743112911 + 0.0im)σˣ₁₋₂σˣ₂₋₂ + -σᶻ₁₋₁ + -σᶻ₂₋₁ + -σᶻ₁₋₂ + -σᶻ₂₋₂\n\n\n\n\n\n","category":"function"},{"location":"#CondensedMatterSOS.energy","page":"Introduction","title":"CondensedMatterSOS.energy","text":"energy(H, maxdegree, solver;\n    cone=NonnegPolyInnerCone{SumOfSquares.COI.HermitianPositiveSemidefiniteConeTriangle}(),\n    sparsity=MonomialSparsity(),\n    non_sparse=SumOfSquares.Certificate.MaxDegree(cone, MonomialBasis, maxdegree),\n    certificate=sparsity isa NoSparsity ? non_sparse : SumOfSquares.Certificate.SparseIdeal(sparsity, non_sparse),\n    kws...\n)\n\nCompute a lower bound to the Ground State Energe (GSE) of the Hamiltonian H using a Sum-of-Squares certificate with monomials of maximum degree maxdegree and solver as a solver. The cone used is cone which is Hermitian PSD matrices by default. The sparsity is exploited to reduce the certificate. The certificate for each reduced part is non_sparse. The overall certificate is certificate. The rest of the keywords are passed to the @constraint macro of the SumOfSquares package.\n\n\n\n\n\n","category":"function"},{"location":"generated/symmetry/","page":"symmetry","title":"symmetry","text":"EditURL = \"https://github.com/blegat/CondensedMatterSOS.jl/blob/master/examples/symmetry.jl\"","category":"page"},{"location":"generated/symmetry/","page":"symmetry","title":"symmetry","text":"using SumOfSquares\ninclude(joinpath(dirname(dirname(pathof(SumOfSquares))), \"examples\", \"symmetry.jl\"))\ninclude(joinpath(dirname(dirname(pathof(SumOfSquares))), \"examples\", \"scaled_perm.jl\"))\n\nstruct KleinElement <: GroupElement\n    id::Int\nend\nBase.:(==)(a::KleinElement, b::KleinElement) = a.id == b.id\n\nPermutationGroups.order(el::KleinElement) = iszero(el.id) ? 1 : 2\nBase.inv(el::KleinElement) = el\n\nfunction Base.:*(a::KleinElement, b::KleinElement)\n    a.id > b.id && return b * a\n    a.id == 0 && return b\n    a.id == b.id && return one(a)\n    a.id == 2 && return KleinElement(1)\n    b.id == 2 && return KleinElement(3)\n    return KleinElement(2)\nend\n\nBase.conj(a::KleinElement, b::KleinElement) = inv(b) * a * b\nBase.:^(a::KleinElement, b::KleinElement) = conj(a, b)\n\nstruct KleinGroup <: Group end\nBase.one(::Union{KleinGroup, KleinElement}) = KleinElement(0)\nPermutationGroups.gens(::KleinGroup) = [KleinElement(1), KleinElement(2)]\nPermutationGroups.order(::Type{T}, G::KleinGroup) where {T} = convert(T, 4)\nfunction Base.iterate(::KleinGroup, prev::KleinElement=KleinElement(-1))\n    id = prev.id + 1\n    if id > 4\n        return nothing\n    else\n        next = KleinElement(id)\n        return next, next\n    end\nend\n\n#SymbolicWedderburn.conjugacy_classes_orbit(KleinGroup())\n\nfunction perm_klein(k::KleinElement, p::Perm)\n    if k.id == 0\n        return k\n    else\n        return KleinElement(k.id^p)\n    end\nend\n\n\"\"\"\nGroup element `k * p = p * k^inv(p)`.\n\"\"\"\nstruct KleinPermElement <: GroupElement\n    p::Perm{Int}\n    k::KleinElement\nend\nBase.isone(el::KleinPermElement) = isone(el.k) && isone(el.p)\n\nfunction Base.hash(el::KleinPermElement, u::UInt64)\n    return hash(el.k, hash(el.p, u))\nend\nfunction Base.:(==)(a::KleinPermElement, b::KleinPermElement)\n    return a.p == b.p && a.k == b.k\nend","category":"page"},{"location":"generated/symmetry/","page":"symmetry","title":"symmetry","text":"k^(inv(p)) * inv(p) * k * p = k^(inv(p)) * k * inv(p) * p = 1","category":"page"},{"location":"generated/symmetry/","page":"symmetry","title":"symmetry","text":"function Base.inv(el::KleinPermElement)\n    inv_p = inv(el.p)\n    KleinPermElement(inv_p, perm_klein(el.k, inv_p))\nend","category":"page"},{"location":"generated/symmetry/","page":"symmetry","title":"symmetry","text":"p * k * q * k' = p * q * k^p * k' k * p * k' * q = k * k'^p * p * q","category":"page"},{"location":"generated/symmetry/","page":"symmetry","title":"symmetry","text":"function Base.:*(a::KleinPermElement, b::KleinPermElement)\n    return KleinPermElement(a.p * b.p, perm_klein(a.k, b.p) * b.k)\nend\nfunction Base.:^(el::KleinPermElement, k::Integer)\n    return Base.power_by_squaring(el, k)\nend\n\nBase.conj(a::KleinPermElement, b::KleinPermElement) = inv(b) * a * b\nBase.:^(a::KleinPermElement, b::KleinPermElement) = conj(a, b)\n\nfunction PermutationGroups.order(el::KleinPermElement)\n    cur = el\n    i = 1\n    while !isone(cur)\n        i += 1\n        cur *= el\n    end\n    return i\nend\n\nstruct KleinPermGroup <: Group end\nBase.one(::Union{KleinPermGroup, KleinPermElement}) = KleinPermElement(Perm(3), one(KleinGroup()))\nfunction PermutationGroups.gens(::KleinPermGroup)\n    els = [KleinPermElement(Perm(3), KleinElement(1))]\n    for g in gens(PermutationGroups.SymmetricGroup(3))\n        push!(els, KleinPermElement(g, one(KleinGroup())))\n    end\n    return els\nend\nPermutationGroups.order(::Type{T}, G::KleinPermGroup) where {T} = convert(T, 6 * 4)\nconst P3 = collect(PermutationGroups.SymmetricGroup(3))\nconst IT = Iterators.product(1:6, 0:3)\nfunction Base.iterate(::KleinPermGroup, args...)\n    el_st = iterate(IT, args...)\n    if el_st === nothing\n        return nothing\n    else\n        el, st = el_st\n        return KleinPermElement(P3[el[1]], KleinElement(el[2])), st\n    end\nend\n\nstruct CyclicElem <: GroupElement\n    n::Int\n    id::Int\nend\nBase.:(==)(a::CyclicElem, b::CyclicElem) = a.n == b.n && a.id == b.id\nBase.inv(el::CyclicElem) = CyclicElem(el.n, (el.n - el.id) % el.n)\n\nfunction Base.:*(a::CyclicElem, b::CyclicElem)\n    return CyclicElem(a.n, (a.id + b.id) % a.n)\nend\nBase.:^(el::CyclicElem, k::Integer) = CyclicElem(el.n, (el.id * k) % el.n)\n\nBase.conj(a::CyclicElem, b::CyclicElem) = inv(b) * a * b\nBase.:^(a::CyclicElem, b::CyclicElem) = conj(a, b)\n\nfunction PermutationGroups.order(el::CyclicElem)\n    return div(el.n, gcd(el.n, el.id))\nend\n\nstruct CyclicGroup <: Group\n    n::Int\nend\nBase.one(c::Union{CyclicGroup, CyclicElem}) = CyclicElem(c.n, 0)\nPermutationGroups.gens(c::CyclicGroup) = [CyclicElem(c.n, 1)]\nPermutationGroups.order(::Type{T}, c::CyclicGroup) where {T} = convert(T, c.n)\nfunction Base.iterate(c::CyclicGroup, prev::CyclicElem=CyclicElem(c.n, -1))\n    id = prev.id + 1\n    if id >= c.n\n        return nothing\n    else\n        next = CyclicElem(c.n, id)\n        return next, next\n    end\nend\n\n#SymbolicWedderburn.conjugacy_classes_orbit(KleinPermGroup())\n\nstruct DirectSum <: GroupElement\n    c::CyclicElem\n    k::KleinPermElement\nend\nfunction PermutationGroups.order(el::DirectSum)\n    return lcm(PermutationGroups.order(el.c), PermutationGroups.order(el.k))\nend\n\nfunction Base.hash(el::DirectSum, u::UInt64)\n    return hash(el.k, hash(el.c, u))\nend\nfunction Base.:(==)(a::DirectSum, b::DirectSum)\n    return a.c == b.c && a.k == b.k\nend\n\nfunction Base.inv(el::DirectSum)\n    DirectSum(inv(el.c), inv(el.k))\nend","category":"page"},{"location":"generated/symmetry/","page":"symmetry","title":"symmetry","text":"k * p * k' * q = k * k'^p * p * q","category":"page"},{"location":"generated/symmetry/","page":"symmetry","title":"symmetry","text":"function Base.:*(a::DirectSum, b::DirectSum)\n    return DirectSum(a.c * b.c, a.k * b.k)\nend\nBase.:^(a::DirectSum, k::Integer) = DirectSum(a.c^k, a.k^k)\n\nBase.conj(a::DirectSum, b::DirectSum) = inv(b) * a * b\nBase.:^(a::DirectSum, b::DirectSum) = conj(a, b)\n\nstruct Lattice1Group <: Group\n    n::Int\nend\nBase.one(el::DirectSum) = DirectSum(one(el.c), one(el.k))\nBase.one(L::Lattice1Group) = DirectSum(CyclicElem(L.n, 0), one(KleinPermGroup()))\nfunction PermutationGroups.gens(L::Lattice1Group)\n    els = DirectSum[]\n    for g in gens(CyclicGroup(L.n))\n        push!(els, DirectSum(g, one(KleinPermGroup())))\n    end\n    for g in gens(KleinPermGroup())\n        push!(els, DirectSum(one(CyclicGroup(L.n)), g))\n    end\n    return els\nend\nPermutationGroups.order(L::Lattice1Group) = PermutationGroups.order(Int, L)\nfunction PermutationGroups.order(::Type{T}, L::Lattice1Group) where {T}\n    return order(T, CyclicGroup(L.n)) * order(T, KleinPermGroup())\nend\nBase.length(L::Lattice1Group) = PermutationGroups.order(L)\nfunction Base.iterate(L::Lattice1Group)\n    el_p, st_p = iterate(CyclicGroup(L.n))\n    el_k, st_k = iterate(KleinPermGroup())\n    return DirectSum(el_p, el_k), ((el_p, st_p), st_k)\nend\nfunction Base.iterate(L::Lattice1Group, st)\n    el_st_k = iterate(KleinPermGroup(), st[2])\n    if el_st_k === nothing\n        el_st_p = iterate(CyclicGroup(L.n), st[1][2])\n        el_st_p === nothing && return nothing\n        el_k, st_k = iterate(KleinPermGroup())\n        return DirectSum(el_st_p[1], el_k), (el_st_p, st_k)\n    end\n    return DirectSum(st[1][1], el_st_k[1]), (st[1], el_st_k[2])\nend","category":"page"},{"location":"generated/symmetry/","page":"symmetry","title":"symmetry","text":"","category":"page"},{"location":"generated/symmetry/","page":"symmetry","title":"symmetry","text":"This page was generated using Literate.jl.","category":"page"}]
}
