import MultivariatePolynomials as MP
using Test
using CondensedMatterSOS
const CMS = CondensedMatterSOS;

@spin sigma[1:4];


A = -im*7*sigmax[2]*sigmaz[3];
B = sigmax[1]*sigmax[2];
C = sigmax[3]*sigmay[4];
D = 3*sigmax[1]*sigmax[2];

@testset "show" begin
    @test sprint(show, A) == "(0 - 7im)sigmaˣ₂sigmaᶻ₃"
    @test sprint(show, MIME"text/print"(), A) == "(0 - 7im)*sigmax[2]*sigmaz[3]"
    @test sprint(show, B) == "sigmaˣ₁sigmaˣ₂"
    @test sprint(show, MIME"text/print"(), B) == "sigmax[1]*sigmax[2]"
    @test sprint(show, C) == "sigmaˣ₃sigmaʸ₄"
    @test sprint(show, MIME"text/print"(), C) == "sigmax[3]*sigmay[4]"
    @test sprint(show, D) == "(3 + 0im)sigmaˣ₁sigmaˣ₂"
    @test sprint(show, MIME"text/print"(), D) == "(3 + 0im)*sigmax[1]*sigmax[2]"
end

# test that in addition to `x == y`, we
# have `isequal(x, y)` and they have equal hash.
# Not that `0.0 == -0.0` but `!isequal(0.0, -0.0)` as they have a different hash.
function test_isequal(x, y)
    @test x == y
    @test isequal(x, y)
    @test hash(x) == hash(y)
end

@testset "equalities" begin
    test_isequal(A, A)
    test_isequal(C, C)
    @test A !== copy(A)
    test_isequal(A, copy(A))
    test_isequal(sigmax[1]*sigmax[2], sigmax[2]*sigmax[1])
    test_isequal(sigmax[3] * B, B * sigmax[3])
end

@testset "single site" begin
    @test sigmax[1] * sigmax[1] == 1
    @test_throws ErrorException("Invalid `index` field for variables, we should have 0 <= `a.index`, `b.index` < 3.") sigmax[1] * CMS.SpinVariable(sigmax[1].id, 3)

    @test sigmax[1]*sigmay[1]==im*sigmaz[1];
    @test sigmay[1]*sigmaz[1]==im*sigmax[1];
    @test sigmaz[1]*sigmax[1]==im*sigmay[1];

    @test sigmax[1]*sigmaz[1]==-im*sigmay[1];
    @test sigmaz[1]*sigmay[1]==-im*sigmax[1];
    @test sigmay[1]*sigmax[1]==-im*sigmaz[1];

    @test sigmax[2]*sigmay[2]==im*sigmaz[2];
    @test sigmay[2]*sigmaz[2]==im*sigmax[2];
    @test sigmaz[2]*sigmax[2]==im*sigmay[2];

    @test sigmax[2]*sigmaz[2]==-im*sigmay[2];
    @test sigmaz[2]*sigmay[2]==-im*sigmax[2];
    @test sigmay[2]*sigmax[2]==-im*sigmaz[2];
end


@testset "unsorted" begin
    @test sigmax[1]*sigmax[2] isa MP.Term{Complex{Int}}
    m = MP.constant_monomial(sigmax[1])
    @test m * m isa MP.Term{Complex{Int}}
    @test collect(MP.variables((sigmax[1]*sigmax[2]))) == [sigmax[1], sigmax[2]]
end

@testset "monomials and terms" begin
    @test sigmax[1]*sigmax[1] == 1
    @test sigmax[2]*sigmax[2] == 1

    @test C*D==3*sigmax[1]*sigmax[2]*sigmax[3]*sigmay[4]
    @test C == sigmax[3]*sigmay[4]
    @test D == 3*sigmax[1]*sigmax[2]

    @test B*D == 3
    @test B == sigmax[1]*sigmax[2]
    @test D == 3*sigmax[1]*sigmax[2]

    @test B*A*C == 7*sigmax[1]*sigmay[3]*sigmay[4]
    @test A == -im*7*sigmax[2]*sigmaz[3]
    @test B == sigmax[1]*sigmax[2]
    @test C == sigmax[3]*sigmay[4]

    @test (3*B)*A*C == 21*sigmax[1]*sigmay[3]*sigmay[4]
    @test A == -im*7*sigmax[2]*sigmaz[3]
    @test B == sigmax[1]*sigmax[2]
    @test C == sigmax[3]*sigmay[4]
end

@testset "Polynomials" begin
    function test_equal(a, b)
        @test iszero(a + (-b))
        @test iszero((-a) + b)
        @test iszero(b + -a)
        @test iszero((-b) + a)
        @test iszero(b - a)
        @test iszero(a - b)
        @test a == b
    end
    test_equal(-sigmax[1] * sigmax[2] - sigmax[2] * sigmax[1], -2sigmax[1] * sigmax[2])
    test_equal(sum(sigmaz), sigmaz[1] + sigmaz[2] + sigmaz[3] + sigmaz[4])
    @test MP.variables(sum(sigmaz)) == sigmaz
    test_equal(sigmax[1] + 1 - sigmax[1], 1)
    test_equal(sigmax[1] + sigmax[2] - sigmax[3], sigmax[1] - sigmax[3] + sigmax[2])
    test_equal(sigmax[1] - sigmax[2] + sigmax[3], sigmax[1] + sigmax[3] - sigmax[2])
end


@testset "isless" begin
    @test   sigmax[2] < sigmax[1]
    @test !(sigmax[1] < sigmax[1])
    @test !(sigmax[1] < sigmay[1])
    @test   sigmay[1] < sigmax[1]
    @test !(sigmax[3] < sigmax[3])
    @test !(sigmax[3] < sigmax[4])
    @test !(sigmax[3] < sigmay[4])
    @test !(sigmay[3] < sigmax[4])
    @test   sigmax[4] < sigmax[3]
    @test   sigmax[4] < sigmay[3]
    @test   sigmay[4] < sigmax[3]
    @test !(sigmax[1] * sigmax[2] < sigmax[2])
    @test !(sigmax[2] * sigmax[4] < sigmax[1])
    @test   sigmax[2] * sigmax[4] < sigmax[1] * sigmax[2]
    @test !(sigmax[1] * sigmax[3] < sigmax[4] * sigmax[1])
    @test   sigmax[1] * sigmax[3] < sigmax[4] * sigmax[1] * sigmax[3]
    @test   sigmay[2] * sigmaz[3] < sigmay[1] * sigmay[3]
    @test   sigmay[2] * sigmaz[3] < sigmay[1] * sigmay[3]
    @test   sigmay[2] * sigmaz[3] < sigmax[1] * sigmay[4]
    @test !(sigmax[1] * sigmax[2] < sigmax[1] * sigmax[2])
    @test !(2 * sigmax[1] < sigmax[2])
    @test sigmax[1] < 2 * sigmax[1]
end

using LinearAlgebra
@testset "matvec with $v" for v in (sigmax[1:2], [1, sigmax[1]], [sigmax[1], sigmay[1]])
    @show v
    exp = v[1]*v[2] + v[2]*v[1] + v[1]^2 + v[2]^2
    @test v' * ones(2, 2) * v == exp
    @test v' * ones(Int, 2, 2) * v == exp
    vv = v * v'
    @test vv isa Matrix{MP.Term{Complex{Int},CMS.SpinMonomial}}
    @test dot(vv, ones(2, 2)) == exp
    @test dot(vv, ones(Int, 2, 2)) == exp
    @test sum(vv .* ones(2, 2)) == exp
    @test sum(vv .* ones(Int, 2, 2)) == exp
end

@testset "monomials" begin
    @test MP.monomial_vector([1, sigmax[1], sigmax[2]]) == [1, sigmax[2], sigmax[1]]
    @test CMS.monomials([sigmax[1], sigmax[2]], 0) == [1]
    for consecutive in [false, true]
        @test MP.monomials([sigmax[1], sigmax[2]], 0, consecutive=consecutive) == [1]
        @test all(MP.monomials([sigmax[1], sigmax[2]], 2, consecutive=consecutive) .== [CMS.SpinVariable(sigmax[1].id, i) * CMS.SpinVariable(sigmax[2].id, j) for i in 0:2 for j in 0:2])
    end
    @test all(MP.monomials([sigmax[1], sigmax[2], sigmax[4]], 2, consecutive=true) .== append!(
        append!(
            [CMS.SpinVariable(sigmax[1].id, i) * CMS.SpinVariable(sigmax[2].id, j) for i in 0:2 for j in 0:2],
            [CMS.SpinVariable(sigmax[2].id, i) * CMS.SpinVariable(sigmax[4].id, j) for i in 0:2 for j in 0:2]
        ),
        [CMS.SpinVariable(sigmax[4].id, i) * CMS.SpinVariable(sigmax[1].id, j) for i in 0:2 for j in 0:2],
    ))
    @test all(CMS.monomials([sigmax[1],sigmax[2]],2) .== [sigmax[1]*sigmax[2]])
    @test all(CMS.monomials([sigmax[1],sigmay[1]],2) .== [im*sigmaz[1],-im*sigmaz[1]])
end
