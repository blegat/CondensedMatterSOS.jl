using MultivariatePolynomials
using Test
using CondensedMatterSOS

const CMS = CondensedMatterSOS;

@spin sigma[1:4];


A = -im*7*sigmax[2]*sigmaz[3];
B = sigmax[1]*sigmax[2];
C = sigmax[3]*sigmay[4];
D = 3*sigmax[1]*sigmax[2];

@testset "equalities" begin
    @test A==A
    @test C==C
    @test sigmax[1]*sigmax[2]==sigmax[2]*sigmax[1]
    @test sigmax[3] * B == B * sigmax[3]
end

@testset "single site" begin
    @test sigmax[1]*sigmax[1]==true;

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
    @test typeof(sigmax[1]*sigmax[2]) == CMS.SpinMonomial
    @test collect(variables((sigmax[1]*sigmax[2]))) == [sigmax[1], sigmax[2]]
end

@testset "monomials and terms" begin
    @test sigmax[1]*sigmax[1]
    @test sigmax[2]*sigmax[2]

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
