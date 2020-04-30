using MultivariatePolynomials
using Test
using CondensedMatterSOS

@spin sigma[1:4];

A = sigmax[1];
B = sigmax[1]*sigmax[2];
C = sigmax[3]*sigmax[4];
D = 3*sigmax[1]*sigmax[2];
E = 4*sigmax[3]*sigmax[4];
@testset "unsorted" begin
    @test A==A
    @test C==C

    @test A = -im*7*sigmax[2]*sigmaz[3];
    @test B = sigmax[1]*sigmax[2];
    @test C = sigmax[3]*sigmay[4];
    @test D = 3*sigmax[1]*sigmax[2];

    @test sigmax[1]*sigmax[1] == true
    @test sigmax[2]*sigmax[2] == true
    @test typeof(sigmax[1]*sigmax[2]) == SpinMonomial
    @test variables((sigmax[1]*sigmax[2])) == [sigmax[1], sigmax[2]]
    @test sigmax[1]*sigmay[1] == im*sigmaz[1]
    @test sigmay[1]*sigmax[1] == -im*sigmaz[1]
    @test C*D==3*sigmax[1]*sigmax[2]*sigmax[3]*sigmay[4]
    @test C == sigmax[3]*sigmay[4]
    @test D == 3*sigmax[1]*sigmax[2]
    @test B*D == 3
    @test B == sigmax[1]*sigmax[2]
    @test D == 3*sigmax[1]*sigmax[2]
    @test B*A*C == -7*sigmax[1]*sigmay[3]*sigmay[4]
    @test A == -im*7*sigmax[2]*sigmaz[3]
    @test B == sigmax[1]*sigmax[2]
    @test C == sigmax[3]*sigmay[4]
    @test (3*B)*A*C == -21*sigmax[1]*sigmay[3]*sigmay[4]
    @test A == -im*7*sigmax[2]*sigmaz[3]
    @test B == sigmax[1]*sigmax[2]
    @test C == sigmax[3]*sigmay[4]
end
