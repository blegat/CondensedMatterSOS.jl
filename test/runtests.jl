using MultivariatePolynomials
using Test
using CondensedMatterSOS

const CMS = CondensedMatterSOS;

@spin sigma[1:4];


A = -im*7*sigmax[2]*sigmaz[3];
B = sigmax[1]*sigmax[2];
C = sigmax[3]*sigmay[4];
D = 3*sigmax[1]*sigmax[2];

@testset "unsorted" begin
    @test A==A
    @test C==C


    @test sigmax[3] * B == B * sigmax[3]

    @test sigmax[1]*sigmax[1]
    @test sigmax[2]*sigmax[2]
    @test typeof(sigmax[1]*sigmax[2]) == CMS.SpinMonomial
    @test collect(variables((sigmax[1]*sigmax[2]))) == [sigmax[1], sigmax[2]]
    @test sigmax[1]*sigmay[1] == im*sigmaz[1]
    @test sigmay[1]*sigmax[1] == -im*sigmaz[1]
    @test C*D==3*sigmax[1]*sigmax[2]*sigmax[3]*sigmay[4]
    @test C == sigmax[3]*sigmay[4]
    @test D == 3*sigmax[1]*sigmax[2]
    @test B*D == 3
    @test B == sigmax[1]*sigmax[2]
    @test D == 3*sigmax[1]*sigmax[2]
    @test_broken B*A*C == -7*sigmax[1]*sigmay[3]*sigmay[4]
    @test A == -im*7*sigmax[2]*sigmaz[3]
    @test B == sigmax[1]*sigmax[2]
    @test C == sigmax[3]*sigmay[4]
    @test_broken (3*B)*A*C == -21*sigmax[1]*sigmay[3]*sigmay[4]
    @test A == -im*7*sigmax[2]*sigmaz[3]
    @test B == sigmax[1]*sigmax[2]
    @test C == sigmax[3]*sigmay[4]
end
