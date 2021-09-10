#===============================================================================

    tests_2d.jl - LR, November 2020

    Some unittests for the HO stuff.

===============================================================================#
using Test
include("../sp_basis_ho1d.jl")

@testset "HO tests" begin
    @testset "1D HO Orbitals" begin
        orb = HOOrbital1D(16)
        @test orb.n == 15 && orb.e == 15.5
        @test orb(0.0) == 0

        orb = HOOrbital1D(1)
        @test orb.n == 0 && orb.e == 0.5
        @test orb(0.0) > 0
    end
end
