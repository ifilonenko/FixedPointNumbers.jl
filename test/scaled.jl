using Base.Test
using FixedPointNumbers

# Checking the struct definition
@testset "Struct Checking" begin
    @test_throws TypeError Scaled{0.5,Int8,Randomized}
    @test_throws TypeError Scaled{0.5,Randomized,Int8}
    @test_throws TypeError Scaled{Int8,0.5,0}
    @test Scaled{Int8,0.5,Randomized} <: Scaled
    @test Scaled{Int8,0.5,NearestNeighbor} <: Scaled
end
