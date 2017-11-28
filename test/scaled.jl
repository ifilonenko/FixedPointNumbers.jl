using Base.Test
using FixedPointNumbers

rows(M::Matrix) = map(x->reshape(getindex(M, x, :), :, size(M)[2]), 1:size(M)[1])
function arr_isapprox(arr,sC_input)
    (row,column) = size(arr)
    s_arr = sC_input(arr)
    for (i,v) in enumerate(rows(s_arr))
        for (i_i,v_i) in enumerate(v)
            if !isapprox(v_i,arr[i,i_i])
                return false
            end
        end
    end
    return true
end

# Checking the struct definition
@testset "Struct Checking" begin
    @test_throws TypeError Scaled{0.5,7,Int8,Randomized}
    @test_throws TypeError Scaled{0.5,7,Randomized,Int8}
    @test_throws TypeError Scaled{Int8,7,0.5,0}
    @test Scaled{Int8,7,0.5,Randomized} <: Scaled
    @test Scaled{Int8,7,0.5,NearestNeighbor} <: Scaled
end

# Checking the Constructor Options
@testset "Constructor Checking" begin
    sC = Scaled{Int8,16,1./128,Randomized}
    @test isapprox(sC(0.05),0.05) # Float check
    sC2 = Scaled{Int16,15,1./12800,Randomized}
    @test isapprox(sC2(2),2) # Int Check
    @test !isapprox(sC2(2),1)
    @test !isapprox(sC2(2),3)
    sC3 = Scaled{Int32,31,1./1280000,Randomized} # Checking scale
    @test !isapprox(sC3(5),4.99999)
    @test isapprox(sC3(5),4.9999999)
    @test isapprox(sC3(0.50),0.500000001)
    @test !isapprox(sC3(0.50),0.501)
    @test arr_isapprox(rand((3,2)),sC3)
    @test arr_isapprox(rand((9,4)),sC3)
end

@testset "Datatype Comparisons" begin
    @test eq(Int32,Int32)
    @test !eq(Int16,Int32)
    @test le(Int32,Int64)
    @test !le(Int16,Int8)
    @test leq(Int16,Int16)
    @test leq(Int8,Int16)
    @test !leq(Int16,Int8)
end

@testset "Upgrading Int size with bit addition" begin
    @test up(Int8,1,6)==Int8
    @test up(Int8,1,7)==Int16
    @test up(Int16,8,7)==Int16
    @test up(Int16,8,16)==Int32
end

@testset "Testing basic operations" begin
    sC = Scaled{Int32,31,1./1280000,Randomized}
    @test isapprox(-sC(0.005),-0.005)
    @test isapprox(sC(-0.005),-0.005)
    @test isapprox(abs(sC(-0.005)),0.005)
    sC2 = Scaled{Int32,30,1./1280000,Randomized}
    arr = rand(2)
    @test isapprox(sC(arr[1])+sC2(arr[2]),sum(arr))
    @test isapprox(sC(arr[1])+sC2(arr[2]),sum(arr))
end
