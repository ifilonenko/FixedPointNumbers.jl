using Base.Test
using FixedPointNumbers

rows(M::Matrix) = map(x->reshape(getindex(M, x, :), :, size(M)[2]), 1:size(M)[1])
num_iters = 100000
rtol=0.8

function arr_isapprox(arr,sC_input)
    (row,column) = size(arr)
    s_arr = sC_input(arr)
    for (i,v) in enumerate(rows(s_arr))
        for (i_i,v_i) in enumerate(v)
            if !isapprox(v_i,arr[i,i_i],num_iters)
                return false
            end
        end
    end
    return true
end

function mult_isapprox(i,j,x::Type{T1}, y::Type{T2}) where {T1,T2 <:Scaled}
    atol = max(eps(x), eps(y))
    diff = atol + rtol*atol
    abs(mean([i for a in 1:num_iters]) - j) <= float(diff)
end


# Checking the struct definition
@testset "Struct Checking" begin
    @test_throws TypeError Scaled{0.5,7,Int8,ExactAndRandomized}
    @test_throws TypeError Scaled{0.5,7,ExactAndRandomized,Int8}
    @test_throws TypeError Scaled{Int8,7,0.5,0}
    @test Scaled{Int8,7,0.5,ExactAndRandomized} <: Scaled
    @test Scaled{Int8,7,0.5,ExactAndNearestNeighbor} <: Scaled
    @test Scaled{Int8,7,0.5,SatAndRandomized} <: Scaled
    @test Scaled{Int8,7,0.5,SatAndNearestNeighbor} <: Scaled
end

# Checking the Constructor Options
@testset "Constructor Checking" begin
    sC = Scaled{Int8,16,1./128,ExactAndRandomized}
    @test isapprox(sC(0.05),0.05,num_iters) # Float check
    sC2 = Scaled{Int16,15,1./12800,ExactAndRandomized}
    @test isapprox(sC2(2),2,num_iters) # Int Check
    @test !isapprox(sC2(2),1,num_iters)
    @test !isapprox(sC2(2),3,num_iters)
    sC3 = Scaled{Int32,31,1./1280000,ExactAndRandomized} # Checking scale
    @test !isapprox(sC3(5),4.99999,num_iters)
    @test isapprox(sC3(5),4.9999999,num_iters)
    @test isapprox(sC3(0.50),0.500000001,num_iters)
    @test !isapprox(sC3(0.50),0.501,num_iters)
    @test arr_isapprox(rand((3,2)),sC3)
    @test arr_isapprox(rand((9,4)),sC3)
    @test_throws InexactError sC(5)
    sC4 = Scaled{Int8,7,1.,SatAndRandomized}
    @test isapprox(sC4(127),127.0,num_iters)
    @test isapprox(sC4(-128),-128.0,num_iters)
    @test isapprox(sC4(127888888),127.0,num_iters)
    @test isapprox(sC4(-12800000),-128.0,num_iters)
    sC5 = Scaled{Int8,7,1.,ExactAndRandomized}
    @test isapprox(sC5(127),127.0,num_iters)
    @test isapprox(sC5(-128),-128.0,num_iters)
    @test_throws InexactError sC5(128)
    @test_throws InexactError sC5(-129)
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
    sC = Scaled{Int16,15,1./12800,ExactAndRandomized}
    @test isapprox(-sC(0.005),-0.005,num_iters)
    @test isapprox(sC(-0.005),-0.005,num_iters)
    @test isapprox(abs(sC(-0.005)),0.005,num_iters)
    sC2 = Scaled{Int32,31,1./12800,ExactAndRandomized}
    sC3 = Scaled{Int16,15,1./12800,ExactAndRandomized}
    for _ in 1:100
        arr = rand(2)
        @test isapprox(sC(arr[1])+sC(arr[2]),sum(arr),num_iters,rtol)
        @test isapprox(sC2(arr[1])+sC3(arr[2]),sum(arr),num_iters,rtol)
        @test isapprox(sC(arr[1])-sC(arr[2]),arr[1]-arr[2],num_iters,rtol)
        @test isapprox(sC2(arr[1])-sC3(arr[2]),arr[1]-arr[2],num_iters,rtol)
        @test mult_isapprox(sC(arr[1])*sC(arr[2]),arr[1]*arr[2],sC,sC)
        @test mult_isapprox(sC2(arr[1])*sC3(arr[2]),arr[1]*arr[2],sC,sC)
    end
end

@testset "Checking add logic" begin
    sC = Scaled{Int32,31,1./12800,ExactAndRandomized}
    sC2 = Scaled{Int16,15,1./12800,ExactAndRandomized}
    arr = rand(2)
    result = sC(arr[1])+sC2(arr[2])
    @test isa(result,Scaled)
    r_type = typeof(result)
    @test get_T(r_type) == Int64
    @test get_f(r_type) == 31+2
    @test get_s(r_type) == 1./12800
    @test get_r(r_type) == ExactAndRandomized
end

@testset "Checking saturated add logic" begin
    sC = Scaled{Int8,7,1./127,SatAndRandomized}
    result = sC(0.5)+sC(2000.5)
    @test isa(result,Scaled)
    r_type = typeof(result)
    @test get_T(r_type) == Int64
    @test get_f(r_type) == 31+2
    @test get_s(r_type) == 1./12800
    @test get_r(r_type) == SatAndRandomized
end

@testset "Checking subtract logic" begin
    sC = Scaled{Int32,31,1./12800,ExactAndRandomized}
    sC2 = Scaled{Int16,15,1./12800,ExactAndRandomized}
    arr = rand(2)
    result = sC(arr[1])-sC2(arr[2])
    @test isa(result,Scaled)
    r_type = typeof(result)
    @test get_T(r_type) == Int64
    @test get_f(r_type) == 31+2
    @test get_s(r_type) == 1./12800
    @test get_r(r_type) == ExactAndRandomized
end


@testset "Checking saturated add logic" begin
    sC = Scaled{Int8,7,1./127,SatAndRandomized}
    result = sC(0.5)-sC(2000.5)
    @test isa(result,Scaled)
    r_type = typeof(result)
    @test get_T(r_type) == Int64
    @test get_f(r_type) == 31+2
    @test get_s(r_type) == 1./12800
    @test get_r(r_type) == SatAndRandomized
end
