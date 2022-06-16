using SpinTrack
using Test
using BenchmarkTools

@testset "SpinTrack.jl" begin
    # Write your tests here.
end

@testset "Performance" begin
    let p = symmetric_hybrid_ring();
        @benchmark get_solution(u0_long, p)
        @test (1.0 > @elapsed get_solution(u0_long, p))
    end
end
