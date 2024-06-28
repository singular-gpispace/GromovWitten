@testset "function.jl" begin
    ve=[(1, 1), (1, 2), (2, 3), (3, 1)]
    G=FeynmanGraph(ve)
    F=FeynmanIntegral(G)
    R,x,q,z=F.S
    @test partition(0,1)== [[0]]
    @test partition(1,0)== [[0]]
    @test  loopterm(z[1],q[1],1,2)==11//192*q[1]^4*z[1]^4 + 3//4*q[1]^4*z[1]^2 + 3*q[1]^4

    @testset "preimg test" begin
        L = [1, 2, 3, 4, 5,3]
        xi = 3
    
        @test preimg(L, xi) == 3
    end

end