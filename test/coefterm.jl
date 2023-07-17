@testset "flip.jl" begin
    ve = [ (1, 2), (2, 3), (3, 1)]
    G=graph(ve)
    R,x,q,z=polynomial_ring(G,"x","q","z")
    a=[0,0,3];
    l=[0,0,0];
    Ω=[1,2,3];
    o=[1,2,3];
    @testset "preimg test" begin
        L = [1, 2, 3, 4, 5,3]
        xi = 3
    
        @test preimg(L, xi) == 3
    end
    @test partition(0,1)== [[0]]
    @test partition(1,0)== [[0]]
    @test  loopterm(z[1],q[1],1,2)==11//192*q[1]^4*z[1]^4 + 3//4*q[1]^4*z[1]^2 + 3*q[1]^4
    @test feynman_integral_branch_type(x,q,z,G,a,aa=1,g=[1,0,0])==115//3*q[3]^6
    @test feynman_integral_branch_type_order(x,q,z,G,a,Ω,aa=1,g=[1,0,0]) == 115//6*q[3]^6
    @test feynman_integral_degree_order(x,q,z,G,o,3,aa=1,g=[1,0,0])==1//24*q[1]^4*q[2]^2 + 1//24*q[1]^4*q[3]^2 + 1//24*q[1]^2*q[2]^4 + 1//12*q[1]^2*q[2]^2*q[3]^2 + 1//24*q[1]^2*q[3]^4 + 1//24*q[2]^4*q[3]^2 + 1//24*q[2]^2*q[3]^4 + 115//6*q[3]^6
    @test feynman_integral_degree(x,q,z,G,3,aa=1,g=[1,0,0])==115//3*q[1]^6 + 1//4*q[1]^4*q[2]^2 + 1//4*q[1]^4*q[3]^2 + 1//4*q[1]^2*q[2]^4 + 1//2*q[1]^2*q[2]^2*q[3]^2 + 1//4*q[1]^2*q[3]^4 + 115//3*q[2]^6 + 1//4*q[2]^4*q[3]^2 + 1//4*q[2]^2*q[3]^4 + 115//3*q[3]^6
    @test feynman_integral_degree_sum_order(x,q,z,G,o,[1,2],aa=1,g=[1,0,0])==1//8*q[1]^2*q[4]^2
    @test feynman_integral_degree_sum(x,q,z,G,3,aa=1,g=[1,0,0])==115//3*q[1]^6 + 1//4*q[1]^4*q[2]^2 + 1//4*q[1]^4*q[3]^2 + 19//4*q[1]^4 + 1//4*q[1]^2*q[2]^4 + 1//2*q[1]^2*q[2]^2*q[3]^2 + 1//4*q[1]^2*q[2]^2 + 1//4*q[1]^2*q[3]^4 + 1//4*q[1]^2*q[3]^2 + 1//12*q[1]^2 + 115//3*q[2]^6 + 1//4*q[2]^4*q[3]^2 + 19//4*q[2]^4 + 1//4*q[2]^2*q[3]^4 + 1//4*q[2]^2*q[3]^2 + 1//12*q[2]^2 + 115//3*q[3]^6 + 19//4*q[3]^4 + 1//12*q[3]^2

end