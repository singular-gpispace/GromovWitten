@testset "flip.jl" begin
    ve = [ (1, 2), (2, 3), (3, 1)]
    G=FeynmanGraph(ve)
    F=FeynmanIntegral(G)
    a=[0,0,3];
    l=[0,0,0];
    Ω=[1,2,3];
    o=[1,2,3];
    g=[1,0,0]
    q=F.S[3]
   # q=Nemo.vars(feynman_integral_degree(F,3,g))
    @test feynman_integral_branch_type(F,a,g)==115//3*q[3]^6
    @test feynman_integral_branch_type_order(F,a,Ω,g) == 115//6*q[3]^6
    @test feynman_integral_degree_order(F,o,3,g)==1//24*q[1]^4*q[2]^2 + 1//24*q[1]^4*q[3]^2 + 1//24*q[1]^2*q[2]^4 + 1//12*q[1]^2*q[2]^2*q[3]^2 + 1//24*q[1]^2*q[3]^4 + 1//24*q[2]^4*q[3]^2 + 1//24*q[2]^2*q[3]^4 + 115//6*q[3]^6
    @test feynman_integral_degree(F,3,g)==115//3*q[1]^6 + 1//4*q[1]^4*q[2]^2 + 1//4*q[1]^4*q[3]^2 + 1//4*q[1]^2*q[2]^4 + 1//2*q[1]^2*q[2]^2*q[3]^2 + 1//4*q[1]^2*q[3]^4 + 115//3*q[2]^6 + 1//4*q[2]^4*q[3]^2 + 1//4*q[2]^2*q[3]^4 + 115//3*q[3]^6
    @test feynman_integral_degree_sum_order(F,o,[1,2],g)==1//24*q[1]^2*q[2]^2 + 1//24*q[1]^2*q[3]^2 + 1//24*q[2]^2*q[3]^2 + 19//8*q[3]^4 + 1//24*q[3]^2
    @test feynman_integral_degree_sum(F,3,g)==115//3*q[1]^6 + 1//4*q[1]^4*q[2]^2 + 1//4*q[1]^4*q[3]^2 + 19//4*q[1]^4 + 1//4*q[1]^2*q[2]^4 + 1//2*q[1]^2*q[2]^2*q[3]^2 + 1//4*q[1]^2*q[2]^2 + 1//4*q[1]^2*q[3]^4 + 1//4*q[1]^2*q[3]^2 + 1//12*q[1]^2 + 115//3*q[2]^6 + 1//4*q[2]^4*q[3]^2 + 19//4*q[2]^4 + 1//4*q[2]^2*q[3]^4 + 1//4*q[2]^2*q[3]^2 + 1//12*q[2]^2 + 115//3*q[3]^6 + 19//4*q[3]^4 + 1//12*q[3]^2

    ve = [ (1, 2), (1, 2), (1, 2)]
    G=FeynmanGraph(ve)
    F=FeynmanIntegral(G)
    o=[1,2];
    g=[1,0]
    q=F.S[3]

    @test feynman_integral_degree_order(F,o,2,g)==5//12*q[1]^4 + 5//12*q[1]^2*q[2]^2 + 5//12*q[1]^2*q[3]^2 + 5//12*q[2]^4 + 5//12*q[2]^2*q[3]^2 + 5//12*q[3]^4
    @test feynman_integral_degree_order(F,o,1,g)==0
    @test feynman_integral_degree(F,2,g)==5//6*q[1]^4 + 5//6*q[1]^2*q[2]^2 + 5//6*q[1]^2*q[3]^2 + 5//6*q[2]^4 + 5//6*q[2]^2*q[3]^2 + 5//6*q[3]^4
    @test feynman_integral_degree(F,1,g)==0
    @test feynman_integral_degree(F,1)==0
    S=polynomial_ring(QQ, :x =>1:nv(G), :q =>1:ne(G), :z => 1:nv(G))
    FF=FeynmanIntegral(G,S)
    Fv=FeynmanIntegral(ve)
end