@testset "flip.jl" begin
    ve = [ (1, 2), (2, 3), (3, 1)]
    G=graph(ve)
    a=[0,0,3];
    l=[0,0,0];
    Ω=[1,2,3];
    o=[1,2,3];
    g=[1,0,0];
   
    q=Nemo.vars(feynman_integral_degree(G ,3))
    @test feynman_integral_branch_type(G,a,g)==115//3*q[3]^6
    @test feynman_integral_degree(G,3,g)==115//3*q[1]^6 + 1//4*q[1]^4*q[2]^2 + 1//4*q[1]^4*q[3]^2 + 1//4*q[1]^2*q[2]^4 + 1//2*q[1]^2*q[2]^2*q[3]^2 + 1//4*q[1]^2*q[3]^4 + 115//3*q[2]^6 + 1//4*q[2]^4*q[3]^2 + 1//4*q[2]^2*q[3]^4 + 115//3*q[3]^6
    @test feynman_integral_degree_sum(G,3,g)==115//3*q[1]^6 + 1//4*q[1]^4*q[2]^2 + 1//4*q[1]^4*q[3]^2 + 19//4*q[1]^4 + 1//4*q[1]^2*q[2]^4 + 1//2*q[1]^2*q[2]^2*q[3]^2 + 1//4*q[1]^2*q[2]^2 + 1//4*q[1]^2*q[3]^4 + 1//4*q[1]^2*q[3]^2 + 1//12*q[1]^2 + 115//3*q[2]^6 + 1//4*q[2]^4*q[3]^2 + 19//4*q[2]^4 + 1//4*q[2]^2*q[3]^4 + 1//4*q[2]^2*q[3]^2 + 1//12*q[2]^2 + 115//3*q[3]^6 + 19//4*q[3]^4 + 1//12*q[3]^2
   

    q=Nemo.vars(feynman_integral_degree_order(G,o,2,g))
    @test feynman_integral_degree_order(G,o,2,g)==1//24*q[1]^2*q[2]^2 + 1//24*q[1]^2*q[3]^2 + 1//24*q[2]^2*q[3]^2 + 19//8*q[3]^4
    @test feynman_integral_branch_type_order(G,a,Ω,g) == 115//6*q[3]^6
    @test feynman_integral_degree_order(G,o,3,g)==1//24*q[1]^4*q[2]^2 + 1//24*q[1]^4*q[3]^2 + 1//24*q[1]^2*q[2]^4 + 1//12*q[1]^2*q[2]^2*q[3]^2 + 1//24*q[1]^2*q[3]^4 + 1//24*q[2]^4*q[3]^2 + 1//24*q[2]^2*q[3]^4 + 115//6*q[3]^6
    @test feynman_integral_degree_sum_order(G,[1,2],o,g)==1//24*q[1]^2*q[2]^2 + 1//24*q[1]^2*q[3]^2 + 1//24*q[2]^2*q[3]^2 + 19//8*q[3]^4 + 1//24*q[3]^2
end