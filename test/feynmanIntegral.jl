@testset "feynmanIntegral.jl" begin
    ve=[(1, 3), (1,2),(1, 2), (2, 4) ,(3, 4),(3,4)]
    G=graph(ve)
   # S,x,q,z=polynomial_ring(G,"x","q","z")
    a=[0,2,1,0,0,1]
    o=[1,3,4,2]
    Ω=[1,3,4,2];
    using Nemo

    q=Nemo.vars(feynman_integral_degree(G ,2))
    @test feynman_integral_degree(G ,2)==8*q[1]^4 + 8*q[2]^2*q[3]^2 + 8*q[4]^4 + 8*q[5]^2*q[6]^2
    @test feynman_integral_branch_type(G,a)==256*q[2]^4*q[3]^2*q[6]^2
    @test sum_of_coeff(feynman_integral_degree(G,4))==25344
    @test substitute(feynman_integral_degree_sum(G,4))==25344*q[1]^8 + 1792*q[1]^6 + 32*q[1]^4
    @test feynman_integral_degree_sum(G,[1,2])==8*q[1]^4 + 8*q[2]^2*q[3]^2 + 8*q[4]^4 + 8*q[5]^2*q[6]^2

    q=Nemo.vars(feynman_integral_degree_order(G,o,4))

    @test feynman_integral_branch_type_order(G,a,o)== 128*q[2]^4*q[3]^2*q[6]^2
    @test sum_of_coeff( feynman_integral_degree_order(G,o,4))==3088
    @test feynman_integral_degree_sum_order(G,[1,3],Ω)==12*q[2]^6 + 76*q[2]^4*q[3]^2 + 4*q[2]^4*q[5]^2 + 4*q[2]^4*q[6]^2 + 76*q[2]^2*q[3]^4 + 16*q[2]^2*q[3]^2*q[5]^2 + 16*q[2]^2*q[3]^2*q[6]^2 + 4*q[2]^2*q[3]^2 + 12*q[3]^6 + 4*q[3]^4*q[5]^2 + 4*q[3]^4*q[6]^2
    @test feynman_integral_degree_sum_order(G,3,o)==12*q[2]^6 + 76*q[2]^4*q[3]^2 + 4*q[2]^4*q[5]^2 + 4*q[2]^4*q[6]^2 + 76*q[2]^2*q[3]^4 + 16*q[2]^2*q[3]^2*q[5]^2 + 16*q[2]^2*q[3]^2*q[6]^2 + 4*q[2]^2*q[3]^2 + 12*q[3]^6 + 4*q[3]^4*q[5]^2 + 4*q[3]^4*q[6]^2
    @test feynman_integral_degree_sum_order(G,[1,2],Ω)==4*q[2]^2*q[3]^2
    #@test  substitute(0)==T(0)
   
end