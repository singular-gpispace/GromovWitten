@testset "feynmanIntegral.jl" begin
    ve=[(1, 3), (1,2),(1, 2), (2, 4) ,(3, 4),(3,4)]
    G=FeynmanGraph(ve)
    F=FeynmanIntegral(G)
    q=F.S[3]
    a=[0,2,1,0,0,1]
    o=[1,3,4,2]
    Ω=[1,3,4,2];
    @test feynman_integral_branch_type_order(F,a,o)==128*q[2]^4*q[3]^2*q[6]^2
    @test feynman_integral_branch_type(F,a)==256*q[2]^4*q[3]^2*q[6]^2
    @test feynman_integral_degree(F ,2)==8*q[1]^4 + 8*q[2]^2*q[3]^2 + 8*q[4]^4 + 8*q[5]^2*q[6]^2
    @test sum_of_coeff( feynman_integral_degree_order(F,o,4))==3088
    @test sum_of_coeff(feynman_integral_degree(F,4))==25344
    @test feynman_integral_degree_sum_order(F,Ω,[1,3])==12*q[2]^6 + 76*q[2]^4*q[3]^2 + 4*q[2]^4*q[5]^2 + 4*q[2]^4*q[6]^2 + 76*q[2]^2*q[3]^4 + 16*q[2]^2*q[3]^2*q[5]^2 + 16*q[2]^2*q[3]^2*q[6]^2 + 4*q[2]^2*q[3]^2 + 12*q[3]^6 + 4*q[3]^4*q[5]^2 + 4*q[3]^4*q[6]^2
    @test feynman_integral_degree_sum_order(F,Ω,3)==12*q[2]^6 + 76*q[2]^4*q[3]^2 + 4*q[2]^4*q[5]^2 + 4*q[2]^4*q[6]^2 + 76*q[2]^2*q[3]^4 + 16*q[2]^2*q[3]^2*q[5]^2 + 16*q[2]^2*q[3]^2*q[6]^2 + 4*q[2]^2*q[3]^2 + 12*q[3]^6 + 4*q[3]^4*q[5]^2 + 4*q[3]^4*q[6]^2
    @test feynman_integral_degree_sum(F,[1,2])==8*q[1]^4 + 8*q[2]^2*q[3]^2 + 8*q[4]^4 + 8*q[5]^2*q[6]^2
    @test feynman_integral_degree_sum_order(F,Ω,[1,2])==4*q[2]^2*q[3]^2
    @test substitute(feynman_integral_degree_sum(F,4))==25344*q[1]^8 + 1792*q[1]^6 + 32*q[1]^4
    @test  substitute(0)==0
    p=6*q[1]^2*q[2]^2*q[3]^4 + 12*q[1]^2*q[2]^2*q[3]^2*q[4]^2 + 6*q[1]^2*q[2]^2*q[4]^4 + 56*q[1]^2*q[3]^6 + 6*q[1]^2*q[3]^4*q[4]^2 + 6*q[1]^2*q[3]^2*q[4]^4 + 56*q[1]^2*q[4]^6 + q[1]
    @test  filter_term(p,q[1],1)==q[1]
end